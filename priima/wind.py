"""
Copyright 2025, Drift+Noise GmbH

This file is part of PRIIMA.
PRIIMA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.
"""

import datetime
import math
from subprocess import call

# from matplotlib.cm import get_cmap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal, ogr, osr

from priima.drift_handler import DriftHandler
from priima.icon import convert_icon_time_to_datetime, process_icon_fields


def compute_drift_from_wind(config, time_range, data_source, order,
                            gcp_list, tide_data):
    """
    Computes the drift vectors
    """
    data_path = config.paths.data_path
    if data_source == "ICON":
        wind_dict, wind_file = read_wind_from_icon(
                time_range, data_path, order, config, tide_data)
    else:
        wind_dict = read_wind(data_path, time_range)

    wind_dict['time'] = convert_icon_time_to_datetime(wind_dict)

    drift_handler = DriftHandler(wind_file, time_range, order,
                                 gcp_list)
    wind_dict = \
        drift_handler.subset_roi_from_meshgrid(wind_dict)
    vice, uice = \
        drift_handler.resample_drift(wind_dict)
    u_drift, v_drift = wind2drift(uice, vice, order)
    # plot_icon_wind(wind_dict, order)
    # test if winddrift is the expected one
    uicewind = wind_dict['uice'].flatten()
    vicewind = wind_dict['vice'].flatten()
    u_winddrift, v_winddrift = wind2drift(uicewind, vicewind, order)
    u_winddrift = np.reshape(u_winddrift, wind_dict['uice'].shape)
    v_winddrift = np.reshape(v_winddrift, wind_dict['vice'].shape)
    winddrift_dict = wind_dict
    winddrift_dict['uice'] = u_winddrift
    winddrift_dict['vice'] = v_winddrift
    if order.plot_wind:
        plot_icon_wind(winddrift_dict, order)

    return u_drift, v_drift


def read_wind_from_icon(time_range, data_path, order, config, tide_data):
    wind_file = process_icon_fields(time_range, data_path, order, config,
                                    tide_data)
    wind_dict = prepare_wind_dict_from_icon(wind_file)

    return wind_dict, wind_file


def prepare_wind_dict_from_icon(path):
    """
    Reads and returns the ICON wind fields given in
    netCDF format. There is no need to select timestemps
    since this is already done from the file selection
    in icon.select_icon_files
    """
    dataset = Dataset(path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    xx = dataset.variables['x'][:]*100000
    yy = dataset.variables['y'][:]*100000
    u_wind = dataset.variables['10u'][:, 0, :, :]
    v_wind = dataset.variables['10v'][:, 0, :, :]
    time = dataset.variables['time'][:]
    time_since = dataset.variables['time'].units

    wind_dict = {"lon": lon, "lat": lat,
                 "xx": xx, "yy": yy,
                 "uice": u_wind, "vice": v_wind,
                 "u_offset": 0,
                 "u_scale": 1,
                 "v_offset": 0,
                 "v_scale": 1,
                 "time": time,
                 "time_since": time_since
                 }

    return wind_dict


def read_wind(path, time_range):
    """
    Reads and returns the ensemble wind fields given in
    netCDF format and dimensions within a defined time slot:
    [time, ensemble, lats, lons]

    :param path: Full path to the netCDF file
    :type path: string

    :param time_range: range between data acquisition and desired forecast
    :type time_range: list of two datetime objects
    """
    dataset = Dataset(path)
    dataset_time = dataset.variables['time'][:]
    time_window_index = select_time_window(dataset_time, time_range,
                                           data_type='ecmwf')

    lon = dataset.variables['longitude'][:]
    lat = dataset.variables['latitude'][:]
    u_wind = dataset.variables['u10'][time_window_index, :, :, :]
    u_wind_scale = dataset.variables['u10'].scale_factor
    u_wind_offset = dataset.variables['u10'].add_offset

    v_wind = dataset.variables['v10'][time_window_index, :, :, :]
    v_wind_scale = dataset.variables['v10'].scale_factor
    v_wind_offset = dataset.variables['v10'].add_offset

    wind_dict = {"longitude": lon, "latitude": lat,
                 "u_wind": u_wind, "v_wind": v_wind,
                 "u_offset": u_wind_offset,
                 "u_scale": u_wind_scale,
                 "v_offset": v_wind_offset,
                 "v_scale": v_wind_scale
                 }

    return wind_dict


def select_time_window(dataset_time, time_range, data_type='topaz'):
    """
    Extracts the relevant time window index between the actual image
    acquisition and the desired forecast.

    :param time: time given as 'hours since 1900-01-01 00:00:0.0'
    :type time: numpy.ndarray

    :param time_range: range between data acquisition and desired forecast
    :type time_range: list of two datetime objects
    """
    start = time_range[0]
    end = time_range[1]
    if data_type == 'topaz':
        epoch = datetime.datetime(1950, 1, 1, 0, 0, 0)
    elif data_type == 'ecmwf':
        epoch = datetime.datetime(1900, 1, 1, 0, 0, 0)
    else:
        print('Unknown data_type string. Should be "topaz" or ecmwf.'
              'Exiting..')
        exit(111)

    start_hours = (start - epoch).days * 24 + start.hour
    end_hours = (end - epoch).days * 24 + end.hour
    time_window_index = np.where(np.logical_and(dataset_time >= start_hours,
                                                dataset_time <= end_hours))[0]

    return time_window_index


def compute_wind(wind_dict):
    """
    Averages across the ensemble members and then applies
    the scaling factor and the offset to compute the
    actuall wind field.
    :param wind_dict: Dictionary containing the wind fields and
                      their attributes.
    :type wind_dict: Dictionary
    """
    u_wind_mean = wind_dict["u_wind"].mean(axis=1)
    u_wind_mean = u_wind_mean*wind_dict["u_scale"] + \
        wind_dict["u_offset"]

    v_wind_mean = wind_dict["v_wind"].mean(axis=1)
    v_wind_mean = v_wind_mean*wind_dict["v_scale"] + \
        wind_dict["v_offset"]

    return [u_wind_mean, v_wind_mean]


def wind2drift(u_wind, v_wind, order):
    """
    Converts wind vectors to drift vectors following the rule of thumb
    that drift velocity is on average 2.5% of the wind velocity with a
    direction 20 degrees right to the wind direction.More info for the rotation
    here: https://www.eol.ucar.edu/content/wind-direction-quick-reference
    For negative latitude we subtract 25deg from 360 to rotate the vector to
    the left as this is the relative motion expected at the Southern
    hemisphere.
    """
    u_drift_all = []
    v_drift_all = []
    if np.isscalar(u_wind):
        u_wind = [u_wind]
        v_wind = [v_wind]
    if order.center[0] > 0:
        rotation_angle = 25
    else:
        rotation_angle = 335

    drift_phi_radians = math.radians(rotation_angle)
    wind_to_drift_scale = 0.025

    for ii in range(len(v_wind)):
        u_drift = \
            math.cos(drift_phi_radians)*u_wind[ii]*wind_to_drift_scale + \
            math.sin(drift_phi_radians)*v_wind[ii]*wind_to_drift_scale
        v_drift = \
            -math.sin(drift_phi_radians)*u_wind[ii]*wind_to_drift_scale + \
            math.cos(drift_phi_radians)*v_wind[ii]*wind_to_drift_scale

        u_drift_all.append(u_drift)
        v_drift_all.append(v_drift)

    return u_drift_all, v_drift_all


def prepare_drift_from_wind_for_qgis(wind_dict):
    driver = gdal.GetDriverByName('GTiff')
    u_wind = wind_dict['uice'].mean(axis=0)  # [0,:,:]
    v_wind = wind_dict['vice'].mean(axis=0)  # [0,:,:]
    u_wind = u_wind[::-1]
    v_wind = v_wind[::-1]

    uout_name = 'uice.tif'
    vout_name = 'vice.tif'
    udataset_out = driver.Create(uout_name, u_wind.shape[1],
                                 u_wind.shape[0], 1, gdal.GDT_Float32)
    vdataset_out = driver.Create(vout_name, v_wind.shape[1],
                                 v_wind.shape[0], 1, gdal.GDT_Float32)

    lonmaxpj, latmaxpj = wind_dict['lon'].max(), wind_dict['lat'].max()
    lonminpj, latminpj = wind_dict['lon'].min(), wind_dict['lat'].min()

    proj = osr.SpatialReference()
    # proj.ImportFromProj4(
    #    "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-45 +k=1 "
    #    "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    # proj.ImportFromEPSG(3995)
    proj.ImportFromEPSG(4326)

    # inSpatialRef = osr.SpatialReference()
    # inSpatialRef.ImportFromEPSG(3857)
    # coordTransform = osr.CoordinateTransformation(inSpatialRef, proj)

    pointmax = ogr.Geometry(ogr.wkbPoint)
    pointmin = ogr.Geometry(ogr.wkbPoint)
    pointmax.AddPoint(lonmaxpj, latmaxpj)
    pointmin.AddPoint(lonminpj, latminpj)
    # pointmax.Transform(coordTransform)
    # pointmin.Transform(coordTransform)

    lonminpj = pointmin.GetX()
    latmaxpj = pointmax.GetY()

    geotransformation = np.array([lonminpj, 0.125, 0, latmaxpj, 0, -0.125])

    udataset_out.SetGeoTransform(geotransformation)
    vdataset_out.SetGeoTransform(geotransformation)

    udataset_out.SetProjection(proj.ExportToWkt())
    vdataset_out.SetProjection(proj.ExportToWkt())

    udataset_out.GetRasterBand(1).WriteArray(u_wind[:, :])
    vdataset_out.GetRasterBand(1).WriteArray(v_wind[:, :])

    udataset_out.FlushCache()
    vdataset_out.FlushCache()
    udataset_out = None
    vdataset_out = None
    uout2 = uout_name.replace(".tif", "projected.tif")
    vout2 = vout_name.replace(".tif", "projected.tif")
    uout = uout_name
    vout = vout_name

    call("gdalwarp -s_srs epsg:4326 "
         "-t_srs epsg:3995 {uin} {uout}".format(uin=uout, uout=uout2),
         shell=True)
    call("gdalwarp -s_srs epsg:4326 "
         "-t_srs epsg:3995 {vin} {vout}".format(vin=vout, vout=vout2),
         shell=True)


def plot_icon_wind(wind_dict, order):

    uwind = wind_dict['uice'].mean(axis=0)[:, :]
    vwind = wind_dict['vice'].mean(axis=0)[:, :]
    lon = wind_dict['lon']
    lat = wind_dict['lat']
    endurance_location = [-52.3, -68.8]

    windspeed = (uwind ** 2 + vwind ** 2) ** 0.5
    plt.figure(figsize=(15, 15))
    lat_ts = order.center[0]
    lon_0 = order.center[1]
    stere_centralized = ccrs.Stereographic(
        central_latitude=lat_ts, central_longitude=lon_0)
    ax = plt.axes(projection=stere_centralized)
    ax.set_extent(
        (-1900000, 1900000, -1900000, 1900000), crs=stere_centralized)

    date_str = wind_dict['time'][0].strftime("%Y/%m/%d %H:%M")
    date_str_for_file = wind_dict['time'][0].strftime("%Y_%m_%d_%H")
    plt.title(date_str)

    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'land', scale='50m', edgecolor='k',
        facecolor=cfeature.COLORS['land']), facecolor='#cc9966')
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'lakes', scale='50m', edgecolor='k',
        facecolor=cfeature.COLORS['water']), facecolor='#99ffff')
    ax.gridlines()

    # plot static point
    x_static, y_static = stere_centralized.transform_point(
        endurance_location[0], endurance_location[1], ccrs.PlateCarree())
    ax.plot(x_static, y_static, "x", markersize=5, markerfacecolor='black')

    yy = np.arange(0, lat.shape[0], 8)
    xx = np.arange(0, lon.shape[0], 8)
    points = np.meshgrid(yy, xx)

    if len(lon.shape) == 1:
        lonm, latm = np.meshgrid(lon, lat)
    else:
        lonm, latm = lon, lat

    transformed = stere_centralized.transform_points(
        ccrs.PlateCarree(), lonm, latm)
    lonmx = transformed[:, :, 0]
    latmy = transformed[:, :, 1]

    ax.contourf(lonmx, latmy, windspeed)
    # normalize
    # uwind = uwind / np.sqrt(uwind ** 2.0 + vwind ** 2.0)
    # vwind = vwind / np.sqrt(uwind ** 2.0 + vwind ** 2.0)

    data_range = np.nanmax(windspeed) - np.nanmin(windspeed)
    full = data_range/2.
    half = full/2.
    flag = full*5

    ind = np.ix_(points[0][0], points[1][:, 0])
#    plt.quiver(lonmx[ind],
#               latmy[ind],
#               uwind[ind],
#               vwind[ind],
#               windspeed[ind],
#               scale=0.001)

    barbs = plt.barbs(lonmx[ind],
                      latmy[ind],
                      uwind[ind],
                      vwind[ind],
                      windspeed[ind],
                      flagcolor='r', barbcolor=['black'],
                      barb_increments=dict(half=half, full=full, flag=flag),
                      length=5)

    # levels = [-10, -5, 0, 5, 10]
    # wspd_contours = plt.contourf(
    #     lonmx[ind], latmy[ind], windspeed[ind], levels=levels,
    #     cmap=get_cmap("rainbow"))
    cbar = plt.colorbar(barbs, ax=ax, orientation="horizontal", pad=.05)

    cbar.set_label('wind in m/s', rotation=0)

    plt.savefig('icon_wind_' + date_str_for_file)
    plt.close()
