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

# plot drift for debugging

import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from pyproj import Proj, transform

# from priima.wind import select_time_window


def create_drift_plot(filename, data_dict, time_range):
    # filename = 'data/ice_drift/'
    #            '20180310_hr-metno-MODEL-topaz4-ARC-b20180309-fv02.0_roi.nc'
    # filename = 'data/ice_drift/'
    #            '20180320_hr-metno-MODEL-topaz4-ARC-b20180319-fv02.0_roi.nc'

    # time_range = time_range = [datetime.datetime(2018, 03, 10, 5, 0),
    #                            datetime.datetime(2018, 03, 10, 8, 0)]
    # time_range = time_range = [datetime.datetime(2018, 03, 20, 6, 0),
    #                            datetime.datetime(2018, 03, 20, 12, 20)]

    # data_dict = read_topaz(filename, time_range)
    lon_drifted, lat_drifted, lon_drifted_mean, lat_drifted_mean = \
        compute_drift(data_dict)
    drift_points = {'lon_0': data_dict['lon'], 'lat_0': data_dict['lat'],
                    'lon_1': lon_drifted, 'lat_1': lat_drifted,
                    'lon_1_mean': lon_drifted_mean,
                    'lat_1_mean': lat_drifted_mean}
    plot_drift_on_map(filename, drift_points, data_dict)


def read_topaz(filename, time_range):
    dataset = Dataset(filename)
    dataset_time = dataset.variables['time'][:]
    time_window_index = select_time_window(
        dataset_time, time_range, data_type='topaz')

    lon = dataset.variables['longitude'][:]
    lat = dataset.variables['latitude'][:]
    uice = dataset.variables['uice'][time_window_index, :, :][::-1]
    vice = dataset.variables['vice'][time_window_index, :, :][::-1]
    # u_scale_factor = dataset.variables['uice'].scale_factor
    # v_scale_factor = dataset.variables['vice'].scale_factor
    # u_offset = dataset.variables['uice'].add_offset
    # v_offset = dataset.variables['vice'].add_offset

    # uice = uice#*u_scale_factor*10000 + u_offset
    # vice = -1*vice#*v_scale_factor*1000 + v_offset

    return {'lon': lon, 'lat': lat, 'uice': uice, 'vice': vice}


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


def compute_drift(data_dict):
    data_dict['uice'] = data_dict['uice'].mean(axis=0)
    data_dict['vice'] = data_dict['vice'].mean(axis=0)
    # compute x from v since it needs to be flipped again
    # (already done in topaz code)
    x_displacement = data_dict['uice'][:, :]*1*60*60
    y_displacement = data_dict['vice'][:, :]*1*60*60

    x_displacement_mean = data_dict['uice'].mean()*1*60*60
    y_displacement_mean = -1*data_dict['vice'].mean()*1*60*60
    print('for plotting: xy x_displacement: ',
          x_displacement_mean, y_displacement_mean)

    proj3995 = Proj("+init=EPSG:3995")
    wgs84 = Proj("+init=EPSG:4326")
    xx, yy = proj3995(data_dict['lon'], data_dict['lat'])
    xx = xx
    yy = yy
    xx = xx + x_displacement
    yy = yy + y_displacement
    lon_drifted, lat_drifted = transform(proj3995, wgs84, xx, yy)
    lon_drifted = np.ma.masked_where(xx.mask, lon_drifted)
    lat_drifted = np.ma.masked_where(yy.mask, lat_drifted)
    print(lon_drifted[0, 0])
    print(data_dict['lon'][0, 0])
    xx_mean, yy_mean = proj3995(data_dict['lon'], data_dict['lat'])
    xx_mean = xx_mean + x_displacement_mean
    yy_mean = yy_mean + y_displacement_mean
    lon_drifted_mean, lat_drifted_mean = \
        transform(proj3995, wgs84, xx_mean, yy_mean)

    return lon_drifted, lat_drifted, lon_drifted_mean, lat_drifted_mean


def plot_drift_on_map(filename, drift_points, data_dict):
    plt.figure(figsize=(12, 12))
    lat_ts = drift_points['lat_0'].mean()
    lon_0 = drift_points['lon_0'].mean()

    stere_centralized = ccrs.Stereographic(
        central_latitude=lat_ts, central_longitude=lon_0)
    ax = plt.axes(projection=stere_centralized)
    ax.set_extent((-600000, 600000, -600000, 600000), crs=stere_centralized)

    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'land', scale='50m', edgecolor='k',
        facecolor=cfeature.COLORS['land']), facecolor='#cc9966')
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'lakes', scale='50m', edgecolor='k',
        facecolor=cfeature.COLORS['water']), facecolor='#99ffff')

    ax.gridlines()

    ziped0 = zip(drift_points['lon_0'].ravel(), drift_points['lat_0'].ravel())
    ziped0 = list(ziped0)
    ziped1 = zip(drift_points['lon_1'].ravel(), drift_points['lat_1'].ravel())
    ziped1 = list(ziped1)
    for ind in range(0, len(ziped1), 1):

        tmpziped0 = ziped0[ind][0]
        tmpziped1 = ziped1[ind][0]

#        xlon0, ylat0 = m(ziped0[ind][0], ziped0[ind][1])
#        xlon1, ylat1 = m(ziped1[ind][0], ziped1[ind][1])
        # print("lat/lon: ",
        #       tmpziped0, ziped0[ind][1], tmpziped1, ziped1[ind][1])
        if not np.ma.is_masked(tmpziped0) or np.ma.is_masked(tmpziped1):
            xlon0, ylat0 = stere_centralized.transform_point(
                tmpziped0, ziped0[ind][1], ccrs.PlateCarree())
            xlon1, ylat1 = stere_centralized.transform_point(
                tmpziped1, ziped1[ind][1], ccrs.PlateCarree())

            xlond = 1.0*xlon1 - xlon0
            ylatd = 1.0*ylat1 - ylat0

        # xlon1 = data_dict["uice"].ravel()[ind]
        # xlat1 = data_dict["vice"].ravel()[ind]
            # print(xlon0, ylat0, xlon1, ylat1)
            if not any([np.ma.is_masked(tmpziped0),
                        np.ma.is_masked(tmpziped1)]):
                plt.arrow(xlon0, ylat0, xlond, ylatd, fc="k", ec="k",
                          linewidth=1, head_width=5000, head_length=12000)

        # m.quiver(xlon0, ylat0, xlon1, ylat1, angles='xy')

    plt.savefig("skata.png")
    # plt.clf()

    xlon0_mean, ylat0_mean = stere_centralized.transform_point(
        drift_points['lon_0'].mean(), drift_points['lat_0'].mean(),
        ccrs.PlateCarree())
    xlon1_mean, ylat1_mean = stere_centralized.transform_point(
        drift_points['lon_1_mean'].mean(),
        drift_points['lat_1_mean'].mean(),
        ccrs.PlateCarree())

    ax.plot(xlon0_mean, ylat0_mean, "x",
            markersize=5, markerfacecolor='black')
    ax.plot(xlon1_mean, ylat1_mean, "x",
            markersize=5, markerfacecolor='black')

    xlon1_mean = 1.0*xlon1_mean - xlon0_mean
    ylat1_mean = 1.0*ylat1_mean - ylat0_mean

    plt.arrow(xlon0_mean, ylat0_mean, xlon1_mean, ylat1_mean, fc="r", ec="r",
              linewidth=2, head_width=30000, head_length=40000)

#    for ind in range(0, len(ziped1)):
#        m.drawgreatcircle(ziped0[ind][0], ziped0[ind][1],
#                          ziped1[ind][0], ziped1[ind][1],
#                          linewidth=3,color='b')
#        xlon, ylat = m(ziped1[ind][0], ziped1[ind][1])
#        nyc = m.plot(xlon, ylat, 'ro')
#        plt.setp(nyc,'markersize',0.5,'markeredgecolor','k')
    outname = filename.replace('.nc', '.png')
    plt.savefig(outname)
    plt.clf()


if __name__ == "__main__":
    # create_drift_plot(filename, time_range)

    # example usage: choose a TOPAZ4 file and
    # time range that is contained in that file
    filename = ('/home/bahlmann/data/L2/IceDrift/TopazHourly/'
                '20220921_hr-metno-MODEL-topaz4-ARC-b20220920-fv02.0.nc')
    time_range = time_range = [datetime.datetime(2022, 9, 21, 5, 0),
                               datetime.datetime(2022, 9, 21, 8, 0)]
    data_dict = read_topaz(filename, time_range)

    create_drift_plot(filename, data_dict, time_range)
