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

import math
import re
from subprocess import call

import fiona
import numpy as np
import shapely.geometry
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta
from osgeo import osr
from pyproj import Proj, transform

from priima.config import Config
from priima.geo_tools import get_half_region_size
from priima.projection import get_projection, get_projection_epsg


def reproject2roi(filename):
    """
    Reproject the images to the relevant polar stereographic projection

    :param filename: Filename including absolute path
    :type filename: string
    """
    warp_command_string = warp_command(filename)
    call(warp_command_string, shell=True)
    print(warp_command_string)

    return filename.with_stem(f"{filename.stem}_roi")


def warp_command(filename):
    epsg_code = get_projection_epsg()
    proj = Proj(epsg_code)
    center_lat, center_lon = Config.instance().center[:]
    center_x_m, center_y_m = proj(center_lon, center_lat)
    half_region_size = get_half_region_size()
    ll_x_m = center_x_m - half_region_size
    ll_y_m = center_y_m - half_region_size
    ur_x_m = center_x_m + half_region_size
    ur_y_m = center_y_m + half_region_size
    output_name = filename.with_stem(f"{filename.stem}_roi")
    resolution = Config.instance().output_resolution

    warp_command = (
        f"gdalwarp -overwrite "
        f"-tr {resolution} {resolution} -r near -multi -srcnodata 0 "
        f"-dstnodata 0 -order 3 -t_srs '+init=epsg:{epsg_code}' "
        f"-te {ll_x_m:.1f} {ll_y_m:.1f} {ur_x_m:.1f} {ur_y_m:.1f} "
        f"{filename} {output_name}"
    )

    return warp_command


def compute_time_range(filename, forecast_duration):
    """
    Computes time range to select appropriate drift data

    :param filename: Filename including path of the image to be forecast
    :type filename: string

    :param forecast_duration: Time step for the forecast in hours
    :type forecast_duration: integer
    """
    match = re.search(
        r"_(?P<start_time>\d{8}T\d{2})\d{4}_", str(filename)
    )
    start_time = parse(match.group('start_time'))

    # time range only reported with hourly resolution: minutes and seconds
    # are ignored in the time range datetime elements.
    time_range = [
        start_time,
        start_time + relativedelta(hours=forecast_duration)
    ]

    return time_range


def compute_displacement(drift, forecast_step):
    """
    Computes the ice displacement by multiplying the drift which is given
    in m s^-1 with the forecast step given in hours
    """
    drift_in_meters = np.array(drift)*forecast_step*60*60

    # simulate a hypothetical Northwards drift
    # drift_in_meters = np.array(drift)

    return drift_in_meters


def update_point_location(gcp_list_dynamic, drift):
    gcp_list_updated = []
    if isinstance(drift[0], list):
        udrift = np.array(drift[0])
        vdrift = np.array(drift[1])
    else:
        udrift = drift[0]
        vdrift = drift[1]
    udrift[np.isnan(udrift)] = 0
    vdrift[np.isnan(vdrift)] = 0
    ind = 0
    outProj = get_projection()
    inProj = Proj(init='epsg:4326')

    for gcp in gcp_list_dynamic:
        drift_point = [udrift[ind], vdrift[ind]]
        ind = ind + 1
        lon = gcp.GCPX
        lat = gcp.GCPY

        lon_projected, lat_projected = gcp.GCPXX, gcp.GCPYY

        point = (lon, lat)
        if gcp.Info == 'is_land':
            gcp_list_updated.append(gcp)
            continue
        else:
            is_inland = is_point_inland(point)

        if is_inland:
            displacement = (0, 0)
            gcp.Info = 'is_land'
        else:
            # drift_point = correct_icon_to_true_north(drift_point)
            displacement = compute_displacement(drift_point, forecast_step=1)

        lon_projected = lon_projected + displacement[0]
        lat_projected = lat_projected + displacement[1]

        lon_drifted, lat_drifted = transform(outProj, inProj,
                                             lon_projected, lat_projected)

        if isinstance(lon_drifted, float):
            gcp.GCPX = lon_drifted
            gcp.GCPY = lat_drifted
        else:
            gcp.GCPX = lon_drifted[0][0]
            gcp.GCPY = lat_drifted[0][0]

        if isinstance(lon_projected, float):
            gcp.GCPXX = lon_projected
            gcp.GCPYY = lat_projected
        else:
            gcp.GCPXX = lon_projected[0][0]
            gcp.GCPYY = lat_projected[0][0]
        gcp_list_updated.append(gcp)

    return gcp_list_updated


def correct_icon_to_true_north(drift):
    l0_grid_central_meridian = 0
    point_lon = -52
    point_lat = -68
    grid_convergence_angle = \
        np.arctan(np.tan(point_lon -
                  l0_grid_central_meridian)*np.sin(point_lat))
    grid_convergence_angle = 360 - math.degrees(grid_convergence_angle)
    u_drift = drift[0]
    v_drift = drift[1]
    grid_convergence_radians = math.radians(grid_convergence_angle)

    u_drift = math.sin(grid_convergence_radians)
    v_drift = math.cos(grid_convergence_radians)

    return np.array([u_drift, v_drift])


def compute_averaged_area_drift(drift):
    udrift = drift[0]
    vdrift = drift[1]
    udrift_aver = udrift.mean()
    vdrift_aver = vdrift.mean()
    drift_magnitude = math.sqrt(udrift_aver**2 + vdrift_aver**2)
    drift_direction = \
        (90 - math.atan2(vdrift_aver, udrift_aver)*180/math.pi) % 360

    return [drift_magnitude, drift_direction]


def is_point_inland(point):
    """
    Determins if a point belongs to the land.
    """
    if point[1] > 0:
        shapefile = "shapefiles/arctic_landmask.shp"
    else:
        shapefile = "shapefiles/antarctic_landmask.shp"

    point = shapely.geometry.Point(point[0], point[1])  # lon, lat

    with fiona.open(shapefile) as fiona_fl:
        for shapefile_record in fiona_fl:
            shape = shapely.geometry.shape(shapefile_record["geometry"])

            if shape.contains(point):
                is_inland = True
                break
            else:
                is_inland = False

    return is_inland


def plot_polygon(poly, points):
    import matplotlib.pyplot as plt
    from descartes.patch import PolygonPatch
    from shapely.geometry import Polygon

    plt.figure()
    ax = plt.axes()
    ax.set_aspect('equal')
    polygon = Polygon(poly)
    patch = PolygonPatch(polygon, facecolor=[0, 0, 0.5], edgecolor=[0, 0, 0],
                         alpha=0.7, zorder=2)
    ax.add_patch(patch)
    plt.plot(points, 'r+')
    plt.xlim(10.5, 12)
    plt.ylim(78.8, 79.5)
    plt.savefig("polygon.png")


def create_transformation_matrix(dataset, gcp_list_initial, gcp_list_ending):
    xpix_drifted = []
    ypix_drifted = []
    lon_projected_list = []
    lat_projected_list = []
    starting_pos_test = []
    ending_pos_test = []

    for point_ind in range(len(gcp_list_initial)):
        gcp_end = gcp_list_ending[point_ind]
        point = [gcp_end.GCPXX, gcp_end.GCPYY]
        gcp = gcp_list_initial[point_ind]

        xpix_tmp, ypix_tmp = compute_pixel_coordinates(dataset, point, gcp)
        xpix_drifted.append(xpix_tmp)
        ypix_drifted.append(ypix_tmp)

        inProj = Proj(init='epsg:4326')
        proj = osr.SpatialReference()
        proj.ImportFromWkt(dataset.GetProjection())

        outProj = Proj(proj.ExportToProj4())
        lon_drifted = gcp_list_ending[point_ind].GCPX
        lat_drifted = gcp_list_ending[point_ind].GCPY

        lon_projected, lat_projected = transform(
            inProj, outProj, lon_drifted, lat_drifted)
        lon_projected_list.append(lon_projected)
        lat_projected_list.append(lat_projected)
        starting_pos_test.append((1*gcp.GCPPixel, gcp.GCPLine))
        ending_pos_test.append((1*gcp_end.GCPPixel, gcp_end.GCPLine))

    xpix_initial = list(
        range(0, dataset.RasterXSize, Config.instance().gcp_separation))
    xpix_initial.append(dataset.RasterXSize)
    ypix_initial = list(
        range(0, dataset.RasterYSize, Config.instance().gcp_separation))
    ypix_initial.append(dataset.RasterXSize)
    starting_pos = [(xi, yi) for xi in xpix_initial for yi in ypix_initial]

    ending_pos = list(zip(xpix_drifted, ypix_drifted))
    ending_pos = np.array(starting_pos) + np.array(ending_pos)

    return starting_pos, ending_pos, lon_projected_list, lat_projected_list


def compute_pixel_coordinates(dataset, point, gcp):
    """
    Dataset loaded from gdal.Open method
    """
    ds = dataset
    image_dimensions = (ds.RasterXSize, ds.RasterYSize)

    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromWkt(ds.GetProjection())
    proj4_string = spatial_ref.ExportToProj4()

    x_min, x_size, _, y_min, _, y_size = ds.GetGeoTransform()

    buffersize = 2000

    x_grid_pos = np.arange(
        0, (image_dimensions[0] + buffersize/4. + 1)*x_size, x_size)
    x_grid_neg = np.arange((buffersize/4. + 1)*x_size*(-1), 0, x_size)
    y_grid_pos = np.arange(
        0, (image_dimensions[0] + buffersize/4. + 1)*abs(y_size), abs(y_size))
    y_grid_neg = np.arange((buffersize/4. + 1)*(y_size), 0, abs(y_size))

    proj = Proj(proj4_string)

    xproj_pos = x_grid_pos + gcp.GCPXX
    xproj_neg = x_grid_neg + gcp.GCPXX
    yproj_pos = y_grid_pos + gcp.GCPYY
    yproj_neg = y_grid_neg + gcp.GCPYY

    lons_pos, lats_pos = proj(xproj_pos, yproj_pos, inverse=True)
    lons_neg, lats_neg = proj(xproj_neg, yproj_neg, inverse=True)
    if xproj_neg.min() <= point[0] <= xproj_neg.max():
        x_index = np.abs(xproj_neg - point[0]).argmin()
        x_index = -1*(len(xproj_neg) - x_index)
    else:
        x_index = 1*np.abs(xproj_pos - point[0]).argmin()

    if yproj_neg.min() <= point[1] <= yproj_neg.max():
        y_index = np.abs(yproj_neg - point[1]).argmin()
        y_index = 1*(len(yproj_neg) - y_index)
    else:
        y_index = -1*np.abs(yproj_pos - point[1]).argmin()

    return x_index, y_index
