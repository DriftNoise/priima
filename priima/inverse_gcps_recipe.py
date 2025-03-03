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

"""
Example script to compute geographical coordinates for image pixels
based on a known projection, image dimension and pixel size

Author: Stefan Hendricks
"""

import csv
from subprocess import call

import numpy as np
import shapefile
from osgeo import gdal, ogr, osr
from pyproj import Proj

from priima.config import Config


def inversion_gcps_recipe(ds, proj, image_order, method='mpl'):
    """
    Dataset loaded from gdal.Open method
    """

    # Lets assume we have a 10k x 10k image
    image_dimensions = (ds.RasterXSize, ds.RasterYSize)

    # in the following projection
    # (for the sake of simplicity I am using a proj4 str here,
    # but this should be easy to obtain from the geotiff tags)
    # proj4_string = "+proj=stere +lat_0=90 +lon_0=-45 +lat_ts=70
    # +ellps=WGS84 +datum=WGS84 +units=m"
    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromWkt(proj)
    proj4_string = spatial_ref.ExportToProj4()

    # further we need to know the extent of the image in projection coordinates
    #
    # The values (in meter) can be taken from the geotiff tags, e.g.:
    #
    #   Corner Coordinates:
    #                   |            |
    #                   v            v
    #   Upper Left  (  -25000.000,   25000.000) ( 19d 4'24.76"W, 81d49'14.99"N)
    #   Lower Left  (  -25000.000,  -25000.000) ( 18d59'32.92"W, 81d22'23.67"N)
    #   Upper Right (   25000.000,   25000.000) ( 15d55'35.24"W, 81d49'14.99"N)
    #   Lower Right (   25000.000,  -25000.000) ( 16d 0'27.08"W, 81d22'23.67"N)
    #                   ^            ^
    #                   |            |
    #
    # Here we chose an example, where the projection center is not the
    # image center, but 500km further south and define the extent as
    # (x_min, y_min, x_max, y_max)
    # extent = (-250000., -750000., 250000., -250000.)
    x_min, x_size, _, y_min, _, y_size = ds.GetGeoTransform()
    extent = (x_min, y_min,
              x_min+x_size*ds.RasterXSize, y_min+y_size*ds.RasterYSize)

    # Step1: Compute the projection coordinates for each pixel

    # 1.a: Start with creating indices (lets say for every 250 pixels)
    # and make sure to have starting and ending pixels
    x_ind = np.arange(0, image_dimensions[0]+1, Config.instance.gcp_separation)
    y_ind = np.arange(0, image_dimensions[1]+1, Config.instance.gcp_separation)

    if x_ind[-1] != image_dimensions[0]:
        x_ind = np.append(x_ind, image_dimensions[0])
        y_ind = np.append(y_ind, image_dimensions[1])

    # 1.b: Convert these to meters according to image extent
    x_coarse = np.array(x_ind)*x_size + extent[0]
    y_coarse = np.array(y_ind)*y_size + extent[1]

    if method == 'optimize':
        polygon_x, polygon_y = get_land_polygons(image_order)
        print('NEW ITERATION')
        print(polygon_x, polygon_y)

        x_ind = np.round((x_coarse - x_coarse[0])/x_size)
        y_ind = np.round((y_coarse - y_coarse[0])/y_size)

        polygon_x = np.array(polygon_x)
        polygon_y = np.array(polygon_y)
        polygon_x_ind = np.round((polygon_x - x_coarse[0])/x_size)
        polygon_y_ind = np.round((polygon_y - y_coarse[0])/y_size)

    # 1.c: So far we only computed the x, y coordinates in 1D, but  we
    #      need them for all pixels (both x, y need to be 2D arrays)
    xx, yy = np.meshgrid(x_coarse, y_coarse)

    # Step 2: Inverse geodetic transformation
    #         (from projection coordinates to lat/lon using projection info)
    p = Proj(proj4_string)
    lons, lats = p(xx, yy, inverse=True)
    polygon_x = polygon_y = polygon_x_ind = polygon_y_ind = 0
    polygon_lons, polygon_lats = p(polygon_x, polygon_y, inverse=True)

    return xx, yy, x_ind, y_ind, lons, lats,\
        polygon_x, polygon_y, polygon_x_ind, polygon_y_ind, \
        polygon_lons, polygon_lats


def get_land_polygons(image_order):
    create_shape(image_order)
    clipped_shape_name = clip_shape(image_order)
    polygon_x, polygon_y = get_polygons_from_shapefile(clipped_shape_name)

    return polygon_x, polygon_y


def create_shape(image_order):
    shape_name = '/'.join(["data", "sar",
                          str(image_order["file"].split("/")[0]),
                          "roi.shp"])
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.CreateDataSource(shape_name)
    layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    feat = ogr.Feature(defn)
    feat.SetField('id', 1)

    box = get_boundary_from_geotiff(image_order)

    ring = ogr.Geometry(ogr.wkbLinearRing)
    source = osr.SpatialReference()
    source.ImportFromEPSG(3995)
    target = osr.SpatialReference()
    target.ImportFromEPSG(3995)
    transform = osr.CoordinateTransformation(source, target)
    for point in box:
        ring.AddPoint(point[0], point[1])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()

    geom = ogr.CreateGeometryFromWkt(poly_wkt)
    geom.Transform(transform)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
    ds = layer = feat = geom = None


def get_boundary_from_geotiff(image_order):
    filename = image_order["file"].replace(".tiff", "_roi.tiff")
    filename = '/'.join(["data", "sar", filename])
    data = gdal.Open(filename, gdal.GA_ReadOnly)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    data = None

    box = [(minx, maxy), (minx, miny), (maxx, miny),
           (maxx, maxy), (minx, maxy)]

    return box


def clip_shape(image_order):
    clipped_name = '/'.join(["data", "sar",
                             str(image_order["file"].split("/")[0]),
                             "clipped.shp"])
    params = dict(
        shapefile="shapefiles/north_shapefiles.shp",
        clipped_shape_name=clipped_name,
        bounding_shape=clipped_name.replace("clipped.shp",
                                            "roi.shp")
        )
    cmd = ("ogr2ogr -clipsrc {shapefile} {clipped_shape_name} "
           "{bounding_shape}").format(**params)

    call(cmd, shell=True)

    return params['clipped_shape_name']


def get_polygons_from_shapefile(clipped_shape_name):
    polygons = shapefile.Reader(clipped_shape_name)
    poly = polygons.shape().points
    with open(('data/sar/case19/points.csv'), 'wb') as ofile:
        writer = csv.writer(ofile, delimiter=',', quoting=csv.QUOTE_ALL)
        writer.writerow(["X", "Y"])
        for point in poly:
            writer.writerow(point)
    np.random.seed(999)
    random_index = np.random.choice(range(0, len(poly)), int(len(poly)*0.05),
                                    replace=False)
    poly_subset = [poly[ii] for ii in random_index]
    poly_subset_x, poly_subset_y = zip(*poly_subset)

    return poly_subset_x, poly_subset_y


if __name__ == '__main__':
    inversion_gcps_recipe()
