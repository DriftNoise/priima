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

Example script to compute geographical coordinates for image pixels
based on a known projection, image dimension and pixel size

Author: Stefan Hendricks
"""

import math

import numpy as np
from osgeo import osr
from pyproj import Proj

from priima.config import Config


def inversion_gcps_recipe(ds, proj):
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

    segments_per_side = math.sqrt(Config.instance().num_gcps) - 1
    x_pix_separation = math.ceil(ds.RasterXSize / segments_per_side)
    y_pix_separation = math.ceil(ds.RasterYSize / segments_per_side)

    # Step1: Compute the projection coordinates for each pixel

    # 1.a: Start with creating indices (lets say for every 250 pixels)
    # and make sure to have starting and ending pixels
    x_ind = np.arange(0, image_dimensions[0]+1, x_pix_separation)
    y_ind = np.arange(0, image_dimensions[1]+1, y_pix_separation)

    if x_ind[-1] != image_dimensions[0]:
        x_ind = np.append(x_ind, image_dimensions[0])
        y_ind = np.append(y_ind, image_dimensions[1])

    # 1.b: Convert these to meters according to image extent
    x_coarse = np.array(x_ind)*x_size + extent[0]
    y_coarse = np.array(y_ind)*y_size + extent[1]

    # 1.c: So far we only computed the x, y coordinates in 1D, but  we
    #      need them for all pixels (both x, y need to be 2D arrays)
    xx, yy = np.meshgrid(x_coarse, y_coarse)

    # Step 2: Inverse geodetic transformation
    #         (from projection coordinates to lat/lon using projection info)
    p = Proj(proj4_string)
    lons, lats = p(xx, yy, inverse=True)
    polygon_x = polygon_y = polygon_x_ind = polygon_y_ind = 0
    polygon_lons, polygon_lats = p(polygon_x, polygon_y, inverse=True)

    return xx, yy, x_ind, y_ind, lons, lats, \
        polygon_x, polygon_y, polygon_x_ind, polygon_y_ind, \
        polygon_lons, polygon_lats


if __name__ == '__main__':
    inversion_gcps_recipe()
