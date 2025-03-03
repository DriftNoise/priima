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

import os
import time

import cv2
import numpy as np
from osgeo import gdal

from priima.fast_cast import fc


def load_image(image_path):
    """
    Loads the geotiff image
    :param image_path: Full path of the image
    :type image_path: string
    """
    dataset = gdal.Open(image_path)
    data = dataset.ReadAsArray()
    projection = dataset.GetProjection()

    return [data, projection, dataset]


def warp_image(forecast_step, image_path, pixels_drift,
               data_source, grid_method, base_path, store_training_data):
    """
    Warps a given image according to the transformation

    :param image_path: Full path of the image
    :type image_path: String

    :param pixels_drift: Dictionary which contains pixel points with their
                         location. Two keys are expected the 'start' and
                         'end' key which represent the starting and ending
                         pixel for a given point.
    :type pixels_drift: Dictionary

    :param base_path: Path and name to save the warped images to
    :type base_path: String
    """

    ds = gdal.Open(image_path)
    image = np.array(ds.GetRasterBand(1).ReadAsArray())
    pts1 = np.float32(pixels_drift['start'])
    pts2 = np.float32(pixels_drift['end'])
    pts1 = pts1.reshape(1, -1, 2)
    pts2 = pts2.reshape(1, -1, 2)
#    pts1 = pts1.reshape(-1, 4, 2)
#    pts2 = pts2.reshape(-1, 4, 2)

    buffersize = 2000

    imbuffer = np.zeros((image.shape[0] + buffersize,
                         image.shape[1] + buffersize), np.float32)

    imbuffer[int(1000):int(1000 + image.shape[0]),
             int(1000):int(1000 + image.shape[1])] = image[:, :]

    pts1 = pts1 + 1000
    pts2 = pts2 + 1000

    if store_training_data:
        np.savetxt(os.path.join(
            base_path, 'pts1_{}.txt'.format(forecast_step)), pts1[0, :, :])
        np.savetxt(os.path.join(
            base_path, 'pts2_{}.txt'.format(forecast_step)), pts2[0, :, :])

    # OpenCV
    if fc is None:
        # XXX: this for loop could probably be just a list comprehension; we
        # need more tests though.
        matches = []
        for ipoint in range(0, pts1.shape[1]):
            matches.append(cv2.DMatch(ipoint, ipoint, 0))

        # print('OPENCV')
        t = time.time()
        tps = cv2.createThinPlateSplineShapeTransformer()
        tps.estimateTransformation(pts2, pts1, matches)
        out_img = tps.warpImage(imbuffer)
        print(f'Warping with OpenCV took {time.time() - t} s')

    # Fast Cast
    else:
        # XXX: this for loop could probably be just a list comprehension; we
        # need more tests though.
        matches = []
        for ipoint in range(0, pts1.shape[1]):
            matches.append([ipoint, ipoint, 0])

        # print('FAST_CAST')
        t = time.time()
        tps = fc.ThinPlateSplineShapeTransformerImpl()
        tps.estimateTransformation(pts2, pts1, matches)
        out_img = tps.warpImage(imbuffer)
        print(f'Warping with FastCast took {time.time() - t} s')

    write_geotiff(forecast_step, image_path, out_img, imbuffer.shape,
                  pixels_drift, grid_method, data_source, base_path)


def write_geotiff(forecast_step, image_path, raster_data, out_shape,
                  pixels_drift, grid_method, data_source, base_path):
    """
    Save image in GTiFF format

    :param dataset: Object containing geospatial information
    :type dataset: 'osgeo.gdal.Dataset'

    :param raster_data: The raster data to be saved
    :type raster_data: numpy.ndarray

    :param base_path: Path and name to save the warped images to
    :type base_path: String
    """

    driver = gdal.GetDriverByName('GTiff')
    _, _, dataset = load_image(image_path)

    base_image_path_name = os.path.join(
        base_path, os.path.basename(image_path))

    out_name = base_image_path_name.replace(
        '.tiff', '_forecast_{0}_{1}_{2}.tiff'.format(
            forecast_step, data_source, grid_method))

    dataset_out = driver.Create(out_name, raster_data.shape[1],
                                raster_data.shape[0], 1, gdal.GDT_UInt16)

    geotransformation = np.array(dataset.GetGeoTransform())
    ind_lonul = [counter for counter, pixel in enumerate(pixels_drift['start'])
                 if pixel == (0, 0)][0]
    indlatlr = [counter for counter, pixel in enumerate(pixels_drift['start'])
                if pixel == (0, 0)][0]
    geotransformation[0] = pixels_drift['lon_projected'][ind_lonul]
    geotransformation[3] = pixels_drift['lat_projected'][indlatlr]
    xpixel_extension = 1000
    ypixel_extension = 1000
    geotransformation[0] = geotransformation[0] \
        - (xpixel_extension + pixels_drift['end'][0][0])*geotransformation[1]
    geotransformation[3] = geotransformation[3] \
        - (ypixel_extension + pixels_drift['end'][0][1])*geotransformation[5]
    dataset_out.SetGeoTransform(geotransformation)
    dataset_out.SetProjection(dataset.GetProjection())

    dataset_out.GetRasterBand(1).WriteArray(raster_data[:, :])

    dataset_out.FlushCache()
    dataset_out = None
