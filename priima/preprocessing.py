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

import time
from pathlib import Path

import cv2
import numpy as np
from osgeo import gdal

from priima.config import Config
from priima.fast_cast import fc


def warp_image(
        *,
        forecast_step: str,
        forecast_duration: int,
        image_path: Path,
        pixels_drift: dict,
        ):
    """
    Warps a given image according to the transformation

    :param forecast_step: string of current forecast step
    :param forecast_duration: overall duration of forecast in hours
    :param image_path: Path to geotiff containing reference GeoTransform
    :pixels_drift: Dictionary which contains pixel points with their
                   location. Two keys are expected the 'start' and
                   'end' key which represent the starting and ending
                   pixel for a given point.
    """

    ds = gdal.Open(str(image_path))
    image = np.array(ds.GetRasterBand(1).ReadAsArray())
    pts1 = np.float32(pixels_drift['start'])
    pts2 = np.float32(pixels_drift['end'])
    pts1 = pts1.reshape(1, -1, 2)
    pts2 = pts2.reshape(1, -1, 2)
#    pts1 = pts1.reshape(-1, 4, 2)
#    pts2 = pts2.reshape(-1, 4, 2)

    _, x_res, _, _, _, _ = ds.GetGeoTransform()
    max_drift_m_per_h = 4000
    max_drift_distance = max_drift_m_per_h * forecast_duration
    # limit the buffersize to the maximum distance the ice could have moved
    buffersize = int(max_drift_distance / x_res)

    imbuffer = np.zeros((image.shape[0] + buffersize * 2,
                         image.shape[1] + buffersize * 2), np.float32)

    imbuffer[int(buffersize):int(buffersize + image.shape[0]),
             int(buffersize):int(buffersize + image.shape[1])] = image[:, :]

    pts1 = pts1 + buffersize
    pts2 = pts2 + buffersize

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

    if out_img.max() == 0:
        err_msg = (
            "TPS point fit invalid. Try decreasing the image resolution "
            "(via --output-resolution)")
        raise ValueError(err_msg)

    write_geotiff(
        forecast_step=forecast_step, image_path=image_path,
        raster_data=out_img, pixels_drift=pixels_drift,
        gtif_out=Config.instance().gtif_out, buffersize=buffersize,
    )


def write_geotiff(
        *,
        forecast_step: str,
        image_path: Path,
        raster_data: np.ndarray,
        pixels_drift: dict,
        gtif_out: bool,
        buffersize: int,
        ):
    """
    Save image in GTiFF format

    :param forecast_step: string of current forecast step
    :param image_path: Path to geotiff containing reference GeoTransform
    :param raster_data: The raster data to be saved
    :param pixels_drift: A dictionary containing the drift information and the
     projected coordinates of the warped GCPs
    :param gtif: If true, a cloud-optimized geotiff is produced
    :param buffersize: Size of no-data frame around the image
    """
    new_image_name = image_path.with_stem(
        f"{image_path.stem}_forecast_{forecast_step}_"
        f"{Config.instance().data_source}").name
    out_name = Config.instance().output_dir / new_image_name

    if gtif_out:
        driver = gdal.GetDriverByName('GTiff')
        dataset_out = driver.Create(str(out_name), raster_data.shape[1],
                                    raster_data.shape[0], 1, gdal.GDT_UInt16)
    else:
        cog_driver = gdal.GetDriverByName('MEM')
        dataset_out = cog_driver.Create(
            '',
            raster_data.shape[1],
            raster_data.shape[0],
            1,
            gdal.GDT_UInt16
        )

    dataset = gdal.Open(str(image_path))
    geotransformation = np.array(dataset.GetGeoTransform())
    ind_lonul = [counter for counter, pixel in enumerate(pixels_drift['start'])
                 if pixel == (0, 0)][0]
    indlatlr = [counter for counter, pixel in enumerate(pixels_drift['start'])
                if pixel == (0, 0)][0]
    geotransformation[0] = pixels_drift['lon_projected'][ind_lonul]
    geotransformation[3] = pixels_drift['lat_projected'][indlatlr]
    xpixel_extension = buffersize
    ypixel_extension = buffersize
    geotransformation[0] = geotransformation[0] \
        - (xpixel_extension + pixels_drift['end'][0][0])*geotransformation[1]
    geotransformation[3] = geotransformation[3] \
        - (ypixel_extension + pixels_drift['end'][0][1])*geotransformation[5]
    dataset_out.SetGeoTransform(geotransformation)
    dataset_out.SetProjection(dataset.GetProjection())

    dataset_out.GetRasterBand(1).WriteArray(raster_data[:, :])

    if not gtif_out:
        dataset_out.BuildOverviews("NEAREST", [2, 4, 8, 16, 32, 64])
        driver = gdal.GetDriverByName('GTiff')
        cog_options = [
            "COPY_SRC_OVERVIEWS=YES",
            "TILED=YES",
            "COMPRESS=LZW"
        ]
        driver.CreateCopy(
            str(out_name),
            dataset_out,
            options=cog_options
        )

    dataset_out.FlushCache()
    dataset_out = None
