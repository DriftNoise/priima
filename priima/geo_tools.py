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

import pyproj
from osgeo import gdal, osr

from priima.config import Config


def get_center_coordinate():
    """Retuns the center coordinate of the image"""
    dataset = gdal.Open(Config.instance().image)
    options = gdal.InfoOptions(format='json')
    center_projected = gdal.Info(
        dataset, options=options)['cornerCoordinates']['center']
    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromWkt(dataset.GetProjection())
    epsg_code = int(spatial_ref.GetAttrValue('AUTHORITY', 1))
    lat, lon = transform_to_WGS84(
        x=center_projected[0], y=center_projected[1], epsg_code=epsg_code)

    return lat, lon


def get_half_region_size():
    """Returns the maximum radius from the center point"""
    dataset = gdal.Open(str(Config.instance().image))
    x_size = dataset.RasterXSize
    y_size = dataset.RasterYSize
    _, x_res, _, _, _, y_res = dataset.GetGeoTransform()
    y_res = y_res if y_res > 0 else y_res * -1
    larger_distance = max(x_res * x_size / 2, y_res * y_size / 2)

    return float(larger_distance)


def transform_to_WGS84(*, x: float, y: float, epsg_code: int):
    """Transforms a projected coordinate to lon / lat"""
    transformer = pyproj.Transformer.from_crs(f"EPSG:{epsg_code}", "EPSG:4326")
    lat, lon = transformer.transform(x, y)
    return lat, lon
