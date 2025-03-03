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

import tempfile
from pathlib import Path
from unittest import TestCase

import numpy as np
from osgeo import gdal, osr

from priima.config import Config
from priima.geo_tools import get_half_region_size, transform_to_WGS84


class TestDriftFromWind(TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        rows, cols = 5, 5
        data = np.zeros((rows, cols), dtype=np.uint8)

        test_file = self.tmpdir / 'test.tif'
        pixel_size = 40.0
        x_origin = 0.0
        y_origin = 0.0
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(str(test_file), cols, rows, 1, gdal.GDT_Byte)
        geotransform = (x_origin, pixel_size, 0, y_origin, 0, -pixel_size)
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3413)
        dataset.SetProjection(srs.ExportToWkt())
        band = dataset.GetRasterBand(1)
        band.WriteArray(data)
        band.SetNoDataValue(0)
        dataset.FlushCache()
        dataset = None

        Config.set_attribute('image', test_file)

    def test_transform_to_polar_stereographic_nh(self):
        x = 0
        y = 0
        epsg_code = 3413
        expected_lat = 90
        expected_lon = -45
        lat, lon = transform_to_WGS84(x=x, y=y, epsg_code=epsg_code)

        self.assertAlmostEqual(lat, expected_lat, places=4)
        self.assertAlmostEqual(lon, expected_lon, places=4)

    def test_get_half_region_size(self):
        expected_half_region_size = 100.0
        half_region_size = get_half_region_size()

        self.assertEqual(expected_half_region_size, half_region_size)
