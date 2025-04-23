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
from datetime import datetime
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import numpy as np
from osgeo import gdal

from priima.config import Config
from priima.shapefile import Shapefile
from priima.warp_forecast import (compute_time_range, update_point_location,
                                  warp_command)


class TestWarpForecast(TestCase):
    def setUp(self):
        gcp = gdal.GCP()
        gcp.GCPX = 14
        gcp.GCPY = 79
        gcp.GCPXX = 1024423.07
        gcp.GCPYY = -615535.48
        gcp.Info = "is_land"
        self.gcp_list = [gcp]
        self.drift = [np.array([0.1]), np.array([-0.2])]

        Config.set_attribute('center', [90.0, 45.0])
        self.landmask = Shapefile(path="shapefiles/arctic_landmask.shp")

    @patch('priima.warp_forecast.is_point_inland')
    def test_update_point_location_breaks_if_point_is_already_in_land(
            self, land_mock):
        update_point_location(self.gcp_list, self.drift, self.landmask)
        self.assertFalse(land_mock.called)

    def test_gcp_list_is_updated_even_if_point_is_in_land(self):
        updated_gcp = update_point_location(
            self.gcp_list, self.drift, self.landmask)

        self.assertEqual(len(self.gcp_list), len(updated_gcp))

    def test_is_point_inland_updates_point_correctly_if_point_is_in_land(
            self):
        self.gcp_list[0].Info = ""

        self.assertEqual(self.gcp_list[0].Info, "")

        updated_gcp = update_point_location(
            self.gcp_list, self.drift, self.landmask)

        self.assertEqual(updated_gcp[0].Info, "is_land")

    def test_update_point_location_returns_expected_gcp_location_for_nh(self):
        # for drift = [0.1, -0.2] m/s  we expect 360m and -720m displacement
        # in 1 hour.
        gcp = gdal.GCP()
        gcp.Info = ""
        gcp.GCPX = 5
        gcp.GCPY = 83
        gcp.GCPXX = 581579.20
        gcp.GCPYY = -488002.89
        self.gcp_list = [gcp]

        gcp_updated = update_point_location(
            self.gcp_list, self.drift, self.landmask)
        displacement_x = math.ceil(gcp_updated[0].GCPXX - 581579.20)
        displacement_y = math.ceil(gcp_updated[0].GCPYY - (-488002.89))

        self.assertEqual(displacement_x, 360)
        self.assertEqual(displacement_y, -720)

    def test_update_point_location_returns_expected_gcp_location_for_sh(self):
        # for drift = [0.1, -0.2] m/s  we expect 360m and -720m displacement
        # in 1 hour.
        gcp = gdal.GCP()
        gcp.Info = ""
        gcp.GCPX = 5
        gcp.GCPY = -65
        gcp.GCPXX = 240419.29
        gcp.GCPYY = 2748005.02
        self.gcp_list = [gcp]

        gcp_updated = update_point_location(
            self.gcp_list, self.drift, self.landmask)
        displacement_x = math.ceil(gcp_updated[0].GCPXX - 240419.29)
        displacement_y = math.ceil(gcp_updated[0].GCPYY - (2748005.02))

        self.assertEqual(displacement_x, 360)
        self.assertEqual(displacement_y, -720)

    def test_compute_time_range_returns_expected_range_as_array(self):
        filename = (
            "/some/path/to/data/forecasts/fast-cast/rifkol-use-case/"
            "S1_EW_HH_20230325T211045_20230325T211345_FC2_Rifkol_30_roi.tiff"
        )
        forecast_duration = 5.0
        time_range = compute_time_range(filename, forecast_duration)

        self.assertEqual(
            time_range,
            [datetime(2023, 3, 25, 21), datetime(2023, 3, 26, 2)]
        )

    def test_warp_command(self):
        filename = Path("blah.tiff")
        with patch(
            'priima.warp_forecast.get_half_region_size',
            return_value=50_000
        ) as get_half_region_size:
            warp_command_string = warp_command(filename)

            self.assertEqual(
                warp_command_string,
                (
                    "gdalwarp -overwrite -tr 100 100 -r near -multi "
                    "-srcnodata 0 -dstnodata 0 -order 3 -t_srs "
                    "'+init=epsg:3413' "
                    "-te -50000.0 -50000.0 50000.0 50000.0 "
                    "blah.tiff blah_roi.tiff"
                )
            )
            get_half_region_size.assert_called_once()
