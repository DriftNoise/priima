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
import os
import tempfile
from unittest import TestCase
from unittest.mock import patch

from priima.config import Config
from priima.image_forecast import netcdf_file_dates, select_drift_files
from priima.order import Order


class TestImageForecast(TestCase):
    def setUp(self):
        self.time_range = [datetime.datetime(2019, 5, 27, 17),
                           datetime.datetime(2019, 5, 30, 12)]
        self.time_range1 = [datetime.datetime(2019, 8, 28),
                            datetime.datetime(2019, 8, 28)]
        Config()

        self.tmpdir = tempfile.TemporaryDirectory()
        self.order_data = dict(
            customer_name="bob_papadopoulos",
            customer_email_address="bob@driftnoise.com",
            order_name="test_order_name",
            sentinel1_order_number=100,
            forecast_duration=12,
            center=[80.36, 58.59],
            region_size=150,
            resolution=30,
            data_source="",
            grid_method="mpl",
            use_tides=False
            )
        self.order = Order(order_data=self.order_data)
        stdout_file = os.path.join(self.tmpdir.name, 'stdout_file')
        stderr_file = os.path.join(self.tmpdir.name, 'stderr_file')
        with open(stdout_file, 'w') as fm:
            fm.write('sentinel_file.tiff')
        with open(stderr_file, 'w') as fm:
            fm.write('')
        self.stdout = open(stdout_file, 'r')
        self.stderr = open(stderr_file, 'r')

    @patch("priima.image_forecast.glob")
    def test_job_worker_error(self, mock_glob):
        self.order.data_source = "TOPAZ"
        filenames_list = \
            ['20190828_hr-metno-MODEL-topaz4-ARC-b20190819-fv02.0.nc',
             '20190828_hr-metno-MODEL-topaz4-ARC-b20190820-fv02.0.nc',
             '20190828_hr-metno-MODEL-topaz4-ARC-b20190821-fv02.0.nc',
             '20190825_hr-metno-MODEL-topaz4-ARC-b20190818-fv02.0.nc',
             '20190827_hr-metno-MODEL-topaz4-ARC-b20190827-fv02.0.nc',
             '20190828_hr-metno-MODEL-topaz4-ARC-b20190827-fv02.0.nc']
        mock_glob.return_value = filenames_list
        forecast_files = select_drift_files(self.time_range1, self.order)

        expected_file = \
            ['20190828_hr-metno-MODEL-topaz4-ARC-b20190827-fv02.0.nc']
        self.assertEqual(forecast_files, expected_file)

    @patch("priima.image_forecast.glob")
    def test_if_multiple_topaz_forecasts_available_only_recent_isreturned(
            self, mock_glob):
        self.order.data_source = "TOPAZ"
        filenames_list = \
            ['20190530_hr-metno-MODEL-topaz4-ARC-b20190529-fv02.0.nc',
             '20190530_hr-metno-MODEL-topaz4-ARC-b20190528-fv02.0.nc',
             '20190530_hr-metno-MODEL-topaz4-ARC-b20190527-fv02.0.nc']
        mock_glob.return_value = filenames_list
        time_range = [datetime.datetime(2019, 5, 30, 12),
                      datetime.datetime(2019, 5, 30, 22)]
        forecast_files = select_drift_files(time_range, self.order)

        expected_file = \
            ['20190530_hr-metno-MODEL-topaz4-ARC-b20190529-fv02.0.nc']
        self.assertEqual(forecast_files, expected_file)

    @patch("priima.image_forecast.glob")
    def test_if_multiple_nextsim_forecasts_available_only_recent_isreturned(
            self, mock_glob):
        self.order.data_source = "NEXTSIM"
        filenames_list = \
            ['20190530_hr-nersc-MODEL-nextsimf-ARC-b20190529-fv00.0.nc',
             '20190530_hr-nersc-MODEL-nextsimf-ARC-b20190528-fv00.0.nc',
             '20190530_hr-nersc-MODEL-nextsimf-ARC-b20190527-fv00.0.nc']
        mock_glob.return_value = filenames_list
        time_range = [datetime.datetime(2019, 5, 30, 12),
                      datetime.datetime(2019, 5, 30, 22)]
        forecast_files = select_drift_files(time_range, self.order)

        expected_file = \
            ['20190530_hr-nersc-MODEL-nextsimf-ARC-b20190529-fv00.0.nc']
        self.assertEqual(forecast_files, expected_file)

    @patch("priima.image_forecast.glob")
    def test_selected_topaz_files_are_within_forecast_range(self, mock_glob):
        self.order.data_source = "TOPAZ"
        filenames_list = \
            ['20190531_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190530_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190529_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190528_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190527_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190526_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc']
        mock_glob.return_value = filenames_list
        forecast_files = select_drift_files(self.time_range, self.order)

        expected_files = \
            ['20190530_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190529_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190528_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc',
             '20190527_hr-metno-MODEL-topaz4-ARC-b20190526-fv02.0.nc']

        forecast_files.sort()
        expected_files.sort()
        self.assertEqual(forecast_files, expected_files)

    @patch("priima.image_forecast.glob")
    def test_selected_nextsim_files_are_within_forecast_range(self, mock_glob):
        self.order.data_source = "NEXTSIM"
        filenames_list = \
            ['20190531_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190530_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190529_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190528_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190527_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190526_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc']
        mock_glob.return_value = filenames_list
        forecast_files = select_drift_files(self.time_range, self.order)

        expected_files = \
            ['20190530_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190529_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190528_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc',
             '20190527_hr-nersc-MODEL-nextsimf-ARC-b20190526-fv02.0.nc']

        forecast_files.sort()
        expected_files.sort()
        self.assertEqual(forecast_files, expected_files)

    def test_netcdf_file_dates_returns_dates_from_single_plain_fname(self):
        netcdf4_file_list = [
            '20230326_hr-metno-MODEL-topaz4-ARC-b20230325-fv02.0.nc'
        ]

        self.assertEqual(
            netcdf_file_dates(netcdf4_file_list), [datetime.date(2023, 3, 26)]
        )

    def test_netcdf_file_dates_returns_dates_from_multiple_plain_fnames(self):
        netcdf4_file_list = [
            '20230326_hr-metno-MODEL-topaz4-ARC-b20230325-fv02.0.nc',
            '20230429_hr-metno-MODEL-topaz4-ARC-b20230429-fv02.0.nc'
        ]

        self.assertEqual(
            netcdf_file_dates(netcdf4_file_list),
            [
                datetime.date(2023, 3, 26),
                datetime.date(2023, 4, 29)
            ]
        )

    def test_netcdf_file_dates_returns_dates_from_multiple_full_paths(self):
        netcdf4_file_list = [
            (
                '/data/L3/ice_drift_forecasts/topaz4/'
                '20230326_hr-metno-MODEL-topaz4-ARC-b20230325-fv02.0.nc'
            ),
            (
                '/data/L3/ice_drift_forecasts/topaz4/'
                '20230429_hr-metno-MODEL-topaz4-ARC-b20230429-fv02.0.nc'
            )
        ]

        self.assertEqual(
            netcdf_file_dates(netcdf4_file_list),
            [
                datetime.date(2023, 3, 26),
                datetime.date(2023, 4, 29)
            ]
        )

    def tearDown(self):
        Config._drop()  # pylint: disable=protected-access
        self.tmpdir.cleanup()
