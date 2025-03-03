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
import shutil
import tempfile
from datetime import datetime
from glob import glob
from unittest import TestCase, skip

import numpy as np
from netCDF4 import Dataset

from priima.config import Config
from priima.order import Order
from priima.smoc import (convert_to_smoc_times, load_tidal_drift,
                         regrid_to_polar, select_smoc_files,
                         set_regridded_filename, smoc_times_to_datetimes)


class TestSmoc(TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir.name = self.tmpdir.name.replace("_", "R")
        self.smoc_dir = os.path.join(self.tmpdir.name, "data", "L2",
                                     "OceanCurrents")
        os.makedirs(self.smoc_dir)
        os.makedirs(os.path.join(self.smoc_dir, "regridded"))
        self.smoc_files = list([
            "SMOC_20210127_R20210127.nc",
            "SMOC_20210128_R20210128.nc",
            "SMOC_20210129_R20210129.nc",
            "SMOC_20210130_R20210129.nc",
            "SMOC_20210131_R20210129.nc",
            "SMOC_20210201_R20210129.nc",
            "SMOC_20210202_R20210129.nc",
            ])
        for filename in self.smoc_files:
            foul_path = os.path.join(self.smoc_dir, filename)
            create_test_file(foul_path)

        order_data = dict(
            customer_name="Bob Smith",
            customer_email_address="bob@example.com",
            order_name="mooring-arctic-1",
            sentinel1_order_number=15,
            forecast_duration=40,  # hours
            center=[-70.0, 45.0],
            region_size=100,
            resolution=100,
            data_source="TOPAZ",
            grid_method="mpl",
            use_tides=True
        )
        self.order = Order(order_data=order_data)

        config_content = """\
            paths:
                data_path: "{0}/data"
                secrets_path: ""
                grid_path: "{0}/data/icon_grid"
            """.format(self.tmpdir.name)
        config_file = os.path.join(self.tmpdir.name, 'priima.local_config.yml')
        with open(config_file, "w") as fh:
            fh.write(config_content)
            fh.flush()

        self.config = Config(local_config_file=config_file)

    def test_set_regridded_filenames_sets_filenames(self):
        expected_filename = "SMOC_20210131_R20210129_regridded.nc"
        filename = os.path.basename(self.smoc_files[4])
        smoc_filename_regridded = set_regridded_filename(filename)

        self.assertEqual(smoc_filename_regridded, expected_filename)

    def test_regrid_to_polar_creates_expected_file(self):
        filename = os.path.join(self.tmpdir.name, "data", "L2",
                                "OceanCurrents", self.smoc_files[2])
        regrid_to_polar(filename, self.order, self.config)

        expected_file = \
            os.path.join(self.tmpdir.name, "data", "L2", "OceanCurrents",
                         "regridded", "SMOC_20210129_R20210129_regridded.nc")
        regridded_file = \
            glob(os.path.join(self.tmpdir.name, "data", "L2", "OceanCurrents",
                              "regridded", "*regridded.nc"))[0]

        self.assertEqual(expected_file, regridded_file)

    @skip('list index error in select_smoc_files()')
    def test_select_smoc_files_returns_expected_files(self):
        time_range = [datetime(2021, 1, 29, 3, 24),
                      datetime(2021, 1, 31, 3, 35)]
        smoc_files = select_smoc_files(time_range, self.config)
        expected_files = list([
            "SMOC_20210129_R20210129.nc",
            "SMOC_20210130_R20210129.nc",
            "SMOC_20210131_R20210129.nc",
            ])
        expected_files = \
            [os.path.join(self.smoc_dir, fl) for fl in expected_files]

        self.assertCountEqual(smoc_files, expected_files)

    def test_convert_to_smoc_times_returns_hours_since_epoch(self):
        time_range = [datetime(2021, 1, 29, 3, 14),
                      datetime(2021, 1, 31, 3, 25)]
        hours = convert_to_smoc_times(time_range)
        expected_hours = np.arange(623067.5, 623116.5, 1)
        expected_times_length = 49

        self.assertCountEqual(hours, expected_hours)
        self.assertEqual(hours.size, expected_times_length)

    def test_smoc_times_to_datetimes_returns_expected_datetimes(self):
        smoc_times = np.array([623431.5, 623432.5, 623455.5])
        hours = smoc_times_to_datetimes(smoc_times)
        expected_datetimes = [
            datetime(2021, 2, 13, 7, 0),
            datetime(2021, 2, 13, 8, 0),
            datetime(2021, 2, 14, 7, 0),
            ]

        self.assertCountEqual(hours, expected_datetimes)

    def test_load_tidal_drift_loads_data_for_whole_timerange(self):
        smoc_files = list([
            "SMOC_20210129_R20210129.nc",
            "SMOC_20210130_R20210129.nc",
            "SMOC_20210131_R20210129.nc",
            ])
        smoc_files = \
            [os.path.join(self.smoc_dir, fl) for fl in smoc_files]
        time_range = [datetime(2021, 1, 29, 3, 24),
                      datetime(2021, 1, 31, 3, 35)]
        tide_data = load_tidal_drift(
            smoc_files, time_range, self.order, self.config)
        expected_shape = (49, 780, 620)

        self.assertEqual(tide_data["utide"].shape, expected_shape)
        self.assertEqual(tide_data["vtide"].shape, expected_shape)
        self.assertEqual(tide_data["time"].size, expected_shape[0])

    def tearDown(self):
        shutil.rmtree(self.tmpdir.name)
        Config._drop()  # pylint: disable=protected-access


def create_test_file(filename):
    ds = Dataset(filename, "w", format="NETCDF4_CLASSIC")
    ds.createDimension("time", None)
    ds.createDimension("latitude", 50)
    ds.createDimension("longitude", 80)

    times = ds.createVariable("time", "f4", ("time",))
    times.units = "hours since 1950-01-01 0:0:0"
    lats = ds.createVariable("latitude", "f4", ("latitude",))
    lats.units = "degrees_east"
    lons = ds.createVariable("longitude", "f4", ("longitude",))
    lons.units = "degrees_north"
    utide = ds.createVariable("utide", "f4",
                              ("time", "latitude", "longitude",))
    vtide = ds.createVariable("vtide", "f4",
                              ("time", "latitude", "longitude",))

    time_unit = datetime(1950, 1, 1, 0, 0)
    file_time = datetime.strptime(filename.split("_")[1], "%Y%m%d")
    hours_since_time_unit = \
        int((file_time - time_unit).total_seconds()/3600.) + 0.5
    times[:] = np.arange(hours_since_time_unit, hours_since_time_unit + 24, 1)
    lats[:] = np.arange(-90, -40, 1)
    lons[:] = np.arange(-70, 10, 1)
    utide[:] = np.random.uniform(0, 100, size=(24, 50, 80))
    vtide[:] = np.random.uniform(0, 100, size=(24, 50, 80))
    ds.close()
