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

import json
import os
import tempfile
from datetime import datetime
from shutil import rmtree
from unittest import TestCase
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset
from pyproj import Proj, transform

from priima.drift_handler import DriftHandler
from priima.order import Order


class TestDriftHandler(TestCase):
    def setUp(self):
        with tempfile.NamedTemporaryFile(mode="w") as temp_fh:
            order_string = dict(
                customer_name="Foufoutos",
                customer_email_address="foufoutos@example.com",
                order_name="order_name",
                sentinel1_order_number=40000,
                forecast_duration=4,
                center=[50, 25],
                region_size=200,
                resolution=30,
                drift_filename="/path/to/drift.nc",
                data_source="TOPAZ",
                grid_method="mpl",
                use_tides="False"
            )
            temp_fh.write(json.dumps(order_string))
            temp_fh.flush()
            temp_fh.seek(0)
            with open(temp_fh.name) as fh:
                order_data = json.load(fh)
                self.order = Order(order_data=order_data)

        self.time_range = [datetime(2022, 12, 10, 20, 0),
                           datetime(2022, 12, 10, 20, 0)]

        self. gcp_list = ["point1", "point2"]
        self.ncfile = "/path/to/drift.nc"

        lon = np.linspace(0, 360, 361)
        lat = np.linspace(0, 90, 91)
        self.drift_dict = dict([
            ("uice", np.random.rand(1, 361, 91)),
            ("vice", np.random.rand(1, 361, 91)),
            ("lon", lon),
            ("lat", lat),
            ("time", 0)
            ])

    def test_the_constructor_sets_attributes(self):
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)

        self.assertEqual(drift_handler.ncfile, self.ncfile)

    def test_subset_roi_from_meshgrid_results_to_expected_roi_for_nh(self):
        # priima assumes 4x larger region size hence effectively the
        # size is 200*4km. 800km around center = 50, 25 is roughly 10 degrees
        # in lon and 6 degrees in lat.
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)

        roi_drift_dict = \
            drift_handler.subset_roi_from_meshgrid(self.drift_dict)
        expected_min_lon = 19
        expected_max_lon = 30
        expected_min_lat = 47
        expected_max_lat = 53

        self.assertEqual(roi_drift_dict["lon"].min(), expected_min_lon)
        self.assertEqual(roi_drift_dict["lon"].max(), expected_max_lon)
        self.assertEqual(roi_drift_dict["lat"].min(), expected_min_lat)
        self.assertEqual(roi_drift_dict["lat"].max(), expected_max_lat)

    def test_subset_roi_from_meshgrid_results_to_expected_roi_for_sh(self):
        # priima assumes 4x larger region size hence effectively the
        # size is 200*4km. 800km around center = 50, 25 is roughly 10 degrees
        # in lon and 6 degrees in lat.
        self.order.center = [-50, 25]
        self.drift_dict["lat"] = np.linspace(0, -90, 91)
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)

        roi_drift_dict = \
            drift_handler.subset_roi_from_meshgrid(self.drift_dict)
        expected_min_lon = 19
        expected_max_lon = 30
        expected_min_lat = -53
        expected_max_lat = -46

        self.assertEqual(roi_drift_dict["lon"].min(), expected_min_lon)
        self.assertEqual(roi_drift_dict["lon"].max(), expected_max_lon)
        self.assertEqual(roi_drift_dict["lat"].min(), expected_min_lat)
        self.assertEqual(roi_drift_dict["lat"].max(), expected_max_lat)

    @patch("priima.drift_handler.DriftHandler._read_drift")
    def test_proper_reader_is_called_when_nextsim_data_are_used(self,
                                                                mock_read):
        self.order.data_source = "NEXTSIM"
        nextsim_path = create_temporary_nextsim_file()
        self.ncfile = nextsim_path
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)
        try:
            drift_handler.compute_drift()
        except Exception as ex:
            print("Test can not read drift data: ", ex.args)
        rmtree(os.path.dirname(nextsim_path))

        self.assertEqual(mock_read.call_count, 1)

    @patch("priima.drift_handler.DriftHandler._read_drift")
    def test_proper_reader_is_called_when_topaz_data_are_used(self,
                                                              mock_read):
        nextsim_path = create_temporary_nextsim_file()
        self.ncfile = nextsim_path
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)
        try:
            drift_handler.compute_drift()
        except Exception as ex:
            print("Test can not read drift data: ", ex.args)
        rmtree(os.path.dirname(nextsim_path))

        self.assertEqual(mock_read.call_count, 1)

    def test_read_drift_loads_nextsim_data(self):
        nextsim_path = create_temporary_nextsim_file()
        self.order.data_source = "NEXTSIM"
        self.ncfile = nextsim_path
        self.time_range = [datetime(2022, 3, 1, 20, 0),
                           datetime(2022, 3, 1, 20, 0)]
        drift_handler = DriftHandler(self.ncfile, self.time_range,
                                     self.order, self.gcp_list)
        drift_dict = drift_handler._read_drift()
        rmtree(os.path.dirname(nextsim_path))

        self.assertEqual(len(drift_dict['lon'].shape), 2)
        self.assertEqual(len(drift_dict['lat'].shape), 2)
        self.assertEqual(len(drift_dict['xx'].shape), 1)
        self.assertEqual(len(drift_dict['yy'].shape), 1)
        self.assertEqual(len(drift_dict['uice'].shape), 3)
        self.assertEqual(len(drift_dict['vice'].shape), 3)
        self.assertEqual(len(drift_dict['time'].shape), 1)

        self.assertAlmostEqual(drift_dict['lon'][0, 0], -84.93638315)
        self.assertAlmostEqual(drift_dict['lat'][0, 0], 41.27724294)


def create_temporary_nextsim_file():
    x_val = np.arange(-3600000, 3801000, 30000)
    y_val = np.arange(-4300000, 2801000, 30000)
    x_dim = len(x_val)
    y_dim = len(y_val)
    xx, yy = np.meshgrid(x_val, y_val)
    inProj = Proj(init="epsg:3413")
    outProj = Proj(init="epsg:4326")
    x_lon, y_lat = transform(inProj, outProj, xx, yy)

    temp_dir = tempfile.mkdtemp()
    filename = "20220301_hr-nersc-MODEL-nextsimf-ARC-b20220223-fv00.0.nc"
    full_path = os.path.join(temp_dir, filename)
    nc = Dataset(full_path, "w", format="NETCDF4")
    nc.createDimension("time", None)
    nc.createDimension("x", x_dim)
    nc.createDimension("y", y_dim)

    times = nc.createVariable("time", "f8", ("time",))
    x = nc.createVariable("x", "f8", ("x",))
    y = nc.createVariable("y", "f8", ("y",))
    lons = nc.createVariable("longitude", "f8", ("y", "x",))
    lats = nc.createVariable("latitude", "f8", ("y", "x",))
    uice = nc.createVariable("vxsi", "f8", ("time", "y", "x",))
    vice = nc.createVariable("vysi", "f8", ("time", "y", "x",))

    x[:] = x_val
    y[:] = y_val
    # field date: 2022-03-01 bulletin date: 2022-02-23
    times[:] = np.arange(44619.0208333333, 44620.0208333334,
                         0.041666666700621136)
    lats[:] = y_lat
    lons[:] = x_lon
    uice[:, :, :] = np.random.rand(24, y_dim, x_dim)
    vice[:, :, :] = np.random.rand(24, y_dim, x_dim)

    nc.close()

    return full_path
