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
import shutil
import tempfile
from pathlib import Path
from unittest import TestCase

from priima.order import Order, load_order


class TestOrder(TestCase):
    def setUp(self):
        self.tempdir = Path(tempfile.mkdtemp())
        self.order_data = dict(
            customer_name="Bob Smith",
            customer_email_address="bob@example.com",
            order_name="mooring-arctic-1",
            sentinel1_order_number=15,
            forecast_duration=40,  # hours
            center=[90.0, 45.0],
            region_size=100,
            resolution=100,
            data_source="TOPAZ",
            grid_method="mpl",
            use_tides="False"
        )

    def test_order_data_arg_is_keyword_only_arg(self):
        error_msg = \
            r"__init__\(\) takes 1 positional argument but 2 were given"
        with self.assertRaisesRegex(TypeError, error_msg):
            Order(self.order_data)  # pylint: disable=missing-kwoa

    def test_attributes_are_set_by_input_data(self):
        order = Order(order_data=self.order_data)

        self.assertEqual(order.customer_name, "Bob Smith")
        self.assertEqual(order.customer_email_address, "bob@example.com")
        self.assertEqual(order.order_name, "mooring-arctic-1")
        self.assertEqual(order.sentinel1_order_number, 15)
        self.assertEqual(order.forecast_duration, 40)
        self.assertEqual(order.center, [90.0, 45.0])
        self.assertEqual(order.region_size, 100)
        self.assertEqual(order.resolution, 100)
        self.assertEqual(order.data_source, "TOPAZ")
        self.assertEqual(order.grid_method, "mpl")
        self.assertEqual(order.use_tides, False)

    def test_plot_wind_default_value_is_false(self):
        order = Order(order_data=self.order_data)

        self.assertFalse(order.plot_wind)

    def test_plot_wind_default_value_is_set_from_order_data(self):
        self.order_data['plot_wind'] = True
        order = Order(order_data=self.order_data)

        self.assertTrue(order.plot_wind)

        self.order_data['plot_wind'] = False
        order = Order(order_data=self.order_data)

        self.assertFalse(order.plot_wind)

    def test_arctic_order_without_explicit_data_source_uses_topaz(self):
        self.order_data['center'] = [78.13, 15.39]  # lat, lon => Arctic
        self.order_data['data_source'] = ""  # intially unset

        order = Order(order_data=self.order_data)

        self.assertEqual(order.data_source, "TOPAZ")

    def test_antarctic_order_without_explicit_data_source_uses_icon(self):
        self.order_data['center'] = [-80.36, 58.59]  # lat, lon => Antarctic
        self.order_data['data_source'] = ""  # intially unset

        order = Order(order_data=self.order_data)

        self.assertEqual(order.data_source, "ICON")

    def test_data_source_defined_in_order_then_order_value_is_used(self):
        antarctic_topaz_order_data = self.order_data
        # lat, lon => Antarctic
        antarctic_topaz_order_data['center'] = [-80.36, 58.59]
        # use TOPAZ in Antarctic (if that seems sensible to you)
        antarctic_topaz_order_data['data_source'] = "TOPAZ"
        antarctic_topaz_order = Order(order_data=antarctic_topaz_order_data)

        arctic_topaz_order_data = self.order_data
        # lat, lon => Antarctic
        arctic_topaz_order_data['center'] = [78.13, 15.39]
        # use TOPAZ in Antarctic (if that seems sensible to you)
        arctic_topaz_order_data['data_source'] = "ICON"
        arctic_topaz_order = Order(order_data=arctic_topaz_order_data)

        self.assertEqual(antarctic_topaz_order.data_source, "TOPAZ")
        self.assertEqual(arctic_topaz_order.data_source, "ICON")

    def test_arctic_order_epsg_is_3413(self):
        self.order_data['center'] = [78.13, 15.39]  # lat, lon => Arctic

        order = Order(order_data=self.order_data)

        self.assertEqual(order.epsg, "+init=epsg:3413")

    def test_antarctic_order_epsg_is_3031(self):
        self.order_data['center'] = [-80.36, 58.59]  # lat, lon => Antarctic

        order = Order(order_data=self.order_data)

        self.assertEqual(order.epsg, "+init=epsg:3031")

    def test_loading_an_order_from_file_sets_order_params(self):
        order_file = Path(self.tempdir, "bob_test_order.json")
        with open(order_file, "w", encoding="utf-8") as fp:
            json.dump(self.order_data, fp)

        order, _ = load_order(order_file=order_file)

        # don't check all fields, just enough to check that we loaded the
        # order data correctly
        self.assertEqual(order.customer_name, "Bob Smith")
        self.assertEqual(order.center, [90, 45])
        self.assertEqual(order.data_source, "TOPAZ")

    def test_loading_a_non_existent_order_file_raises_filenotfound_error(self):
        order_file = Path(self.tempdir, "i-dont-exist")

        error_msg = f"{order_file} does not exist!"
        with self.assertRaisesRegex(FileNotFoundError, error_msg):
            load_order(order_file=order_file)

    def test_loading_an_order_without_order_file_arg_raises_error(self):
        error_msg = (
            "missing 1 required keyword-only argument: 'order_file'"
        )
        with self.assertRaisesRegex(TypeError, error_msg):
            load_order()  # pylint: disable=missing-kwoa

    # TODO: test entering invalid input data

    def tearDown(self):
        if self.tempdir.exists():
            shutil.rmtree(self.tempdir)
