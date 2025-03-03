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
from unittest import TestCase, skip

from priima.order import Order
from priima.wind import prepare_wind_dict_from_icon, wind2drift


class TestIconLibrary(TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.icon_file_regrided = ("tests/data/ICON_iko_single_level_elements"
                                   "_world_combined_10M_U_U_regridded.nc")

        order_data = dict(
            customer_name="Bob Smith",
            customer_email_address="bob@example.com",
            order_name="priima_order",
            sentinel1_order_number=15,
            forecast_duration=40,  # hours
            center=[70.0, 45.0],
            region_size=100,
            resolution=100,
            data_source="TOPAZ",
            grid_method="mpl",
            use_tides=False
        )
        self.order = Order(order_data=order_data)

    @skip('test data file missing')
    def test_prepare_wind_dict_from_icon_returns_expected_keys(self):
        wind_dict = prepare_wind_dict_from_icon(self.icon_file_regrided)

        wind_dict_keys = wind_dict.keys()

        self.assertTrue('lon' in wind_dict_keys)
        self.assertTrue('lat' in wind_dict_keys)
        self.assertTrue('xx' in wind_dict_keys)
        self.assertTrue('yy' in wind_dict_keys)
        self.assertTrue('uice' in wind_dict_keys)
        self.assertTrue('vice' in wind_dict_keys)
        self.assertTrue('u_offset' in wind_dict_keys)
        self.assertTrue('v_offset' in wind_dict_keys)
        self.assertTrue('u_scale' in wind_dict_keys)
        self.assertTrue('v_scale' in wind_dict_keys)

    def test_wind2drift_returns_expected_components_for_nh_rotation(self):

        u_component_initial = 0
        v_component_initial = 4
        u_component_rotated, v_component_rotated = \
            wind2drift(u_component_initial, v_component_initial, self.order)
        u_component_expected = 0.042
        v_component_expected = 0.091

        self.assertEqual([round(u_component_rotated[0], 3),
                          round(v_component_rotated[0], 3)],
                         [u_component_expected, v_component_expected])

    def test_wind2drift_returns_expected_components_for_sh_rotation(self):
        self.order.center = [-70, 45]
        u_component_initial = 0
        v_component_initial = 4
        u_component_rotated, v_component_rotated = \
            wind2drift(u_component_initial, v_component_initial, self.order)
        u_component_expected = -0.042
        v_component_expected = 0.091

        self.assertEqual([round(u_component_rotated[0], 3),
                          round(v_component_rotated[0], 3)],
                         [u_component_expected, v_component_expected])

#    def tearDown(self):
#        self.tmpdir.cleanup()
