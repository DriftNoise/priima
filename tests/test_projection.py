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
import tempfile
from unittest import TestCase

from priima.order import Order
from priima.projection import get_projection


class TestTopaz(TestCase):
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

    def test_get_projection_returns_expected_proj4_for_nh(self):
        projection = get_projection(self.order)

        expected_projection = (
           "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 "
           "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

        self.assertEqual(projection.srs, expected_projection)

    def test_get_projection_returns_expected_proj4_for_sh(self):
        self.order.center[0] = -50
        projection = get_projection(self.order)

        expected_projection = (
           "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 "
           "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

        self.assertEqual(projection.srs, expected_projection)
