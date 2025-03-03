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

from unittest import TestCase

from priima.config import Config
from priima.projection import get_projection, get_projection_epsg


class TestTopaz(TestCase):

    def test_get_projection_returns_expected_proj4_for_nh(self):
        Config.set_attribute('center', [50, 0])
        projection = get_projection()

        expected_projection = (
           "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 "
           "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

        self.assertEqual(projection.srs, expected_projection)

    def test_get_projection_returns_expected_proj4_for_sh(self):
        Config.set_attribute('center', [-50, 0])
        projection = get_projection()

        expected_projection = (
           "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 "
           "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

        self.assertEqual(projection.srs, expected_projection)

    def test_arctic_epsg_is_3413(self):
        Config.set_attribute('center', [78.13, 15.39])

        self.assertEqual(get_projection_epsg(), 3413)

    def test_antarctic_epsg_is_3031(self):
        Config.set_attribute('center', [-80.36, 58.59])

        self.assertEqual(get_projection_epsg(), 3031)
