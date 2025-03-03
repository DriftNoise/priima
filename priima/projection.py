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

from pyproj import Proj


def get_projection(order):
    if order.center[0] > 0:
        # epsg 3413
        projection = Proj("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 "
                          "+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    else:
        # epsg 3031
        projection = Proj("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 "
                          "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

    return projection
