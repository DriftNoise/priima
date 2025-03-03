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


class Order:
    """
    Data structure to hold information needed to run the forecast for this
    customer's order.

    The input data is JSON-formatted; for example:
    ```
    {
        "customer_name": "Bob Smith",
        "customer_email_address": "bob@example.com",
        "order_name": "mooring station 1",
        "sentinel1_order_number": 15,
        "forecast_duration": 40,
        "center: [lat, lon]",
        "region_size": 60,
        "resolution": 50,
        "data_source": "TOPAZ",   # or ICON
        "grid_method": "mpl",     # or optimize
        "plot_wind": false,       # or true; only runs with ICON
        "use_tides": "False"
    }
    ```
    """
    def __init__(self, *, order_data):
        self.customer_name = order_data['customer_name']
        self.customer_email_address = order_data['customer_email_address']
        self.order_name = order_data['order_name']
        self.center = order_data['center']
        self.region_size = order_data['region_size']
        self.resolution = order_data['resolution']
        self.sentinel1_order_number = order_data['sentinel1_order_number']
        self.forecast_duration = order_data['forecast_duration']
        if order_data['data_source'] == "" and self.center[0] > 0:
            self.data_source = "TOPAZ"
        elif order_data['data_source'] == "" and self.center[0] < 0:
            self.data_source = "ICON"
        else:
            self.data_source = order_data['data_source']
        self.grid_method = order_data['grid_method']
        if order_data['use_tides'] == "False":
            self.use_tides = False
        elif order_data['use_tides'] == "True":
            self.use_tides = True
        else:
            print("Invalid key used for 'use_tides'. It should be True/False")

        if 'plot_wind' in order_data:
            self.plot_wind = order_data['plot_wind']
        else:
            self.plot_wind = False

    @property
    def epsg(self):
        # If the centre latitude (the first element of self.center) is
        # greater than zero, then we're in the Arctic (epsg:3413); otherwise
        # we're in the Antarctic (epsg:3031).
        epsg_code = 3413 if self.center[0] > 0 else 3031
        proj_str = f'+init=epsg:{epsg_code}'

        return proj_str


def load_order(*, order_file):
    if not os.path.exists(order_file):
        raise FileNotFoundError("{} does not exist!".format(order_file))

    with open(order_file) as fp:
        order_data = json.load(fp)

    order = Order(order_data=order_data)

    return order, order_data
