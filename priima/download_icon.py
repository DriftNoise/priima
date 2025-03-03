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

import requests
from bs4 import BeautifulSoup

url00u = 'https://opendata.dwd.de/weather/nwp/icon/grib/00/u_10m/'
url00v = 'https://opendata.dwd.de/weather/nwp/icon/grib/00/v_10m/'
url06u = 'https://opendata.dwd.de/weather/nwp/icon/grib/06/u_10m/'
url06v = 'https://opendata.dwd.de/weather/nwp/icon/grib/06/v_10m/'
url12u = 'https://opendata.dwd.de/weather/nwp/icon/grib/12/u_10m/'
url12v = 'https://opendata.dwd.de/weather/nwp/icon/grib/12/v_10m/'
url18u = 'https://opendata.dwd.de/weather/nwp/icon/grib/18/u_10m/'
url18v = 'https://opendata.dwd.de/weather/nwp/icon/grib/18/v_10m/'

url_list = [url00u, url00v, url06u, url06v, url12u, url12v, url18u, url18v]


def main():
    for url in url_list:

        request = requests.get(url, allow_redirects=True)
        soup = BeautifulSoup(request.content, 'html.parser')
        links = soup.findAll('a')
        links = [url + link['href'] for link in links
                 if link['href'].startswith('icon')]

        for link in links:
            if int(link.split('_')[-3]) > 100:
                continue
            file_name = os.path.join('/data/L2/Wind/ICON',
                                     link.split('/')[-1])

            if not os.path.exists(file_name):
                print("Downloading file: ", file_name)
                req = requests.get(link, stream=True)
                with open(file_name, 'wb') as fl:
                    for chunk in req.iter_content(chunk_size=1024*1024):
                        if chunk:
                            fl.write(chunk)
                req = None


if __name__ == '__main__':
    main()
