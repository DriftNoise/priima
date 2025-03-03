```
Copyright 2025, Drift+Noise GmbH

PRIIMA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.
```

# PRIIMA

PRIIMA stands for PRedicted Ice IMAges and combines (1) a static Sentinel-1 SAR
image of sea ice with (2) a modelled sea ice drift forecast to produce (3) a series of
manipulated Sentinel-1 images visualizing potential future positions of sea-ice.
Merging all the produced forecast Sentinel-1 images into an animated
gif or a video file creates animations of sea-ice movement. You
can think of it as something similar to the moving rain radar images
in a modern weather app, showing you a predicted rain pattern of the
next hours. An example PRIIMA animation is shown
[here](https://driftnoise.com/priima.html).

Beside the actual Sentinel-1 image PRIIMA needs either
Arctic sea-ice drift forecasts from the Arctic
[TOPAZ4 system](https://www.researchgate.net/publication/258687774_TOPAZ4_An_ocean-sea_ice_data_assimilation_system_for_the_North_Atlantic_and_Arcticasystem)
available via the [Copernicus Marine Service](https://marine.copernicus.eu/)
or for both polar regions - a wind forecast from the global
[ICON weather model](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)
provided by the DWD.

PRIIMA was initially built during an ESA kick start initiative with the same name.

More information on the involved algorithms and methods are explained in the
[algorithm explanation and basic principles documentation](docs/algorithm-explanation.md).

## Usage

From a local development environment:

```
$ source venv/bin/activate
$ PYTHONPATH=$PWD python priima/image_forecast.py --order-file orders/myorder
```

or equivalently from a Docker container:

```
$ docker-compose run --volume /path/to/local/data/dir:/data --volume /path/to/order/dir:/orders \
    priima python priima/image_forecast.py --order-file /orders/<order-file>.json --no-send
```

To run PRIIMA for a specific Sentinel-1 file, put the TIFF into the order's
output data path (`<data-path>/forecasts/<customer-name>/<order-name>`) and
then specify the file to run via the `--sentinel-1-file` command line
option:

```
$ PYTHONPATH=$PWD python priima/image_forecast.py --order-file orders/myorder \
    --sentinel-1-file <sentinel-1-filename>
```

## [Local development setup instructions](docs/dev-setup-instructions.md)
## [ICON usage notes](docs/icon-notes.md)

Where `<data-path>` is configured in the `priima.local_config.yml` via the
`data_path` key.

## [Algorithm explanation and basic principles](docs/algorithm-explanation.md)
