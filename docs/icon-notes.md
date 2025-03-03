<!-- 
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
-->

# Notes on Using ICON

## In the Arctic

Some notes on getting PRIIMA to work with ICON data in the Arctic.

### ICON Data

* PRIIMA needs the wind data at the stated path (`<data-path>/L2/Wind/ICON/`)
* do not unzip the files! PRIIMA needs zipped files as input, because it
  cleans the directory of previously unzipped files on startup

### Current Grid Workarounds

* download `ICON_GLOBAL2WORLD_0125_EASY.tar.bz2` from
  [https://opendata.dwd.de/weather/lib/cdo/](https://opendata.dwd.de/weather/lib/cdo/)

```
$ mkdir -p /path/to/data/icon_grid  # ensure the icon_grid dir exists
$ wget -O /path/to/data/data/icon_grid/ICON_GLOBAL2WORLD_0125_EASY.tar.bz2 https://opendata.dwd.de/weather/lib/cdo/ICON_GLOBAL2WORLD_0125_EASY.tar.bz2
```

* extract the tar archive to the directory you set as your `<grid_path>` in
  the `priima.local_config.yml` config file:

```
$ cd /path/to/data/icon_grid
$ tar -xvjf ICON_GLOBAL2WORLD_0125_EASY.tar.bz2
```

* the path to the files `target_grid_world_0125.txt` and
  `weights_icogl2world_0125.nc` should now be `<grid_path>/ICON_GLOBAL2WORLD_0125_EASY/<file>`
* additionally we need a TOPAZ grid file for the Arctic that is hardcoded as
  `20180426_hr-metno-MODEL-topaz4-ARC-b20180424-fv02.0.nc` and needs to be
  placed in `<grid_path>/POLAR_GRID/`

```
$ mkdir -p /path/to/data/icon_grid/POLAR_GRID/
$ cp /path/to/data/L3/ice_drift_forecasts/topaz4/<YYYYMMDD>_hr-metno-MODEL-topaz4-ARC-b<YYYYMMDD>-fv02.0.nc \
     /path/to/data/icon_grid/POLAR_GRID/20180426_hr-metno-MODEL-topaz4-ARC-b20180424-fv02.0.nc
```

* this file has the typical TOPAZ4 file naming convention, so we copied a
  recent TOPAZ4 file to `<grid_path>/POLAR_GRID` and renamed it to
  `20180426_hr-metno-MODEL-topaz4-ARC-b20180424-fv02.0.nc` as expected by the
  code
* This is not completely validated yet but works until a better solution is found
