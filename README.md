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

PRIIMA stands for **PR**edicted **I**ce **IMA**ges and combines a static SAR
or otherwise remotely sensed image of sea ice with a modelled sea ice drift
forecast to produce a series of manipulated images visualizing potential
future positions of sea-ice.
Merging all the produced forecast satellite images into an animated gif or a
video file creates animations of sea-ice movement. You can think of it as
something similar to the moving rain radar images in a modern weather app,
showing you a predicted rain pattern of the next hours. An example PRIIMA
animation is shown on the
[PRIIMA project overview](https://driftnoise.com/priima.html). The development
of PRIIMA was funded by the [ESA kickstart program](https://business.esa.int/news/kick-start-activities-new-funding-opportunity-for-innovative-applications-ideas).

Possible model inputs are either sea-ice drift forecasts from the pan-Arctic
[NeXtSIM system](https://data.marine.copernicus.eu/product/ARCTIC_ANALYSISFORECAST_PHY_ICE_002_011/description)
available via the Copernicus Marine Service or, the globally available wind
forecast from the
[ICON weather model](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)
provided by the German Weather Service DWD. While NeXtSIM forecasts the ice
drift over the next nine days, ICON only forecasts wind for about 3 days into
the future. NeXtSIM ice drift can be used as-is, while the wind values are
scaled to guesstimate the expected ice drift from them by applying the Nansen
rule (2.5% of wind speed and wind direction +/- 25° on Northern/Southern
hemisphere).

## Usage

Pull the repository and build a `priima` docker image by executing

```
docker compose build
```

### Preparations

Get **SAR imagery** in GRD format from your source of trust (for example the
[Copernicus Dataspace](https://browser.dataspace.copernicus.eu) holding
Sentinel-1 imagery or the [EODMS Portal](https://www.eodms-sgdot.nrcan-rncan.gc.ca/))
and (optional but recommended) transform the images into a suitable projection
(polar stereographic) by running a geometric correction on them. This has the
effect that the images will be correctly converted from radar geometry (GRD)
into a projected map format using a DEM. GDAL can also do the projection but
will not take into account SAR geometry artifacts and, more severely, GDAL is
flawed when it comes to the poles or to the antimeridian (180° E/W). The
`priima` program operates with polar stereographic projections and will use
GDAL to reproject incoming images. Lastly, it is important that the used
image file contains the acquisition date in the format `_<YYYYMMDD>T<hhmmss>`.

**Forecast ice drift or wind data** is also required. It needs to cover the
timespan of the image forecast, starting with the timestamp of the image and
going into the future for as many hours as needed. All forecast files need
to be present in a single directory. For NEXTSIM, files of the product
`cmems_mod_arc_phy_anfc_nextsim_hm` are required. ICON wind forecast files
can be found [here](https://opendata.dwd.de/weather/nwp/icon/grib/).

### Running a Forecast

A `priima` docker container can then be run like so:

```
docker compose run --rm --volume /input/model/data/:/model_data \
--volume /out/data/:/out_dir --volume /image/data:/image_data priima \
--image /image_data/s1_9572_arctic_20241115T060000.tif \
--forecast-duration 12 --data-source NEXTSIM
```

You'll see that there are three mounted **volumes**:
* the input model data directory needs to be mapped to `/model_data`. Data for
  the time in question needs to be already present and is not downloaded
  automatically.
* the output directory needs to be mapped to `/out_dir`. An appropriate output
  subdirectory will be created there.
* the directory containing the input image (specified with `--image`) needs to
  be mounted in the container

Followed by the docker image name `priima` (or a custom name you may have
given during build).

Furthermore, the following parameters are passed:
* `--image` the path of the input image, starting from the mount point of the
  volume
* `--forecast-duration` the forecast duration in hours
* `--data-source` either `NEXTSIM` or `ICON`, the model with which to forecast
  the ice drift

There are some other parameters possible for fine tuning results:
* `--gcp-separation` the distance between the ground control points to warp the
  image, in pixels. The default is `1000`.
* `--output-resolution` the image resolution of the output image, default is
  `100`.
* `--gtif-out` to output images in GeoTiff format instead of cloud-optimized
  GeoTiffs (default)
* `--video` will output a `.avi` video by combining all `.tif` images found in
  the output directory

## Further Documentation

* [Algorithm explanation and basic principles](docs/algorithm-explanation.md)
