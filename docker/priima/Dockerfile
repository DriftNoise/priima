# Copyright 2025, Drift+Noise GmbH

# This file is part of PRIIMA.
# PRIIMA is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.

FROM python:3.10-bookworm AS base

RUN apt-get update && apt-get install -y gdal-bin=3.6.2+dfsg-1+b2 \
 libgdal-dev=3.6.2+dfsg-1+b2 nco=5.1.4-1+deb12u1 cdo=2.1.1-1+deb12u1 \
 python3-dev=3.11.2-1+b1 python3-venv=3.11.2-1+b1 libffi-dev=3.4.4-1 \
 libopencv-dev=4.6.0+dfsg-12 cmake=3.25.1-1 libboost-all-dev=1.74.0.3 \
 --no-install-recommends \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN adduser --disabled-password --gecos "" ossi

ENV APP_ROOT=/home/ossi/priima
RUN mkdir -p /data && chown -R ossi:ossi /data
WORKDIR $APP_ROOT
RUN mkdir logs && chown -R ossi:ossi .

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

COPY matplotlibrc $APP_ROOT/matplotlibrc
COPY devops $APP_ROOT/devops
COPY base_requirements.txt $APP_ROOT/base_requirements.txt
COPY requirements.txt $APP_ROOT/requirements.txt
COPY Makefile $APP_ROOT/Makefile
COPY fast_cast $APP_ROOT/fast_cast
RUN chown -R ossi:ossi $APP_ROOT/

USER ossi
RUN ./devops/bin/install-deps \
 && make fast_cast

FROM base AS create_landmasks

RUN pip install geopandas==1.0.1 \
 && mkdir $APP_ROOT/shapefiles \
 && wget -q -O $APP_ROOT/shapefiles/land-polygons-split-4326.zip \
 https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip \
 && unzip -o -q $APP_ROOT/shapefiles/land-polygons-split-4326.zip \
 -d $APP_ROOT/shapefiles/
RUN python3 devops/bin/create_landmask_shapefiles.py \
 --data-path $APP_ROOT/shapefiles

FROM base

COPY --from=create_landmasks $APP_ROOT/shapefiles/arctic_landmask.* \
 $APP_ROOT/shapefiles/
COPY --from=create_landmasks $APP_ROOT/shapefiles/antarctic_landmask.*\
 $APP_ROOT/shapefiles/

RUN mkdir $APP_ROOT/icon_grid \
 && wget -q -O $APP_ROOT/icon_grid/ICON_GLOBAL2WORLD_0125_EASY.tar.bz2 \
 https://opendata.dwd.de/weather/lib/cdo/ICON_GLOBAL2WORLD_0125_EASY.tar.bz2 \
 && tar -xjf $APP_ROOT/icon_grid/ICON_GLOBAL2WORLD_0125_EASY.tar.bz2 \
 -C $APP_ROOT/icon_grid/
RUN chown -R ossi:ossi $APP_ROOT/icon_grid/

COPY --chown=ossi:ossi priima $APP_ROOT/priima

ENV PATH=$APP_ROOT/venv/bin:$PATH
ENV PYTHONPATH=$APP_ROOT

ENTRYPOINT ["python3", "priima/image_forecast.py"]
