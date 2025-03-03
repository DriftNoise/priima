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

FROM python:3.7-buster

RUN apt-get update && apt-get install -y build-essential libpython3-dev virtualenv libgdal-dev gdal-bin nco cdo libffi-dev cmake libopencv-dev libboost-all-dev mailutils
RUN adduser --disabled-password --gecos "" priimauser

ENV APP_ROOT=/home/priimauser/priima
RUN mkdir -p $APP_ROOT && chown -R priimauser:priimauser $APP_ROOT
RUN mkdir -p /data && chown -R priimauser:priimauser /data
RUN mkdir -p $APP_ROOT/logs && chown -R priimauser:priimauser $APP_ROOT/logs
RUN mkdir -p $APP_ROOT/orders && chown -R priimauser:priimauser $APP_ROOT/orders
WORKDIR $APP_ROOT

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

COPY priima $APP_ROOT/priima
COPY bin $APP_ROOT/bin
COPY priima.config.yml $APP_ROOT/priima.config.yml
COPY matplotlibrc $APP_ROOT/matplotlibrc
COPY devops $APP_ROOT/devops
COPY base_requirements.txt $APP_ROOT/base_requirements.txt
COPY requirements.txt $APP_ROOT/requirements.txt
COPY Makefile $APP_ROOT/Makefile
COPY shapefiles $APP_ROOT/shapefiles
# The secrets file isn't copied because it is not currently intended to send
# data to users via email; this is the only current use of secrets
# information in PRIIMA.
# COPY docker/priima/secrets.json $APP_ROOT/secrets.json
COPY icon_grid $APP_ROOT/icon_grid
COPY fast_cast $APP_ROOT/fast_cast
# need to remove any local pre-existing builds, otherwise FastCast won't
# build inside Docker
RUN rm -rf $APP_ROOT/fast_cast/build
RUN chown -R priimauser:priimauser $APP_ROOT/fast_cast/

USER priimauser
RUN ./devops/bin/install-deps
RUN make fast_cast

ENV PATH $APP_ROOT/venv/bin:$PATH
ENV PYTHONPATH $APP_ROOT
