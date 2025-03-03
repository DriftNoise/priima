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

import errno
import os
from datetime import datetime, timedelta
from glob import glob
from subprocess import call

import numpy as np
from netCDF4 import Dataset

from icon_grid.POLAR_GRID.antarctic_grid_generator import make_antarctic_grid


def set_regridded_filename(filename):
    smoc_filename_regridded = filename.replace(".nc", "_regridded.nc")

    return smoc_filename_regridded


def regrid_to_polar(filename, order, config):
    """
    Makes external system call to run CDO and to convert
    the crs from a NetCDF file, from WGS84 to polar grid.
    The desired polar grid is defined by providing an example file
    from TOPAZ data. Then it renames using nco the longitude=>lon
    and latitude=>lat.
    """
    output_name = os.path.basename(set_regridded_filename(filename))
    base_dir = os.path.dirname(filename)
    out_dir = os.path.join(base_dir, "regridded")
    out_path = os.path.join(out_dir, output_name)
    if order.center[0] > 0:
        grid = os.path.join(
            config.paths.grid_path,
            "POLAR_GRID",
            "20180426_hr-metno-MODEL-topaz4-ARC-b20180424-fv02.0.nc")
        if not os.path.exists(grid):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                (f'{grid}. You are missing a needed grid file. '
                  'Refer to docs/icon-notes.md for help on obtaining it.'))
    else:
        grid = os.path.join(
            config.paths.grid_path,
            "POLAR_GRID",
            "antarctic_grid.nc")
        if not os.path.exists(os.path.dirname(grid)):
            os.makedirs(os.path.dirname(grid))
        if not os.path.exists(grid):
            make_antarctic_grid(grid)
        if not os.path.exists(grid):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                f'{grid}. The creation of the grid file went wrong.')

    call('cdo -r -remapbil,{grid} {filename} {out_path}'.format(
        grid=grid, filename=filename, out_path=out_path),
        shell=True)

    call('ncrename -v longitude,lon {out_path}'''.format(
        out_path=out_path), shell=True)
    call('ncrename -v latitude,lat {out_path}'''.format(
        out_path=out_path), shell=True)

    call('ncks -A -v x,y {grid} {out_path}'''.format(
        grid=grid, out_path=out_path), shell=True)

    return out_path


def select_smoc_files(time_range, config):
    search_path = os.path.join(
            config.paths.data_path, "L2", "OceanCurrents", "*.nc")
    smoc_files = glob(search_path)
    smoc_files.sort()
    start_time = int(time_range[0].strftime("%Y%m%d"))
    end_time = int(time_range[1].strftime("%Y%m%d"))
    all_days = [str(tm) for tm in range(start_time, end_time + 1, 1)]
    files_in_timerange = []
    for day in all_days:
        files_in_timerange.append(
            [fl for fl in smoc_files if day in os.path.basename(fl).split(
                "_")[1]][0])

    return files_in_timerange


def load_tidal_drift(smoc_files, time_range, order, config):
    hours = convert_to_smoc_times(time_range)

    first_time = True
    for sf in smoc_files:
        sf_regridded = regrid_to_polar(sf, order, config)
        ds = Dataset(sf_regridded)
        timetmp = ds.variables["time"][:]
        timetmp_list = timetmp.data.tolist()
        indexes = \
            [ind for ind, val in enumerate(timetmp_list) if val in set(hours)]

        if not indexes:
            break

        utmp = ds.variables["utide"][indexes, :, :]
        vtmp = ds.variables["vtide"][indexes, :, :]
        timetmp = timetmp[indexes]

        if first_time:
            utide = utmp.copy()
            vtide = vtmp.copy()
            time = timetmp.copy()
        else:
            utide = np.concatenate((utide, utmp), axis=0)
            vtide = np.concatenate((vtide, vtmp), axis=0)
            time = np.concatenate((time, timetmp), axis=0)
        first_time = False
    tide_data = {"utide": utide, "vtide": vtide, "time": time}

    return tide_data


def convert_to_smoc_times(time_range):
    """
    Gets a timerange list of datetime objects and converts it
    into hours since 1950 01 Jan 00:00:00
    """
    time_epoch = datetime(1950, 1, 1, 0, 0)
    start_time = time_range[0]
    end_time = time_range[1]
    start_hour = \
        int((start_time - time_epoch).total_seconds()/3600.) + 0.5
    end_hour = \
        int((end_time - time_epoch).total_seconds()/3600.) + 1.5
    hours = np.arange(start_hour, end_hour, 1)

    return hours


def smoc_times_to_datetimes(smoc_times):
    """
    Gets time object from SMOC data and converts it
    into datetime objects. Round half hours to the next lower hour.
    All SMOC times are at XX:30. The best is to round them, in order to be
    able to compare them later with the ICON times
    """
    epoch = datetime(1950, 1, 1, 0, 0)
    smoc_datetimes = [epoch + timedelta(hours=tm - 0.5) for tm in smoc_times]

    return smoc_datetimes
