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

import bz2
import datetime
import errno
import glob
import os
from pathlib import Path
from subprocess import call

import numpy as np
from netCDF4 import Dataset

from priima.config import Config
from priima.grid_generator import make_antarctic_grid, make_arctic_grid


def process_icon_fields(time_step):

    converted_icon_path = os.path.join("/var", "tmp", "converted")
    Path(converted_icon_path).mkdir(parents=True, exist_ok=True)

    target_grid = os.path.join(
        Config.instance().grid_path,
        "ICON_GLOBAL2WORLD_0125_EASY",
        "target_grid_world_0125.txt")
    if not os.path.exists(target_grid):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT),
            (f'{target_grid}. You are missing a needed grid file. '
              'Refer to docs/icon-notes.md for help on obtaining it.'))
    weights = os.path.join(
        Config.instance().grid_path,
        "ICON_GLOBAL2WORLD_0125_EASY",
        "weights_icogl2world_0125.nc")
    if not os.path.exists(weights):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT),
            (f'{weights}. You are missing a needed grid file. '
              'Refer to docs/icon-notes.md for help on obtaining it.'))

    grib_files = select_icon_files(time_step)
    print("Forecast will use following ICON files: {}".format(grib_files))

    decompress_tmp_dir = Path("/var") / "tmp" / "decompressed"
    decompress_tmp_dir.mkdir(exist_ok=True, parents=True)
    out_list = []
    for ifile in grib_files:
        decompressed_file = decompress_bz2(
            input_file=ifile, output_path=decompress_tmp_dir)
        outfile = decompressed_file.replace('grib2', 'nc')
        outfile = os.path.join(converted_icon_path, os.path.basename(outfile))

        cmd = ('cdo -f nc remap,{target_grid},{weights} {ifile} '
               '{outfile}'.format(target_grid=target_grid,
                                  weights=weights, ifile=decompressed_file,
                                  outfile=outfile))
        call(cmd, shell=True)
        out_list.append(outfile)

    combined_file = combine_ncfiles(out_list)
    regridded_file = regrid_regular2polar(
        combined_file)
    rotate_velocity(regridded_file)

    return regridded_file


def combine_ncfiles(out_list):
    u_files = [fl for fl in out_list if '_U_' in fl]
    v_files = [fl for fl in out_list if '_V_' in fl]
    start_time = u_files[0].rsplit('_')[-2]
    end_time = u_files[len(u_files) - 1].rsplit('_')[-2]
    output_name_u = \
        os.path.join(os.path.dirname(
            u_files[0]), 'ICON_iko_single_level_elements_world_U_10M_' +
                          start_time + '_' + end_time + '.nc')
    output_name_v = \
        os.path.join(os.path.dirname(
            v_files[0]), 'ICON_iko_single_level_elements_world_V_10M_' +
                         start_time + '_' + end_time + '.nc')

    u_files = ' '.join(u_files)
    v_files = ' '.join(v_files)

    call("rm {0}".format(output_name_u), shell=True)
    call("rm {0}".format(output_name_v), shell=True)
    call('cdo mergetime {u_files} {output_name_u}'.format(
        u_files=u_files, output_name_u=output_name_u),
         shell=True)
    call('cdo mergetime {v_files} {output_name_v}'.format(
        v_files=v_files, output_name_v=output_name_v),
         shell=True)
    output_name_all = output_name_v.replace('_V_', '_combined_')
    call("rm {0}".format(output_name_all), shell=True)
    call('cdo merge {output_name_u} {output_name_v} {output_name_all}'
         ''.format(output_name_u=output_name_u, output_name_v=output_name_v,
                   output_name_all=output_name_all), shell=True)

    return output_name_all


def regrid_regular2polar(filename):
    """
    Makes external system call to run CDO and to convert
    the crs from a NetCDF file, from regular to polar grid.
    The desired polar grid is defined by providing an example file
    from TOPAZ data. Then it renames using nco the longitude=>lon
    and latitude=>lat.
    """
    output_name = filename.replace('.nc', '_regridded.nc')
    if Config.instance().center[0] > 0:
        grid = os.path.join(
            Config.instance().grid_path,
            "arctic_grid.nc")
        if not os.path.exists(grid):
            make_arctic_grid(grid)
        if not os.path.exists(grid):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                f'{grid}. The creation of the grid file went wrong.')
    else:
        grid = os.path.join(
            Config.instance().grid_path,
            "antarctic_grid.nc")
        if not os.path.exists(grid):
            make_antarctic_grid(grid)
        if not os.path.exists(grid):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                f'{grid}. The creation of the grid file went wrong.')

    call('cdo -r -remapbil,{grid} {filename} {output_name}'.format(
        grid=grid, filename=filename, output_name=output_name),
        shell=True)

    call('ncrename -v longitude,lon {output_name}'''.format(
        output_name=output_name), shell=True)
    call('ncrename -v latitude,lat {output_name}'''.format(
        output_name=output_name), shell=True)

    call('ncks -A -v x,y {grid} {output_name}'''.format(
        grid=grid, output_name=output_name), shell=True)

    return output_name


def rotate_velocity(path):
    """
    Rotate winds from earth relative to polar grid relative
    some information is here:
    https://www-k12.atmos.washington.edu/~ovens/wrfwinds.html
    https://journals.ametsoc.org/view/journals/apme/47/11/2008jamc1746.1.xml
    https://apps.ecmwf.int/codes/grib/format/grib1/grids/5/
    """
    from priima.wind import prepare_wind_dict_from_icon

    wind_dict = prepare_wind_dict_from_icon(path)
    ds = Dataset(path, 'r+')
    for timeind in range(len(wind_dict["time"])):
        uin = wind_dict["uice"][timeind, ::]
        vin = wind_dict["vice"][timeind, ::]
        earth_lons = wind_dict["lon"]

        # The LOV value (e.g. - -95.0) (single value in degrees)
        # "Lov = orientation of the grid; i.e. the east longitude value of
        # the meridian which is parallel to the Y-axis (or columns of
        # the grid) along which latitude increases as the Y-coordinate
        # increases (the orientation longitude may or may not appear
        # on a particular grid). In our case, this is -45 for the Arctic and
        # 0 for the Antarctic
        if Config.instance().center[0] > 0:
            lov_lon = -45
        else:
            lov_lon = 0

        if lov_lon > 0.:
            lov_lon = lov_lon-360.
        dtr = np.pi/180.0             # Degrees to radians

        # Compute rotation constant which is also
        # known as the Lambert cone constant.  In the case
        # of a polar stereographic projection, this is one.
        # (positive one for Arctic, negative one for Antarctic)
        if Config.instance().center[0] > 0:
            rotcon_p = 1.0
        else:
            rotcon_p = -1.0

        angles = rotcon_p*(earth_lons-lov_lon)*dtr
        sinx2 = np.sin(angles)
        cosx2 = np.cos(angles)

        # Return the grid relative winds
        uout = cosx2*uin-sinx2*vin
        vout = sinx2*uin+cosx2*vin
        ds['10u'][timeind, ::] = uout
        ds['10v'][timeind, ::] = vout
    ds.close()


def select_icon_files(time_step):
    """
    Selects relevant files which are closer in time to the model run.
    The model runs 4 times per day at 00, 06, 12, 18 hours.
    """
    pattern = os.path.join(Config.instance().data_path, 'icon_global_*bz2')
    files = glob.glob(pattern)
    timestamp = [fl.split('/')[-1].split('_')[-4] for fl in files]
    forchours = [fl.split('/')[-1].split('_')[-3].split('.')[0]
                 for fl in files]
    component = [fl.split('/')[-1].split('_')[-2].split('.')[0]
                 for fl in files]
    timestamp = zip(timestamp, forchours, component)

    # sort everything
    timestamp = sorted(timestamp, key=lambda x: (x[0], x[1]))
    component = [cmpn[2] for cmpn in timestamp]
    hours = [int(day[0][8:10]) for day in timestamp]
    days = [int(day[0][6:8]) for day in timestamp]
    months = [int(month[0][4:6]) for month in timestamp]
    years = [int(year[0][0:4]) for year in timestamp]
    forchours = [datetime.timedelta(hours=int(hrs[1])) for hrs in timestamp]

    timestamp_forecasted = []
    for tm in range(len(days)):
        timestamp_forecasted.append(
            datetime.datetime(
                years[tm], months[tm], days[tm], hours[tm]) + forchours[tm])

    timestamp_filtered = []
    idd = [ind for ind, value in enumerate(timestamp_forecasted)
           if value == time_step]
    if len(idd) == 0:
        raise IndexError("No ICON forecast available, PRIIMA will exit!")

    timestamp_filtered_tmp = [timestamp[ii] for ii in idd]
    hrs = [int(hr) for hr in np.asarray(timestamp_filtered_tmp)[:, 1]]
    closest2start_index = np.where(np.asarray(hrs).min() == hrs)[0]

    for ii in range(len(closest2start_index)):
        timestamp_filtered.append(
            timestamp_filtered_tmp[closest2start_index[ii]])

    files_filtered = []
    for tm, hr, cm in timestamp_filtered:
        date_suffix = "{tm}_{hr}".format(tm=tm, hr=hr)
        ind = [ind for ind, value in enumerate(files)
               if date_suffix in value and cm in value][0]
        files_filtered.append(files[ind])

    return files_filtered


def convert_icon_time_to_datetime(icon_dict):
    time = []
    time_since = icon_dict['time_since']
    date_since = time_since.split(' ')[2]
    hms_time_since = time_since.split(' ')[3]
    time_since_concatenated = date_since + " " + hms_time_since
    epoch = datetime.datetime.strptime(time_since_concatenated,
                                       "%Y-%m-%d %H:%M:%S")
    for tm in icon_dict['time']:
        time.append(epoch + datetime.timedelta(minutes=tm))

    return time


def decompress_bz2(*, input_file, output_path):
    decompressed_file = output_path / Path(input_file).with_suffix('').name
    with open(str(decompressed_file), 'wb') as new_file:
        decompressor = bz2.BZ2Decompressor()
        compr_file = open(input_file, 'rb').read()
        new_file.write(decompressor.decompress(compr_file))

    return str(decompressed_file)
