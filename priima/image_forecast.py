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

import csv
import datetime
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np
from osgeo import gdal

from priima.config import Config
from priima.drift_handler import DriftHandler
from priima.gcp import create_gcp_from_ul_lr, export_gcp_2csv
from priima.geo_tools import get_center_coordinate
from priima.preprocessing import warp_image
from priima.shapefile import base_shapefile_name, convert_csv_to_shp
from priima.video import create_video
from priima.warp_forecast import (compute_averaged_area_drift,
                                  compute_time_range,
                                  create_transformation_matrix, reproject2roi,
                                  update_point_location)
from priima.wind import compute_drift_from_wind


def main():
    center_coords = get_center_coordinate()
    Config.set_attribute('center', [center_coords[0], center_coords[1]])

    output_dir = create_output_directory(
        image_fname=Path(Config.instance().image),
        data_source=Config.instance().data_source,
        forecast_duration=Config.instance().forecast_duration
    )
    Config.set_attribute("output_dir", output_dir)

    initial_sar_file = reproject2roi(Path(Config.instance().image))
    ds = gdal.Open(str(initial_sar_file))
    proj = ds.GetProjection()
    gcp_list = create_gcp_from_ul_lr(ds, proj)

    export_gcp_2csv(gcp_list)

    time_range = compute_time_range(
        initial_sar_file, Config.instance().forecast_duration)

    # for local point warp and time varying
    dt = (time_range[1] - time_range[0]).total_seconds()/60./60.
    dt_hours = int(dt // 1)
    forecast_steps = range(dt_hours)  # + [dt_decimal]
    gcp_list_dynamic = gcp_list[:]
    drifted_coord = []

    if Config.instance().data_source in ('TOPAZ', 'NEXTSIM'):

        avrg_drift = []
        for gld in gcp_list_dynamic:
            drifted_coord.append([gld.GCPX, gld.GCPY, time_range[0]])

        ncfile_list = select_drift_files(time_range)
        file_dates = netcdf_file_dates(ncfile_list)
        for it in forecast_steps:
            it = it + 1
            progress_msg = (
                "Forecast step: {0} of {1}"
                ).format(it, len(forecast_steps))
            print(progress_msg)
            time_range_step = time_range[0] + datetime.timedelta(hours=it)
            file_index = file_dates.index(time_range_step.date())
            ncfile = ncfile_list[file_index]
            time_range_step = [time_range_step, time_range_step]
            try:
                drift_handler = DriftHandler(
                    ncfile, time_range_step, gcp_list_dynamic
                )
                drift = \
                    drift_handler.compute_drift()
                avrg_drifti = compute_averaged_area_drift(drift)
                avrg_drifti.insert(
                    0, time_range_step[0].strftime('%Y%m%dT%H'))
                avrg_drift.append(avrg_drifti)

            except ValueError:
                continue
            gcp_list_dynamic = update_point_location(gcp_list_dynamic, drift)

            for gld in gcp_list_dynamic:
                drifted_coord.append([gld.GCPX, gld.GCPY, time_range_step[0]])

            gcp_list = create_gcp_from_ul_lr(ds, proj)

            starting_pos, ending_pos, \
                lon_projected_list, lat_projected_list = \
                create_transformation_matrix(ds, gcp_list, gcp_list_dynamic)

            # print("Drift is: ", drift)
            pixels_drift = {'start': starting_pos, 'end': ending_pos,
                            'lon_projected': lon_projected_list,
                            'lat_projected': lat_projected_list}
            forecast_stepi = '{0:d}_{1:02d}'.format(
                Config.instance().forecast_duration, it
            )
            warp_image(
                forecast_step=forecast_stepi,
                forecast_duration=Config.instance().forecast_duration,
                image_path=initial_sar_file,
                pixels_drift=pixels_drift,
            )
        shapename = \
            '{}.csv'.format(base_shapefile_name(initial_sar_file))
        with open(Config.instance().output_dir / shapename, mode='w') as fh:
            writer = csv.writer(fh, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['x', 'y', 'timestep'])
            for line in drifted_coord:
                writer.writerow(line)

        save_averaged_drift(avrg_drift)

    elif Config.instance().data_source == 'ICON':
        for gld in gcp_list_dynamic:
            drifted_coord.append([gld.GCPX, gld.GCPY, time_range[0]])

        for it in forecast_steps:
            it = it + 1
            print(f"\nForecast step {it} of {dt_hours}")
            time_step = time_range[0] + datetime.timedelta(hours=it)
            drift = compute_drift_from_wind(time_step, gcp_list_dynamic)

            gcp_list_dynamic = update_point_location(gcp_list_dynamic, drift)
            for gld in gcp_list_dynamic:
                drifted_coord.append([gld.GCPX, gld.GCPY, time_step])

            gcp_list = create_gcp_from_ul_lr(ds, proj)
            starting_pos, ending_pos, lon_projected_list, \
                lat_projected_list = \
                create_transformation_matrix(ds, gcp_list, gcp_list_dynamic)

            pixels_drift = {'start': starting_pos, 'end': ending_pos,
                            'lon_projected': lon_projected_list,
                            'lat_projected': lat_projected_list}

            forecast_stepi = '{0:d}_{1:02d}'.format(
                Config.instance().forecast_duration, it
            )
            warp_image(
                forecast_step=forecast_stepi,
                forecast_duration=Config.instance().forecast_duration,
                image_path=initial_sar_file,
                pixels_drift=pixels_drift,
            )
            shapename = \
                '{}.csv'.format(base_shapefile_name(initial_sar_file))
        with open(Config.instance().output_dir / shapename, mode='w') as fh:
            writer = csv.writer(fh, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['x', 'y', 'timestep'])
            for line in drifted_coord:
                writer.writerow(line)
    else:
        raise ValueError(
            "Unknown data source {}".format(Config.instance().data_source)
        )

    ds = None
    # Create Shapefiles for easy viewing
    convert_csv_to_shp(initial_sar_file)
    # Make a video if requested
    if Config.instance().video:
        video_fpath = create_video(image_dir=Config.instance().output_dir)
        log_msg = "Created video output here: %s"
        print(log_msg % video_fpath)


def create_output_directory(
    *,
    image_fname: Path,
    data_source: str,
    forecast_duration: int
):
    """Builds the path to and creates an output directory"""
    match = re.search(
        r"_(?P<start_time>\d{8}T\d{6})", str(image_fname)
    )
    try:
        start_time_string = match.group('start_time')
    except AttributeError as exc:
        err_msg = (
            "No date of format _<YYYYMMDD>T<hhmmss> found in filename %s"
        )
        raise ValueError(err_msg % image_fname) from exc
    sat = "sat"
    if "s1" in Path(Config.instance().image).name.lower():
        sat = "s1"
    elif "rcm" in Path(Config.instance().image).name.lower():
        sat = "rcm"
    output_dir_name = (
        f"{sat}_{start_time_string}_priima_{data_source}_"
        f"{forecast_duration}h"
    )
    output_dir = Config.instance().out_path / output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def select_drift_files(time_range):
    """
    Selects appropriate drift NetCDF files to run the forecast.
    It always selects the forecast days that are closer to the conducted date.

    :param time_range:
        The time range for which the relevant topaz files will be selected
    :type time_range:
        list of datetime.datetime objects
    """
    start_date = time_range[0].replace(hour=0, minute=0)
    end_date = time_range[1].replace(hour=0, minute=0)
    date_range = [start_date + datetime.timedelta(days=num_days)
                  for num_days in range(0, (end_date - start_date).days + 1)]
    if Config.instance().data_source in ("TOPAZ", "NEXTSIM"):
        ice_drift_file_glob = str(Config.instance().data_path / '*.nc')
    else:
        print("Unknown data_source, should be TOPAZ/NEXTSIM")

    all_ice_drift_files = [Path(file) for file in glob(ice_drift_file_glob)]
    # filter for only TOPAZ4 files; these files start with a datestamp in
    # the format YYYYMMDD, hence we filter for files which only start with
    # 8 digits.  The aggregated TOPAZ4 drift files all start with the word
    # `aggregated`, hence they are filtered out of the list of all ice drift
    # files found in the given path.
    topaz_list = [
        fname for fname in all_ice_drift_files
        if re.match(r'^\d{8}', fname.name)]

    topaz_list_forecast_day = \
        [topaz_file.name.split("_")[0]
         for topaz_file in topaz_list]
    topaz_list_day_of_forecast_run = \
        [topaz_file.name.split("-")[-2]
         for topaz_file in topaz_list]

    matched_files = []
    for date in date_range:
        indexes = [ind for ind, val in enumerate(topaz_list_forecast_day)
                   if val == date.strftime('%Y%m%d')]
        if len(indexes) == 1:
            matched_files.append(topaz_list[indexes[0]])
        elif len(indexes) > 1:
            # select all data that they match for a given forecast day and
            # create a subset from the full list
            topaz_list_subset = \
                [topaz_list[index] for index in indexes]
            day_of_forecast_run_indexes = \
                [topaz_list_day_of_forecast_run[index] for index in indexes]
            max_day_of_forecast_run_index = \
                day_of_forecast_run_indexes.index(
                        max(day_of_forecast_run_indexes))
            # from the subset select the most recent forecast
            matched_files.append(topaz_list_subset
                                 [max_day_of_forecast_run_index])
        else:
            continue

    if len(matched_files) == 0:
        error_msg = (
            "No drift files found matching date in image to be projected")
        raise ValueError(error_msg)

    print("Forecast will use following NetCDF files: {}".format(matched_files))
    return matched_files


def file_validation(filenames):
    for fl in filenames:
        try:
            ds = gdal.Open(str(fl), gdal.GA_ReadOnly)
            ds_stats = ds.GetRasterBand(1).GetStatistics(0, 1)
            if np.allclose(ds_stats[0:2], 0.0):
                print("File contains only 0s, PRIIMA will STOP!!")
                sys.exit(10)
        except Exception as ex:
            print("Error in {0}. File seems to be not valid."
                  "PRIIMA will STOP!!!!".format(ex))
            sys.exit(20)


def save_averaged_drift(avrg_drift):
    with open(
        Config.instance().output_dir / 'averaged_drift.csv', mode='w'
            ) as drfile:
        dr_writer = csv.writer(
            drfile, delimiter=',',
            quotechar='"', quoting=csv.QUOTE_MINIMAL)
        dr_writer.writerow(
            ['Timestamp', 'Drift Magnitude', 'Drift Direction'])
        for di in avrg_drift:
            dr_writer.writerow(di)


def netcdf_file_dates(netcdf_file_list):
    regex = re.compile(r'^(?P<datestamp>\d{8})_.*\.nc$')
    file_dates = []
    for fname in netcdf_file_list:
        match = regex.search(fname.name)
        datestamp = match.group('datestamp')
        date = datetime.datetime.strptime(datestamp, "%Y%m%d").date()
        file_dates.append(date)

    return file_dates


if __name__ == "__main__":
    main()
