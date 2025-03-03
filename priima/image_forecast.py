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

import argparse
import csv
import datetime
import os
import re
import sys
from glob import glob
from pathlib import Path
from subprocess import call

import numpy as np
from osgeo import gdal
from slugify import slugify

from priima.config import Config
from priima.drift_handler import DriftHandler
from priima.gcp import create_gcp_from_ul_lr, export_gcp_2csv
from priima.icon import set_up_icon_subpath
from priima.order import load_order
from priima.preprocessing import warp_image
from priima.shapefile import base_shapefile_name, convert_csv_to_shp
from priima.smoc import load_tidal_drift, select_smoc_files
from priima.warp_forecast import (compute_averaged_area_drift,
                                  compute_time_range,
                                  create_transformation_matrix, reproject2roi,
                                  update_point_location)
from priima.wind import compute_drift_from_wind


def main():
    args = parse_arguments()

    order, order_data = load_order(order_file=args.order_file)

    if args.forecast_duration is not None:
        forecast_duration = float(args.forecast_duration)
    else:
        forecast_duration = order.forecast_duration

    data_source = order.data_source
    grid_method = order.grid_method

    # create working directory for this order's forecast
    customer_name_slug = slugify(order.customer_name)
    order_name_slug = slugify(order.order_name)
    config = Config()
    order_data_path = os.path.join(
        config.paths.data_path,
        'forecasts', customer_name_slug, order_name_slug)
    if not os.path.exists(order_data_path):
        os.makedirs(order_data_path)

    raw_scene_full_path = Path(
        order_data_path, args.sentinel_1_file
    ).as_posix()
    initial_sar_file = reproject2roi(raw_scene_full_path, order)
    ds = gdal.Open(initial_sar_file)
    proj = ds.GetProjection()
    gcp_list = create_gcp_from_ul_lr(ds, proj, grid_method, order_data)

    export_gcp_2csv(gcp_list)

    time_range = compute_time_range(initial_sar_file, forecast_duration)

    # for local point warp and time varying
    dt = (time_range[1] - time_range[0]).total_seconds()/60./60.
    dt_hours = int(dt // 1)
    forecast_steps = range(dt_hours)  # + [dt_decimal]
    gcp_list_dynamic = gcp_list[:]
    drifted_coord = []

    if data_source in ('TOPAZ', 'NEXTSIM'):

        avrg_drift = []
        for gld in gcp_list_dynamic:
            drifted_coord.append([gld.GCPX, gld.GCPY, time_range[0]])

        ncfile_list = select_drift_files(time_range, order)
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
                drift_handler = DriftHandler(ncfile, time_range_step, order,
                                             gcp_list_dynamic)
                drift = \
                    drift_handler.compute_drift()
                avrg_drifti = compute_averaged_area_drift(drift)
                avrg_drifti.insert(
                    0, time_range_step[0].strftime('%Y%m%dT%H'))
                avrg_drift.append(avrg_drifti)

            except ValueError:
                continue
            gcp_list_dynamic = update_point_location(
                order, gcp_list_dynamic, drift)

            for gld in gcp_list_dynamic:
                drifted_coord.append([gld.GCPX, gld.GCPY, time_range_step[0]])

            gcp_list = create_gcp_from_ul_lr(
                ds, proj, grid_method, order_data)

            starting_pos, ending_pos, \
                lon_projected_list, lat_projected_list = \
                create_transformation_matrix(
                    ds, gcp_list, gcp_list_dynamic, grid_method)

            # print("Drift is: ", drift)
            pixels_drift = {'start': starting_pos, 'end': ending_pos,
                            'lon_projected': lon_projected_list,
                            'lat_projected': lat_projected_list}
            forecast_stepi = '{0:d}_{1:02d}'.format(int(forecast_duration),
                                                    it)
            warp_image(forecast_stepi, initial_sar_file, pixels_drift,
                       data_source, grid_method, order_data_path,
                       args.store_training_data)
        shapename = \
            '{}.csv'.format(base_shapefile_name(initial_sar_file))
        with open(os.path.join(order_data_path,
                               shapename), mode='w') as fh:
            writer = csv.writer(fh, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['x', 'y', 'timestep'])
            for line in drifted_coord:
                writer.writerow(line)

        save_averaged_drift(order_data_path, avrg_drift)

    elif data_source == 'ICON':
        # clean up ICON folders
        set_up_icon_subpath(config.paths.icon_path)

        call(("rm {}/converted/*".format(config.paths.icon_path)), shell=True)
        call(("rm {}/*.grib2".format(config.paths.icon_path)), shell=True)
        tide_data = 0
        if order.use_tides:
            smoc_files = select_smoc_files(time_range, config)
            tide_data = load_tidal_drift(
                smoc_files, time_range, order, config)
        for gld in gcp_list_dynamic:
            drifted_coord.append([gld.GCPX, gld.GCPY, time_range[0]])

        for it in forecast_steps:
            it = it + 1
            print(f"\nForecast step {it} of {dt_hours}")
            time_step = time_range[0] + datetime.timedelta(hours=it)
            drift = compute_drift_from_wind(
                config, time_step, data_source, order,
                gcp_list_dynamic, tide_data)

            gcp_list_dynamic = update_point_location(
                order, gcp_list_dynamic, drift)
            for gld in gcp_list_dynamic:
                drifted_coord.append([gld.GCPX, gld.GCPY, time_step])

            gcp_list = create_gcp_from_ul_lr(ds, proj, grid_method, order_data)
            starting_pos, ending_pos, lon_projected_list, \
                lat_projected_list = \
                create_transformation_matrix(
                        ds, gcp_list, gcp_list_dynamic, grid_method)

            pixels_drift = {'start': starting_pos, 'end': ending_pos,
                            'lon_projected': lon_projected_list,
                            'lat_projected': lat_projected_list}

            forecast_stepi = '{0:d}_{1:02d}'.format(int(forecast_duration),
                                                    it)
            warp_image(forecast_stepi, initial_sar_file, pixels_drift,
                       data_source, grid_method, order_data_path,
                       args.store_training_data)
            shapename = \
                '{}.csv'.format(base_shapefile_name(initial_sar_file))
        with open(os.path.join(order_data_path,
                               shapename), mode='w') as fh:
            writer = csv.writer(fh, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['x', 'y', 'timestep'])
            for line in drifted_coord:
                writer.writerow(line)
    else:
        raise ValueError("Unknown data source {}".format(data_source))

    ds = None
    # Create Shapefiles for easy viewing
    convert_csv_to_shp(order_data_path, initial_sar_file)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate the ice image of tomorrow')
    parser.add_argument(
        '--forecast-duration',
        help='Number of hours in the future to run forecast')
    parser.add_argument(
        '--order-file',
        help='File describing parameters of an order to process')

    parser.add_argument(
        '--store-training-data', default=False, action='store_true',
        help='Store the pts*.txt files needed \
            for neural network model training')
    parser.add_argument(
        '--sentinel-1-file',
        help='Name of the Sentinel-1 file to process')
    args = parser.parse_args()

    return args


def select_drift_files(time_range, order):
    """
    Selects appropriate drift NetCDF files to run the forecast.
    It always selects the forecast days that are closer to the conducted date.

    :param time_range:
        The time range for which the relevant topaz files will be selected
    :type time_range:
        list of datetime.datetime objects
    :param order:
        An order object
    :type order:
        priima.order.Order object
    """
    start_date = time_range[0].replace(hour=0, minute=0)
    end_date = time_range[1].replace(hour=0, minute=0)
    date_range = [start_date + datetime.timedelta(days=num_days)
                  for num_days in range(0, (end_date - start_date).days + 1)]
    if order.data_source == "TOPAZ":
        ice_drift_file_glob = os.path.join(
            Config.instance.paths.data_path,
            "L3", "ice_drift_forecasts", "topaz4", '*.nc')
    elif order.data_source == "NEXTSIM":
        ice_drift_file_glob = os.path.join(
            Config.instance.paths.data_path,
            "L3", "ice_drift_forecasts", "nextsim", '*.nc')
    else:
        print("Unknown data_source, should be TOPAZ/NEXTSIM")

    all_ice_drift_files = glob(ice_drift_file_glob)
    # filter for only TOPAZ4 files; these files start with a datestamp in
    # the format YYYYMMDD, hence we filter for files which only start with
    # 8 digits.  The aggregated TOPAZ4 drift files all start with the word
    # `aggregated`, hence they are filtered out of the list of all ice drift
    # files found in the given path.
    topaz_list = [
        fname for fname in all_ice_drift_files
        if re.match(r'^\d{8}', os.path.basename(fname))]

    topaz_list_forecast_day = \
        [os.path.basename(topaz_file).split("_")[0]
         for topaz_file in topaz_list]
    topaz_list_day_of_forecast_run = \
        [os.path.basename(topaz_file).split("-")[-2]
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
            ds = gdal.Open(fl, gdal.GA_ReadOnly)
            ds_stats = ds.GetRasterBand(1).GetStatistics(0, 1)
            if np.allclose(ds_stats[0:2], 0.0):
                print("File contains only 0s, PRIIMA will STOP!!")
                sys.exit(10)
        except Exception as ex:
            print("Error in {0}. File seems to be not valid."
                  "PRIIMA will STOP!!!!".format(ex))
            sys.exit(20)


def save_averaged_drift(order_data_path, avrg_drift):
    with open(os.path.join(order_data_path,
                           'averaged_drift.csv'), mode='w') as drfile:
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
        match = regex.search(os.path.basename(fname))
        datestamp = match.group('datestamp')
        date = datetime.datetime.strptime(datestamp, "%Y%m%d").date()
        file_dates.append(date)

    return file_dates


if __name__ == "__main__":
    main()
