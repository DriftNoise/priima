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
import re
import zipfile
from pathlib import Path
from subprocess import call

from priima.config import Config


def compress_shapefiles(shapefiles, sentinel1_path):
    shapefiles_dir = base_shapefile_name(sentinel1_path)
    output_zip_file = os.path.join(
        Config.instance().output_dir, f"{shapefiles_dir}.zip"
    )
    with zipfile.ZipFile(output_zip_file, 'w') as fp:
        for shapefile in shapefiles:
            fp.write(
                shapefile,
                arcname=os.path.basename(shapefile),
                compress_type=zipfile.ZIP_DEFLATED
            )

    return output_zip_file


def base_shapefile_name(sentinel1_path):
    sentinel1_fname = Path(sentinel1_path).name
    match = re.search(
        r'_(?P<datestamp>\d{8})T(?P<timestamp>\d{6})_', sentinel1_fname
    )
    datestamp = match.group('datestamp')
    timestamp = match.group('timestamp')
    shapename = f'trajectories_{datestamp}_{timestamp}'

    return shapename


def convert_csv_to_shp(initial_sar_file):
    trajectory_filename = base_shapefile_name(initial_sar_file)
    trajectories_path = Config.instance().output_dir / trajectory_filename

    trajectories_path.mkdir(exist_ok=True)
    base_ogr_to_ogr_cmd = [
        'ogr2ogr -s_srs EPSG:4326 -t_srs EPSG:4326',
        '-oo X_POSSIBLE_NAMES=x* -oo Y_POSSIBLE_NAMES=y*',
        '-f "ESRI Shapefile"',
        '{0}/{1}.shp'.format(trajectories_path, trajectory_filename),
        '{0}.csv'.format(trajectories_path),
        ]

    base_ogr_to_ogr_cmd = " ".join(base_ogr_to_ogr_cmd)
    try:
        call(base_ogr_to_ogr_cmd, shell=True)
        files_to_zip = trajectories_path.glob('*')
        trajectories_filename = compress_shapefiles(
            files_to_zip, initial_sar_file
        )
    except Exception as error:
        print(error)

    return trajectories_filename
