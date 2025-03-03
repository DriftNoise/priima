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

import shutil
import tempfile
from pathlib import Path
from unittest import TestCase
from zipfile import ZipFile

from priima.shapefile import base_shapefile_name, compress_shapefiles


class TestShapefile(TestCase):
    def setUp(self):
        self.tempdir = Path(tempfile.mkdtemp())

    def test_base_shapefile_name_converts_sar_fname_to_shapefile_name(self):
        sentinel1_path = (
            "/some/path/S1_EW_HH_sub_20230325T211045_20230325T211345_"
            "order_name_30_roi.tiff"
        )

        self.assertEqual(
            "trajectories_20230325_211045",
            base_shapefile_name(sentinel1_path)
        )

    def test_compress_shapefiles_creates_expected_zip_file(self):
        sentinel1_path = (
            "/some/path/S1_EW_HH_sub_20230325T211045_20230325T211345_"
            "order_name_30_roi.tiff"
        )
        data_path = self.tempdir

        # create a set of files to compress in a very similar manner to how
        # this is used in a production setting: i.e. to bundle shapefiles.
        shapefiles_to_compress = [
            Path(self.tempdir, "trajectories.shx"),
            Path(self.tempdir, "trajectories.dbf"),
            Path(self.tempdir, "trajectories.prj"),
            Path(self.tempdir, "trajectories.shp"),
        ]
        for path in shapefiles_to_compress:
            path.touch()

        zip_file = compress_shapefiles(
            shapefiles_to_compress, data_path, sentinel1_path
        )

        self.assertEqual(
            Path(zip_file),
            Path(self.tempdir, "trajectories_20230325_211045.zip")
        )
        self.assertTrue(Path(zip_file).exists())

        files_in_archive = ZipFile(zip_file).namelist()

        self.assertEqual(
            sorted(files_in_archive),
            [
                "trajectories.dbf",
                "trajectories.prj",
                "trajectories.shp",
                "trajectories.shx",
            ]
        )

    def tearDown(self):
        if self.tempdir.exists():
            shutil.rmtree(self.tempdir)
