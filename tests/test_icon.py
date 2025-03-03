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
import tempfile
from datetime import datetime
from pathlib import Path
from unittest import TestCase

from priima.config import Config
from priima.icon import select_icon_files


class TestIconLibrary(TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.config = Config.instance()
        data_path = Path(self.tmpdir.name)
        Config.set_attribute("data_path", data_path)

        self.icon_files = list([
            "icon_global_icosahedral_single-level_"
            "2019012900_001_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012900_002_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012900_003_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012900_001_V_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012900_002_V_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012900_003_V_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012912_001_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012912_002_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012912_001_V_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012912_002_V_10M.grib2.bz2",
            ])
        for filename in self.icon_files:
            filename = os.path.join(data_path, filename)
            with open(filename, 'w') as fh:
                fh.write('test')

    def test_select_icon_files_selects_file_closest_to_model_run(self):
        time_step = datetime(2019, 1, 29, 14)
        icon_files = select_icon_files(time_step)
        icon_files = [os.path.basename(icf) for icf in icon_files]
        expected_icon_files = list([
            "icon_global_icosahedral_single-level_"
            "2019012912_002_U_10M.grib2.bz2",
            "icon_global_icosahedral_single-level_"
            "2019012912_002_V_10M.grib2.bz2"
            ])

        self.assertCountEqual(icon_files, expected_icon_files)

    def test_if_no_relevant_forecasts_exist_index_error_is_raising(self):
        time_step = datetime(2019, 1, 25, 14)
        with self.assertRaises(IndexError) as error:
            select_icon_files(time_step)

        self.assertTrue("No ICON forecast available, PRIIMA will exit!"
                        in error.exception.args)
