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

import os.path
import tempfile
import unittest

from priima.config import Config


class TestConfig(unittest.TestCase):
    """
    Tests of the Config() class
    """
    def setUp(self):
        self.temp_files = []

    def test_config_file_path_default_value(self):
        config = Config()
        self.assertEqual(config.config_file, "priima.config.yml")

    def test_loading_empty_config_file_handles_nonexistent_fields(self):
        config_file = self.create_config_file("paths:\n")
        config = Config(config_file=config_file, local_config_file="")

        self.assertEqual(config.paths.data_path, "data")
        self.assertEqual(config.paths.secrets_path, "secrets")
        self.assertEqual(config.paths.grid_path, "data/icon_grid")
        self.assertEqual(config.paths.icon_path, "data/L2/Wind/ICON")
        self.assertEqual(config.gcp_separation, 1200)

    def test_raise_error_when_config_file_nonexistent(self):
        with self.assertRaises(IOError):
            Config(config_file="nonexistent")

    def test_raise_error_when_config_file_cant_be_parsed(self):
        unparseable = self.create_config_file('[non-yaml')
        with self.assertRaises(ValueError):
            Config(config_file=unparseable)

    def test_params_set_to_configured_values(self):
        config_content = """\
        paths:
            data_path: "/my/test/data/path"
            secrets_path: "/my/test/secrets/path"
            grid_path: "/my/test/grid/path"
        gcp_separation: 1200
        """
        config_file = self.create_config_file(config_content)
        config = Config(config_file=config_file, local_config_file="")
        self.assertEqual(config.paths.data_path,
                         "/my/test/data/path")
        self.assertEqual(config.paths.secrets_path,
                         "/my/test/secrets/path")
        self.assertEqual(config.paths.grid_path,
                         "/my/test/grid/path")
        self.assertEqual(config.gcp_separation, 1200)

    def test_local_config_overrides_default_settings(self):
        base_config_text = """\
        paths:
            data_path: /my/test/data/path"
            secrets_path: /my/test/secrets/path"
            grid_path: /my/test/grid/path"
        gcp_separation: 1200
        """
        base_config = self.create_config_file(base_config_text)

        local_config_text = """\
        paths:
            data_path: "/my/local/data/path"
            secrets_path: "/my/local/secrets/path"
            grid_path: "/my/local/grid/path"
        gcp_separation: 800
        """
        local_config = self.create_config_file(local_config_text)

        config = Config(
            config_file=base_config, local_config_file=local_config)
        self.assertEqual(config.paths.data_path,
                         "/my/local/data/path")
        self.assertEqual(config.paths.secrets_path,
                         "/my/local/secrets/path")
        self.assertEqual(config.paths.grid_path,
                         "/my/local/grid/path")
        self.assertEqual(config.gcp_separation, 800)

    def test_only_one_config_instance_can_be_instantiated(self):
        config_content = """\
        paths:
            data_path: "/my/test/data/path"
        """
        config_file = self.create_config_file(config_content)
        config = Config(config_file=config_file, local_config_file="")

        self.assertEqual(config.paths.data_path, "/my/test/data/path")

        with self.assertRaises(Exception):
            Config(config_file=config_file)

    def test_instance_property_returns_current_instance(self):
        config_content = """\
        paths:
            data_path: "/my/test/data/path"
        """
        config_file = self.create_config_file(config_content)
        config = Config(config_file=config_file, local_config_file="")

        self.assertEqual(Config.instance, config)

    def create_config_file(self, content):
        temp_fh = tempfile.NamedTemporaryFile(mode="w", delete=False)
        temp_fh.write(content)
        temp_fh.flush()

        self.temp_files.append(temp_fh.name)

        return temp_fh.name

    def tearDown(self):
        Config._drop()  # pylint: disable=protected-access
        for fname in self.temp_files:
            if os.path.exists(fname):
                os.remove(fname)


if __name__ == "__main__":
    unittest.main()
