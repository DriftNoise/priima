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

"""
Config: a module containing the Config class
"""

import logging
import os.path

import yaml

LOG = logging.getLogger('priima.config')


class PathConfig:
    """
    Container to hold the path configuration data
    """
    def __init__(self, path_config=None):
        self.data_path = "data"
        self.secrets_path = "secrets"
        self.grid_path = "data/icon_grid"
        if path_config is not None:
            if "data_path" in path_config:
                self.data_path = path_config["data_path"]
            if "secrets_path" in path_config:
                self.secrets_path = path_config["secrets_path"]
            if "grid_path" in path_config:
                self.grid_path = path_config["grid_path"]
            else:
                self.grid_path = os.path.join(self.data_path, "icon_grid")

        self.icon_path = os.path.join(self.data_path, "L2", "Wind", "ICON")


class MetaConfig(type):
    """
    Metaclass for construction of the Config as a Singleton

    See the text on metaclasses on this site
    https://www.python-course.eu/python3_metaclasses.php
    for more information about how to use a metaclass to create a singleton.

    By using a metaclass here we can expose the `instance()` method as a
    property rather than it having to be a method that one has to call.
    After all, it's a read-only value we want to have access to, hence a
    property makes sense in this case.
    """
    __instance = None

    def __call__(cls, *args, **kwargs):
        if cls.__instance is not None:
            raise Exception("Cannot instantiate Config more than once")
        else:
            cls.__instance = super(MetaConfig, cls).__call__(*args, **kwargs)
        return cls.__instance

    @property
    def instance(cls):
        return cls.__instance

    def _drop(cls):
        """
        Drop the instance (for testing purposes).

        Solution found on StackOverflow:
        https://stackoverflow.com/a/1578800/10874800
        """
        cls.__instance = None


class Config(metaclass=MetaConfig):
    """
    Store global configuration state

    Since this configuration is intended to be read once at program start,
    we use the Singleton pattern to ensure that only one instance of this
    object is created per program run.

    Configuration file format (YAML):

    .. code-block:: guess

        paths:
            data_path: "/data/"

    **Usage:**

    .. code-block:: python

        # at program start
        from priima.config import Config

        config = Config()  # use default config file path

        config = Config(
            config_file=os.path.join("path", "to", "priima.config.yml")

    Local configuration overrides can be placed in "priima.local_config.yml"
    """
    def __init__(self, config_file="priima.config.yml",
                 local_config_file="priima.local_config.yml"):
        self.logger = logging.getLogger('priima.config.Config')
        self.config_file = config_file
        self.local_config_file = local_config_file
        self.paths = PathConfig()
        self.gcp_separation = 1200

        self._load()

    def _load(self):
        """
        Load the config into the Config object.
        """
        if not os.path.exists(self.config_file):
            error_msg = \
                "Config file '{0}' can't be found".format(self.config_file)
            raise IOError(error_msg)

        config_yaml = _slurp(self.config_file)
        config_data = self._parse_config(config_yaml)

        if config_data["paths"] is not None:
            self.paths = PathConfig(config_data["paths"])

        if "gcp_separation" in config_data:
            if config_data["gcp_separation"] is not None:
                self.gcp_separation = config_data["gcp_separation"]

        if os.path.exists(self.local_config_file):
            local_config_yaml = _slurp(self.local_config_file)
            local_config_data = self._parse_config(local_config_yaml)
            if local_config_data["paths"] is not None:
                self.paths = PathConfig(local_config_data["paths"])
            if "gcp_separation" in local_config_data:
                if local_config_data["gcp_separation"] is not None:
                    self.gcp_separation = local_config_data["gcp_separation"]

    def _parse_config(self, config_yaml):
        """
        Parse the input yaml configuration data and return the configuration
        data structure.
        """
        try:
            config_data = yaml.safe_load(config_yaml)
        except yaml.YAMLError as error:
            error_msg = \
                f"Unable to parse config file '{self.config_file}': {error}"
            raise ValueError(error_msg) from error

        return config_data


def _slurp(filename):
    """
    Read the contents of the given file and return as a string
    """
    with open(filename) as fh:
        contents = "".join(fh.readlines())

    return contents
