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

from pathlib import Path
from typing import Optional

from configargparse import ArgParser


class ConfigInstance:
    forecast_duration: int
    store_training_data: bool
    image: str


class Config:
    _instance: Optional[ConfigInstance] = None

    @classmethod
    def instance(cls) -> ConfigInstance:
        if cls._instance is None:
            parser = ArgParser(
                description=(
                    "Generate the ice image of tomorrow"
                )
            )
            parser.add_argument(
                '--image',
                help='Path of the image file to process')
            parser.add_argument(
                '--forecast-duration', type=int,
                help='Number of hours in the future to run forecast')
            parser.add_argument(
                '--data-source', choices=['ICON', 'NEXTSIM'],
                help='Model from which to derive the ice drift'
            )
            parser.add_argument(
                '--gcp-separation', default=5000, type=int,
                help='Space between ground control points in pixels'
            )
            parser.add_argument(
                '--output-resolution', default=100, type=int,
                help="Resolution of the output images"
            )
            parser.add_argument(
                '--video', action="store_true",
                help="Create a video of the warped images"
            )
            parser.add_argument(
                '--gtif-out', default=False, action='store_true',
                help="Disables the output of a cloud-optimized geotiff"
            )
            parser.add_argument(
                '--plot-wind', action='store_true', default=False,
                help="Whether or not to create plots of wind data"
            )
            args = parser.parse_known_args()[0]
            args.data_path = Path("/model_data/")
            args.out_path = Path("/out_dir")
            args.grid_path = Path("/home/ossi/priima/icon_grid/")
            cls._instance = args

        return cls._instance

    @classmethod
    def set_attribute(cls, key: str, value: any) -> None:
        """Sets attribute in the config instance"""
        if cls._instance is None:
            cls.instance()
        setattr(cls._instance, key, value)
