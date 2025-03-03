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
from pathlib import Path
from typing import Literal

import geopandas as gpd


def filter_shapefile_by_latitude(
        *,
        input_shapefile: Path,
        output_shapefile: Path,
        latitude_threshold: int,
        polar_region: Literal['arctic', 'antarctic'],
):
    gdf = gpd.read_file(input_shapefile)
    if not gdf.crs.is_geographic:
        err_msg = "input shapefile %s is not of geographic CRS"
        raise ValueError(err_msg % input_shapefile)
    
    if polar_region == 'arctic':
        gdf_filtered = gdf[gdf.geometry.centroid.y > latitude_threshold]
    else:
        gdf_filtered = gdf[gdf.geometry.centroid.y < latitude_threshold]

    gdf_filtered.to_file(output_shapefile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-path", type=Path, help="Path to read from and write to")
    args = parser.parse_args()

    input_shapefile = (
        args.data_path / "land-polygons-split-4326" / "land_polygons.shp"
    )
    # create landmask for NH
    output_arctic = args.data_path / "arctic_landmask.shp"  
    filter_shapefile_by_latitude(
        input_shapefile=input_shapefile,
        output_shapefile=output_arctic,
        latitude_threshold=50,
        polar_region='arctic',
    )
    # create landmask for SH
    output_antarctic = args.data_path / "antarctic_landmask.shp"
    filter_shapefile_by_latitude(
        input_shapefile=input_shapefile,
        output_shapefile=output_antarctic,
        latitude_threshold=-60,
        polar_region='antarctic',
    )
