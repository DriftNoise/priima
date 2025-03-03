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

from osgeo import gdal

from priima.inverse_gcps_recipe import inversion_gcps_recipe


def export_gcp_2csv(gcp_list):
    with open('gcp_reprojected.csv', "w") as ofile:
        writer = csv.writer(ofile, delimiter=',', quoting=csv.QUOTE_ALL)
        writer.writerow(["X", "Y", "XX", "YY"])
        for gcpi in gcp_list:
            line = [gcpi.GCPX, gcpi.GCPY, gcpi.GCPXX, gcpi.GCPYY]
            writer.writerow(line)


def create_gcp_from_ul_lr(dst, proj):
    gcp_list = []
    ii = 1
    ilon = -1
    xx, yy, x_ind, y_ind, lons, lats, \
        polygon_x, polygon_y, polygon_x_ind, polygon_y_ind, \
        polygon_lons, polygon_lats = inversion_gcps_recipe(dst, proj)

    for ix in x_ind:
        ilon = ilon + 1
        ilat = 0
        for iy in y_ind:
            gcpi = gdal.GCP()

            gcpi.GCPLine = int(ix)
            gcpi.GCPPixel = int(iy)

            gcpi.GCPX = lons[ilat, ilon]
            gcpi.GCPY = lats[ilat, ilon]
            gcpi.GCPXX = xx[ilat, ilon]
            gcpi.GCPYY = yy[ilat, ilon]

            gcpi.Id = str(ii)
            ii = ii + 1
            ilat = ilat + 1
            gcp_list.append(gcpi)
            gcpi = None
    if False:
        for ind, ix in enumerate(polygon_x_ind):
            gcpi = gdal.GCP()

            gcpi.GCPLine = ix
            gcpi.GCPPixel = polygon_y_ind[ind]

            gcpi.GCPX = polygon_lons[ind]
            gcpi.GCPY = polygon_lats[ind]
            gcpi.GCPXX = polygon_x[ind]
            gcpi.GCPYY = polygon_y[ind]
            gcpi.Info = 'is_land'

            gcpi.Id = str(ii)
            ii = ii + 1
            gcp_list.append(gcpi)
            gcpi = None

    return gcp_list
