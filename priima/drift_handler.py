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

import datetime
import math

import numpy as np
# import matplotlib.pyplot as plt
import scipy.interpolate
from matplotlib import path
from netCDF4 import Dataset
from pyproj import Proj, transform

from priima.projection import get_projection


class DriftHandler:
    def __init__(self, ncfile, time_range, order, gcp_list):
        self.ncfile = ncfile
        self.time_range = time_range
        self.order = order
        self.gcp_list = gcp_list

    def compute_drift(self):
        drift_dict = self._read_drift()
        drift_dict = self.subset_roi_from_meshgrid(drift_dict)
        vice, uice = self.resample_drift(drift_dict)
        vice = 1*np.array(vice)
        uice = np.array(uice)

        return [uice, vice]

    def subset_roi_from_meshgrid(self, drift_dict):
        half_region_size = self.order.region_size*1000.0/2.0*4
        cx, cy = self._project_coord(invert=False)
        cxul = cx - half_region_size
        cxlr = cx + half_region_size
        cyul = cy + half_region_size
        cylr = cy - half_region_size

        bbox = [cxul, cxlr, cylr, cyul]
        bbox = np.array(bbox)
        mypath = np.array([bbox[[0, 1, 1, 0]], bbox[[2, 2, 3, 3]]]).T
        p = path.Path(mypath)

        if 'xx' not in drift_dict:
            projection = get_projection(self.order)
            wgs84 = Proj("+init=EPSG:4326")

            if len(drift_dict['lon'].shape) == 1:
                drift_dict['lon'], drift_dict['lat'] = \
                    np.meshgrid(drift_dict['lon'], drift_dict['lat'])

            drift_dict['xx'], drift_dict['yy'] = \
                transform(wgs84, projection, drift_dict['lon'],
                          drift_dict['lat'])

            xxm, yym = drift_dict['xx'], drift_dict['yy']

        else:
            xxm, yym = np.meshgrid(drift_dict['xx'], drift_dict['yy'])

        points = np.vstack((xxm.flatten(), yym.flatten())).T
        n, m = np.shape(drift_dict['lon'])
        inside = p.contains_points(points).reshape(n, m)
        ii, jj = np.meshgrid(range(m), range(n))

        i0, i1 = min(ii[inside]), max(ii[inside])
        j0, j1 = min(jj[inside]), max(jj[inside])
        lon = drift_dict['lon'][j0:j1, i0:i1]
        lat = drift_dict['lat'][j0:j1, i0:i1]

        xx = xxm[j0:j1, i0:i1]
        yy = yym[j0:j1, i0:i1]

        uice = 1*drift_dict['uice'][:, j0:j1, i0:i1]
        vice = 1*drift_dict['vice'][:, j0:j1, i0:i1]
        time = drift_dict["time"]

        drift_dict = {'lon': lon, 'lat': lat, 'xx': xx, 'yy': yy,
                      'uice': uice, 'vice': vice, 'time': time}

        return drift_dict

    def _project_coord(self, invert=False):
        projection = get_projection(self.order)
        if invert:
            inProj = projection
            outProj = Proj(init='epsg:4326')
        else:
            inProj = Proj(init='epsg:4326')
            outProj = projection

        lat, lon = self.order.center
        xx, yy = transform(inProj, outProj, lon, lat)

        return xx, yy

    def resample_drift(self, drift_dict):
        uice = drift_dict['uice'].mean(axis=0)
        vice = drift_dict['vice'].mean(axis=0)
        uice[uice == 0] = np.nan
        vice[vice == 0] = np.nan
        uice[uice == -32767] = np.nan
        vice[vice == -32767] = np.nan

        uice = np.nanmean(uice, axis=0)
        vice = np.nanmean(vice, axis=0)

        ugcp = []
        vgcp = []
        length = drift_dict['uice'].shape[1]*drift_dict['uice'].shape[2]
        for gcp in self.gcp_list:
            if len(drift_dict["lon"].shape) == 1:
                drift_dict['lon'], drift_dict['lat'] = \
                    np.meshgrid(drift_dict['lon'], drift_dict['lat'])
            lon = drift_dict['lon'].reshape(length, 1)
            lat = drift_dict['lat'].reshape(length, 1)
            ufield = drift_dict['uice']  # .mean(axis=0)
            vfield = drift_dict['vice']  # .mean(axis=0)
            ufield = ufield.reshape(length, 1)
            vfield = vfield.reshape(length, 1)
            upoints = []
            vpoints = []
            uvalues = [ufield[x] for x in ufield.nonzero()[0]]
            vvalues = [vfield[x] for x in vfield.nonzero()[0]]
            for z in ufield.nonzero():
                upoints.append(np.hstack((lon[z], lat[z])))
            for z in vfield.nonzero():
                vpoints.append(np.hstack((lon[z], lat[z])))

            ugcptmp = scipy.interpolate.griddata(
                upoints[0], uvalues, [gcp.GCPX, gcp.GCPY])
            vgcptmp = scipy.interpolate.griddata(
                vpoints[0], vvalues, [gcp.GCPX, gcp.GCPY])

            ugcp.append(ugcptmp)
            vgcp.append(vgcptmp)

        return vgcp, ugcp

    def _read_drift(self):
        """

        """
        ind = 0
        if not isinstance(self.ncfile, list):
            ncfiles = [self.ncfile]

        if self.order.data_source == "TOPAZ":
            uice_var_name = "uice"
            vice_var_name = "vice"
            data_type = "topaz"
            xy_scale_factor = 100000
        # else: implicitly meant NEXTSIM
        else:
            uice_var_name = "vxsi"
            vice_var_name = "vysi"
            data_type = "nextsim"
            xy_scale_factor = 1

        for filename in ncfiles:
            dataset = Dataset(filename)
            dataset_time = dataset.variables['time'][:]
            time_window_index = self._select_time_window(
                dataset_time, data_type=data_type)

            if ind == 0:
                lon = dataset.variables['longitude'][:]
                lat = dataset.variables['latitude'][:]
                xx = dataset.variables['x'][:] * xy_scale_factor
                yy = dataset.variables['y'][:] * xy_scale_factor
                uice = \
                    dataset.variables[uice_var_name][time_window_index, :, :]
                vice = \
                    dataset.variables[vice_var_name][time_window_index, :, :]
                ind = 1
                time = dataset.variables['time']
                continue

            uice = np.concatenate(
                    (uice, dataset.variables[uice_var_name][
                        time_window_index, :, :]), axis=0)
            vice = np.concatenate(
                    (vice, dataset.variables[vice_var_name][
                        time_window_index, :, :]), axis=0)

        return {'lon': lon, 'lat': lat, 'xx': xx, 'yy': yy,
                'uice': uice, 'vice': vice, 'time': time}

    def _select_time_window(self, dataset_time, data_type='topaz'):
        """
        Extracts the relevant time window index between the actual image
        acquisition and the desired forecast.

        :param time: time given as 'hours since 1900-01-01 00:00:0.0'
        :type time: numpy.ndarray
        """
        start = self.time_range[0]
        end = self.time_range[1]
        if data_type == 'topaz':
            epoch = datetime.datetime(1950, 1, 1, 0, 0, 0)
            hours_coefficient = 24
        elif data_type == 'nextsim':
            epoch = datetime.datetime(1900, 1, 1, 0, 0, 0)
            hours_coefficient = 24
            dataset_time = dataset_time * hours_coefficient
            dataset_time = np.array(
                    [math.floor(item) for item in dataset_time.data])
        else:
            print('Unknown data_type string. Should be "topaz" or ecmwf.'
                  'Exiting..')
            exit(111)

        start_hours = (start - epoch).days * hours_coefficient + start.hour
        end_hours = (end - epoch).days * hours_coefficient + end.hour
        time_window_index = \
            np.where(np.logical_and(dataset_time >= start_hours,
                                    dataset_time <= end_hours))[0]

        return time_window_index
