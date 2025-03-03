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

import numpy as np
from netCDF4 import Dataset
from pyproj import Proj, transform


def make_antarctic_grid(filepath):
    xcells = 780
    ycells = 620

    # units are in 100 km
    xinc = 0.125
    yinc = 0.125
    x0 = -48
    y0 = -38

    xx = np.arange(x0, x0 + xinc*xcells, xinc)
    yy = np.arange(y0, y0 + yinc*ycells, yinc)
    tt = np.arange(0, 24, 1)
    dummy_data = np.random.rand(len(tt), ycells, xcells)

    outProj = Proj(init='epsg:3031')
    inProj = Proj(init='epsg:4326')

    xx_mesh, yy_mesh = np.meshgrid(xx, yy)
    lon, lat = transform(outProj, inProj, xx_mesh*100000, yy_mesh*100000)

    nc = Dataset(filepath, 'w', clobber=False)
    nc.createDimension('x', xcells)
    nc.createDimension('y', ycells)
    nc.createDimension('time', None)

    xo = nc.createVariable('x', 'f', ('x'))
    xo.units = 'm'
    xo.standard_name = 'projection_x_coordinate'
    xo.axis = 'X'

    yo = nc.createVariable('y', 'f', ('y'))
    yo.units = 'm'
    yo.standard_name = 'projection_y_coordinate'
    yo.axis = 'Y'

    lono = nc.createVariable('longitude', 'f', ('y', 'x'))
    lono.units = 'degrees_east'
    lono.standard_name = 'longitude'
    lono.long_name = 'longitude'
    lono._CoordinateAxisType = "Lon"

    lato = nc.createVariable('latitude', 'f', ('y', 'x'))
    lato.units = 'degrees_north'
    lato.standard_name = 'latitude'
    lato.long_name = 'latitude'
    lato._CoordinateAxisType = "Lat"

    to = nc.createVariable('time', 'double', ('time'))
    to.units = 'hours since 1950-01-01 00:00:00'
    to.standard_name = 'time'

    dummyo = nc.createVariable('dummy', 'f4', ('time', 'y', 'x'))
    dummyo.units = 'm'
    dummyo.coordinates = 'longitude latitude'
    dummyo.missing_value = -32767.
    dummyo.standard_name = 'dummy_data'
    dummyo.grid_mapping = 'stereographic'
    dummyo.cell_methods = "area:mean"

    crso = nc.createVariable('stereographic', 'i4')
    crso.grid_mapping_name = "polar_stereographic"
    crso.latitude_of_projection_origin = -90.
    crso.longitude_of_projection_origin = 0.
    crso.scale_factor_at_projection_origin = 1.
    crso.straight_vertical_longitude_from_pole = 0.
    crso.false_easting = 0.
    crso.false_northing = 0.

    nc.Conventions = 'CF-1.4'

    xo[:] = xx
    yo[:] = yy
    lono[:] = lon
    lato[:] = lat
    to[:] = tt
    dummyo[:] = dummy_data

    nc.close()


if __name__ == "__main__":
    make_antarctic_grid()
