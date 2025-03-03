<!-- 
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
-->

# Algorithm explanation and basic principles

## GEOLOCATION and GCP

PRIIMA geolocates and resamples forecast and image data. The image and the
TOPAZ4 data are in NSIDC Polar Stereographic projection. ICON wind fields
have an unstructured grid in GRIB2 format.  The algorithm is capable to
perform a temporal alignment, and selects only the relevant times starting
from the SAR image timestamp included in the filename, and expanding to the
`forecast_duration` which we set in a given order.

PRIIMA sets an initial grid of Ground Control Points (GCP) based on the SAR
image which will use as input.  For each GCP computes the horizontal u, v
drift components.

Sometimes it is necessary to exactly control the number of GCPs set, e.g. to
reproduce the training data used for the neural network produced in the
Fast-Cast project. For this, refer to
[Hacks to control the number of GCPs](./docs/priima-hacks-to-set-GCPs.md).

**relevant scripts/functions**
1. inverse_gcps_recipe.py
    1. inversion_gcps_recipe(): It uses the image dimensions and projection information to workout and
       return a grid of GCP. Usually we use that in "mpl" mode which means multiple point locations.
       The method "optimized" can also be used. In this case a dynamic grid is being created, with more
       points around the shorelines, so that PRIIMA can better constrain the ice motion there. This
       is totally experimental, and needs further testing and development.


## INPUT DATA

The relevant region information (central coordinates and extension) are
extracted from a given order request, and the satellite data are then
resampled to match the desired ROI. To resample the drift information a
rectangle box is being formed surrounding the ROI and reprojected to NSIDC
Polar Stereographic projection.  TOPAZ4 data are mapped onto a mesh grid
which is then compared to the ROI box. The outcome of this comparison is a
binary mask which is true for all the points located within the ROI box. The
mask is then applied to the TOPAZ4 data to extract the desired ROI

ICON data are converted from GRIB2 to a NetCDF format and from the
unstructured icosehdral grid to a global regular latitude/longitude grid
using the Climate Data Operators (CDO). Then the NetCDF files are
reprojected to match the TOPAZ4 projection. The next processing steps are
the same as in the TOPAZ4  case.

**relevant scripts/functions**
1. drift_handler.py:
    1. compute_drift(): Higher level function to compute the drift. Calls several other functions.
    2. subset_roi_from_meshgrid(): Extracts only the relevant region from the TOPAZ4/NEXTSIM data.
    3. resample_drift(): Samples the drift for the specified GCP. It performs a two dimensional interpolation
       to match the drift to the GCP locations.

2. icon.py:
    1. process_icon_fields(): Reads the target grid and weights, selects the relevant files, reprojects them
       to a regular lat/lon grid, and saves them as NetCDF file under the converted/ folder. Then combines
       the u, v files into 1, and reprojects again to polar projection.

    2. regrid_regular2polar(): Uses CDO and an example NetCDF file (TOPAZ4 data) to remap to polar projection.
       Also renames variables as necessary.

    3. rotate_velocity(): This is tricky and also not totally clear to me. It seems that the wind field is
       earth relative and we need to rotate the vectors to be polar grid relative. The problem is that when
       using the wind fields PRIIMA was morphing the images to a wrong direction. This function seems to solve
       that problem and it is configured to work for the Southern Hemisphere. For the Northern Hemisphere it
       may be that we need to change the sign of the rotcon_p constant. **DO VALIDATE** results when changing
       Hemisphere! Check what PRIIMA is producing, a) trajectories b) images. See how ice is moving for every image.
       Compare with subsequent SAR image, and with online wind models where you can see how the wind is blowing e.g.
       use this site: https://earth.nullschool.net/#2020/02/02/0600Z/wind/surface/level/stereographic=-47.75,82.37,3000/loc=17.557,82.319 .
       More information on comparisons here:
       https://www.dropbox.com/home/Projects/FASTCAST/Panos_Notes?preview=project_notes.odf

    4. couple_tide_and_wind(): Activate with use_tides: True key in a given order. It will linearly add the tidal
       component to the wind. We use an INFLATION = 40 factor to compensate for the physical difference between
       wind and ice drift. In a later step, wind needs to be scaled down to 2.5% of its magnitude and will
       have a 25 degrees rotation as well to create the "ice drift equivalent". In order to add the tidal
       component we need to scale that, so that we create a "wind equivalent" from the tide. 2.5% --> 100/2.5 = 40.

3. wind.py
    1. wind2drift(): Converts the wind into ice drift equivalent by scaling the magnitude down to 2.5%, and rotating
       the vectors 25 degrees right to the wind direction for Northern, and 25 degrees left for the Southern Hemisphere.

## DRIFT VECTOR CALCULATION

Drift vectors are computed for each of the initial grid points produced from
the `inversion_gcps_recipe()` function, by linearly interpolating the point
positions to the drift field. The drift information is translated into
displacement and the new point locations are calculated. This step is
repeated for every hourly drift file until the algorithm reaches the end of
the forecast window.

The displacement is calculated by propagating the time component (hourly
steps hence t=1h and converting to International System gives 60min * 60sec)
to the drift vector for each direction (x, y) using the following equation
system:

x_displacement = uice * t * 60 * 60                                                                                                 (1)
y_displacement = vice * t * 60 * 60                                                                                                 (2)

## WARPING PROCESS

We use the openCV library to perform a Thin Plate Warping on the image. With
this method we take a sparse set of control points with their corresponding
displacements, and a mapping function f is computed from the algorithm,
which maps the control points of the initial image to their new locations.
The warping process needs the displacement information to be in the pixel
space, whilst the displacement of the points using the drift field, is
of course in the geographical space (expressed in meters, and locations in
geographical coordinates). PRIIMA is performing the translation between
pixel and geographical space.

Given an image, a set of initial control point locations, and their
corresponding new locations (in pixel space), the warping process will
create a new image which its control points will be placed in their new
locations. There is a tricky part here, when for example the average move of
the image is on one direction. Imagine that the image is within a frame. If
you move the image right, then you will lose those pixels that are moved
outside of the frame. To overcome this, PRIIMA is creating a buffer zone
around the image filled with 0s. That way the image is allowed to move
within the buffer, without losing any information.

**relevant scripts/functions**
1. warp_forecast.py:
    1. update_point_location(): It updates the location of the control points by adding the corresponding displacement. It makes
       also a test if a point is in land or not. If yes it will turn the gcp.Info attribute to True, and the point will not
       move any more for the rest time steps.

    2. compute_pixel_coordinates(): It converts geodetic displacement information into pixel displacement. For more information
       look at: https://www.dropbox.com/home/Projects/ESA_PRIIMA/Reports/deliverables/final?preview=SSA_priima_final_ver_0_1.docx

    3. create_transformation_matrix(): Creates the matrix with the starting and ending positions for all control points, in both
       pixel and geodetic coordinates. This is the input for the warping function.

2. preprocessing.py:
    1. warp_image(): Warps the image by using the transformation matrix. It makes use of a pixel buffer, where it places the
       initial image inside that buffer. That ensures that if the image is forced to move on one side due to the drift,
       we will still recover the whole image, and it will not be trimmed. But this trick, needs to be considered later, when
       saving the geotiff.

    2. write_geotiff(): Saves the image as GeoTIFF. The buffer inserted in the previous step, needs to be considered, when
       constructing the geotransformation vector to ensure a proper georeference. Hence i.e. the upper left corner x coordinate
       is defined as the initial upper left x coordinate (from the initial SAR image) minus the outcome of adding the buffer pixels
       with the relative pixel displacement due to the drift and then multiplying it with the grid resolution to get it in meters.



