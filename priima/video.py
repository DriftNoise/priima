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

import cv2
import numpy
from osgeo import gdal
from PIL import Image


def create_video(image_dir: Path):
    """Creates a video out of the forecasted satellite images"""
    jpeg_dir = Path("/tmp") / "jpegs"
    jpeg_dir.mkdir(exist_ok=True, parents=True)

    tif_files = image_dir.glob('*.tif')

    for tif in tif_files:
        jpeg_path = jpeg_dir / tif.with_suffix('.jpg').name
        convert_tif_to_jpeg(tif_path=tif, jpeg_path=jpeg_path)

    jpeg_files = sorted(jpeg_dir.glob('*.jpg'))

    frame = cv2.imread(jpeg_files[0])
    height, width, _ = frame.shape

    video_name = str(image_dir / f"{image_dir.name}.avi")
    video = cv2.VideoWriter(video_name, 0, 3, (width, height))

    for image in jpeg_files:
        video.write(cv2.imread(image))

    cv2.destroyAllWindows()
    video.release()

    return video_name


def convert_tif_to_jpeg(tif_path, jpeg_path):
    """Creates a tif file into a jpeg file of a defined size"""
    input_array = gtif_scaling(input_path=tif_path)
    image = Image.fromarray(input_array.astype('uint8'))
    square_frame_size = 1000
    img = image.resize((square_frame_size, square_frame_size))
    img.save(jpeg_path, "JPEG", quality=95)


def gtif_scaling(input_path: Path):
    """
    Employs the SAR image scaling to output a scaled numpy array.
    Source is the SAR-processor repository.
    """
    input_dataset = gdal.Open(str(input_path), gdal.GA_ReadOnly)
    input_array = input_dataset.GetRasterBand(1).ReadAsArray()
    input_dataset = None

    # scale the float data up to max uint8/16 so that the conversion to
    # uint8/16 doesn't lose too much information, after which we can take
    # the log10 of the data to convert to decibels (as the image enhancement
    # code had done) which we can then linearly stretch to fit the uint8/16
    # range.
    max_value = 2**8 - 1
    numpy_type = numpy.uint8

    if input_array.max() == 0:
        return None

    scale_factor = float(max_value)/float(input_array.max())
    input_array = scale_factor*input_array
    nonzero_mask = numpy.nonzero(input_array)
    input_array[nonzero_mask] = numpy.log10(input_array[nonzero_mask])

    min_val = numpy.percentile(input_array[nonzero_mask], 1)
    max_val = numpy.percentile(input_array[nonzero_mask], 99)

    # filter integer overflow artifacts
    input_array[numpy.where(input_array > max_val)] = max_val
    input_array[
        numpy.where((input_array < min_val) & (input_array != 0))] = min_val

    # see https://en.wikipedia.org/wiki/Normalization_(image_processing)
    new_min = 1
    new_max = max_value

    input_array[nonzero_mask] = \
        (input_array[nonzero_mask] - min_val) * \
        (new_max - new_min)/(max_val - min_val) + new_min
    input_array = numpy.rint(input_array).astype(numpy_type)

    return input_array
