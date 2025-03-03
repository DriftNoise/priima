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

from unittest import TestCase

import cv2
import numpy as np
import pytest

from priima.fast_cast import fc


class TestFastCast(TestCase):

    @pytest.mark.skipif(fc is None, reason="Fast Cast Library not found.")
    def test_thin_plate_spline_shape_transformer(self):

        image = np.random.rand(300, 300)

        pts1 = np.float32([[[100., 100.], [100., 160.],
                            [100., 220.], [100., 280.]]])
        pts2 = np.float32([[[103., 101.], [103., 161.],
                            [103., 221.], [103., 282.]]])

        match1 = [1, 2, 0.0]
        match2 = [5, 16, 00.0]

        # run OpenCV
        matches = [
            cv2.DMatch(match1[0], match1[1], match1[2]),
            cv2.DMatch(match2[0], match2[1], match2[2])
        ]
        tps = cv2.createThinPlateSplineShapeTransformer()
        tps.estimateTransformation(pts2, pts1, matches)
        out_img = tps.warpImage(image)

        # run FastCast
        matches2 = [match1, match2]
        tps2 = fc.ThinPlateSplineShapeTransformerImpl()
        tps2.estimateTransformation(pts2, pts1, matches2)
        out_img2 = tps2.warpImage(image)

        # assert that result is the same
        self.assertTrue(np.all(out_img == out_img2))
        # make sure that we do not compare zero matrices here
        self.assertTrue(np.any(out_img2 != 0))
        # make sure that something actually changed
        self.assertTrue(np.any(out_img2 != image))
