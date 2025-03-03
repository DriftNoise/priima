// This file contains both BSD-3-Clause and GPL-3.0 licensed code.
// The entire modified work is licensed under GPL-3.0 due to GPL compatibility rules

// Copyright 2021, Kai Wah Chan, Mitja Echim, Andreas Folkers

// This file is part of PRIIMA.
// PRIIMA is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
// PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with
// PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.
// 
// /*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

// #define __OPENCV_BUILD

// #include "precomp.hpp"
#include <vector>
#include <cmath>
#include <iostream>

#include "opencv2/calib3d.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/shape.hpp"

#include "opencv2/core/utility.hpp"
// #include "opencv2/core/private.hpp"

#include "opencv2/opencv_modules.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL pbcvt_ARRAY_API
#include <boost/python.hpp>
#include <pyboostcvconverter/pyboostcvconverter.hpp>

namespace fast_cast
{

class ThinPlateSplineShapeTransformerImpl /*CV_FINAL*/ : public cv::ThinPlateSplineShapeTransformer
{
public:

    /* Constructors */
    ThinPlateSplineShapeTransformerImpl()
    {
        regularizationParameter=0;
        name_ = "ShapeTransformer.TPS";
        tpsComputed=false;
        transformCost = 0;
    }

    ThinPlateSplineShapeTransformerImpl(double _regularizationParameter)
    {
        regularizationParameter=_regularizationParameter;
        name_ = "ShapeTransformer.TPS";
        tpsComputed=false;
        transformCost = 0;
    }

    /* Destructor */
    ~ThinPlateSplineShapeTransformerImpl() override /*CV_OVERRIDE*/
    {
    }

    //! the main operators
    virtual void estimateTransformation(cv::InputArray _pts1, cv::InputArray _pts2, std::vector<cv::DMatch> &_matches) override /*CV_OVERRIDE*/;
    virtual float applyTransformation(cv::InputArray inPts, cv::OutputArray output=cv::noArray()) override /*CV_OVERRIDE*/;
    virtual void warpImage(cv::InputArray transformingImage, cv::OutputArray output,
                           int flags=cv::INTER_LINEAR, int borderMode=cv::BORDER_CONSTANT, const cv::Scalar& borderValue=cv::Scalar()) const override /*CV_OVERRIDE*/;
              
    void estimateTransformationWrapper(cv::Mat _pts1, cv::Mat _pts2, boost::python::list _matches);                
    cv::Mat warpImageWrapper(cv::Mat transformingImage) const;


    //! Setters/Getters
    virtual void setRegularizationParameter(double _regularizationParameter) override /*CV_OVERRIDE*/ { regularizationParameter=_regularizationParameter; }
    virtual double getRegularizationParameter() const override /*CV_OVERRIDE*/ { return regularizationParameter; }

    static float distance(cv::Point2f p, cv::Point2f q);
    static float distance(float px_, float py_, float qx_, float qy_);
    static cv::Point2f _applyTransformation(const cv::Mat &shapeRef, const cv::Point2f point, const cv::Mat &tpsParameters);
    // static cv::Point2f _applyTransformation(float *shapeRef[], const cv::Point2f point, const cv::Mat &tpsParameters);

private:
    bool tpsComputed;
    double regularizationParameter;
    float transformCost;
    cv::Mat tpsParameters;
    cv::Mat shapeReference;

protected:
    cv::String name_;
};

cv::Ptr <cv::ThinPlateSplineShapeTransformer> createThinPlateSplineShapeTransformer(double regularizationParameter)
{
    return cv::Ptr<cv::ThinPlateSplineShapeTransformer>( new ThinPlateSplineShapeTransformerImpl(regularizationParameter) );
}

} // fast_cast