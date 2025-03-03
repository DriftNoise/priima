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

/*M///////////////////////////////////////////////////////////////////////////////////////
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

#include "tps_trans.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/python/stl_iterator.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <iomanip>

#define PORTABLE (false)

namespace fast_cast
{

using namespace boost::python;
using namespace cv;

//https://stackoverflow.com/questions/39821367/very-fast-approximate-logarithm-natural-log-function-in-c
float __int_as_float(int32_t a) { 
    float r; 
    memcpy (&r, &a, sizeof(r)); 
    return r;
}

int32_t __float_as_int(float a) { 
    int32_t r; 
    memcpy (&r, &a, sizeof(r)); 
    return r;
}

/* compute natural logarithm, maximum error 0.85756 ulps */
float fast_logf(float a)
{
    float m, r, s, t, i, f;
    int32_t e;

    if ((a > 0.0f) && (a <= 3.40282347e+38f)) { // 0x1.fffffep+127
#if PORTABLE
        m = frexpf (a, &e);
        if (m < 0.666666667f) {
            m = m + m;
            e = e - 1;
        }
        i = (float)e;
#else // PORTABLE
        i = 0.0f;
        /* fix up denormal inputs */
        if (a < 1.175494351e-38f){ // 0x1.0p-126
            a = a * 8388608.0f; // 0x1.0p+23
            i = -23.0f;
        }
        e = (__float_as_int (a) - 0x3f2aaaab) & 0xff800000;
        m = __int_as_float (__float_as_int (a) - e);
        i = fmaf ((float)e, 1.19209290e-7f, i); // 0x1.0p-23
#endif // PORTABLE
        /* m in [2/3, 4/3] */
        f = m - 1.0f;
        s = f * f;
        /* Compute log1p(f) for f in [-1/3, 1/3] */
        r = fmaf (-0.130187988f, f, 0.140889585f); // -0x1.0aa000p-3, 0x1.208ab8p-3
        t = fmaf (-0.121489584f, f, 0.139809534f); // -0x1.f19f10p-4, 0x1.1e5476p-3
        r = fmaf (r, s, t);
        r = fmaf (r, f, -0.166845024f); // -0x1.55b2d8p-3
        r = fmaf (r, f,  0.200121149f); //  0x1.99d91ep-3
        r = fmaf (r, f, -0.249996364f); // -0x1.fffe18p-3
        r = fmaf (r, f,  0.333331943f); //  0x1.5554f8p-2
        r = fmaf (r, f, -0.500000000f); // -0x1.000000p-1
        r = fmaf (r, s, f);
        r = fmaf (i, 0.693147182f, r); //   0x1.62e430p-1 // log(2) 
    } else {
        r = a + a;  // silence NaNs if necessary
        if (a  < 0.0f) r =  0.0f / 0.0f; //  NaN
        if (a == 0.0f) r = -1.0f / 0.0f; // -Inf
    }
    return r;
}

/* natural log on [0x1.f7a5ecp-127, 0x1.fffffep127]. Maximum relative error 9.4529e-5 */
float faster_logf(float a)
{
    float m, r, s, t, i, f;
    int32_t e;

    e = (__float_as_int (a) - 0x3f2aaaab) & 0xff800000;
    m = __int_as_float (__float_as_int (a) - e);
    i = (float)e * 1.19209290e-7f; // 0x1.0p-23
    /* m in [2/3, 4/3] */
    f = m - 1.0f;
    s = f * f;
    /* Compute log1p(f) for f in [-1/3, 1/3] */
    r = fmaf (0.230836749f, f, -0.279208571f); // 0x1.d8c0f0p-3, -0x1.1de8dap-2
    t = fmaf (0.331826031f, f, -0.498910338f); // 0x1.53ca34p-2, -0x1.fee25ap-2
    r = fmaf (r, s, t);
    r = fmaf (r, s, f);
    r = fmaf (i, 0.693147182f, r); // 0x1.62e430p-1 // log(2) 
    return r;
}

inline float logarithm(float a)
{
    // if (fabs(logf(a) - fast_logf(a)) > 1e-6) {
    //     std::cout << std::setprecision(20);
    //     std::cout << a << std::endl;
    //     std::cout << logf(a) << std::endl;
    //     std::cout << faster_logf(a) << std::endl;
    //     exit(1);
    // }
    // return std::log(a);
    return logf(a);
    // return fast_logf(a);
    // return faster_logf(a);
}

float ThinPlateSplineShapeTransformerImpl::distance(Point2f p, Point2f q)
{
    Point2f diff = p - q;
    float norma = diff.x*diff.x + diff.y*diff.y;// - 2*diff.x*diff.y;
    assert(norma >= 0.0);
    // if (norma<0) norma=0; // ??????
    //else norma = std::sqrt(norma);
    return norma*logarithm(norma+FLT_EPSILON);
}

float ThinPlateSplineShapeTransformerImpl::distance(float px_, float py_, float qx_, float qy_)
{
	float norma = (px_-qx_)*(px_-qx_) + (py_-qy_)*(py_-qy_);
    assert(norma >= 0.0);
	return norma*logarithm(norma+FLT_EPSILON);
}

Point2f ThinPlateSplineShapeTransformerImpl::_applyTransformation(const Mat &shapeRef, const Point2f point, const Mat &tpsParameters)
// Point2f ThinPlateSplineShapeTransformerImpl::_applyTransformation(float *shapeRef[], const Point2f point, const Mat &tpsParameters)
{
    Point2f out;

    const float a1_x=tpsParameters.at<float>(tpsParameters.rows-3,0);
    const float ax_x=tpsParameters.at<float>(tpsParameters.rows-2,0);
    const float ay_x=tpsParameters.at<float>(tpsParameters.rows-1,0);
    const float affine_x=a1_x+ax_x*point.x+ay_x*point.y;
    
    const float a1_y=tpsParameters.at<float>(tpsParameters.rows-3,1);
    const float ax_y=tpsParameters.at<float>(tpsParameters.rows-2,1);
    const float ay_y=tpsParameters.at<float>(tpsParameters.rows-1,1);
    const float affine_y=a1_y+ax_y*point.x+ay_y*point.y;

    float nonrigid_x = 0.0;
    float nonrigid_y = 0.0;

    for (int j=0; j<shapeRef.rows; j++) {

        const float dist = distance(shapeRef.at<float>(j,0),shapeRef.at<float>(j,1), point.x, point.y);

        nonrigid_x += tpsParameters.at<float>(j,0) * dist;
        nonrigid_y += tpsParameters.at<float>(j,1) * dist;
    }

    out.x = affine_x + nonrigid_x;
    out.y = affine_y + nonrigid_y;

    return out;
}

/* public methods */

Mat ThinPlateSplineShapeTransformerImpl::warpImageWrapper(Mat transformingImage) const {

    Mat output;
    warpImage(transformingImage, output);
    return output;
}

void ThinPlateSplineShapeTransformerImpl::warpImage(InputArray transformingImage, OutputArray output,
                                      int flags, int borderMode, const Scalar& borderValue) const
{
    // CV_INSTRUMENT_REGION();

    CV_Assert(tpsComputed==true);

    Mat theinput = transformingImage.getMat();
    Mat mapX(theinput.rows, theinput.cols, CV_32FC1);
    Mat mapY(theinput.rows, theinput.cols, CV_32FC1);

    #pragma omp parallel for
    for (int row = 0; row < theinput.rows; row++)
    {
        for (int col = 0; col < theinput.cols; col++)
        {
            Point2f pt = _applyTransformation(shapeReference, Point2f(float(col), float(row)), tpsParameters);
            mapX.at<float>(row, col) = pt.x;
            mapY.at<float>(row, col) = pt.y;
        }
    }
    remap(transformingImage, output, mapX, mapY, flags, borderMode, borderValue);
}

float ThinPlateSplineShapeTransformerImpl::applyTransformation(InputArray inPts, OutputArray outPts)
{
    // CV_INSTRUMENT_REGION();

    CV_Assert(tpsComputed);
    Mat pts1 = inPts.getMat();
    CV_Assert((pts1.channels()==2) && (pts1.cols>0));

    //Apply transformation in the complete set of points
    // Ensambling output //
    if (outPts.needed())
    {
        outPts.create(1,pts1.cols, CV_32FC2);
        Mat outMat = outPts.getMat();
        for (int i=0; i<pts1.cols; i++)
        {
            Point2f pt=pts1.at<Point2f>(0,i);
            outMat.at<Point2f>(0,i)=_applyTransformation(shapeReference, pt, tpsParameters);
        }
    }

    return transformCost;
}

void ThinPlateSplineShapeTransformerImpl::estimateTransformationWrapper(cv::Mat _pts1, cv::Mat _pts2, boost::python::list _matches) {

    // convert:
    // https://stackoverflow.com/questions/3761391/boostpython-python-list-to-stdvector
    std::vector<DMatch> matches; // = std::vector<DMatch>( boost::python::stl_input_iterator<DMatch>( _matches ), boost::python::stl_input_iterator<DMatch>( ) );
    for (int i = 0; i < len(_matches); ++i)
    {
        if (len(_matches[i]) == 3) {
            matches.emplace_back(boost::python::extract<int>(_matches[i][0]), 
                                 boost::python::extract<int>(_matches[i][1]),
                                 boost::python::extract<double>(_matches[i][2]));
        }
        else if (len(_matches[i]) == 4) {
            matches.emplace_back(boost::python::extract<int>(_matches[i][0]), 
                                 boost::python::extract<int>(_matches[i][1]),
                                 boost::python::extract<int>(_matches[i][2]),
                                 boost::python::extract<int>(_matches[i][3]));
        }
        else {
            throw std::invalid_argument("A match must have length 3 or 4 but was given as " + boost::lexical_cast<std::string>(len(_matches[i])));
        }
        //matches.push_back(boost::python::extract<Match>(_matches[i]));
        //std::cout << boost::python::extract<DMatch>(_matches[i]) << std::endl;
        // matches.emplace_back();
    } 
    estimateTransformation(_pts1, _pts2, matches);
}

void ThinPlateSplineShapeTransformerImpl::estimateTransformation(InputArray _pts1, InputArray _pts2, std::vector<DMatch>& _matches )
{
    // CV_INSTRUMENT_REGION();
    Mat pts1 = _pts1.getMat();
    Mat pts2 = _pts2.getMat();
    CV_Assert((pts1.channels()==2) && (pts1.cols>0) && (pts2.channels()==2) && (pts2.cols>0));
    CV_Assert(_matches.size()>1);

    if (pts1.type() != CV_32F)
        pts1.convertTo(pts1, CV_32F);
    if (pts2.type() != CV_32F)
        pts2.convertTo(pts2, CV_32F);

    // Use only valid matchings //
    std::vector<DMatch> matches;
    for (size_t i=0; i<_matches.size(); i++)
    {
        if (_matches[i].queryIdx<pts1.cols &&
            _matches[i].trainIdx<pts2.cols)
        {
            matches.push_back(_matches[i]);
        }
    }

    // Organizing the correspondent points in matrix style //
    Mat shape1((int)matches.size(),2,CV_32F); // transforming shape
    Mat shape2((int)matches.size(),2,CV_32F); // target shape
    for (int i=0, end = (int)matches.size(); i<end; i++)
    {
        Point2f pt1=pts1.at<Point2f>(0,matches[i].queryIdx);
        shape1.at<float>(i,0) = pt1.x;
        shape1.at<float>(i,1) = pt1.y;

        Point2f pt2=pts2.at<Point2f>(0,matches[i].trainIdx);
        shape2.at<float>(i,0) = pt2.x;
        shape2.at<float>(i,1) = pt2.y;
    }
    shape1.copyTo(shapeReference);

    // Building the matrices for solving the L*(w|a)=(v|0) problem with L={[K|P];[P'|0]}

    //Building K and P (Needed to build L)
    Mat matK((int)matches.size(),(int)matches.size(),CV_32F);
    Mat matP((int)matches.size(),3,CV_32F);
    for (int i=0, end=(int)matches.size(); i<end; i++)
    {
        for (int j=0; j<end; j++)
        {
            if (i==j)
            {
                matK.at<float>(i,j)=float(regularizationParameter);
            }
            else
            {
                matK.at<float>(i,j) = distance(Point2f(shape1.at<float>(i,0),shape1.at<float>(i,1)),
                                               Point2f(shape1.at<float>(j,0),shape1.at<float>(j,1)));
            }
        }
        matP.at<float>(i,0) = 1;
        matP.at<float>(i,1) = shape1.at<float>(i,0);
        matP.at<float>(i,2) = shape1.at<float>(i,1);
    }

    //Building L
    Mat matL=Mat::zeros((int)matches.size()+3,(int)matches.size()+3,CV_32F);
    Mat matLroi(matL, Rect(0,0,(int)matches.size(),(int)matches.size())); //roi for K
    matK.copyTo(matLroi);
    matLroi = Mat(matL,Rect((int)matches.size(),0,3,(int)matches.size())); //roi for P
    matP.copyTo(matLroi);
    Mat matPt;
    transpose(matP,matPt);
    matLroi = Mat(matL,Rect(0,(int)matches.size(),(int)matches.size(),3)); //roi for P'
    matPt.copyTo(matLroi);

    //Building B (v|0)
    Mat matB = Mat::zeros((int)matches.size()+3,2,CV_32F);
    for (int i=0, end = (int)matches.size(); i<end; i++)
    {
        matB.at<float>(i,0) = shape2.at<float>(i,0); //x's
        matB.at<float>(i,1) = shape2.at<float>(i,1); //y's
    }

    //Obtaining transformation params (w|a)
    solve(matL, matB, tpsParameters, DECOMP_LU);
    //tpsParameters = matL.inv()*matB;

    //Setting transform Cost and Shape reference
    Mat w(tpsParameters, Rect(0,0,2,tpsParameters.rows-3));
    Mat Q=w.t()*matK*w;
    transformCost=fabs(Q.at<float>(0,0)*Q.at<float>(1,1));//fabs(mean(Q.diag(0))[0]);//std::max(Q.at<float>(0,0),Q.at<float>(1,1));
    tpsComputed=true;
}

#if (PY_VERSION_HEX >= 0x03000000)

    static void *init_ar() {
#else
    static void init_ar(){
#endif
    Py_Initialize();

    import_array();
    return NUMPY_IMPORT_ARRAY_RETVAL;
    }

BOOST_PYTHON_MODULE (FastCast) {

	//using namespace XM;
	init_ar();

	//initialize converters
	to_python_converter<cv::Mat,pbcvt::matToNDArrayBoostConverter>();
	pbcvt::matFromNDArrayBoostConverter();

	//expose module-level functions
    def("createThinPlateSplineShapeTransformer", createThinPlateSplineShapeTransformer);
    class_<ThinPlateSplineShapeTransformerImpl>("ThinPlateSplineShapeTransformerImpl")
        .def("estimateTransformation", &ThinPlateSplineShapeTransformerImpl::estimateTransformationWrapper)
        // .def("applyTransformation", &ThinPlateSplineShapeTransformerImpl::applyTransformation)
        .def("warpImage", &ThinPlateSplineShapeTransformerImpl::warpImageWrapper)
    ;

}

} // fast_cast
