#include "CacheWeights.h"
#include "WeightBase.h"
#include "WeightSpline.h"

#include <algorithm>
#include <iostream>
#include <exception>
#include <limits>
#include <cmath>

#include <hemi/hemi_error.h>
#include <hemi/launch.h>
#include <hemi/grid_stride_range.h>

// The constructor
Cache::Weight::Spline::Spline(
    Cache::Weights::Results& weights,
    Cache::Parameters::Values& parameters,
    Cache::Parameters::Clamps& lowerClamps,
    Cache::Parameters::Clamps& upperClamps,
    std::size_t splines, std::size_t knots)
    : Cache::Weight::Base("spline",weights,parameters),
      fLowerClamp(lowerClamps), fUpperClamp(upperClamps),
      fSplinesReserved(splines), fSplinesUsed(0),
      fSplineKnotsReserved(knots), fSplineKnotsUsed(0) {

    std::cout << "Reserved " << GetName() << " Splines: "
            << GetSplinesReserved() << std::endl;
    if (GetSplinesReserved() < 1) return;

    fTotalBytes += GetSplinesReserved()*sizeof(int);      // fSplineResult
    fTotalBytes += GetSplinesReserved()*sizeof(short);    // fSplineParameter
    fTotalBytes += (1+GetSplinesReserved())*sizeof(int);  // fSplineIndex

    fSplineKnotsReserved = 2*fSplinesReserved + 3*fSplineKnotsReserved;

    std::cout << "Reserved " << GetName()
            << " Spline Knots: " << GetSplineKnotsReserved()
            << std::endl;
    fTotalBytes += GetSplineKnotsReserved()*sizeof(WEIGHT_BUFFER_FLOAT);  // fSpineKnots

    std::cout << "Approximate Memory Size for " << GetName()
            << ": " << fTotalBytes/1E+9
            << " GB" << std::endl;

    try {
        // Get the CPU/GPU memory for the spline index tables.  These are
        // copied once during initialization so do not pin the CPU memory into
        // the page set.
        fSplineResult.reset(new hemi::Array<int>(GetSplinesReserved(),false));
        fSplineParameter.reset(
            new hemi::Array<short>(GetSplinesReserved(),false));
        fSplineIndex.reset(new hemi::Array<int>(1+GetSplinesReserved(),false));

        // Get the CPU/GPU memory for the spline knots.  This is copied once
        // during initialization so do not pin the CPU memory into the page
        // set.
        fSplineKnots.reset(
            new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetSplineKnotsReserved(),false));
    }
    catch (std::bad_alloc&) {
        std::cerr << "Failed to allocate memory, so stopping" << std::endl;
        throw std::runtime_error("Not enough memory available");
    }

    // Initialize the caches.  Don't try to zero everything since the
    // caches can be huge.
    fSplineIndex->hostPtr()[0] = 0;
}

// The destructor
Cache::Weight::Spline::~Spline() {}

int Cache::Weight::Spline::FindPoints(const TSpline3* s) {
    return s->GetNp();
}

void Cache::Weight::Spline::AddSpline(int resultIndex,
                                             int parIndex,
                                             TGraph* graph) {
    int NP = graph->GetN();
    double xMin = graph->GetX()[0];
    double xMax = graph->GetX()[NP-1];
    int sIndex = ReserveSpline(resultIndex,parIndex,xMin,xMax,NP);
    for (int i=0; i<NP; ++i) {
        double x;
        double y;
        graph->GetPoint(i,x,y);
        double m = 0; // for now only do linear interpolation, so slope is not used
        if (!std::isfinite(m)) std::runtime_error("invalid spline slope");
        SetSplineKnot(sIndex,i,x,y,m);
    }
}

// Reserve space in the internal structures for spline with uniform knots.
// The knot values will still need to be set using set spline knot.
int Cache::Weight::Spline::ReserveSpline(
    int resIndex, int parIndex, double low, double high, int points) {
    if (resIndex < 0) {
        std::cerr << "Invalid result index"
               << std::endl;
        throw std::runtime_error("Negative result index");
    }
    if (fWeights.size() <= resIndex) {
        std::cerr << "Invalid result index"
               << std::endl;
        throw std::runtime_error("Result index out of bounds");
    }
    if (parIndex < 0) {
        std::cerr << "Invalid parameter index"
               << std::endl;
        throw std::runtime_error("Negative parameter index");
    }
    if (fParameters.size() <= parIndex) {
        std::cerr << "Invalid parameter index"
               << std::endl;
        throw std::runtime_error("Parameter index out of bounds");
    }
    if (high <= low) {
        std::cerr << "Invalid spline bounds"
               << std::endl;
        throw std::runtime_error("Invalid spline bounds");
    }
    // if (points < 3) {
    //     std::cerr << "Insufficient points in spline"
    //            << std::endl;
    //     throw std::runtime_error("Invalid number of spline points");
    // }
    int newIndex = fSplinesUsed++;
    if (fSplinesUsed > fSplinesReserved) {
        std::cerr << "Not enough space reserved for splines"
                  << std::endl;
        throw std::runtime_error("Not enough space reserved for splines");
    }
    fSplineResult->hostPtr()[newIndex] = resIndex;
    fSplineParameter->hostPtr()[newIndex] = parIndex;
    if (fSplineIndex->hostPtr()[newIndex] != fSplineKnotsUsed) {
        std::cerr << "Last spline knot index should be at old end of splines"
                  << std::endl;
        throw std::runtime_error("Problem with control indices");
    }
    int knotIndex = fSplineKnotsUsed;
    fSplineKnotsUsed += 2; // Space for the upper and lower bound
    fSplineKnotsUsed += 3*points; // Space for the knots.
    if (fSplineKnotsUsed > fSplineKnotsReserved) {
        std::cerr << "Not enough space reserved for spline knots"
               << std::endl;
        throw std::runtime_error("Not enough space reserved for spline knots");
    }
    fSplineIndex->hostPtr()[newIndex+1] = fSplineKnotsUsed;
    // Save values needed to calculate the spline offset index.  If the input
    // value is x, the index is
    // v = (x-CD[dataIndex])*CD[dataIndex+1].
    // i = v;
    // v = v - i;
    double invStep =  1.0*(points-1.0)/(high-low);
    fSplineKnots->hostPtr()[knotIndex] = low;
    fSplineKnots->hostPtr()[knotIndex+1] = invStep;

    return newIndex;
}

void Cache::Weight::Spline::SetSplineKnot(
    int sIndex, int kIndex, double place, double value, double slope) {
    SetSplineKnotValue(sIndex,kIndex,value);
    SetSplineKnotSlope(sIndex,kIndex,slope);
    SetSplineKnotPlace(sIndex,kIndex,place);
}

void Cache::Weight::Spline::SetSplineKnotValue(
    int sIndex, int kIndex, double value) {
    if (sIndex < 0) {
        std::cerr << "Requested spline index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (GetSplinesUsed() <= sIndex) {
        std::cerr << "Requested spline index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (kIndex < 0) {
        std::cerr << "Requested control point index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    int knotIndex = fSplineIndex->hostPtr()[sIndex] + 2 + 3*kIndex;
    if (fSplineIndex->hostPtr()[sIndex+1] <= knotIndex) {
        std::cerr << "Requested control point index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    fSplineKnots->hostPtr()[knotIndex] = value;
}

void Cache::Weight::Spline::SetSplineKnotSlope(
    int sIndex, int kIndex, double slope) {
    if (sIndex < 0) {
        std::cerr << "Requested spline index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (GetSplinesUsed() <= sIndex) {
        std::cerr << "Requested spline index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (kIndex < 0) {
        std::cerr << "Requested control point index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    int knotIndex = fSplineIndex->hostPtr()[sIndex] + 2 + 3*kIndex;
    if (fSplineIndex->hostPtr()[sIndex+1]+1 <= knotIndex) {
        std::cerr << "Requested control point index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    fSplineKnots->hostPtr()[knotIndex+1] = slope;
}

void Cache::Weight::Spline::SetSplineKnotPlace(
    int sIndex, int kIndex, double place) {
    if (sIndex < 0) {
        std::cerr << "Requested spline index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (GetSplinesUsed() <= sIndex) {
        std::cerr << "Requested spline index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    if (kIndex < 0) {
        std::cerr << "Requested control point index is negative"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    int knotIndex = fSplineIndex->hostPtr()[sIndex] + 2 + 3*kIndex;
    if (fSplineIndex->hostPtr()[sIndex+1] <= knotIndex+2) {
        std::cerr << "Requested control point index is to large"
                  << std::endl;
        std::runtime_error("Invalid control point being set");
    }
    fSplineKnots->hostPtr()[knotIndex+2] = place;
}

int Cache::Weight::Spline::GetSplineParameterIndex(int sIndex) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    return fSplineParameter->hostPtr()[sIndex];
}

double Cache::Weight::Spline::GetSplineParameter(int sIndex) {
    int i = GetSplineParameterIndex(sIndex);
    if (i<0) {
        throw std::runtime_error("Spine parameter index out of bounds");
    }
    if (fParameters.size() <= i) {
        throw std::runtime_error("Spine parameter index out of bounds");
    }
    return fParameters.hostPtr()[i];
}

int Cache::Weight::Spline::GetSplineKnotCount(int sIndex) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int k = fSplineIndex->hostPtr()[sIndex+1]-fSplineIndex->hostPtr()[sIndex]-2;
    return k/3; // should be divided by 3 instead of 2?
}

double Cache::Weight::Spline::GetSplineLowerBound(int sIndex) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int knotsIndex = fSplineIndex->hostPtr()[sIndex];
    return fSplineKnots->hostPtr()[knotsIndex];
}

double Cache::Weight::Spline::GetSplineUpperBound(int sIndex) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int knotCount = GetSplineKnotCount(sIndex);
    double lower = GetSplineLowerBound(sIndex);
    int knotsIndex = fSplineIndex->hostPtr()[sIndex];
    double step = fSplineKnots->hostPtr()[knotsIndex+1];
    return lower + (knotCount-1)/step;
}

double Cache::Weight::Spline::GetSplineLowerClamp(int sIndex) {
    int i = GetSplineParameterIndex(sIndex);
    if (i<0) {
        throw std::runtime_error("Spine lower clamp index out of bounds");
    }
    if (fLowerClamp.size() <= i) {
        throw std::runtime_error("Spine lower clamp index out of bounds");
    }
    return fLowerClamp.hostPtr()[i];
}

double Cache::Weight::Spline::GetSplineUpperClamp(int sIndex) {
    int i = GetSplineParameterIndex(sIndex);
    if (i<0) {
        throw std::runtime_error("Spine upper clamp index out of bounds");
    }
    if (fUpperClamp.size() <= i) {
        throw std::runtime_error("Spine upper clamp index out of bounds");
    }
    return fUpperClamp.hostPtr()[i];
}

double Cache::Weight::Spline::GetSplineKnotValue(int sIndex, int knot) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int knotsIndex = fSplineIndex->hostPtr()[sIndex];
    int count = GetSplineKnotCount(sIndex);
    if (knot < 0) {
        throw std::runtime_error("Knot index invalid");
    }
    if (count <= knot) {
        throw std::runtime_error("Knot index invalid");
    }
    return fSplineKnots->hostPtr()[knotsIndex+2+3*knot];
}

double Cache::Weight::Spline::GetSplineKnotSlope(int sIndex, int knot) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int knotsIndex = fSplineIndex->hostPtr()[sIndex];
    int count = GetSplineKnotCount(sIndex);
    if (knot < 0) {
        throw std::runtime_error("Knot index invalid");
    }
    if (count <= knot) {
        throw std::runtime_error("Knot index invalid");
    }
    return fSplineKnots->hostPtr()[knotsIndex+2+3*knot+1];
}

double Cache::Weight::Spline::GetSplineKnotPlace(int sIndex, int knot) {
    if (sIndex < 0) {
        throw std::runtime_error("Spline index invalid");
    }
    if (GetSplinesUsed() <= sIndex) {
        throw std::runtime_error("Spline index invalid");
    }
    int knotsIndex = fSplineIndex->hostPtr()[sIndex];
    int count = GetSplineKnotCount(sIndex);
    if (knot < 0) {
        throw std::runtime_error("Knot index invalid");
    }
    if (count <= knot) {
        throw std::runtime_error("Knot index invalid");
    }
    return fSplineKnots->hostPtr()[knotsIndex+2+3*knot+2];
}

#include "CacheAtomicMult.h"

namespace {
    // Interpolate one point.  This is the only place that changes when the
    // interpolation method changes.  This accepts a normalized "x" value, and
    // an array of control points with "dim" entries..  The control points
    // will be at (0, 1.0, 2.0, ... , dim-1).  The input variable "x" must be
    // a "floating point" index. If the index "x" is out of range, then this
    // turns into a linear extrapolation of the boundary points (try to avoid
    // that).
    //
    // Example: If the control points have dim of 5, the index "x" must be
    // greater than zero, and less than 5.  Assuming linear interpolation, an
    // input value of 2.1 results in the linear interpolation between element
    // [2] and element [3], or "(1.0-0.1)*p[2] + 0.1*p[3])".
    HEMI_DEV_CALLABLE_INLINE
    double HEMIInterp(int ix, double x, const WEIGHT_BUFFER_FLOAT* points, int dim) {
        double x1 = points[3*ix+2];
        double x2 = points[3*(ix+1)+2];
        double step = x2-x1;

        double fx = (x - x1)/step;
        double fxx = fx*fx;
        double fxxx = fx*fxx;

        double p1 = points[3*ix];
        double m1 = points[3*ix+1]*step;
        double p2 = points[3*(ix+1)];
        double m2 = points[3*(ix+1)+1]*step;

        // Cubic spline with the points and slopes.
        // double v = (p1*(2.0*fxxx-3.0*fxx+1.0) + m1*(fxxx-2.0*fxx+fx)
        //             + p2*(3.0*fxx-2.0*fxxx) + m2*(fxxx-fxx));

        // Linear interpolation
        double v = p1 + fx*(p2-p1);

        return v;
    }

    // A function to be used as the kernel on either the CPU or GPU.  This
    // must be valid CUDA coda.
    HEMI_KERNEL_FUNCTION(HEMISplinesKernel,
                         double* results,
                         const double* params,
                         const double* lowerClamp,
                         const double* upperClamp,
                         const WEIGHT_BUFFER_FLOAT* knots,
                         const int* rIndex,
                         const short* pIndex,
                         const int* sIndex,
                         const int NP) {
        for (int i : hemi::grid_stride_range(0,NP)) {
            const int id0 = sIndex[i];
            const int id1 = sIndex[i+1];
            const int dim = (id1-id0-2)/3;
            const double x = params[pIndex[i]];
#ifndef HEMI_DEV_CODE
            if (dim>15) std::runtime_error("To many bins in spline");
#endif
            int ix = 0;
            // Check to find a point that is less than x.
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 1
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 2
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 3
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 4
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 5
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 6
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 7
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 8
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 9
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 10
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 11
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 12
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 13
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 14
            if (x > knots[id0+2+3*(ix+1)+2] && ix < dim-2) ++ix; // 15

            const double s = 1.0/knots[id0+1];
            double v = HEMIInterp(ix, x, &knots[id0+2], dim);
#ifndef HEMI_DEV_CODE
#ifdef CACHE_DEBUG
            if (i < PRINT_STEP) {
                std::cout << "Splines kernel " << i
                          << " iEvt " << rIndex[i]
                          << " iPar " << pIndex[i]
                          << " = " << x
                          << " m " << knots[id0] << " d "  << knots[id0+1]
                          << " (" << x << "," << ix << ")"
                          << " --> " << v
                          << " s: " << s
                          << " d: " << dim
                          << std::endl;
                for (int k = 0; k < dim; ++k) {
                    std::cout << "        " << k
                              << " x: " << knots[id0+2+3*k+2]
                              << " y: " << knots[id0+2+3*k]
                              << " m: " << knots[id0+2+3*k+1]
                              << std::endl;
                }
            }
#endif
#endif
            const double lc = lowerClamp[pIndex[i]];
            if (v < lc) v = lc;
            const double uc = upperClamp[pIndex[i]];
            if (v > uc) v = uc;

            CacheAtomicMult(&results[rIndex[i]], v);
        }
    }
}

bool Cache::Weight::Spline::Apply() {
    if (GetSplinesUsed() < 1) return false;

    HEMISplinesKernel splinesKernel;
    hemi::launch(splinesKernel,
                 fWeights.writeOnlyPtr(),
                 fParameters.readOnlyPtr(),
                 fLowerClamp.readOnlyPtr(),
                 fUpperClamp.readOnlyPtr(),
                 fSplineKnots->readOnlyPtr(),
                 fSplineResult->readOnlyPtr(),
                 fSplineParameter->readOnlyPtr(),
                 fSplineIndex->readOnlyPtr(),
                 GetSplinesUsed()
        );

    return true;
}

// An MIT Style License

// Copyright (c) 2022 Clark McGrew

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/cmake/gundam-build.sh"
// End:
