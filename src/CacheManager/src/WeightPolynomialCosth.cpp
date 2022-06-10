#include "WeightPolynomialCosth.h"

#include <algorithm>
#include <iostream>
#include <exception>
#include <limits>
#include <cmath>

#include <hemi/hemi_error.h>
#include <hemi/launch.h>
#include <hemi/grid_stride_range.h>

// #include "Logger.h"
// LoggerInit([](){
//   Logger::setUserHeaderStr("[Cache]");
// })

// The constructor
Cache::Weight::PolynomialCosth::PolynomialCosth(
    Cache::Weights::Results& weights,
    Cache::Parameters::Values& parameters,
    std::size_t polys)
    : Cache::Weight::Base("polynomialCosth",weights,parameters),
      fPolysReserved(polys), fPolysUsed(0) {

    std::cout << "Cached Weights: reserved PolynomialCosths: "
           << GetPolysReserved()
           << std::endl;
    if (GetPolysReserved() < 1) return;

    fTotalBytes += GetPolysReserved()*sizeof(int);   // fPolyResult
    fTotalBytes += GetPolysReserved()*sizeof(short); // fPolyParameter
    fTotalBytes += GetPolysReserved()*sizeof(short);   // fPolyOrder
    fTotalBytes += GetPolysReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fPolyArg

    std::cout << "Approximate Memory Size: " << fTotalBytes/1E+9
           << " GB" << std::endl;

    try {
        // Get the CPU/GPU memory for the PolynomialCosth index tables.  These
        // are copied once during initialization so do not pin the CPU memory
        // into the page set.
        fPolyResult.reset(new hemi::Array<int>(GetPolysReserved(),false));
        fPolyParameter.reset(new hemi::Array<short>(GetPolysReserved(),false));
        fPolyOrder.reset(new hemi::Array<short>(GetPolysReserved(),false));
        fPolyArg.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetPolysReserved(),false));
    }
    catch (std::bad_alloc&) {
        std::cerr << "Failed to allocate memory, so stopping" << std::endl;
        throw std::runtime_error("Not enough memory available");
    }

}

// The destructor
Cache::Weight::PolynomialCosth::~PolynomialCosth() {}

// void Cache::Weight::SetPolynominal(int idx, std::vector<int> orders, std::vector<double> ranges)
// {
//     int newIndex = fNPolysUsed++;
//     if (fNPolysUsed > fNPolysReserved) {
//         std::cerr << "Not enough space reserved for NPolys"
//                   << std::endl;
//         throw std::runtime_error("Not enough space reserved for NPolys");
//     }
//     fStartIdx->hostPtr()[newIndex] = idx;

//     for (int i=0;i<orders.size();i++)
//     {
// 	    newIndex = fOrdersUsed++;
//     	if (fOrdersUsed > fOrdersReserved) 
//         {
//             std::cerr << "Not enough space reserved for Orders"
//                       << std::endl;
//             throw std::runtime_error("Not enough space reserved for Orders");
//         }
//         fOrder->hostPtr()[newIndex] = orders[i];
//     }

//     for (int i=0;i<ranges.size();i++)
//     {
// 	    newIndex = fRangesUsed++;
//     	if (fRangesUsed > fRangesReserved) 
//         {
//             std::cerr << "Not enough space reserved for Ranges"
//                       << std::endl;
//             throw std::runtime_error("Not enough space reserved for Ranges");
//         }
//         fRange->hostPtr()[newIndex] = ranges[i];
//     }

// }

// Reserve space for another PolynomialCosth parameter.
int Cache::Weight::PolynomialCosth::ReservePoly(int resIndex, int parIndex, int order, double polyArg) {
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
    int newIndex = fPolysUsed++;
    if (fPolysUsed > fPolysReserved) {
        std::cerr << "Not enough space reserved for Polys"
                  << std::endl;
        throw std::runtime_error("Not enough space reserved for results");
    }
    fPolyResult->hostPtr()[newIndex] = resIndex;
    fPolyParameter->hostPtr()[newIndex] = parIndex;
    fPolyOrder->hostPtr()[newIndex] = order;
    fPolyArg->hostPtr()[newIndex] = polyArg;
    return newIndex;
}

#include "CacheAtomicMult.h"

namespace {
    // A function to be used as the kernen on a CPU or GPU.  This must be
    // valid CUDA.  This accumulates the PolynomialCosth parameters into the
    // results.
    HEMI_KERNEL_FUNCTION(HEMIPolysKernel,
                         double* results,
                         const double* params,
                         const short* orders,
                         const WEIGHT_BUFFER_FLOAT* polyArgs,
                         const int* rIndex,
                         const short* pIndex,
                         const int NP) {
        for (int i : hemi::grid_stride_range(0,NP)) {
            double val = 0.;
            for (int j = orders[i] ; j>=0 ; j-- )
                val = val*polyArgs[i] + params[pIndex[i]+j];
            CacheAtomicMult(&results[rIndex[i]], val);
#ifndef HEMI_DEV_CODE
#ifdef CACHE_DEBUG
            if (rIndex[i] < PRINT_STEP) {
                std::cout << "Polys kernel " << i
                       << " iEvt " << rIndex[i]
                       << " iPar " << pIndex[i]
                       << " = " << params[pIndex[i]]
                       << std::endl;
            }
#endif
#endif
        }
    }
}

bool Cache::Weight::PolynomialCosth::Apply() {
    if (GetPolysUsed() < 1) return false;

    HEMIPolysKernel polysKernel;
    hemi::launch(polysKernel,
                 fWeights.writeOnlyPtr(),
                 fParameters.readOnlyPtr(),
                 fPolyOrder->readOnlyPtr(),
                 fPolyArg->readOnlyPtr(),
                 fPolyResult->readOnlyPtr(),
                 fPolyParameter->readOnlyPtr(),
                 GetPolysUsed());

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
