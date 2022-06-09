#include "WeightSourcePhiVar.h"

#include <algorithm>
#include <iostream>
#include <exception>
#include <limits>
#include <cmath>

#include <hemi/hemi_error.h>
#include <hemi/launch.h>
#include <hemi/grid_stride_range.h>

// The constructor
Cache::Weight::SourcePhiVar::SourcePhiVar(
    Cache::Weights::Results& weights,
    Cache::Parameters::Values& parameters,
    std::size_t sphis)
    : Cache::Weight::Base("SourcePhiVar",weights,parameters),
      fSPhisReserved(sphis), fSPhisUsed(0) {

    std::cout << "Cached Weights: reserved SourcePhiVars: "
           << GetSPhisReserved()
           << std::endl;
    if (GetSPhisReserved() < 1) return;

    fTotalBytes += GetSPhisReserved()*sizeof(int);   // fSPhiResult
    fTotalBytes += GetSPhisReserved()*sizeof(short); // fSPhiParameter
    fTotalBytes += GetSPhisReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fCosPhis

    std::cout << "Approximate Memory Size: " << fTotalBytes/1E+9
           << " GB" << std::endl;

    try {
        // Get the CPU/GPU memory for the normalization index tables.  These
        // are copied once during initialization so do not pin the CPU memory
        // into the page set.
        fSPhiResult.reset(new hemi::Array<int>(GetSPhisReserved(),false));
        fSPhiParameter.reset(new hemi::Array<short>(GetSPhisReserved(),false));
        fCosPhis.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetSPhisReserved(),false));
    }
    catch (std::bad_alloc&) {
        std::cerr << "Failed to allocate memory, so stopping" << std::endl;
        throw std::runtime_error("Not enough memory available");
    }

}

// The destructor
Cache::Weight::SourcePhiVar::~SourcePhiVar() {}

// Reserve space for another normalization parameter.
int Cache::Weight::SourcePhiVar::ReserveSPhi(int resIndex, int parIndex, double cosphis) {
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
    int newIndex = fSPhisUsed++;
    if (fSPhisUsed > fSPhisReserved) {
        std::cerr << "Not enough space reserved for SPhis"
                  << std::endl;
        throw std::runtime_error("Not enough space reserved for results");
    }
    fSPhiResult->hostPtr()[newIndex] = resIndex;
    fSPhiParameter->hostPtr()[newIndex] = parIndex;
    fCosPhis->hostPtr()[newIndex] = cosphis;
    return newIndex;
}

#include "CacheAtomicMult.h"

namespace {
    // A function to be used as the kernen on a CPU or GPU.  This must be
    // valid CUDA.  This accumulates the normalization parameters into the
    // results.
    HEMI_KERNEL_FUNCTION(HEMISPhisKernel,
                         double* results,
                         const double* params,
                         const WEIGHT_BUFFER_FLOAT* cosphis,
                         const int* rIndex,
                         const short* pIndex,
                         const int NP) {
        for (int i : hemi::grid_stride_range(0,NP)) {
            CacheAtomicMult(&results[rIndex[i]], 1+params[pIndex[i]]*cosphis[i]);
#ifndef HEMI_DEV_CODE
#ifdef CACHE_DEBUG
            if (rIndex[i] < PRINT_STEP) {
                std::cout << "SPhis kernel " << i
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

bool Cache::Weight::SourcePhiVar::Apply() {
    if (GetSPhisUsed() < 1) return false;

    HEMISPhisKernel sphisKernel;
    hemi::launch(sphisKernel,
                 fWeights.writeOnlyPtr(),
                 fParameters.readOnlyPtr(),
                 fCosPhis->readOnlyPtr(),
                 fSPhiResult->readOnlyPtr(),
                 fSPhiParameter->readOnlyPtr(),
                 GetSPhisUsed());

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
