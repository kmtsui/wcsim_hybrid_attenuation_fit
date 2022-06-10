#include "WeightAttenuationZ.h"

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
Cache::Weight::AttenuationZ::AttenuationZ(
    Cache::Weights::Results& weights,
    Cache::Parameters::Values& parameters,
    std::size_t attenzs)
    : Cache::Weight::Base("attenuationz",weights,parameters),
      fAttenzsReserved(attenzs), fAttenzsUsed(0) {

    std::cout << "Cached Weights: reserved AttenuationZs: "
           << GetAttenzsReserved()
           << std::endl;
    if (GetAttenzsReserved() < 1) return;

    fTotalBytes += GetAttenzsReserved()*sizeof(int);   // fNormResult
    fTotalBytes += GetAttenzsReserved()*sizeof(short); // fNormParameter
    fTotalBytes += GetAttenzsReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fR
    fTotalBytes += GetAttenzsReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fOmega
    fTotalBytes += GetAttenzsReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fZ0
    fTotalBytes += GetAttenzsReserved()*sizeof(WEIGHT_BUFFER_FLOAT); // fDz

    std::cout << "Approximate Memory Size: " << fTotalBytes/1E+9
           << " GB" << std::endl;

    try {
        // Get the CPU/GPU memory for the normalization index tables.  These
        // are copied once during initialization so do not pin the CPU memory
        // into the page set.
        fAttenzResult.reset(new hemi::Array<int>(GetAttenzsReserved(),false));
        fAttenzParameter.reset(new hemi::Array<short>(GetAttenzsReserved(),false));
        fR.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetAttenzsReserved(),false));
        fOmega.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetAttenzsReserved(),false));
        fZ0.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetAttenzsReserved(),false));
        fDz.reset(new hemi::Array<WEIGHT_BUFFER_FLOAT>(GetAttenzsReserved(),false));
    }
    catch (std::bad_alloc&) {
        std::cerr << "Failed to allocate memory, so stopping" << std::endl;
        throw std::runtime_error("Not enough memory available");
    }

}

// The destructor
Cache::Weight::AttenuationZ::~AttenuationZ() {}

// Reserve space for another normalization parameter.
int Cache::Weight::AttenuationZ::ReserveAttenz(int resIndex, int parIndex, double R, double omega, double z0, double dz) {
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
    int newIndex = fAttenzsUsed++;
    if (fAttenzsUsed > fAttenzsReserved) {
        std::cerr << "Not enough space reserved for Attenzs"
                  << std::endl;
        throw std::runtime_error("Not enough space reserved for results");
    }
    fAttenzResult->hostPtr()[newIndex] = resIndex;
    fAttenzParameter->hostPtr()[newIndex] = parIndex;
    fR->hostPtr()[newIndex] = R;
    fOmega->hostPtr()[newIndex] = omega;
    fZ0->hostPtr()[newIndex] = z0;
    fDz->hostPtr()[newIndex] = dz;
    return newIndex;
}

#include "CacheAtomicMult.h"

namespace {
    // A function to be used as the kernen on a CPU or GPU.  This must be
    // valid CUDA.  This accumulates the normalization parameters into the
    // results.
    HEMI_KERNEL_FUNCTION(HEMIAttenzsKernel,
                         double* results,
                         const double* params,
                         const WEIGHT_BUFFER_FLOAT* R,
                         const WEIGHT_BUFFER_FLOAT* omega,
                         const WEIGHT_BUFFER_FLOAT* z0,
                         const WEIGHT_BUFFER_FLOAT* dz,
                         const int* rIndex,
                         const short* pIndex,
                         const int NP) {
        for (int i : hemi::grid_stride_range(0,NP)) {
            const double alpha_z0 = params[pIndex[i]] + params[pIndex[i]+1]*z0[i];
            const double da = params[pIndex[i]+1]*dz[i];
            double val = 1.;
            if (fabs(da)>1.e-9)
            {
                val = pow(1+da/alpha_z0,-R[i]/da)*omega[i];
            }
            else
            {
                val = exp(-R[i]/alpha_z0)*omega[i];
            }
            CacheAtomicMult(&results[rIndex[i]], val);
#ifndef HEMI_DEV_CODE
#ifdef CACHE_DEBUG
            if (rIndex[i] < PRINT_STEP) {
                std::cout << "Attenzs kernel " << i
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

bool Cache::Weight::AttenuationZ::Apply() {
    if (GetAttenzsUsed() < 1) return false;

    HEMIAttenzsKernel attenzsKernel;
    hemi::launch(attenzsKernel,
                 fWeights.writeOnlyPtr(),
                 fParameters.readOnlyPtr(),
                 fR->readOnlyPtr(),
                 fOmega->readOnlyPtr(),
                 fZ0->readOnlyPtr(),
                 fDz->readOnlyPtr(),
                 fAttenzResult->readOnlyPtr(),
                 fAttenzParameter->readOnlyPtr(),
                 GetAttenzsUsed());

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
