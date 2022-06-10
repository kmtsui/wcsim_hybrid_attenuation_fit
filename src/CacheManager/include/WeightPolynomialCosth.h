#ifndef WeightPolynomialCosth_hxx_seen
#define WeightPolynomialCosth_hxx_seen

#include "WeightBase.h"

#include "hemi/array.h"

#include <cstdint>
#include <memory>

namespace Cache {
    namespace Weight {
        class Base;
        class PolynomialCosth;
    }
}

/// A class apply a PolynomialCosth parameter to the cached event weights.  This
/// will be used in Cache::Weights to run the GPU for this type of
/// reweighting.
class Cache::Weight::PolynomialCosth: public Cache::Weight::Base {
private:

    ///////////////////////////////////////////////////////////////////////
    /// An array of indices into the results that go for each PolynomialCosth.
    /// This is copied from the CPU to the GPU once, and is then constant.
    std::size_t fPolysReserved;
    std::size_t fPolysUsed;
    std::unique_ptr<hemi::Array<int>> fPolyResult;

    /// An array of indices into the parameters that go for each
    /// PolynomialCosth.  This is copied from the CPU to the GPU once, and is
    /// then constant.
    std::unique_ptr<hemi::Array<short>> fPolyParameter;
    std::unique_ptr<hemi::Array<short>> fPolyOrder;
    std::unique_ptr<hemi::Array<WEIGHT_BUFFER_FLOAT>> fPolyArg;

    // std::size_t fNPolysReserved;
    // std::size_t fNPolysUsed;
    // std::unique_ptr<hemi::Array<short>> fStartIdx;
    // std::size_t fOrdersReserved;
    // std::size_t fOrdersUsed;
    // std::unique_ptr<hemi::Array<short>> fOrder;
    // std::size_t fRangesReserved;
    // std::size_t fRangesUsed;
    // std::unique_ptr<hemi::Array<short>> fRange;
    // std::size_t fCoeffsReserved;
    // std::unique_ptr<hemi::Array<WEIGHT_BUFFER_FLOAT>> fCoeff;

public:

    // Construct the class.  This should allocate all the memory on the host
    // and on the GPU.  The PolynomialCosths are applied to the event weights
    // which are managed by the Weights class.
    PolynomialCosth(Cache::Weights::Results& weights,
                  Cache::Parameters::Values& parameters,
                  std::size_t polys);

    // Deconstruct the class.  This should deallocate all the memory
    // everyplace.
    virtual ~PolynomialCosth();

    // void SetPolynominal(int idx, std::vector<int> orders, std::vector<double> ranges);

    /// Return the number of PolynomialCosth parameters that are reserved
    std::size_t GetPolysReserved() {return fPolysReserved;}

    /// Return the number of PolynomialCosth parameters that are used.
    std::size_t GetPolysUsed() {return fPolysUsed;}

    // std::size_t GetNPolysReserved() {return fNPolysReserved;}
    // std::size_t GetOrdersReserved() {return fOrdersReserved;}
    // std::size_t GetRangesReserved() {return fRangesReserved;}

    // /// Return the number of PolynomialCosth coefficient parameters that are reserved
    // std::size_t GetCoeffsReserved() {return fCoeffsReserved;}

    /// Add a PolynomialCosth parameter. This takes a result index, and a
    /// parameter index as inputs.
    int ReservePoly(int resIndex, int parIndex, int order, double polyArg);

    /// Apply the PolynomialCosths to the event weight cache.  This will run a
    /// HEMI kernel to modify the weights cache.
    bool Apply();
};

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
#endif
