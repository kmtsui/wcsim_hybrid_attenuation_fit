#ifndef TOYTHROWER_HH
#define TOYTHROWER_HH

#include <algorithm>
#include <iostream>

#include "TDecompChol.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TRandom3.h"
#include "TVectorT.h"

using TMatrixD = TMatrixT<double>;
using TMatrixDSym = TMatrixTSym<double>;
using TVectorD = TVectorT<double>;

class ToyThrower
{
    private:
        TMatrixD* covmat;
        TMatrixD* L_matrix;
        TVectorD* R_vector;

        unsigned int npar;
        unsigned int force_limit;

        ToyThrower(int nrows);

    public:
        ToyThrower(const TMatrixD &cov, bool do_setup = true, double decomp_tol = 0xCAFEBABE);
        ToyThrower(const TMatrixDSym &cov, bool do_setup = true, double decomp_tol = 0xCAFEBABE);
        ~ToyThrower();

        void SetMatrix(const TMatrixD& cov);
        void SetMatrix(const TMatrixDSym& cov);
        bool SetupDecomp(double decomp_tol = 0xCAFEBABE);
        bool ForcePosDef(double val, double decomp_tol = 0xCAFEBABE);

        bool CholDecomp();
        bool IncompCholDecomp(const double tol = 1.0E-3, bool modify_cov = true);

        void Throw(TVectorD& toy);
        void Throw(std::vector<double>& toy);

        double ThrowSinglePar(double nom, double err) const;

};

#endif
