#ifndef __AnaFitParameters_hh__
#define __AnaFitParameters_hh__

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <TDecompLU.h>
#include <TMatrixTSym.h>
using TMatrixDSym = TMatrixTSym<double>;

#include "AnaSample.hh"
#include "BinManager.hh"
#include "ParameterFunction.hh"
#include "EigenDecomp.hh"
#include "ColorOutput.hh"
#include "ToyThrower.hh"

// some error codes
const int PASSEVENT = -1;
const int BADBIN    = -2;

class AnaFitParameters
{
public:
    AnaFitParameters(const std::string& par_name, const int pmttype);
    ~AnaFitParameters();

    void InitParameters(std::vector<std::string> names, std::vector<double> priors, std::vector<double> steps, 
                        std::vector<double> lows, std::vector<double> highs, std::vector<bool> fixed);
    void InitEventMap(std::vector<AnaSample*>& sample);
    void ApplyParameters(std::vector<double>& params);
    void ReWeight(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params);
    double GetWeight(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params);
    int GetParBin(int nsample, int nevent) const { return m_evmap[nsample][nevent]; }

    std::string GetName() const { return m_name; }

    void GetParNames(std::vector<std::string>& vec) const { vec = pars_name; }
    void GetParPriors(std::vector<double>& vec) const { vec = pars_prior; }
    void GetParValues(std::vector<double>& vec) const { vec = pars_value; }
    void GetParOriginal(std::vector<double>& vec) const { vec = pars_original; }
    double GetParOriginal(int i) const { return pars_original.at(i); }
    double GetParPrior(int i) const { return pars_prior.at(i); }
    double GetParValue(int i) const { return pars_value.at(i); }
    void GetParSteps(std::vector<double>& vec) const { vec = pars_step; }
    void GetParFixed(std::vector<bool>& vec) const { vec = pars_fixed; }
    void GetParLimits(std::vector<double>& vec1, std::vector<double>& vec2) const
    {
        vec1 = pars_limlow;
        vec2 = pars_limhigh;
    }

    void SetParNames(std::vector<std::string>& vec) { pars_name = vec; }
    void SetParPriors(std::vector<double>& vec) { pars_prior = vec; }
    void SetParValues(std::vector<double>& vec) { pars_value = vec; }
    void SetParSteps(std::vector<double>& vec) { pars_step = vec; }
    void SetParLimits(std::vector<double>& vec1, std::vector<double>& vec2)
    {
        pars_limlow  = vec1;
        pars_limhigh = vec2;
    }
    void SetParFixed(std::vector<bool>& vec) { pars_fixed = vec; }

    int GetNpar() const { return Npar; }
    void SetInfoFrac(double frac) { m_info_frac = frac; }
    double GetInfoFrac() const { return m_info_frac; }
    bool DoRNGstart() const { return m_rng_start; }
    void SetRNGstart(bool flag = true) { m_rng_start = flag; }
    void SetWeightCap(double cap, bool flag = true) { m_do_cap_weights = flag; m_weight_cap = cap; }
    void SetThrowPars(bool flag = true) { m_pars_throw = flag; }
    void ThrowPars();

    inline void SetBinning(const std::string& binning) { m_bm.SetBinning(binning); }

    inline const std::vector<std::string>& GetBinVar() const { return m_binvar; }
    inline void SetBinVar(std::vector<std::string> vec) { m_binvar = vec; }

    void SetParameterFunction(const std::string& func_name);
    inline int GetParameterFunctionType() const { return m_func_type; }
    ParameterFunction* GetParameterFunction() const { return m_func ; }

    inline void SetPMTType(const int val) { m_pmttype = val; }
    inline int GetPMTType() const { return m_pmttype; }

    void SetCovarianceMatrix(const TMatrixDSym& covmat, bool decompose = false);
    TMatrixDSym* GetCovMat() const { return covariance; }
    bool HasCovMat() const { return covariance != nullptr; }
    double GetChi2(const std::vector<double>& params) const;

    bool IsDecomposed() const { return m_decompose; }
    TMatrixDSym* GetOriginalCovMat() const { return original_cov; }
    TMatrixDSym GetOriginalCovMat(const TMatrixDSym& cov, unsigned int start_idx) const
    {
        return eigen_decomp->GetOriginalCovMat(cov, start_idx);
    }
    std::vector<double> GetOriginalParameters(const std::vector<double>& param) const
    {
        return eigen_decomp->GetOriginalParameters(param);
    }
    std::vector<double> GetOriginalParameters(const std::vector<double>& param,
                                              unsigned int start_idx) const
    {
        return eigen_decomp->GetOriginalParameters(param, start_idx);
    }

    void SetSpline(const std::vector<std::string> file_name, const std::vector<std::string> spline_name);
    void LoadSpline(std::vector<AnaSample*>& sample);
    void ReWeightSpline(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params);


protected:
    bool CheckDims(const std::vector<double>& params) const;

    std::size_t Npar;
    std::string m_name;
    std::vector<std::string> pars_name;
    std::vector<double> pars_original;
    std::vector<double> pars_prior;
    std::vector<double> pars_value;
    std::vector<double> pars_step;
    std::vector<double> pars_limlow;
    std::vector<double> pars_limhigh;
    std::vector<bool> pars_fixed;

    std::vector<std::vector<int>> m_evmap;
    bool m_rng_start;
    bool m_do_cap_weights;
    double m_weight_cap;
    double m_info_frac;
    bool m_pars_throw;

    int m_pmttype;
    std::vector<std::string> m_binvar;
    BinManager m_bm;

    ParameterFunction* m_func;
    int m_func_type;

    EigenDecomp* eigen_decomp;
    TMatrixDSym* covariance;
    TMatrixDSym* covarianceI;
    TMatrixDSym* original_cov;
    bool m_decompose;

    std::vector<int> pol_orders; // order of polynomial in each piece 
    std::vector<double> pol_range; // applicable range for each polynomial

    // spline reweight for indirect photon
    bool m_spline;
    std::vector<std::string> m_spline_file_name;
    std::vector<std::string> m_spline_name;
    std::vector<std::vector<std::vector<TGraph*>>> spline;
    //std::vector<TGraph> spline;

    const std::string TAG = color::GREEN_STR + "[AnaFitParameters]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + "[AnaFitParameters ERROR]: " + color::RESET_STR;
    const std::string WAR = color::RED_STR + "[AnaFitParameters WARNING]: " + color::RESET_STR;
};

#endif
