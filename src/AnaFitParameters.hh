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
    void ReWeight(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params);

    std::string GetName() const { return m_name; }

    void GetParNames(std::vector<std::string>& vec) const { vec = pars_name; }
    void GetParPriors(std::vector<double>& vec) const { vec = pars_prior; }
    void GetParOriginal(std::vector<double>& vec) const { vec = pars_original; }
    double GetParOriginal(int i) const { return pars_original.at(i); }
    double GetParPrior(int i) const { return pars_prior.at(i); }
    void GetParSteps(std::vector<double>& vec) const { vec = pars_step; }
    void GetParFixed(std::vector<bool>& vec) const { vec = pars_fixed; }
    void GetParLimits(std::vector<double>& vec1, std::vector<double>& vec2) const
    {
        vec1 = pars_limlow;
        vec2 = pars_limhigh;
    }

    void SetParNames(std::vector<std::string>& vec) { pars_name = vec; }
    void SetParPriors(std::vector<double>& vec) { pars_prior = vec; }
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

    inline void SetBinning(const std::string& binning) { m_bm.SetBinning(binning); }

    inline const std::vector<std::string>& GetBinVar() const { return m_binvar; }
    inline void SetBinVar(std::vector<std::string> vec) { m_binvar = vec; }

    void SetParameterFunction(const std::string& func_name);

    inline void SetPMTType(const int val) { m_pmttype = val; }
    inline int GetPMTType() const { return m_pmttype; }

protected:
    bool CheckDims(const std::vector<double>& params) const;

    std::size_t Npar;
    std::string m_name;
    std::vector<std::string> pars_name;
    std::vector<double> pars_original;
    std::vector<double> pars_prior;
    std::vector<double> pars_step;
    std::vector<double> pars_limlow;
    std::vector<double> pars_limhigh;
    std::vector<bool> pars_fixed;

    std::vector<std::vector<int>> m_evmap;
    bool m_rng_start;
    bool m_do_cap_weights;
    double m_weight_cap;
    double m_info_frac;

    int m_pmttype;
    std::vector<std::string> m_binvar;
    BinManager m_bm;

    ParameterFunction* m_func;
};

#endif
