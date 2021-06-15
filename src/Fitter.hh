#ifndef __Fitter_hh__
#define __Fitter_hh__

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <TFile.h>
#include <TGraph.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TVectorT.h>
#include <TMath.h>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "AnaSample.hh"
#include "AnaFitParameters.hh"

struct MinSettings
{
    std::string minimizer;
    std::string algorithm;
    std::string likelihood;
    int print_level;
    int strategy;
    double tolerance;
    double max_iter;
    double max_fcn;
};


class Fitter
{
public:
    Fitter(TDirectory* dirout, const int seed);
    Fitter(TDirectory* dirout, const int seed, const int num_threads);
    ~Fitter();
    double CalcLikelihood(const double* par);
    void InitFitter(std::vector<AnaFitParameters*>& fitpara);

    void FixParameter(const std::string& par_name, const double& value);
    bool Fit(const std::vector<AnaSample*>& samples, bool stat_fluc=false);
    void ParameterScans(const std::vector<int>& param_list, unsigned int nsteps);

    void SetMinSettings(const MinSettings& ms);
    void SetSeed(int seed);
    void SetZeroSyst(bool flag) { m_zerosyst = flag; }
    void SetNumThreads(const unsigned int num) { m_threads = num; }
    void SetSaveFreq(int freq, bool flag = true)
    {
        m_freq = freq;
        m_save = flag;
    }
    void SetSaveEvents(bool flag = true) { m_save_events = flag; };

    // Declaration of leaf types
    int sampleId;
    double nPE;
    double R;
    double costh;
    double cosths;
    double costhm;
    double phim;
    double omega;
    int PMT_id;
    int mPMT_id;
    std::vector<double> weight;

    void InitOutputTree()
    {
        m_outtree->Branch("sampleId", &sampleId, "sampleId/I");
        m_outtree->Branch("nPE", &nPE, "nPE/D");
        m_outtree->Branch("R", &R, "R/D");
        m_outtree->Branch("costh", &costh, "costh/D");
        m_outtree->Branch("cosths", &cosths, "cosths/D");
        m_outtree->Branch("costhm", &costhm, "costhm/D");
        m_outtree->Branch("phim", &phim, "phim/D");
        m_outtree->Branch("omega", &omega, "omega/D");
        m_outtree->Branch("PMT_id", &PMT_id, "PMT_id/I");
        m_outtree->Branch("mPMT_id", &mPMT_id, "mPMT_id/I");
        m_outtree->Branch("weight", &weight);
    }

private:
    double FillSamples(std::vector<std::vector<double>>& new_pars);
    void SaveParams(const std::vector<std::vector<double>>& new_pars);
    void SaveEventHist(bool is_final = false);
    void SaveEventTree(std::vector<std::vector<double>>& par_results);
    void SaveChi2();
    void SaveResults(const std::vector<std::vector<double>>& parresults,
                     const std::vector<std::vector<double>>& parerrors);

    ROOT::Math::Minimizer* m_fitter;
    ROOT::Math::Functor* m_fcn;

    TTree* m_outtree;
    TRandom3* rng;
    TDirectory* m_dir;
    bool m_save;
    bool m_save_events;
    bool m_zerosyst;
    int m_threads;
    int m_npar, m_calls, m_freq;
    std::vector<std::string> par_names;
    std::vector<double> par_prefit;
    std::vector<double> par_postfit;
    std::vector<int> par_type;
    std::vector<int> par_pmttype;
    std::vector<std::string> par_var;
    std::vector<double> par_var_low;
    std::vector<double> par_var_high;
    std::vector<double> vec_chi2_stat;
    std::vector<double> vec_chi2_sys;
    std::vector<double> vec_chi2_reg;
    std::vector<AnaFitParameters*> m_fitpara;
    std::vector<AnaSample*> m_samples;

    MinSettings min_settings;
};
#endif
