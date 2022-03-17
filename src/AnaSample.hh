#ifndef __AnaSample_hh__
#define __AnaSample_hh__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TGraph.h>

#include "AnaEvent.hh"
#include "AnaTree.hh"
#include "BinManager.hh"
#include "Likelihoods.hh"
#include "ColorOutput.hh"

class AnaSample
{
public:
    AnaSample(int sample_id, const std::string& name, const std::string& binning, int pmt_type);
    ~AnaSample();

    void InitEventMap();

    inline int GetNPMTs() const { return m_pmts.size(); }
    inline void ClearPMTs() { m_pmts.clear(); }
    inline void AddPMT(const AnaEvent& pmt_entry) { m_pmts.push_back(pmt_entry); }

    AnaEvent* GetPMT(const unsigned int evnum);

    void LoadEventsFromFile(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name);

    void PrintStats() const;
    void MakeHistos();
    void SetNorm(const double val) { m_norm = val; }

    void SetLLHFunction(const std::string& func_name);
    double CalcLLH() const;

    void FillEventHist(bool reset_weights = false);
    void FillDataHist(bool stat_fluc = false);

    void WriteEventHist(TDirectory* dirout, const std::string& bsname);
    void WriteDataHist(TDirectory* dirout, const std::string& bsname);

    inline double GetNorm() const { return m_norm; }
    inline int GetSampleID() const { return m_sample_id; }
    inline std::string GetName() const { return m_name; }
    inline int GetPMTType() const { return m_pmttype; }

    inline void SetCut(const std::string& var_name, const double cut_low, const double cut_high) { m_cutvar.push_back(var_name); m_cutlow.push_back(cut_low); m_cuthigh.push_back(cut_high);}
    inline void ResetCut() { m_cutvar.clear(); m_cutlow.clear(); m_cuthigh.clear(); }

    inline void MaskPMT(const int nPMT) { m_pmtmask = nPMT; }
    inline void SetnPMTpermPMT(const int nPMTpermPMT) { m_nPMTpermPMT = nPMTpermPMT; }
    inline void MaskmPMT(std::vector<int> vec) { m_mPMTmask = vec; }

    inline const std::vector<std::string>& GetBinVar() const { return m_binvar; }
    inline void SetBinVar(std::vector<std::string> vec) { m_binvar = vec; }

    inline void SetStatFluc(bool val) { m_stat_fluc = val; }

    inline void SetScatter(double time1, double time2, double time3, double factor)
    {
        m_scatter = true; m_scatter_time1 = time1; m_scatter_time2 = time2; m_scatter_time3 = time3; m_scatter_factor = factor; 
    }

    void SetScatterMap(double time1, double time2, double time3, const TH1D& hist);

    void SetTemplate(const TH2D& hist, double offset, bool combine, bool template_only);
    inline bool UseTemplate() { return m_template; } 
    TH2D* GetTemplate() { return m_htimetof_pmt_pred;}

    void SetPMTEff(const TH1D& hist);
    inline void SetPMTEffVar(double val) { m_eff_var  = true; m_eff_sig = val; }

    inline void UnsetScatter() { m_scatter = false; m_scatter_map = false; }

    inline void SetTimeOffset(bool offset, double width) { m_time_offset = offset; m_time_offset_width = width; }

    inline void SetTimeSmear(bool smear, double mean, double width) { m_time_smear = smear; m_time_smear_mean = mean; m_time_smear_width = width; }

    inline void SetZ0(double val) { m_z0 = val; }

    void InitToy();

protected:
    int m_sample_id;
    int m_nbins;
    double m_norm;
    int m_pmttype;
    int m_pmtmask;
    int m_nPMTpermPMT;
    std::vector<int> m_mPMTmask;
    bool m_stat_fluc;

    bool m_scatter;
    double m_scatter_time1;
    double m_scatter_time2;
    double m_scatter_time3;
    double m_scatter_factor;

    bool m_scatter_map;
    TH1D* m_h_scatter_map;

    bool m_time_offset;
    double m_time_offset_width;

    bool m_time_smear;
    double m_time_smear_mean;
    double m_time_smear_width;

    double m_z0;

    bool m_use_eff;
    bool m_eff_var;
    double m_eff_sig;

    AnaTree* selTree;

    std::vector<std::string> m_binvar;

    std::vector<std::string> m_cutvar;
    std::vector<double> m_cutlow;
    std::vector<double> m_cuthigh;

    std::string m_name;
    std::string m_binning;
    std::vector<AnaEvent> m_pmts;  // Geometry info etc. for each PMT

    BinManager m_bm;
    CalcLLHFunc* m_llh;

    TH1D* m_hpred; // direct PE prediction
    //TH1D* m_hpred_indirect; // indirect PE prediction
    TH1D* m_hpred_err2; // MC stat error in PE prediction
    TH1D* m_hdata; // data histogram
    TH1D* m_hdata_control; // data histogram in control region
    TH1D* m_hdata_pmt; // data histogram per PMT 
    TH1D* m_hdata_pmt_control; // data histogram in control region per PMT
    TH1D* m_pmt_eff;

    // template fit of timetof distribution for indirect PE
    TH2D* m_htimetof_pred; // timetof distribution for fit
    TH2D* m_htimetof_pred_w2; 
    TH2D* m_htimetof_data;
    TH2D* m_htimetof_pmt_pred; // timetof distribution for each PMT
    TH2D* m_htimetof_pmt_data;
    bool m_template, m_template_combine, m_template_only;
    double m_timetof_offset;

    const std::string TAG = color::GREEN_STR + "[AnaSample]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + "[AnaSample ERROR]: " + color::RESET_STR;
    const std::string WAR = color::RED_STR + "[AnaSample WARNING]: " + color::RESET_STR;

};

#endif
