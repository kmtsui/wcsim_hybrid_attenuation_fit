#ifndef ANATREE_HH
#define ANATREE_HH

#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH2F.h>

#include "AnaEvent.hh"
#include "ColorOutput.hh"

class AnaTree
{
private:
    TChain* fChain; //! pointer to the analyzed TTree or TChain
    TFile* f_pmt; 
    TTree* t_pmt; // pointer to PMT tree

    std::vector<int> pmt_mask;
    std::vector<int> mpmt_mask;
    bool m_maskpmt;
    bool m_maskmpmt;

    // Declaration of leaf types
    double nHits;
    double nPE;
    double R;
    double costh;
    double cosths;
    double phis;
    double costhm;
    double phim;
    double omega;
    double dz;
    double xpos, ypos, zpos;
    double timetof;
    int PMT_id;
    int mPMT_id;
    double weight;

    // histogram that stores all the hits, to speed-up loading process
    TH2F* data_hist;
    bool use_hist;

    const std::string TAG = color::GREEN_STR + "[AnaTree]: " + color::RESET_STR;
    const std::string ERR = color::RED_STR + "[AnaTree ERROR]: " + color::RESET_STR;
    const std::string WAR = color::RED_STR + "[AnaTree WARNING]: " + color::RESET_STR;

public:
    AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name, const std::string& hist_name="");
    ~AnaTree();

    void MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT = 19);
    void MaskmPMT(std::vector<int> vec, int nPMTpermPMT = 19);

    long int GetEntry(long int entry) const;
    void SetDataBranches();
    void SetPMTBranches();
    std::vector<AnaEvent> GetPMTs();
    bool GetDataEntry(unsigned long entry, double& time, double& charge, int& pmtID);
    unsigned long GetDataEntries() const { return fChain->GetEntries(); }
    int GetPMTEntries() const { return t_pmt->GetEntries(); }

    inline bool UseDataHist() const { return use_hist; }
    TH2F* GetDataHist() { return data_hist; }

    double GetEventVar(const std::string& var) const
    {
        if(var == "R")
            return R;
        else if(var == "costh")
            return costh;
        else if(var == "costhm")
            return costhm;
        else if(var == "cosths")
            return cosths;
        else if(var == "phis")
            return phis;
        else if(var == "phim")
            return phim;
        else if(var == "timetof")
            return timetof;
        else if(var == "nPE")
            return nPE;
        else if(var == "omega")
            return omega;
        else if(var == "PMT_id")
            return PMT_id;
        else if(var == "mPMT_id")
            return mPMT_id;
        else
        {
            std::cout << ERR <<" Error! Variable "<<var<<" not available in AnaTree"<<std::endl;
            return -1;
        }
    }

};

#endif
