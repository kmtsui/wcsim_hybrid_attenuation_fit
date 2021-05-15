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

#include "AnaEvent.hh"

class AnaTree
{
private:
    TChain* fChain; //! pointer to the analyzed TTree or TChain
    TFile* f_pmt; // pointer to PMT tree
    TTree* t_pmt;

    std::vector<int> pmt_mask;
    bool m_maskpmt;

    // Declaration of leaf types
    double nHits;
    double nPE;
    double R;
    double costh;
    double cosths;
    double omega;
    double timetof;
    int PMT_id;
    double weight;

    void SetLeafs();

    std::vector<std::string> m_bins;
    std::vector<std::string> m_cuts;
    std::vector<TLeaf*> leafs_bins;
    std::vector<TLeaf*> leafs_cuts;
    std::vector<double> vars_bin;
    std::vector<double> vars_cut;

public:
    AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name);
    ~AnaTree();

    void MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT = 19);

    long int GetEntry(long int entry) const;
    void SetDataBranches(std::vector<std::string> binvar, std::vector<std::string> cutvar);
    void SetPMTBranches();
    std::vector<AnaEvent> GetPMTs();
    void GetData(std::vector<std::vector<double>>& data_vec, std::vector<std::vector<double>>& cut_vec, std::vector<double>& weight_vec);
    bool GetDataEntry(unsigned long entry, std::vector<double>& data_vec, std::vector<double>& cut_vec, double& weight);
    unsigned long GetDataEntries() const { return fChain->GetEntries(); }

    double GetEventVar(const std::string& var) const
    {
        if(var == "R")
            return R;
        else if(var == "costh")
            return costh;
        else if(var == "cosths")
            return cosths;
        else if(var == "timetof")
            return timetof;
        else if(var == "nPE")
            return nPE;
        else if(var == "omega")
            return omega;
        else if(var == "PMT_id")
            return PMT_id;
        else
        {
            std::cout<<" Error! Variable "<<var<<" not available in AnaTree"<<std::endl;
            return -1;
        }
    }

};

#endif
