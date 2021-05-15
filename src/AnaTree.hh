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
    double dist;
    double costh;
    double cosths;
    double omega;
    double timetof;
    int PMT_id;
    double weight;


public:
    AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name);
    ~AnaTree();

    void MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT = 19);

    long int GetEntry(long int entry) const;
    void SetBranches();
    void SetPMTBranches();
    std::vector<std::vector<AnaEvent>> GetEvents();

};

#endif
