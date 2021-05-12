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
#include "AnaSample.hh"

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
    double timetof;
    int PMT_id;
    float weight;

    std::vector<double> m_timetof_range;
    std::vector<double> m_cosths_range;

    inline void SetTimetofCut(std::vector<double> vec) { m_timetof_range = vec; }
    inline const std::vector<double>& GetTimetofCut() const { return m_timetof_range; }
    inline void SetCosthsCut(std::vector<double> vec) { m_cosths_range = vec; }
    inline const std::vector<double>& GetCosthsCut() const { return m_cosths_range; }

public:
    AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name);
    ~AnaTree();

    void MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT = 19);

    long int GetEntry(long int entry) const;
    void SetBranches();
    void SetPMTBranches();
    void GetEvents(AnaSample* ana_sample);
};

#endif
