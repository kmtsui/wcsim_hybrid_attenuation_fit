#include "AnaTree.hh"

AnaTree::AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name)
{
    fChain = new TChain(tree_name.c_str());
    fChain->Add(file_name.c_str());

    std::cout << TAG<<"Loading data files: "<<file_name.c_str()<<std::endl;

    std::string single_file_name = fChain->GetFile()->GetName();

    f_pmt = new TFile(single_file_name.c_str());
    t_pmt = (TTree*)f_pmt->Get(pmt_tree_name.c_str());

    std::cout << TAG<<"Loading PMT tree from : "<<f_pmt->GetName()<<std::endl;

    m_maskpmt = false;

}

AnaTree::~AnaTree()
{
    if(fChain != nullptr)
        delete fChain->GetCurrentFile();

}

void AnaTree::MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT)
{
    pmt_mask.clear();
    m_maskpmt = false;

    if (nPMT<=0) {
        std::cout << TAG << "In AnaTree::MaskPMT(), nPMT = " << nPMT << " <=0\n"
                  << "No PMT is masked"<<std::endl;

        return;
    }

    int nPMT_total = t_pmt->GetEntries();
    if (mPMT) nPMT_total = nPMT_total/nPMTpermPMT;

    if (nPMT>=nPMT_total) {
        std::cout << TAG << "In AnaTree::MaskPMT(), nPMT = " << nPMT << " >= total nPMT = " << nPMT_total <<"\n"
                  << "No PMT is masked"<<std::endl;

        return;
    }

    std::cout << TAG<< "In AnaTree::MaskPMT(), enabling "<< nPMT << " out of "<< nPMT_total << " PMTs" << std::endl;
    double PMT_frac = (nPMT+0.)/(nPMT_total);
    int PMT_count = 0;
    for (int i=0;i<nPMT_total;i++)
    {
        if ((PMT_count+0.)/(i+1.)<PMT_frac && PMT_count<nPMT) 
        {
            if (!mPMT) pmt_mask.emplace_back(0);
            else 
            {    
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    pmt_mask.emplace_back(0);
                }
            }
            PMT_count++;
        } 
        else 
        {
            if (!mPMT) pmt_mask.emplace_back(1);
            else
            {
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    pmt_mask.emplace_back(1);
                }
            }
        }
    }

    m_maskpmt = true;
}

long int AnaTree::GetEntry(long int entry) const
{
    // Read contents of entry.
    if(fChain == nullptr)
        return -1;
    else
        return fChain->GetEntry(entry);
}

void AnaTree::SetDataBranches()
{    
    // Set branch addresses and branch pointers

    fChain->SetBranchAddress("nPE", &nPE);
    fChain->SetBranchAddress("timetof", &timetof);
    fChain->SetBranchAddress("PMT_id", &PMT_id);

}

void AnaTree::SetPMTBranches()
{
    // Set branch addresses and branch pointers

    t_pmt->SetBranchAddress("R", &R);
    t_pmt->SetBranchAddress("costh", &costh);
    t_pmt->SetBranchAddress("costhm", &costhm);
    t_pmt->SetBranchAddress("cosths", &cosths);
    t_pmt->SetBranchAddress("phis", &phis);
    t_pmt->SetBranchAddress("phim", &phim);
    t_pmt->SetBranchAddress("omega", &omega);
    t_pmt->SetBranchAddress("PMT_id", &PMT_id);
    t_pmt->SetBranchAddress("mPMT_id", &mPMT_id);
    
}

std::vector<AnaEvent> AnaTree::GetPMTs()
{
    std::vector<AnaEvent> pmt_vec;

    if(t_pmt == nullptr)
    {
        std::cout << TAG<<"[Error] Reading no PMTs in AnaTree::GetPMTs()"<<std::endl;
        return pmt_vec;
    }
    unsigned long nentries = fChain->GetEntries();

    SetPMTBranches();
    for(long int jentry = 0; jentry < t_pmt->GetEntries(); jentry++)
    {
        t_pmt->GetEntry(jentry);

        if (m_maskpmt) if (pmt_mask.at(PMT_id)) continue;

        AnaEvent ev(jentry);
        ev.SetR(R);
        ev.SetCosth(costh);
        ev.SetCosthm(costhm);
        ev.SetCosths(cosths);
        ev.SetPhis(phis);
        ev.SetPhim(phim);
        ev.SetOmega(omega); 
        ev.SetPMTID(PMT_id);
        ev.SetmPMTID(mPMT_id);
        
        std::vector<double> reco_var;
        reco_var.emplace_back(costh); reco_var.emplace_back(R);
        ev.SetRecoVar(reco_var);

        pmt_vec.push_back(ev);

    }

    return pmt_vec;
}

bool AnaTree::GetDataEntry(unsigned long entry, double& time, double& charge, int& pmtID)
{
    fChain->GetEntry(entry);

    if (m_maskpmt) if (pmt_mask.at(PMT_id)) return false;

    time = timetof;

    charge = nPE;

    pmtID = PMT_id;

    return true;
}