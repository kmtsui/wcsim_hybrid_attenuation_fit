#include "AnaTree.hh"

AnaTree::AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name)
{
    fChain = new TChain(tree_name.c_str());
    fChain->Add(file_name.c_str());

    f_pmt = fChain->GetFile();
    t_pmt = (TTree*)f_pmt->Get(pmt_tree_name.c_str());

    m_maskpmt = false;

}

AnaTree::~AnaTree()
{
    if(fChain != nullptr)
        delete fChain->GetCurrentFile();
    if (f_pmt != nullptr)
        delete f_pmt;
}

void AnaTree::MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT)
{
    pmt_mask.clear();
    m_maskpmt = false;

    if (nPMT<=0) {
        std::cout << "In AnaTree::MaskPMT(), nPMT = " << nPMT << " <=0\n"
                  << "No PMT is masked"<<std::endl;
    }

    int nPMT_total = t_pmt->GetEntries();
    if (mPMT) nPMT_total = nPMT_total/nPMTpermPMT;

    if (nPMT>=nPMT_total) {
        std::cout << "In AnaTree::MaskPMT(), nPMT = " << nPMT << " >= total nPMT = " << nPMT_total <<"\n"
                  << "No PMT is masked"<<std::endl;
    }

    std::cout<< "Masking "<< nPMT << " out of "<< nPMT_total << "PMTs" << std::endl;
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

void AnaTree::SetBranches()
{
    // Set branch addresses and branch pointers

    fChain->SetBranchAddress("nPE", &nPE);
    fChain->SetBranchAddress("dist", &dist);
    fChain->SetBranchAddress("costh", &costh);
    fChain->SetBranchAddress("cosths", &cosths);
    fChain->SetBranchAddress("timetof", &timetof);
    fChain->SetBranchAddress("PMT_id", &PMT_id);

}

void AnaTree::SetPMTBranches()
{
    // Set branch addresses and branch pointers

    t_pmt->SetBranchAddress("dist", &dist);
    t_pmt->SetBranchAddress("costh", &costh);
    t_pmt->SetBranchAddress("cosths", &cosths);
    t_pmt->SetBranchAddress("PMT_id", &PMT_id);
    
}
void AnaTree::GetEvents(AnaSample* ana_sample)
{
    if(fChain == nullptr || ana_sample == nullptr)
        return;
    
    long int nentries = fChain->GetEntries();
    long int nbytes   = 0;

    std::cout  << "Reading events for "<<ana_sample->GetName()<<"...\n";
    SetBranches();
    for(long int jentry = 0; jentry < nentries; jentry++)
    {
        nbytes += fChain->GetEntry(jentry); 

        if (timetof<m_timetof_range.at(0) || timetof>m_timetof_range.at(1)) continue;
        if (cosths<m_cosths_range.at(0) || cosths>m_cosths_range.at(1)) continue;
        if (m_maskpmt) if (pmt_mask.at(PMT_id)) continue;

        AnaEvent ev(jentry);
        ev.SetPE(nPE);
        ev.SetR(dist);
        ev.SetCosth(-costh);
        ev.SetCosths(cosths); 
        ev.SetTimetof(timetof); 
        ev.SetPMTID(PMT_id);
        
        std::vector<double> reco_var;
        reco_var.emplace_back(-costh); reco_var.emplace_back(dist);
        ev.SetRecoVar(reco_var);

        ana_sample->AddEvent(ev);

    }

    std::cout  << "Reading PMT geometry for "<<ana_sample->GetName()<<"...\n";
    SetPMTBranches();
    for(long int jentry = 0; jentry < t_pmt->GetEntries(); jentry++)
    {
        t_pmt->GetEntry(jentry);

        if (cosths<m_cosths_range.at(0) || cosths>m_cosths_range.at(1)) continue;
        if (m_maskpmt) if (pmt_mask.at(PMT_id)) continue;

        AnaEvent ev(jentry);
        ev.SetR(dist);
        ev.SetCosth(-costh);
        ev.SetCosths(cosths);
        ev.SetPMTID(PMT_id);
        
        std::vector<double> reco_var;
        reco_var.emplace_back(-costh); reco_var.emplace_back(dist);
        ev.SetRecoVar(reco_var);

        ana_sample->AddPMT(ev);

    }

    ana_sample->PrintStats();
}
