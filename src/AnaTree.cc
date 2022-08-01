#include "AnaTree.hh"

AnaTree::AnaTree(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name, const std::string& hist_name)
    : data_hist(nullptr)
{
    fChain = new TChain(tree_name.c_str());
    fChain->Add(file_name.c_str());

    std::cout << TAG<<"Loading data files: "<<file_name.c_str()<<std::endl;

    // only need the first file to read PMT tree
    std::string single_file_name = fChain->GetFile()->GetName();

    f_pmt = new TFile(single_file_name.c_str());
    t_pmt = (TTree*)f_pmt->Get(pmt_tree_name.c_str());

    std::cout << TAG<<"Loading PMT tree from : "<<f_pmt->GetName()<<std::endl;

    m_maskpmt = false;
    m_maskmpmt = false;

    use_hist = false;
    if(data_hist != nullptr)
        delete data_hist;
    if (hist_name.size()>0)
    {
        TH2F* hist = (TH2F*)f_pmt->Get(hist_name.c_str());
        if (!hist) return;

        data_hist = (TH2F*)hist->Clone();
        data_hist->SetDirectory(0);
        use_hist = true;

        for (int i=1; i<fChain->GetListOfFiles()->GetEntries(); i++ )
        {
            TFile* f = new TFile(fChain->GetListOfFiles()->At(i)->GetTitle());
            if (f && f->Get(hist_name.c_str()))
            {
                TH2F* h = (TH2F*)f->Get(hist_name.c_str());
                if (h) data_hist->Add(h);
            }
            else std::cout << WAR << "Cannot find TH2F " << hist_name.c_str() << " in " << f->GetName() << std::endl;
            f->Close();
        }

        std::cout << TAG << "Loaded " << hist_name << " as data histogram " << std::endl;
    }
}

AnaTree::~AnaTree()
{
    if(fChain != nullptr)
        delete fChain->GetCurrentFile();

    if(data_hist != nullptr)
        delete data_hist;
}

void AnaTree::MaskPMT(int nPMT, bool mPMT, int nPMTpermPMT)
{
    // mask the whole PMT
    pmt_mask.clear();
    m_maskpmt = false;

    if (nPMT<=0) {
        std::cout << TAG << "In MaskPMT(), nPMT = " << nPMT << " <=0\n"
                  << "No PMT is masked"<<std::endl;

        return;
    }

    int nPMT_total = t_pmt->GetEntries();
    if (mPMT) nPMT_total = nPMT_total/nPMTpermPMT;

    if (nPMT>=nPMT_total) {
        std::cout << WAR << "In MaskPMT(), nPMT = " << nPMT << " >= total nPMT = " << nPMT_total <<"\n"
                  << "No PMT is masked"<<std::endl;

        return;
    }

    std::cout << TAG<< "In MaskPMT(), enabling "<< nPMT << " out of "<< nPMT_total << " PMTs" << std::endl;
    double PMT_frac = (nPMT+0.)/(nPMT_total);
    int PMT_count = 0;
    for (int i=0;i<nPMT_total;i++)
    {
        if ((PMT_count+0.)/(i+1.)<PMT_frac && PMT_count<nPMT) 
        {
            if (!mPMT) pmt_mask.emplace_back(0);
            else 
            {   // this assumes the PMTs are properly ordered, the whole mPMT is turned on/off
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

void AnaTree::MaskmPMT(std::vector<int> vec, int nPMTpermPMT)
{
    // mask the small PMT in mPMT
    mpmt_mask.clear();
    m_maskmpmt = false;

    mpmt_mask.resize(nPMTpermPMT,0);
    bool mask = false;
    for (auto m : vec) 
    {
        if ( m<0 || m>=nPMTpermPMT )
            std::cout << ERR << "In MaskmPMT(), mPMT_id = " << m << " is invalid, this PMT is not masked" << std::endl;
        else
        {
            std::cout << TAG << "Masking mPMT_id = " << m << std::endl;
            mpmt_mask[m] = 1;
            mask = true;
        }
    }

    if (mask) m_maskmpmt = true;
}

void AnaTree::MaskPMTid(std::vector<int> vec)
{
    // mask PMT with specific id
    bool mask = false;

    if (!m_maskpmt)
    {
        pmt_mask.clear();
        pmt_mask.resize(t_pmt->GetEntries(),0);
    }

    for (auto id : vec) 
    {
        if ( id<0 || id>=pmt_mask.size() )
            std::cout << ERR << "In MaskPMTid(), PMT_id = " << id << " is invalid, this PMT is not masked" << std::endl;
        else
        {
            std::cout << TAG << "Masking PMT_id = " << id << std::endl;
            pmt_mask[id] = 1;
            mask = true;
        }
    }

    if (!m_maskpmt && mask) m_maskpmt = true;

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
    // Set branch addresses and branch pointers for PMT hits

    fChain->SetBranchAddress("nPE", &nPE);
    fChain->SetBranchAddress("timetof", &timetof);
    fChain->SetBranchAddress("PMT_id", &PMT_id);

}

void AnaTree::SetPMTBranches()
{
    // Set branch addresses and branch pointers for PMT geometry

    t_pmt->SetBranchAddress("R", &R);
    t_pmt->SetBranchAddress("costh", &costh);
    t_pmt->SetBranchAddress("costhm", &costhm);
    t_pmt->SetBranchAddress("cosths", &cosths);
    t_pmt->SetBranchAddress("phis", &phis);
    t_pmt->SetBranchAddress("phim", &phim);
    t_pmt->SetBranchAddress("omega", &omega);
    t_pmt->SetBranchAddress("dz", &dz);
    t_pmt->SetBranchAddress("xpos", &xpos);
    t_pmt->SetBranchAddress("ypos", &ypos);
    t_pmt->SetBranchAddress("zpos", &zpos);
    t_pmt->SetBranchAddress("PMT_id", &PMT_id);
    t_pmt->SetBranchAddress("mPMT_id", &mPMT_id);
    
}

std::vector<AnaEvent> AnaTree::GetPMTs()
{
    std::vector<AnaEvent> pmt_vec;

    if(t_pmt == nullptr)
    {
        std::cout << ERR <<"Reading no PMTs in GetPMTs()"<<std::endl;
        return pmt_vec;
    }

    // enforce maskpmt if maskmpmt is true
    if (m_maskmpmt)
        if (!m_maskpmt)
        {
            m_maskpmt = true;
            pmt_mask.clear();
            pmt_mask.resize(t_pmt->GetEntries(),0);
        }

    SetPMTBranches();
    for(long int jentry = 0; jentry < t_pmt->GetEntries(); jentry++)
    {
        t_pmt->GetEntry(jentry);

        if (m_maskpmt) if (pmt_mask.at(PMT_id)) continue;
        // Mask small PMT at run time to avoid possible ordering problem
        if (m_maskmpmt) 
            if (mpmt_mask.at(mPMT_id))
            {
                pmt_mask.at(PMT_id) = 1;
                continue;
            }

        AnaEvent ev(jentry);
        ev.SetR(R);
        ev.SetCosth(costh);
        ev.SetCosthm(costhm);
        ev.SetCosths(cosths);
        ev.SetPhis(phis);
        ev.SetPhim(phim);
        ev.SetOmega(omega); 
        ev.SetDz(dz); 
        ev.SetPos(std::vector<double>{xpos,ypos,zpos}); 
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

    // ignore PMT hits from masked PMT
    if (m_maskpmt) if (pmt_mask.at(PMT_id)) return false;

    time = timetof;

    charge = nPE;

    pmtID = PMT_id;

    return true;
}