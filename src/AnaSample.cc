#include "AnaSample.hh"

AnaSample::AnaSample(int sample_id, const std::string& name, const std::string& binning, int pmt_type)
    : m_sample_id(sample_id)
    , m_norm(1.0)
    , m_pmttype(pmt_type)
    , m_pmtmask(0)
    , m_nPMTpermPMT(19)
    , m_stat_fluc(false)
    , m_name(name)
    , m_binning(binning)
    , m_hpred(nullptr)
    , m_hpred_err2(nullptr)
    , m_hdata(nullptr)
    , m_hdata_control(nullptr)
    , m_hdata_pmt(nullptr)
    , m_hdata_pmt_control(nullptr)
    , m_scatter(false)
    , m_scatter_map(false)
    , m_h_scatter_map(nullptr)
    , m_time_offset(false)
    , m_time_smear(false)
    , m_z0(0.)
    , m_htimetof_data(nullptr)
    , m_htimetof_pred(nullptr)
    , m_htimetof_pred_w2(nullptr)
    , m_htimetof_pmt_pred(nullptr)
    , m_htimetof_pmt_data(nullptr)
    , m_pmt_eff(nullptr)
    , m_use_eff(false)
    , m_eff_var(false)
    , m_eff_sig(0.0)
    , m_template(false)
    , selTree(nullptr)
{
    TH1::SetDefaultSumw2(true);

    std::cout << TAG  << "Sample: " << m_name << " (ID: " << m_sample_id << ")" << std::endl;

    m_bm.SetBinning(m_binning);
    //m_bm.Print();
    m_nbins = m_bm.GetNbins();

    m_llh = new PoissonLLH;

    MakeHistos(); // with default binning

    std::cout << TAG << "MakeHistos called." << std::endl;
}

AnaSample::~AnaSample()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    
    if(m_hpred_err2 != nullptr)
        delete m_hpred_err2;

    if(m_hdata != nullptr)
        delete m_hdata;

    if(m_hdata_control != nullptr)
        delete m_hdata_control;

    if(m_hdata_pmt != nullptr)
        delete m_hdata_pmt;

    if(m_hdata_pmt_control != nullptr)
        delete m_hdata_pmt_control;

    if(m_h_scatter_map != nullptr)
        delete m_h_scatter_map;

    if(m_htimetof_data != nullptr)
        delete m_htimetof_data;
    
    if(m_htimetof_pred != nullptr)
        delete m_htimetof_pred;

    if(m_htimetof_pred_w2 != nullptr)
        delete m_htimetof_pred_w2;

    if(m_htimetof_pmt_data != nullptr)
        delete m_htimetof_pmt_data;
    
    if(m_htimetof_pmt_pred != nullptr)
        delete m_htimetof_pmt_pred;

    if(m_pmt_eff != nullptr)
        delete m_pmt_eff;

}

void AnaSample::LoadEventsFromFile(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name)
{
    selTree = new AnaTree(file_name, tree_name, pmt_tree_name);
    if (m_pmttype==1 && m_mPMTmask.size()>0) selTree->MaskmPMT(m_mPMTmask, m_nPMTpermPMT);
    if (m_pmtmask>0) selTree->MaskPMT(m_pmtmask, m_pmttype, m_nPMTpermPMT);

    std::cout << TAG  << "Reading events for from "<<file_name<<"...\n";

    std::cout << TAG  << "Reading PMT geometry...\n";
    std::vector<AnaEvent> pmt_vec = selTree->GetPMTs();

    m_pmts.clear();
    // Add PMT geometry
    for (long int i=0;i<pmt_vec.size();i++) {
        bool skip = false;
        for (int j=0;j<m_cutvar.size();j++) {
            if (m_cutvar[j]=="nPE"||m_cutvar[j]=="timetof") continue;
            double val = pmt_vec[i].GetEventVar(m_cutvar[j]); // check if PMT passes the cuts
            if (val<m_cutlow[j] || val>m_cuthigh[j]) {
                skip=true;
                break;
            }
        }
        pmt_vec[i].SetSampleType(m_sample_id);
        double eff = 1.0;
        if (m_use_eff) eff = m_pmt_eff->GetBinContent(pmt_vec[i].GetPMTID()+1);
        if (m_eff_var) eff *= gRandom->Gaus(1,m_eff_sig);
        pmt_vec[i].SetEff(eff);

        pmt_vec[i].SetZ0(m_z0);
        if (!skip) AddPMT(pmt_vec[i]);
    }

    double timetof, nPE;
    selTree->SetDataBranches();

    unsigned long nDataEntries = selTree->GetDataEntries();

    if(m_hdata_pmt != nullptr)
        delete m_hdata_pmt;

    int nPMTs = selTree->GetPMTEntries();
    m_hdata_pmt = new TH1D("","",nPMTs,0,nPMTs);

    // determine the random offset for each PMT
    std::vector<double> timetof_shift;
    if (m_time_offset)
    {
        for (int i=0;i<nPMTs;i++)
            timetof_shift.push_back( gRandom->Gaus(0,m_time_offset_width) );
    }

    // determine the random smearing for each PMT
    std::vector<double> time_resolution;
    if (m_time_smear)
    {
        for (int i=0;i<nPMTs;i++)
        {
            double resol = -1;
            while (resol<0)
                resol = gRandom->Gaus(m_time_smear_mean,m_time_smear_width);
            time_resolution.push_back(resol);   
        }
    }

    // Setup control region if the scattering option is on
    if (m_scatter || m_scatter_map)
    {
        if(m_hdata_pmt_control != nullptr)
            delete m_hdata_pmt_control;
        m_hdata_pmt_control = new TH1D("","",nPMTs,0,nPMTs);
    }

    int pmtID;

    std::cout << TAG<<"Reading PMT hit data..."<<std::endl;
    for (unsigned long i=0;i<nDataEntries;i++)
    {
        if (!selTree->GetDataEntry(i,timetof,nPE,pmtID)) continue;

        if (m_time_offset) timetof += timetof_shift[pmtID];
        if (m_time_smear) timetof += gRandom->Gaus(0,time_resolution[pmtID]);

        bool skip = false;

        for (int j=0;j<m_cutvar.size();j++) {
            if (m_cutvar[j]=="timetof")
            {
                if (timetof<m_cutlow[j] || timetof>m_cuthigh[j]) {
                    skip = true;   
                    break;
                }
            }
            else if (m_cutvar[j]=="nPE")
            {
                if (nPE<m_cutlow[j] || nPE>m_cuthigh[j]) {
                    skip = true;   
                    break;
                }
            }
        }

        if (m_template)
        {
            m_htimetof_pmt_data->Fill(pmtID+0.5,timetof+m_timetof_offset,nPE);
        }

        if (skip) continue;
        if (m_scatter || m_scatter_map)
        {
            if (timetof>=m_scatter_time1 && timetof<m_scatter_time2) m_hdata_pmt->Fill(pmtID+0.5, nPE);
            else if (timetof>=m_scatter_time2 && timetof<m_scatter_time3) m_hdata_pmt_control->Fill(pmtID+0.5, nPE);
        }
        else m_hdata_pmt->Fill(pmtID+0.5, nPE);
    }

    PrintStats();
}

AnaEvent* AnaSample::GetPMT(const unsigned int evnum)
{
#ifndef NDEBUG
    if(m_pmts.empty())
    {
        std::cerr << ERR << " In AnaSample::GetPMT()" << std::endl;
        std::cerr << ERR << " No PMTs are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_pmts.size())
    {
        std::cerr << ERR << " In AnaSample::GetPMT()" << std::endl;
        std::cerr << ERR << " PMT number out of bounds in " << m_name << " sample." << std::endl;
        return nullptr;
    }
#endif

    return &m_pmts[evnum];
}

void AnaSample::PrintStats() const
{
    const double mem_kb = sizeof(m_pmts) * m_pmts.size() / 1000.0;
    std::cout << TAG << "Sample " << m_name << " ID = " << m_sample_id << std::endl;
    std::cout << TAG << "Num of PMTs   = " << m_pmts.size() << std::endl;
    std::cout << TAG << "Memory used   = " << mem_kb << " kB." << std::endl;
}

void AnaSample::MakeHistos()
{
    if(m_hpred != nullptr)
        delete m_hpred;
    m_hpred = new TH1D(Form("%s_pred", m_name.c_str()), Form("%s_pred", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred->SetDirectory(0);

    if(m_hpred_err2 != nullptr)
        delete m_hpred_err2;
    m_hpred_err2 = new TH1D(Form("%s_pred_err2", m_name.c_str()), Form("%s_pred_err2", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred_err2->SetDirectory(0);

    if(m_hdata != nullptr)
        delete m_hdata;
    m_hdata = new TH1D(Form("%s_data", m_name.c_str()), Form("%s_data", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hdata->SetDirectory(0);

    if(m_hdata_control != nullptr)
        delete m_hdata_control;
    m_hdata_control = new TH1D(Form("%s_data_control", m_name.c_str()), Form("%s_data_tail", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hdata_control->SetDirectory(0);
}

void AnaSample::InitEventMap()
{
    if (m_binvar[0]=="unbinned")
    {
        m_nbins = m_pmts.size();
        MakeHistos();
        std::cout << TAG<<"Using unbinned histogram for fit"<<std::endl;

        int count = 0;
        for(auto& e : m_pmts)
        {
            e.SetSampleBin(count);
            count++;
        }
    }
    else for(auto& e : m_pmts)
    {
        std::vector<double> binvar;
        for (auto t : m_binvar)
            binvar.push_back(e.GetEventVar(t));
        const int b = m_bm.GetBinIndex(binvar);
#ifndef NDEBUG
        if(b < 0)
        {
            std::cout << TAG << "In InitEventMap()\n"
                      << "No bin for current PMT." << std::endl;
            std::cout << TAG << "PMT Var: " << std::endl;
            for(const auto val : e.GetRecoVar())
                std::cout << TAG << "\t" << val << std::endl;
        }
#endif
        e.SetSampleBin(b);
    }

    if (m_template)
    {
        for(auto& e : m_pmts)
        {
            const int pmtID = e.GetPMTID();
            std::vector<double> timetof_nom;
            std::vector<double> timetof_nom_sig2;
            for (int i=1;i<=m_htimetof_pmt_pred->GetNbinsY();i++)
            {
                timetof_nom.push_back(m_htimetof_pmt_pred->GetBinContent(pmtID+1,i));
                double sig2 = m_htimetof_pmt_pred->GetBinError(pmtID+1,i);
                if (sig2>0) sig2 = sig2*sig2/timetof_nom[i-1]/timetof_nom[i-1];
                timetof_nom_sig2.push_back(sig2);
            }
            e.SetTimetofNom(timetof_nom);
            e.SetTimetofNomSig2(timetof_nom_sig2);
            e.SetTimetofPred(timetof_nom);
        }

        int nx = m_template_combine ? 1 : m_nbins ;
        int ny = m_htimetof_pmt_pred->GetNbinsY();
        if(m_htimetof_pred != nullptr)
            delete m_htimetof_pred;
        if(m_htimetof_pred_w2 != nullptr)
            delete m_htimetof_pred_w2;
        if(m_htimetof_data != nullptr)
            delete m_htimetof_data;
        m_htimetof_pred = new TH2D(Form("%s_timetof_pred", m_name.c_str()), Form("%s_timetof_pred", m_name.c_str()), nx, 0, nx, ny, 0, ny );
        m_htimetof_pred->SetDirectory(0);
        m_htimetof_pred_w2 = new TH2D(Form("%s_timetof_pred_w2", m_name.c_str()), Form("%s_timetof_pred_w2", m_name.c_str()), nx, 0, nx, ny, 0, ny );
        m_htimetof_pred_w2->SetDirectory(0);
        m_htimetof_data = new TH2D(Form("%s_timetof_data", m_name.c_str()), Form("%s_timetof_data", m_name.c_str()), nx, 0, nx, ny, 0, ny );
        m_htimetof_data->SetDirectory(0);
    }
}

void AnaSample::FillEventHist(bool reset_weights)
{
#ifndef NDEBUG
    if(m_hpred == nullptr)
    {
        std::cerr << ERR << "In FillEventHist() h_pred is a nullptr!"
                  << "Returning from function." << std::endl;
        return;
    }
#endif
    m_hpred->Reset();
    if (m_scatter||m_scatter_map) m_hpred_err2->Reset();

    if (m_template) 
    {
        m_htimetof_pred->Reset();
        m_htimetof_pred_w2->Reset();
    }

    for(const auto& e : m_pmts)
    {
        const double weight = reset_weights ? e.GetEvWghtMC() : e.GetEvWght();
        const int reco_bin  = e.GetSampleBin();
        m_hpred->Fill(reco_bin + 0.5, weight); // direct PE prediction

        if (m_scatter||m_scatter_map)
        {
            m_hpred->Fill(reco_bin + 0.5, e.GetPEIndirect()); // indirect PE prediction
            m_hpred_err2->Fill(reco_bin + 0.5, e.GetPEIndirectErr()); // MC stat error of indirect PE prediction
        }

        if (m_template)
        {
            std::vector<double> timetof_pred = e.GetTimetofPred();
            std::vector<double> timetof_nom_sig2 = e.GetTimetofNomSig2();
            for (int i=1;i<=m_htimetof_pred->GetNbinsY();i++)
            {
                if (!m_template_combine) 
                {
                    m_htimetof_pred->Fill(reco_bin + 0.5,i-0.5,timetof_pred[i-1]);
                    m_htimetof_pred_w2->Fill(reco_bin + 0.5,i-0.5,timetof_nom_sig2[i-1]*timetof_pred[i-1]*timetof_pred[i-1]);
                }
                else
                {
                    m_htimetof_pred->Fill(0.5,i-0.5,timetof_pred[i-1]);
                    m_htimetof_pred_w2->Fill(0.5,i-0.5,timetof_nom_sig2[i-1]*timetof_pred[i-1]*timetof_pred[i-1]);
                }
            }
        }
    }

    return;
}

void AnaSample::FillDataHist(bool stat_fluc)
{
#ifndef NDEBUG
    if(m_hdata == nullptr)
    {
        std::cerr << ERR << "In FillDataHist() m_hdata is a nullptr!"
                  << "Returning from function." << std::endl;
        return;
    }
#endif
    m_hdata->Reset();
    if (m_scatter || m_scatter_map) m_hdata_control->Reset();

    if(stat_fluc) 
        std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

    for(auto& e : m_pmts)
    {
        const int pmtID = e.GetPMTID();
        double weight = m_hdata_pmt->GetBinContent(pmtID+1)*m_norm;

        if (stat_fluc)
        {
            weight = gRandom->Poisson(weight);
#ifndef NDEBUG
            if(weight <= 0.0)
            {
                std::cout << ERR << "In FillDataHist()\n"
                            << "In Sample " <<  m_name << ", PMT " << pmtID
                            << " has 0 (or negative) entries. This may cause a problem with chi2 computations."
                            << std::endl;
            }
#endif
        }

        e.SetPE(weight);
        const int reco_bin  = e.GetSampleBin();
        m_hdata->Fill(reco_bin + 0.5, weight);

        if (m_scatter || m_scatter_map) 
        {
            double weight_control = 0;

            if (m_hdata_pmt_control->GetBinContent(pmtID+1)>0)
            {
                weight_control = m_hdata_pmt_control->GetBinContent(pmtID+1)*m_norm;
                if (stat_fluc)
                {
                    weight_control = gRandom->Poisson(weight_control);
#ifndef NDEBUG
                    if(weight_control <= 0.0)
                    {
                        std::cout << ERR << "In FillDataHist()\n"
                                    << "In Sample " <<  m_name << ", PMT " << pmtID
                                    << " has 0 (or negative) entries. This may cause a problem with chi2 computations."
                                    << std::endl;
                    }
#endif
                }
                double indirect_pred = weight_control;
                // m_scatter_factor: constant scale factor for indirect PE prediction
                // m_scatter_map: PMT specific scale factor for indirect PE prediction
                indirect_pred *= m_scatter ? m_scatter_factor : m_h_scatter_map->GetBinContent(pmtID+1);
                // stat error from measured PE in control region
                double indirect_err2 = 1./weight_control;
                // MC stat error from scale factor estimation
                if (m_scatter_map)
                    indirect_err2 += m_h_scatter_map->GetBinError(pmtID+1)*m_h_scatter_map->GetBinError(pmtID+1);
                indirect_err2 *= indirect_pred*indirect_pred;

                e.SetPEIndirect(indirect_pred);
                e.SetPEIndirectErr(indirect_err2);
            }

            m_hdata_control->Fill(reco_bin + 0.5, weight_control);
        }

        if (m_template)
        {
            for (int i=1;i<=m_htimetof_data->GetNbinsY();i++)
            {
                if (!m_template_combine) 
                {
                    m_htimetof_data->Fill(reco_bin + 0.5,i-0.5,m_htimetof_pmt_data->GetBinContent(pmtID+1,i));
                }
                else
                {
                    m_htimetof_data->Fill(0.5,i-0.5,m_htimetof_pmt_data->GetBinContent(pmtID+1,i));
                }
            }
        }
    }

    return;
}

void AnaSample::SetLLHFunction(const std::string& func_name)
{
    if(m_llh != nullptr)
        delete m_llh;

    if(func_name.empty())
    {
        std::cout << TAG << "Likelihood function name empty. Setting to Poisson by default." << std::endl;
        m_llh = new PoissonLLH;
    }
    else if(func_name == "Poisson")
    {
        std::cout << TAG << "Setting likelihood function to Poisson." << std::endl;
        m_llh = new PoissonLLH;
    }
    else if(func_name == "Effective")
    {
        std::cout << TAG << "Setting likelihood function to Tianlu's effective likelihood." << std::endl;
        m_llh = new EffLLH;
    }
    else if(func_name == "Barlow")
    {
        std::cout << TAG << "Setting likelihood function to Barlow-Beeston." << std::endl;
        m_llh = new BarlowLLH;
    }
    else
    {
        std::cout << TAG << "Likelihood function not defined. Setting to Poisson by default." << std::endl;
        m_llh = new PoissonLLH;
    }
}

double AnaSample::CalcLLH() const
{
    const unsigned int nbins = m_hpred->GetNbinsX();
    double* exp_w  = m_hpred->GetArray();
    //double* exp_w2 = m_hpred->GetSumw2()->GetArray();
    double* exp_w2 = m_hpred_err2->GetArray();
    double* data   = m_hdata->GetArray();

    double chi2 = 0.0;
    for(unsigned int i = 1; i <= nbins; ++i)
        chi2 += (*m_llh)(exp_w[i], exp_w2[i], data[i]);

    if (m_template)
    {
        if (m_template_only) chi2 = 0.0;
        for(unsigned int i=1;i<=m_htimetof_pred->GetNbinsX();i++)
        {
            for (int j=1;j<=m_htimetof_pred->GetNbinsY();j++)
            {
                chi2 += (*m_llh)(m_htimetof_pred->GetBinContent(i,j), m_htimetof_pred_w2->GetBinContent(i,j), m_htimetof_data->GetBinContent(i,j));
            }
        }
    }

    return chi2;
}

void AnaSample::WriteEventHist(TDirectory* dirout, const std::string& bsname)
{
    dirout->cd();
    if(m_hpred != nullptr)
        m_hpred->Write(Form("evhist_sam%d_pred%s", m_sample_id, bsname.c_str()));

    if (m_scatter || m_scatter_map)
    {
        if (m_scatter_map)
            if (m_hpred_err2 != nullptr)
                m_hpred_err2->Write(Form("evhist_err2_sam%d_pred%s", m_sample_id, bsname.c_str()));
    }

    if (m_template)
    {
        if (m_htimetof_pred != nullptr)
            m_htimetof_pred->Write(Form("evhist_timetof_sam%d_pred%s", m_sample_id, bsname.c_str()));
    }
}

void AnaSample::WriteDataHist(TDirectory* dirout, const std::string& bsname)
{
    dirout->cd();
    if(m_hdata != nullptr)
        m_hdata->Write(Form("evhist_sam%d_data%s", m_sample_id, bsname.c_str()));

    if (m_scatter || m_scatter_map)
        if(m_hdata_control != nullptr)
            m_hdata_control->Write(Form("evhist_control_sam%d_data%s", m_sample_id, bsname.c_str())); 

    if (m_template)
    {
        if (m_htimetof_data != nullptr)
            m_htimetof_data->Write(Form("evhist_timetof_sam%d_data%s", m_sample_id, bsname.c_str()));
    }        
}

void AnaSample::SetScatterMap(double time1, double time2, double time3, const TH1D& hist)
{
    // read the pre-calculated scale factor for indirect PE estimation
    m_scatter_map = true; m_scatter_time1 = time1; m_scatter_time2 = time2; m_scatter_time3 = time3;
    if(m_h_scatter_map != nullptr)
        delete m_h_scatter_map;
    m_h_scatter_map = new TH1D(hist);
    m_h_scatter_map->SetDirectory(0);
}

void AnaSample::SetTemplate(const TH2D& hist, double offset, bool combine, bool template_only)
{
    m_template = true;
    if(m_htimetof_pmt_pred != nullptr)
        delete m_htimetof_pmt_pred;
    if(m_htimetof_pmt_data != nullptr)
        delete m_htimetof_pmt_data;
    m_htimetof_pmt_pred = new TH2D(hist);
    m_htimetof_pmt_pred->SetDirectory(0);
    m_htimetof_pmt_data = (TH2D*)m_htimetof_pmt_pred->Clone();
    m_htimetof_pmt_data->Reset();
    m_htimetof_pmt_data->SetDirectory(0);

    m_timetof_offset = offset;
    m_template_combine = combine;
    m_template_only = template_only;
}

void AnaSample::SetPMTEff(const TH1D& hist)
{
    m_use_eff = true;

    if(m_pmt_eff != nullptr)
        delete m_pmt_eff;
    m_pmt_eff = new TH1D(hist);
    m_pmt_eff->SetDirectory(0);
}

void AnaSample::InitToy()
{
    for(auto& e : m_pmts)
    {
        double eff = 1.0;
        if (m_use_eff) eff = m_pmt_eff->GetBinContent(e.GetPMTID()+1);
        if (m_eff_var) eff *= gRandom->Gaus(1,m_eff_sig);
        e.SetEff(eff);
    }

    if ( m_time_offset || m_time_smear ) 
    {
        int nPMTs = selTree->GetPMTEntries();
        // determine the random offset for each PMT
        std::vector<double> timetof_shift;
        if (m_time_offset)
        {
            for (int i=0;i<nPMTs;i++)
                timetof_shift.push_back( gRandom->Gaus(0,m_time_offset_width) );
        }

        // determine the random smearing for each PMT
        std::vector<double> time_resolution;
        if (m_time_smear)
        {
            for (int i=0;i<nPMTs;i++)
            {
                double resol = -1;
                while (resol<0)
                    resol = gRandom->Gaus(m_time_smear_mean,m_time_smear_width);
                time_resolution.push_back(resol);   
            }
        }

        unsigned long nDataEntries = selTree->GetDataEntries();
        double timetof, nPE;
        int pmtID;
        if (m_template) m_htimetof_pmt_data->Reset();
        m_hdata_pmt->Reset();
        m_hdata_pmt_control->Reset();
        for (unsigned long i=0;i<nDataEntries;i++)
        {
            if (!selTree->GetDataEntry(i,timetof,nPE,pmtID)) continue;

            if (m_time_offset) timetof += timetof_shift[pmtID];
            if (m_time_smear) timetof += gRandom->Gaus(0,time_resolution[pmtID]);

            bool skip = false;

            for (int j=0;j<m_cutvar.size();j++) {
                if (m_cutvar[j]=="timetof")
                {
                    if (timetof<m_cutlow[j] || timetof>m_cuthigh[j]) {
                        skip = true;   
                        break;
                    }
                }
                else if (m_cutvar[j]=="nPE")
                {
                    if (nPE<m_cutlow[j] || nPE>m_cuthigh[j]) {
                        skip = true;   
                        break;
                    }
                }
            }

            if (m_template)
            {
                m_htimetof_pmt_data->Fill(pmtID+0.5,timetof+m_timetof_offset,nPE);
            }

            if (skip) continue;
            if (m_scatter || m_scatter_map)
            {
                if (timetof>=m_scatter_time1 && timetof<m_scatter_time2) m_hdata_pmt->Fill(pmtID+0.5, nPE);
                else if (timetof>=m_scatter_time2 && timetof<m_scatter_time3) m_hdata_pmt_control->Fill(pmtID+0.5, nPE);
            }
            else m_hdata_pmt->Fill(pmtID+0.5, nPE);
        }
    }
}