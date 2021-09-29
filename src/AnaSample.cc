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
    , m_hpred_indirect(nullptr)
    , m_hdata(nullptr)
    , m_hdata_control(nullptr)
    , m_hdata_unbinned(nullptr)
    , m_hdata_unbinned_control(nullptr)
    , m_scatter(false)
    , m_scatter_map(false)
    , m_h_scatter_map(nullptr)
    , m_time_offset(false)
    , m_time_smear(false)
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

    if(m_hpred_indirect != nullptr)
        delete m_hpred_indirect;

    if(m_hdata != nullptr)
        delete m_hdata;

    if(m_hdata_control != nullptr)
        delete m_hdata_control;

    if(m_hdata_unbinned != nullptr)
        delete m_hdata_unbinned;

    if(m_hdata_unbinned_control != nullptr)
        delete m_hdata_unbinned_control;

    if(m_h_scatter_map != nullptr)
        delete m_h_scatter_map;

}

void AnaSample::LoadEventsFromFile(const std::string& file_name, const std::string& tree_name, const std::string& pmt_tree_name)
{
    AnaTree selTree(file_name, tree_name, pmt_tree_name);
    if (m_pmttype==1 && m_mPMTmask.size()>0) selTree.MaskmPMT(m_mPMTmask, m_nPMTpermPMT);
    if (m_pmtmask>0) selTree.MaskPMT(m_pmtmask, m_pmttype, m_nPMTpermPMT);

    std::cout << TAG  << "Reading events for from "<<file_name<<"...\n";

    std::cout << TAG  << "Reading PMT geometry...\n";
    std::vector<AnaEvent> pmt_vec = selTree.GetPMTs();

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
        if (!skip) AddPMT(pmt_vec[i]);
    }

    double timetof, nPE;
    selTree.SetDataBranches();

    unsigned long nDataEntries = selTree.GetDataEntries();

    if(m_hdata_unbinned != nullptr)
        delete m_hdata_unbinned;

    int nPMTs = selTree.GetPMTEntries();
    m_hdata_unbinned = new TH1D("","",nPMTs,0,nPMTs);

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
        if(m_hdata_unbinned_control != nullptr)
            delete m_hdata_unbinned_control;
        m_hdata_unbinned_control = new TH1D("","",nPMTs,0,nPMTs);
    }

    int pmtID;

    std::cout << TAG<<"Reading PMT hit data..."<<std::endl;
    for (unsigned long i=0;i<nDataEntries;i++)
    {
        if (!selTree.GetDataEntry(i,timetof,nPE,pmtID)) continue;

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

        if (skip) continue;
        if (m_scatter || m_scatter_map)
        {
            if (timetof>=m_scatter_time1 && timetof<m_scatter_time2) m_hdata_unbinned->Fill(pmtID+0.5, nPE);
            else if (timetof>=m_scatter_time2 && timetof<m_scatter_time3) m_hdata_unbinned_control->Fill(pmtID+0.5, nPE);
        }
        else m_hdata_unbinned->Fill(pmtID+0.5, nPE);
    }

    PrintStats();
}

AnaEvent* AnaSample::GetEvent(const unsigned int evnum) // not used anymore
{
#ifndef NDEBUG
    if(m_events.empty())
    {
        std::cerr << " In AnaSample::GetEvent()" << std::endl;
        std::cerr << " No events are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_events.size())
    {
        std::cerr << " In AnaSample::GetEvent()" << std::endl;
        std::cerr << " Event number out of bounds in " << m_name << " sample." << std::endl;
        return nullptr;
    }
#endif

    return &m_events[evnum];
}

AnaEvent* AnaSample::GetPMT(const unsigned int evnum)
{
#ifndef NDEBUG
    if(m_pmts.empty())
    {
        std::cerr << " In AnaSample::GetPMT()" << std::endl;
        std::cerr << " No PMTs are found in " << m_name << " sample." << std::endl;
        return nullptr;
    }
    else if(evnum >= m_pmts.size())
    {
        std::cerr << " In AnaSample::GetPMT()" << std::endl;
        std::cerr << " PMT number out of bounds in " << m_name << " sample." << std::endl;
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

    if(m_hpred_indirect != nullptr)
        delete m_hpred_indirect;
    m_hpred_indirect = new TH1D(Form("%s_pred_indirect", m_name.c_str()), Form("%s_pred_tail", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hpred_indirect->SetDirectory(0);

    if(m_hdata != nullptr)
        delete m_hdata;
    m_hdata = new TH1D(Form("%s_data", m_name.c_str()), Form("%s_data", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hdata->SetDirectory(0);

    if(m_hdata_control != nullptr)
        delete m_hdata_control;
    m_hdata_control = new TH1D(Form("%s_data_control", m_name.c_str()), Form("%s_data_tail", m_name.c_str()), m_nbins, 0, m_nbins);
    m_hdata_control->SetDirectory(0);

    // if(m_h_scatter_map != nullptr)
    //     delete m_h_scatter_map;
    // m_h_scatter_map = new TH1D(Form("%s_scatter_map", m_name.c_str()), Form("%s_scatter_map", m_name.c_str()), m_nbins, 0, m_nbins);
    // m_h_scatter_map->SetDirectory(0);
}

void AnaSample::InitEventMap()
{
//     for(auto& e : m_events)
//     {
//         std::vector<double> binvar;
//         for (auto t : m_binvar)
//             binvar.push_back(e.GetEventVar(t));
//         const int b = m_bm.GetBinIndex(binvar);
// #ifndef NDEBUG
//         if(b < 0)
//         {
//             std::cout << TAG << "In AnaSample::InitEventMap()\n"
//                       << "No bin for current event." << std::endl;
//             std::cout << TAG << "Event kinematics: " << std::endl;
//             for(const auto val : e.GetRecoVar())
//                 std::cout << TAG << "\t" << val << std::endl;
//         }
// #endif
//         e.SetSampleBin(b);
//     }

//     std::sort(m_events.begin(), m_events.end(), [](const AnaEvent a, const AnaEvent b){ return a.GetSampleBin() < b.GetSampleBin(); });

    if (m_binvar[0]=="unbinned")
    {
        m_nbins = m_hdata_unbinned->GetNbinsX();
        MakeHistos();
        std::cout << TAG<<"Using unbinned histogram for fit"<<std::endl;

        for(auto& e : m_pmts)
        {
            const int b = e.GetPMTID();
            e.SetSampleBin(b);
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
            std::cout << TAG << "In AnaSample::InitEventMap()\n"
                      << "No bin for current PMT." << std::endl;
            std::cout << TAG << "PMT Var: " << std::endl;
            for(const auto val : e.GetRecoVar())
                std::cout << TAG << "\t" << val << std::endl;
        }
#endif
        e.SetSampleBin(b);
    }
}

void AnaSample::FillEventHist(bool reset_weights)
{
#ifndef NDEBUG
    if(m_hpred == nullptr)
    {
        std::cerr << "In AnaSample::FillEventHist() h_pred is a nullptr!"
                  << "Returning from function." << std::endl;
        return;
    }
#endif
    m_hpred->Reset();
    //if (m_scatter) m_hpred_indirect->Reset();
    if (m_scatter||m_scatter_map) m_hpred_err2->Reset();

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

        // if (m_scatter)
        // {
        //     // e.GetTailPE() gives the (fractional) number of PE in the scatter control region
        //     // m_scatter_factor is the scale factor from tail scatter pe to peak scatter pe
        //     const double scatter_pe_at_peak = weight*e.GetTailPE()*m_scatter_factor;
        //     m_hpred->Fill(reco_bin + 0.5, scatter_pe_at_peak);
        //     const double tailpe = (weight+scatter_pe_at_peak)*e.GetTailPE();
        //     m_hpred_indirect->Fill(reco_bin + 0.5, tailpe);
        // }
    }

    return;
}

void AnaSample::FillDataHist(bool stat_fluc)
{
#ifndef NDEBUG
    if(m_hdata == nullptr)
    {
        std::cerr << "In AnaSample::FillDataHist() m_hdata is a nullptr!"
                  << "Returning from function." << std::endl;
        return;
    }
#endif
    m_hdata->Reset();
    if (m_scatter || m_scatter_map) m_hdata_control->Reset();

    for(auto& e : m_pmts)
    {
        const int pmtID = e.GetPMTID();
        double weight = m_hdata_unbinned->GetBinContent(pmtID+1);
        if (m_scatter) 
        {
            double val = m_hdata_unbinned_control->GetBinContent(pmtID+1)*m_scatter_factor; // constant scale factor for indirect PE prediction
            double err = 1./m_hdata_unbinned_control->GetBinContent(pmtID+1)*val*val; // stat error from measured PE in control region
            e.SetPEIndirect(val);
            e.SetPEIndirectErr(err);
            // weight -= m_hdata_unbinned_control->GetBinContent(pmtID+1)*m_scatter_factor;
            // if (weight<0) weight=0;
        }
        else if (m_scatter_map)
        {
            //weight -= m_hdata_unbinned_control->GetBinContent(pmtID+1)*m_h_scatter_map->GetBinContent(pmtID+1);
            //std::cout << TAG<<pmtID<<" "<<m_hdata_unbinned->GetBinContent(pmtID+1)<<" "<<m_hdata_unbinned_control->GetBinContent(pmtID+1)<<" "<<m_h_scatter_map->GetBinContent(pmtID+1)<<std::endl;
            //if (weight<0) weight=0;
            double val = m_hdata_unbinned_control->GetBinContent(pmtID+1)*m_h_scatter_map->GetBinContent(pmtID+1); // predefined scale factor for indirect PE prediction
            double err = m_h_scatter_map->GetBinError(pmtID+1); // error from scale factor estimation
            err = (err*err + 1./m_hdata_unbinned_control->GetBinContent(pmtID+1))*val*val; 
            // double val = m_h_scatter_map->GetBinContent(pmtID+1);
            // double err = m_h_scatter_map->GetBinError(pmtID+1);
            // err = err*err*val*val;
            e.SetPEIndirect(val);
            e.SetPEIndirectErr(err);
        }
        e.SetPE(weight);
        const int reco_bin  = e.GetSampleBin();
        m_hdata->Fill(reco_bin + 0.5, weight);
        if (m_scatter || m_scatter_map) m_hdata_control->Fill(reco_bin + 0.5, m_hdata_unbinned_control->GetBinContent(pmtID+1));
    }

    // if (m_scatter)
    // {
    //     for (int j = 1; j <= m_hdata->GetNbinsX(); ++j)
    //     {
    //         double val = m_hdata->GetBinContent(j);
    //         if (val<=0) continue;
    //         double ratio = m_hdata_control->GetBinContent(j)/val;
    //         double correction = 1.-ratio*m_scatter_factor;
    //         m_hdata->SetBinContent(j, val*correction);
    //     }
    // }


    m_hdata->Scale(m_norm);
    if (m_scatter || m_scatter_map) m_hdata_control->Scale(m_norm);

    if(stat_fluc) 
    {
        std::cout << TAG << "Applying statistical fluctuations..." << std::endl;

        for(int j = 1; j <= m_hdata->GetNbinsX(); ++j)
        {
            double val = m_hdata->GetBinContent(j);
            val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val <= 0.0)
            {
                std::cout << TAG   << "In AnaSample::FillDataHist()\n"
                            << "In Sample " <<  m_name << ", bin " << j
                            << " has 0 (or negative) entries. This may cause a problem with chi2 computations."
                            << std::endl;
            }
#endif
            m_hdata->SetBinContent(j, val);
        }

        if (m_scatter || m_scatter_map) for(int j = 1; j <= m_hdata_control->GetNbinsX(); ++j) // for now ignore the inconsistency of the indirect PE prediction
        {
            double val = m_hdata_control->GetBinContent(j);
            val = gRandom->Poisson(val);
#ifndef NDEBUG
            if(val <= 0.0)
            {
                std::cout << TAG   << "In AnaSample::FillDataHist()\n"
                            << "In Sample " <<  m_name << ", bin " << j
                            << " has 0 (or negative) entries at tail. This may cause a problem with chi2 computations."
                            << std::endl;
            }
#endif
            m_hdata_control->SetBinContent(j, val);
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

    // if (m_scatter)
    // {
    //     exp_w  = m_hpred_indirect->GetArray();
    //     exp_w2 = m_hpred_indirect->GetSumw2()->GetArray();
    //     data   = m_hdata_control->GetArray();
    //     for(unsigned int i = 1; i <= nbins; ++i)
    //         chi2 += (*m_llh)(exp_w[i], exp_w2[i], data[i]);
    // }

    return chi2;
}

double AnaSample::CalcChi2() const // deprecated
{
    if(m_hdata == nullptr)
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Need to define data histogram." << std::endl;
        return 0.0;
    }

    int nbins = m_hpred->GetNbinsX();
    if(nbins != m_hdata->GetNbinsX())
    {
        std::cerr << "[ERROR]: In AnaSample::CalcChi2()\n"
                  << "[ERROR]: Binning mismatch between data and mc.\n"
                  << "[ERROR]: MC bins: " << nbins << ", Data bins: " << m_hdata->GetNbinsX()
                  << std::endl;
        return 0.0;
    }

    double chi2 = 0.0;
    for(int j = 1; j <= nbins; ++j)
    {
        double obs = m_hdata->GetBinContent(j);
        double exp = m_hpred->GetBinContent(j);
        if(exp > 0.0)
        {
            // added when external fake datasets (you cannot reweight when simply 0)
            // this didn't happen when all from same MC since if exp=0 then obs =0

            chi2 += 2 * (exp - obs);
            if(obs > 0.0)
                chi2 += 2 * obs * TMath::Log(obs / exp);

            if(chi2 < 0.0)
            {
#ifndef NDEBUG
                std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                          << "[WARNING]: Stat chi2 is less than 0: " << chi2 << ", setting to 0."
                          << std::endl;
                std::cerr << "[WARNING]: exp and obs is: " << exp << " and " << obs << "."
                          << std::endl;
#endif
                chi2 = 0.0;
            }
        }
    }

    if(chi2 != chi2)
    {
        std::cerr << "[WARNING]: In AnaSample::CalcChi2()\n"
                  << "[WARNING]: Stat chi2 is nan, setting to 0." << std::endl;
        chi2 = 0.0;
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
        if(m_hpred_indirect != nullptr)
            m_hpred_indirect->Write(Form("evhist_indirect_sam%d_pred%s", m_sample_id, bsname.c_str()));

        if (m_scatter_map)
            if (m_hpred_err2 != nullptr)
                m_hpred_err2->Write(Form("evhist_err2_sam%d_pred%s", m_sample_id, bsname.c_str()));
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
}

void AnaSample::SetScatterMap(double time1, double time2, double time3, const TH1D& hist)
{
    // read the prefined scale factor for indirect PE estimation
    m_scatter_map = true; m_scatter_time1 = time1; m_scatter_time2 = time2; m_scatter_time3 = time3;
    if(m_h_scatter_map != nullptr)
        delete m_h_scatter_map;
    m_h_scatter_map = new TH1D(hist);
    m_h_scatter_map->SetDirectory(0);
}