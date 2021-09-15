#include "Fitter.hh"

Fitter::Fitter(TDirectory* dirout, const int seed, const int num_threads)
    : rng(new TRandom3(seed))
    , m_fitter(nullptr)
    , m_fcn(nullptr)
    , m_dir(dirout)
    , m_save(false)
    , m_save_events(true)
    , m_zerosyst(false)
    , m_freq(10000)
    , m_threads(num_threads)
    , m_npar(0)
    , m_calls(0)
{
    gRandom = rng;

    min_settings.minimizer = "Minuit2";
    min_settings.algorithm = "Migrad";
    min_settings.print_level = 2;
    min_settings.strategy  = 1;
    min_settings.tolerance = 1E-2;
    min_settings.max_iter  = 1E6;
    min_settings.max_fcn   = 1E9;
}

Fitter::Fitter(TDirectory* dirout, const int seed)
    : Fitter(dirout, seed, 1)
{
}

Fitter::~Fitter()
{
    m_dir = nullptr;
    if(rng != nullptr)
        delete rng;
    if(m_fitter != nullptr)
        delete m_fitter;
    if(m_fcn != nullptr)
        delete m_fcn;
}

void Fitter::SetSeed(int seed)
{
    if(rng == nullptr)
    {
        rng     = new TRandom3(seed);
        gRandom = rng; // Global pointer
    }
    else
        rng->SetSeed(seed);
}

void Fitter::FixParameter(const std::string& par_name, const double& value)
{
    auto iter = std::find(par_names.begin(), par_names.end(), par_name);
    if(iter != par_names.end())
    {
        const int i = std::distance(par_names.begin(), iter);
        m_fitter->SetVariable(i, par_names.at(i).c_str(), value, 0);
        m_fitter->FixVariable(i);
        std::cout << "Fixing parameter " << par_names.at(i) << " to value " << value
                  << std::endl;
    }
    else
    {
        std::cerr  << "In function Fitter::FixParameter()\n"
                   << "Parameter " << par_name << " not found!" << std::endl;
    }
}

void Fitter::SetMinSettings(const MinSettings& ms)
{
    min_settings = ms;
    if(m_fitter != nullptr)
    {
        m_fitter->SetStrategy(min_settings.strategy);
        m_fitter->SetPrintLevel(min_settings.print_level);
        m_fitter->SetTolerance(min_settings.tolerance);
        m_fitter->SetMaxIterations(min_settings.max_iter);
        m_fitter->SetMaxFunctionCalls(min_settings.max_fcn);
    }
}

void Fitter::InitFitter(std::vector<AnaFitParameters*>& fitpara)
{
    m_fitpara = fitpara;
    std::vector<double> par_step, par_low, par_high;
    std::vector<bool> par_fixed;

    //TRandom3 rng(0);
    for(std::size_t i = 0; i < m_fitpara.size(); i++)
    {
        m_npar += m_fitpara[i]->GetNpar();

        std::vector<std::string> vec0;
        m_fitpara[i]->GetParNames(vec0);
        par_names.insert(par_names.end(), vec0.begin(), vec0.end());

        std::vector<double> vec1, vec2;
        m_fitpara[i]->GetParPriors(vec1);
        if(m_fitpara[i]->DoRNGstart())
        {
            std::cout << "Randomizing start point for " << m_fitpara[i]->GetName() << std::endl;
            for(auto& p : vec1)
                p += (p * rng->Gaus(0.0, 0.1));
        }
        par_prefit.insert(par_prefit.end(), vec1.begin(), vec1.end());

        m_fitpara[i]->GetParSteps(vec1);
        par_step.insert(par_step.end(), vec1.begin(), vec1.end());

        m_fitpara[i]->GetParLimits(vec1, vec2);
        par_low.insert(par_low.end(), vec1.begin(), vec1.end());
        par_high.insert(par_high.end(), vec2.begin(), vec2.end());

        std::vector<bool> vec3;
        m_fitpara[i]->GetParFixed(vec3);
        par_fixed.insert(par_fixed.end(), vec3.begin(), vec3.end());
    }

    if(m_npar == 0)
    {
        std::cerr << "No fit parameters were defined." << std::endl;
        return;
    }

    std::cout << "===========================================" << std::endl;
    std::cout << "           Initilizing fitter              " << std::endl;
    std::cout << "===========================================" << std::endl;

    std::cout << "Minimizer settings..." << std::endl
              << "Minimizer: " << min_settings.minimizer << std::endl
              << "Algorithm: " << min_settings.algorithm << std::endl
              << "Likelihood: " << min_settings.likelihood << std::endl
              << "Strategy : " << min_settings.strategy << std::endl
              << "Print Lvl: " << min_settings.print_level << std::endl
              << "Tolerance: " << min_settings.tolerance << std::endl
              << "Max Iterations: " << min_settings.max_iter << std::endl
              << "Max Fcn Calls : " << min_settings.max_fcn << std::endl;

    m_fitter = ROOT::Math::Factory::CreateMinimizer(min_settings.minimizer.c_str(), min_settings.algorithm.c_str());
    m_fcn    = new ROOT::Math::Functor(this, &Fitter::CalcLikelihood, m_npar);

    m_fitter->SetFunction(*m_fcn);
    m_fitter->SetStrategy(min_settings.strategy);
    m_fitter->SetPrintLevel(min_settings.print_level);
    m_fitter->SetTolerance(min_settings.tolerance);
    m_fitter->SetMaxIterations(min_settings.max_iter);
    m_fitter->SetMaxFunctionCalls(min_settings.max_fcn);
    std::cout<<"m_npar = "<<m_npar<<std::endl;
    for(int i = 0; i < m_npar; ++i)
    {
        m_fitter->SetVariable(i, par_names[i], par_prefit[i], par_step[i]);
        m_fitter->SetVariableLimits(i, par_low[i], par_high[i]);

        if(par_fixed[i] == true)
            m_fitter->FixVariable(i);
    }
    par_var_fixed = par_fixed;

    std::cout << "Number of defined parameters: " << m_fitter->NDim() << std::endl
              << "Number of free parameters   : " << m_fitter->NFree() << std::endl
              << "Number of fixed parameters  : " << m_fitter->NDim() - m_fitter->NFree()
              << std::endl;

    TH1D h_prefit("hist_prefit_par_all", "hist_prefit_par_all", m_npar, 0, m_npar);
    TVectorD v_prefit_original(m_npar);
    TVectorD v_prefit_start(m_npar, par_prefit.data());

    int num_par = 1;
    for(int i = 0; i < m_fitpara.size(); ++i)
    {
        for(int j = 0; j < m_fitpara[i]->GetNpar(); ++j)
        {
            h_prefit.SetBinContent(num_par, m_fitpara[i]->GetParPrior(j));
            if(m_fitpara[i]->HasCovMat())
            {
                TMatrixDSym* cov_mat = m_fitpara[i]->GetCovMat();
                h_prefit.SetBinError(num_par, std::sqrt((*cov_mat)[j][j]));
            }
            else
                h_prefit.SetBinError(num_par, 0);

            v_prefit_original[num_par-1] = m_fitpara[i]->GetParOriginal(j);
            num_par++;
        }
    }


    m_dir->cd();
    h_prefit.Write();
    v_prefit_original.Write("vec_prefit_original");
    v_prefit_start.Write("vec_prefit_start");
}

bool Fitter::Fit(const std::vector<AnaSample*>& samples, bool stat_fluc)
{
    std::cout << "Starting to fit." << std::endl;
    m_samples = samples;

    if(m_fitter == nullptr)
    {
        std::cerr  << "In Fitter::Fit()\n"
                   << "Fitter has not been initialized." << std::endl;
        return false;
    }

    for(const auto& s : m_samples)
    {
        s->FillEventHist();
        s->FillDataHist(stat_fluc);
        s->SetLLHFunction(min_settings.likelihood);
    }

    SaveEventHist();

    bool did_converge = false;
    std::cout << "Fit prepared." << std::endl;
    std::cout << "Calling Minimize, running " << min_settings.algorithm << std::endl;
    did_converge = m_fitter->Minimize();

    if(!did_converge)
    {
        std::cout  << "Fit did not converge while running " << min_settings.algorithm
                   << std::endl;
        std::cout  << "Failed with status code: " << m_fitter->Status() << std::endl;
    }
    else
    {
        std::cout << "Fit converged." << std::endl
                  << "Status code: " << m_fitter->Status() << std::endl;

        std::cout << "Calling HESSE." << std::endl;
        did_converge = m_fitter->Hesse();
    }

    if(!did_converge)
    {
        std::cout  << "Hesse did not converge." << std::endl;
        std::cout  << "Failed with status code: " << m_fitter->Status() << std::endl;
    }
    else
    {
        std::cout << "Hesse converged." << std::endl
                  << "Status code: " << m_fitter->Status() << std::endl;
    }

    if(m_dir)
        SaveChi2();

    const int ndim        = m_fitter->NDim();
    const int nfree       = m_fitter->NFree();
    const double* par_val = m_fitter->X();
    const double* par_err = m_fitter->Errors();
    double cov_array[ndim * ndim];
    m_fitter->GetCovMatrix(cov_array);

    std::vector<double> par_val_vec(par_val, par_val + ndim);
    std::vector<double> par_err_vec(par_err, par_err + ndim);

    unsigned int par_offset = 0;
    TMatrixDSym cov_matrix(ndim, cov_array);
    for(const auto& fit_param : m_fitpara)
    {
        if(fit_param->IsDecomposed())
        {
            cov_matrix  = fit_param->GetOriginalCovMat(cov_matrix, par_offset);
            par_val_vec = fit_param->GetOriginalParameters(par_val_vec, par_offset);
        }
        par_offset += fit_param->GetNpar();
    }

    par_postfit = par_val_vec;

    TMatrixDSym cor_matrix(ndim);
    for(int r = 0; r < ndim; ++r)
    {
        for(int c = 0; c < ndim; ++c)
        {
            cor_matrix[r][c] = cov_matrix[r][c] / std::sqrt(cov_matrix[r][r] * cov_matrix[c][c]);
            if(std::isnan(cor_matrix[r][c]))
                cor_matrix[r][c] = 0;
        }
    }

    TVectorD postfit_globalcc(ndim);
    for(int i = 0; i < ndim; ++i)
        postfit_globalcc[i] = m_fitter->GlobalCC(i);

    TVectorD postfit_param(ndim, &par_val_vec[0]);
    std::vector<std::vector<double>> res_pars;
    std::vector<std::vector<double>> err_pars;
    int k = 0;
    for(int i = 0; i < m_fitpara.size(); i++)
    {
        const unsigned int npar = m_fitpara[i]->GetNpar();
        std::vector<double> vec_res;
        std::vector<double> vec_err;

        for(int j = 0; j < npar; j++)
        {
            vec_res.push_back(par_val_vec[k]);
            vec_err.push_back(std::sqrt(cov_matrix[k][k]));
            k++;
        }

        res_pars.emplace_back(vec_res);
        err_pars.emplace_back(vec_err);
    }

    for(int i = 0; i < m_fitpara.size(); i++)
    {
        res_pars.emplace_back(par_val_vec[i]);
        err_pars.emplace_back(std::sqrt(cov_matrix[i][i]));
    }

    if(k != ndim)
    {
        std::cout << "Number of parameters does not match." << std::endl;
        return false;
    }

    m_dir->cd();
    cov_matrix.Write("res_cov_matrix");
    cor_matrix.Write("res_cor_matrix");
    postfit_param.Write("res_vector");
    postfit_globalcc.Write("res_globalcc");

    SaveResults(res_pars, err_pars);
    SaveEventHist(true);

    if(m_save_events)
        SaveEventTree(res_pars);

    if(!did_converge)
        std::cout  << "Not valid fit result." << std::endl;
    std::cout << "Fit routine finished. Results saved." << std::endl;

    return did_converge;
}


double Fitter::FillSamples(std::vector<std::vector<double>>& new_pars)
{
    double chi2      = 0.0;
    bool output_chi2 = false;
    if((m_calls < 1001 && (m_calls % 100 == 0 || m_calls < 20))
       || (m_calls > 1001 && m_calls % 1000 == 0))
        output_chi2 = true;

    unsigned int par_offset = 0;
    for(int i = 0; i < m_fitpara.size(); ++i)
    {
        if(m_fitpara[i]->IsDecomposed())
        {
            new_pars[i] = m_fitpara[i]->GetOriginalParameters(new_pars[i]);
        }
        par_offset += m_fitpara[i]->GetNpar();

        m_fitpara[i]->ApplyParameters(new_pars[i]);
    }

    for(int s = 0; s < m_samples.size(); ++s)
    {
        const unsigned int num_pmts = m_samples[s]->GetNPMTs();
        const int pmttype = m_samples[s]->GetPMTType();
#pragma omp parallel for num_threads(m_threads)
        for(unsigned int i = 0; i < num_pmts; ++i)
        {
            AnaEvent* ev = m_samples[s]->GetPMT(i);
            ev->ResetEvWght();
            for(int j = 0; j < m_fitpara.size(); ++j)
            {
                if (m_fitpara[j]->GetPMTType()>=0 && m_fitpara[j]->GetPMTType() != pmttype) continue;
                m_fitpara[j]->ReWeight(ev, pmttype, s, i, new_pars[j]);
            }
        }

        m_samples[s]->FillEventHist();
        double sample_chi2 = m_samples[s]->CalcLLH();
        chi2 += sample_chi2;

        if(output_chi2)
        {
            std::cout << "Chi2 for sample " << m_samples[s]->GetName() << " is "
                      << sample_chi2 << std::endl;
        }
    }

    return chi2;
}

double Fitter::CalcLikelihood(const double* par)
{
    m_calls++;

    bool output_chi2 = false;
    if((m_calls < 1001 && (m_calls % 100 == 0 || m_calls < 20))
       || (m_calls > 1001 && m_calls % 1000 == 0))
        output_chi2 = true;

    int k           = 0;
    double chi2_sys = 0.0;
    double chi2_reg = 0.0;
    std::vector<std::vector<double>> new_pars;
    for(int i = 0; i < m_fitpara.size(); ++i)
    {
        const unsigned int npar = m_fitpara[i]->GetNpar();
        std::vector<double> vec;
        for(int j = 0; j < npar; ++j)
        {
            vec.push_back(par[k++]);
        }

        chi2_sys += m_fitpara[i]->GetChi2(vec);

        new_pars.push_back(vec);

        if(output_chi2)
        {
            std::cout << "Chi2 contribution from " << m_fitpara[i]->GetName() << " is "
                      << m_fitpara[i]->GetChi2(vec) << std::endl;
        }

    }


    double chi2_stat = FillSamples(new_pars);
    vec_chi2_stat.push_back(chi2_stat);
    vec_chi2_sys.push_back(chi2_sys);
    vec_chi2_reg.push_back(chi2_reg);

    if(m_calls % m_freq == 0 && m_save)
    {
        SaveParams(new_pars);
        SaveEventHist();
    }

    if(output_chi2)
    {
        std::cout << "Func Calls: " << m_calls << std::endl;
        std::cout << "Chi2 total: " << chi2_stat + chi2_sys + chi2_reg << std::endl;
        std::cout << "Chi2 stat : " << chi2_stat << std::endl
                  << "Chi2 syst : " << chi2_sys  << std::endl
                  << "Chi2 reg  : " << chi2_reg  << std::endl;
    }

    return chi2_stat + chi2_sys + chi2_reg;
}

void Fitter::SaveEventHist(bool is_final)
{
    for(auto& sample : m_samples)
    {
        std::stringstream ss;
        if(!is_final)
            sample->WriteDataHist(m_dir, ss.str());

        if(is_final)
            ss << "_final";
        else
            ss << "_iter" << m_calls;

        sample->WriteEventHist(m_dir, ss.str());
    }
}


void Fitter::SaveParams(const std::vector<std::vector<double>>& new_pars)
{
    std::vector<double> temp_vec;
    for(size_t i = 0; i < m_fitpara.size(); i++)
    {
        const unsigned int npar = m_fitpara[i]->GetNpar();
        const std::string name  = m_fitpara[i]->GetName();
        std::stringstream ss;

        ss << "hist_" << name << "_iter" << m_calls;
        TH1D h_par(ss.str().c_str(), ss.str().c_str(), npar, 0, npar);

        std::vector<std::string> vec_names;
        m_fitpara[i]->GetParNames(vec_names);
        for(int j = 0; j < npar; j++)
        {
            h_par.GetXaxis()->SetBinLabel(j + 1, vec_names[j].c_str());
            h_par.SetBinContent(j + 1, new_pars[i][j]);
            temp_vec.emplace_back(new_pars[i][j]);
        }
        m_dir->cd();
        h_par.Write();
    }

    TVectorD root_vec(temp_vec.size(), &temp_vec[0]);
    root_vec.Write(Form("vec_par_all_iter%d", m_calls));
}

void Fitter::SaveChi2()
{
    TH1D h_chi2stat("chi2_stat_periter", "chi2_stat_periter", m_calls + 1, 0, m_calls + 1);
    TH1D h_chi2sys("chi2_syst_periter", "chi2_syst_periter", m_calls + 1, 0, m_calls + 1);
    TH1D h_chi2reg("chi2_reg_periter", "chi2_reg_periter", m_calls + 1, 0, m_calls + 1);
    TH1D h_chi2tot("chi2_total_periter", "chi2_total_periter", m_calls + 1, 0, m_calls + 1);

    for(size_t i = 0; i < vec_chi2_stat.size(); i++)
    {
        h_chi2stat.SetBinContent(i + 1, vec_chi2_stat[i]);
        h_chi2sys.SetBinContent(i + 1, vec_chi2_sys[i]);
        h_chi2reg.SetBinContent(i + 1, vec_chi2_reg[i]);
        h_chi2tot.SetBinContent(i + 1, vec_chi2_sys[i] + vec_chi2_stat[i] + vec_chi2_reg[i]);
    }

    std::vector<double> v_chi2_prefit = {vec_chi2_stat.front(), vec_chi2_sys.front(), vec_chi2_reg.front()};
    v_chi2_prefit.push_back(std::accumulate(v_chi2_prefit.begin(), v_chi2_prefit.end(), 0.0, std::plus<double>()));

    std::vector<double> v_chi2_pstfit = {vec_chi2_stat.back(), vec_chi2_sys.back(), vec_chi2_reg.back()};
    v_chi2_pstfit.push_back(std::accumulate(v_chi2_pstfit.begin(), v_chi2_pstfit.end(), 0.0, std::plus<double>()));

    TVectorD chi2_prefit(v_chi2_prefit.size(), v_chi2_prefit.data());
    TVectorD chi2_pstfit(v_chi2_pstfit.size(), v_chi2_pstfit.data());

    m_dir->cd();
    h_chi2stat.Write();
    h_chi2sys.Write();
    h_chi2reg.Write();
    h_chi2tot.Write();

    chi2_prefit.Write("chi2_tuple_prefit");
    chi2_pstfit.Write("chi2_tuple_postfit");
}

void Fitter::SaveResults(const std::vector<std::vector<double>>& par_results,
                         const std::vector<std::vector<double>>& par_errors)
{
    for(std::size_t i = 0; i < m_fitpara.size(); i++)
    {
        const unsigned int npar = m_fitpara[i]->GetNpar();
        const std::string name  = m_fitpara[i]->GetName();
        std::vector<double> par_original;
        m_fitpara[i]->GetParOriginal(par_original);

        TMatrixDSym* cov_mat = m_fitpara[i]->GetOriginalCovMat();

        std::stringstream ss;

        ss << "hist_" << name << "_result";
        TH1D h_par_final(ss.str().c_str(), ss.str().c_str(), npar, 0, npar);

        ss.str("");
        ss << "hist_" << name << "_prior";
        TH1D h_par_prior(ss.str().c_str(), ss.str().c_str(), npar, 0, npar);

        ss.str("");
        ss << "hist_" << name << "_error_final";
        TH1D h_err_final(ss.str().c_str(), ss.str().c_str(), npar, 0, npar);

        ss.str("");
        ss << "hist_" << name << "_error_prior";
        TH1D h_err_prior(ss.str().c_str(), ss.str().c_str(), npar, 0, npar);

        std::vector<std::string> vec_names;
        m_fitpara[i]->GetParNames(vec_names);
        for(int j = 0; j < npar; j++)
        {
            h_par_final.GetXaxis()->SetBinLabel(j + 1, vec_names[j].c_str());
            h_par_final.SetBinContent(j + 1, par_results[i][j]);
            h_par_prior.GetXaxis()->SetBinLabel(j + 1, vec_names[j].c_str());
            h_par_prior.SetBinContent(j + 1, par_original[j]);
            h_err_final.GetXaxis()->SetBinLabel(j + 1, vec_names[j].c_str());
            h_err_final.SetBinContent(j + 1, par_errors[i][j]);

            double err_prior = 0.0;
            if(cov_mat != nullptr)
                err_prior = std::sqrt((*cov_mat)(j,j));

            h_err_prior.GetXaxis()->SetBinLabel(j + 1, vec_names[j].c_str());
            h_err_prior.SetBinContent(j + 1, err_prior);
        }

        m_dir->cd();
        h_par_final.Write();
        h_par_prior.Write();
        h_err_final.Write();
        h_err_prior.Write();
    }

}

void Fitter::ParameterScans(const std::vector<int>& param_list, unsigned int nsteps)
{
    std::cout << "Performing parameter scans..." << std::endl;

    //Internally Scan performs steps-1, so add one to actually get the number of steps
    //we ask for.
    unsigned int adj_steps = nsteps+1;
    double* x = new double[adj_steps] {};
    double* y = new double[adj_steps] {};

    for(const auto& p : param_list)
    {
        std::cout << "Scanning parameter " << p
                  << " (" << m_fitter->VariableName(p) << ")." << std::endl;

        bool success = m_fitter->Scan(p, adj_steps, x, y);

        TGraph scan_graph(nsteps, x, y);
        m_dir->cd();

        std::stringstream ss;
        ss << "par_scan_" << std::to_string(p);
        scan_graph.Write(ss.str().c_str());
    }

    delete[] x;
    delete[] y;
}

void Fitter::SaveEventTree(std::vector<std::vector<double>>& res_params)
{
    m_outtree = new TTree("PMTTree", "PMTTree");
    weight.resize(m_fitpara.size());
    InitOutputTree();

    for(size_t s = 0; s < m_samples.size(); s++)
    {
        sampleId = s;
        const unsigned int num_pmts = m_samples[s]->GetNPMTs();
        const int pmttype = m_samples[s]->GetPMTType();
        for(int i = 0; i < num_pmts; i++)
        {
            AnaEvent* ev = m_samples[s]->GetPMT(i);
            for(size_t j = 0; j < m_fitpara.size(); j++)
            {
                if (m_fitpara[j]->GetPMTType()>=0 && m_fitpara[j]->GetPMTType() != pmttype) weight[j] = 1.;
                else weight[j] = m_fitpara[j]->GetWeight(ev, pmttype, s, i, res_params[j]);
            }

            nPE     = ev->GetPE();
            R       = ev->GetR();
            costh   = ev->GetCosth();
            cosths  = ev->GetCosths();
            costhm  = ev->GetCosthm();
            phim    = ev->GetPhim();
            omega   = ev->GetOmega();
            PMT_id  = ev->GetPMTID();
            mPMT_id = ev->GetmPMTID();
            m_outtree->Fill();
        }
    }
    m_dir->cd();
    m_outtree->Write();
}

void Fitter::RunMCMCScan(int step, double stepsize, bool do_force_posdef, double force_padd, bool do_incompl_chol, double dropout_tol)
{
    const int ndim        = m_fitter->NDim();
    double cov_array[ndim * ndim];
    m_fitter->GetCovMatrix(cov_array);
    TMatrixDSym cov_matrix(ndim, cov_array);

    ToyThrower* toy_thrower = new ToyThrower(cov_matrix, false, 1E-48);
    if(do_force_posdef)
    {
        if(!toy_thrower->ForcePosDef(force_padd, 1E-48))
        {
            std::cout << "Covariance matrix could not be made positive definite.\n"
                      << "Exiting." << std::endl;
            return;
        }
    }

    // LU decomposition
    // Use the lower-triangular L matrix to generate correlated variables
    if(do_incompl_chol)
    {
        std::cout << "Performing incomplete Cholesky decomposition." << std::endl;
        toy_thrower->IncompCholDecomp(dropout_tol, true);
    }
    else
    {
        std::cout << "Performing ROOT Cholesky decomposition." << std::endl;
        toy_thrower->SetupDecomp(1E-48);
    }

    InitMCMCOutputTree();
    // First MCMC step is the best-fit point
    par_mcmc = par_postfit;
    m_chi2 = CalcLikelihood(par_mcmc.data());
    m_jump = 0;
    m_mcmctree->Fill();

    std::vector<double> toy(ndim, 0.0);
    for (int i=0; i<step; i++)
    {
        m_jump = 0;

        toy_thrower->Throw(toy);
        for (int j=0;j<ndim;j++)
        {
            toy[j] = par_var_fixed[j] ? 0: stepsize*toy[j] ; // do not move the fixed variables
            toy[j] += par_mcmc[j];
        }

        double chi2_next = CalcLikelihood(toy.data());
        // Metropolis-Hastings Algorithm
        if ( rng->Uniform() < exp(-chi2_next/2.+m_chi2/2.) )
        {
            m_chi2 = chi2_next;
            par_mcmc = toy;
            m_jump = 1;
        }
        m_mcmctree->Fill();
    }

    m_dir->cd();
    m_mcmctree->Write();
}