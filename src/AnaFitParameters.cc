#include "AnaFitParameters.hh"

AnaFitParameters::AnaFitParameters(const std::string& par_name, const int pmttype)
    : m_name(par_name)
    , m_pmttype(pmttype)
    , Npar(0)
    , m_rng_start(false)
    , m_do_cap_weights(false)
    , m_weight_cap(1000)
    , m_info_frac(1.00)
    , m_decompose(false)
    , eigen_decomp(nullptr)
    , covariance(nullptr)
    , covarianceI(nullptr)
    , original_cov(nullptr)

{
    m_func = new Identity;
    m_func_type = kIdentity;
    std::cout<<"Setting up parameter "<<m_name<<std::endl;
}

AnaFitParameters::~AnaFitParameters()
{
    if(eigen_decomp != nullptr)
        delete eigen_decomp;
    if(covariance != nullptr)
        delete covariance;
    if(covarianceI != nullptr)
        delete covarianceI;
    if(original_cov != nullptr)
        delete original_cov;

}


bool AnaFitParameters::CheckDims(const std::vector<double>& params) const
{
    bool vector_size = false;

    if(params.size() == pars_prior.size())
    {
        vector_size = true;
    }

    else
    {
        std::cerr << "Dimension of parameter vector does not match priors.\n"
                  << "Prams size is: " << params.size() << std::endl
                  << "Prior size is: " << pars_prior.size() << std::endl;
        vector_size = false;
    }

    return vector_size;
}


void AnaFitParameters::SetParameterFunction(const std::string& func_name)
{
    if(m_func != nullptr)
        delete m_func;

    if(func_name.empty())
    {
        std::cout << "Parameter function name empty. Setting to identity by default." << std::endl;
        m_func = new Identity;
        m_func_type = kIdentity;
    }
    else if(func_name == "Identity")
    {
        std::cout << "Setting function to Identity." << std::endl;
        m_func = new Identity;
        m_func_type = kIdentity;
    }
    else if(func_name == "Attenuation")
    {
        std::cout << "Setting function to Attenuation." << std::endl;
        m_func = new Attenuation;
        m_func_type = kAttenuation;
    }
    // else if(func_name == "Scatter")
    // {
    //     std::cout << "Setting function to Scatter." << std::endl;
    //     m_func = new Scatter;
    //     m_func_type = kScatter;
    // }
    else
    {
        std::cout << "Invalid function name. Setting to identity by default." << std::endl;
        m_func = new Identity;
        m_func_type = kIdentity;
    }
}

void AnaFitParameters::InitParameters(std::vector<std::string> names, std::vector<double> priors, std::vector<double> steps, 
                                      std::vector<double> lows, std::vector<double> highs, std::vector<bool> fixed)
{
    SetParNames(names);
    SetParPriors(priors);
    SetParSteps(steps);
    SetParLimits(lows,highs);
    SetParFixed(fixed);

    Npar = pars_name.size();
    pars_original = pars_prior;

    std::cout<<"Number of parameters = "<<Npar<<std::endl;

    if(m_decompose)
    {
        pars_prior = eigen_decomp -> GetDecompParameters(pars_prior);
        pars_limlow = std::vector<double>(Npar, -100);
        pars_limhigh = std::vector<double>(Npar, 100);

        const int idx = eigen_decomp -> GetInfoFraction(m_info_frac);
        for(int i = idx; i < Npar; ++i)
            pars_fixed[i] = true;

        std::cout << "Decomposed parameters.\n"
                  << "Keeping the " << idx << " largest eigen values.\n"
                  << "Corresponds to " << m_info_frac * 100.0
                  << "\% total variance.\n";
    }

}

void AnaFitParameters::InitEventMap(std::vector<AnaSample*> &sample)
{
    m_evmap.clear();

    std::vector<bool> params_used(Npar,false);

    for(std::size_t s=0; s < sample.size(); s++)
    {
        std::vector<int> sample_map;
        for(int i=0; i < sample[s] -> GetNPMTs(); i++)
        {
            AnaEvent* ev = sample[s] -> GetPMT(i);

            std::vector<double> binvar;
            for (auto t : m_binvar)
                binvar.push_back(ev->GetEventVar(t));

            const int bin = m_bm.GetBinIndex(binvar);

            if (m_pmttype == -1 || sample[s]->GetPMTType() == m_pmttype)
            {
                sample_map.push_back(bin);
                params_used[bin]=true;
            }
            else sample_map.push_back(PASSEVENT);
        }
        std::cout<<"In AnaFitParameters::InitEventMap, built event map for sample "<< sample[s]->GetName() << " of total "<< sample[s] -> GetNPMTs() << "PMTs"<<std::endl;
        m_evmap.push_back(sample_map);
    }

    for (int i=0;i<Npar;i++)
        if (!params_used[i]) pars_fixed[i]=true;
}

void AnaFitParameters::ReWeight(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params)
{
#ifndef NDEBUG
    if(m_evmap.empty()) //need to build an event map first
    {
        std::cerr  << "In AnaFitParameters::ReWeight()\n"
                   << "Need to build event map index for " << m_name << std::endl;
        return;
    }
#endif

    if (m_pmttype >=0 && pmttype != m_pmttype) return;

    const int bin = m_evmap[nsample][nevent];
    if(bin == PASSEVENT || bin == BADBIN)
        return;
    else
    {
#ifndef NDEBUG
        if(bin > params.size())
        {
            std::cout  << "In AnaFitParameters::ReWeight()\n"
                        << "Number of bins in " << m_name << " does not match num of parameters.\n"
                        << "Setting event weight to zero." << std::endl;
            event -> AddEvWght(0.0);
        }
#endif
        double wgt = (*m_func)(params[bin],*event);

        // if (m_func_type==kScatter) event -> SetTailPE(wgt);
        // else event -> AddEvWght(wgt);
        event -> AddEvWght(wgt);
    }
}

double AnaFitParameters::GetWeight(AnaEvent* event, int pmttype, int nsample, int nevent, std::vector<double>& params)
{
#ifndef NDEBUG
    if(m_evmap.empty()) //need to build an event map first
    {
        std::cerr  << "In AnaFitParameters::ReWeight()\n"
                   << "Need to build event map index for " << m_name << std::endl;
        return 1.;
    }
#endif

    if (m_pmttype >=0 && pmttype != m_pmttype) return 1.;

    const int bin = m_evmap[nsample][nevent];
    if(bin == PASSEVENT || bin == BADBIN)
        return 1.;
    else
    {
#ifndef NDEBUG
        if(bin > params.size())
        {
            std::cout  << "In AnaFitParameters::GetWeight()\n"
                        << "Number of bins in " << m_name << " does not match num of parameters.\n"
                        << "Setting event weight to zero." << std::endl;
            return 1.;
        }
#endif
        double wgt = (*m_func)(params[bin],*event);
        return wgt;
    }
}

void AnaFitParameters::SetCovarianceMatrix(const TMatrixDSym& covmat, bool decompose)
{
    if(covariance != nullptr)
        delete covariance;
    if(covarianceI != nullptr)
        delete covarianceI;
    if(original_cov != nullptr)
        delete original_cov;
    if(eigen_decomp != nullptr)
        delete eigen_decomp;

    if(decompose)
    {
        m_decompose  = true;
        eigen_decomp = new EigenDecomp(covmat);
        original_cov = new TMatrixDSym(covmat);
        covariance   = new TMatrixDSym(eigen_decomp->GetEigenCovMat());
        covarianceI  = new TMatrixDSym(eigen_decomp->GetEigenCovMat());
    }
    else
    {
        original_cov = new TMatrixDSym(covmat);
        covariance  = new TMatrixDSym(covmat);
        covarianceI = new TMatrixDSym(covmat);
    }

    double det = 0;
    TDecompLU inv_test;
    TMatrixD inv_matrix(*covariance);
    if(inv_test.InvertLU(inv_matrix, 1E-48, &det))
    {
        covarianceI->SetMatrixArray(inv_matrix.GetMatrixArray());
        std::cout << "Covariance matrix inverted successfully." << std::endl;
    }
    else
    {
        std::cerr << "In AnaFitParameters::SetCovarianceMatrix():\n"
                  << "Covariance matrix is non invertable. Determinant is " << det
                  << std::endl;
        return;
    }

    std::cout << "Covariance matrix size: " << covariance->GetNrows()
              << " x " << covariance->GetNrows() << " for " << this->m_name << std::endl;
}

double AnaFitParameters::GetChi2(const std::vector<double>& params) const
{
    if(covariance == nullptr)
        return 0.0;

    if(CheckDims(params) == false)
    {
        std::cout << "In AnaFitParameters::GetChi2()\n"
                  << "Dimension check failed. Returning zero." << std::endl;
        return 0.0;
    }

    double chi2 = 0;
    for(int i = 0; i < covarianceI->GetNrows(); i++)
    {
        for(int j = 0; j < covarianceI->GetNrows(); j++)
        {
            chi2
                += (params[i] - pars_prior[i]) * (params[j] - pars_prior[j]) * (*covarianceI)(i, j);
        }
    }

    return chi2;
}

