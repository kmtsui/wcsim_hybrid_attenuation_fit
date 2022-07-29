std::vector<int> pol_orders; // order of polynomial in each piece 
std::vector<double> pol_range; // applicable range for each polynomial
std::vector<double> CalcPol(const double* par, std::vector<double> costh_array)
{
    // Derive the coefficients for each polynomial
    std::vector<std::vector<double>> pol_coeff;
    std::vector<double> pol_p0, pol_p1;
    int par_index = 0;
    for (int i=0;i<pol_orders.size();i++)
    {
        std::vector<double> coeff;
        if (i==0) // for the first polynomial, the coefficients are unconstrained
        {
            for (int j=0;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[j]);
                par_index++;
            } 
        }
        else // for others, we need to match the 0-th and 1-st order derivatives 
        {
            coeff.push_back(pol_p0[i-1]);
            coeff.push_back(pol_p1[i-1]);
            for (int j=2;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[par_index]);
                par_index++;
            } 
        }
        pol_coeff.emplace_back(coeff);
        // std::cout<<"pol "<<i<<": ";
        // for (int j=0;j<pol_coeff[i].size();j++)
        //     std::cout<<pol_coeff[i][j]<<" ";
        // std::cout<<std::endl;
        double p0 = 0;
        double p1 = 0;
        for (int j=0;j<=pol_orders[i];j++) // store the 0-th and 1-st order derivatives at end-point as boundary conditions
        {
            p0 += coeff[j]*TMath::Power(pol_range[i+1]-pol_range[i],j);
            p1 += coeff[j]*j*TMath::Power(pol_range[i+1]-pol_range[i],j-1);
        } 
        pol_p0.push_back(p0);
        pol_p1.push_back(p1);
    }

    std::vector<double> angular_response;
    for (int k=0;k<costh_array.size();k++) // actually calculate the angular response
    {
        double costh = costh_array[k];
        double val = 0;
        for (int i=0;i<pol_orders.size();i++)
        {
            if (costh>=pol_range[i] && costh<pol_range[i+1])
            {
                //std::cout<<"Using pol"<<i<<std::endl;
                for (int j=0;j<=pol_orders[i];j++)
                {
                    val += pol_coeff[i][j]*TMath::Power(costh-pol_range[i],j);
                }
            }
        }
        //std::cout<<"CalcPol:: costh = "<<costh<<", val = "<<val<<std::endl;
        angular_response.push_back(val);
    }
        

    return angular_response;
}

std::vector<double> data_val; 
std::vector<double> costh_val;
TMatrixDSym cov_mat;
double CalcLikelihood(const double* par)
{
    std::vector<double> fit_val = CalcPol(par,costh_val);

    double chi2 = 0;
    for (int i=0;i<fit_val.size();i++)
        for (int j=0;j<fit_val.size();j++)
    {
        double x = fit_val[i]-data_val[i];
        double y = fit_val[j]-data_val[j];
        chi2 += x*y*cov_mat(i,j);
    }

    return chi2;
}

void FitAngularResponsePol()
{
    TFile f("/bundle/data/T2K/users/kmtsui/LI/fitter/TN_update/diffuser4_400nm_nominal_combined_simplest.root");
    //TFile f("/bundle/data/T2K/users/kmtsui/LI/fitter/TN/diffusr4_400nm_nominal_mPMT.root");
    TVectorD* res_vector = (TVectorD*)f.Get("res_vector");
    TMatrixDSym* res_cov_matrix = (TMatrixDSym*)f.Get("res_cov_matrix");
    int startingIndex = 42; // starting index of the angular response parameters
    int nParameters = 40;   // number of angular response parameters
    double costh_min = 0.0;
    double costh_max = 1.0;
    TH1D* hist_postfit = new TH1D("","",nParameters,costh_min,costh_max);

    // Setup the polynomials to fit
    // BnL example here consists of a 3rd order polynomial in [0.0,0.7], and a 3rd order polynomial in [0.7,1.0]
    // int orders[] = {3,3};
    // double ranges[] = {0,0.7,1.0};

    // mPMT requires more polynomial 
    int orders[] = {3,3,3,4};
    double ranges[] = {0.0,0.3,0.6,0.75,1.0};

    bool plot_extra_pol = false; // use this with par_pol to plot an extra polynomials for comparison
    double par_pol[] = {0.00118237,0.0225071,0.388372,0.243631,0.0605638,0.386826,-0.971796,9.61525,-5.79266,28.192,-61.1739};

    int ndof = 0;
    std::vector<int> index_array;
    for (int i=0;i<nParameters;i++)
    {
        int idx = i+startingIndex;
        double val = (*res_vector)[idx];
        double err = sqrt((*res_cov_matrix)[idx][idx]);
        if (err>0) // remove the fixed variables which are not fit in optical_fit
        {
            if (ndof==0) hist_postfit->SetMinimum(val*0.9);
            hist_postfit->SetBinContent(i+1,val);
            hist_postfit->SetBinError(i+1,err);
            ndof++;
            index_array.push_back(idx);
            data_val.push_back(val);
            double costh = hist_postfit->GetBinCenter(i+1);
            costh_val.push_back(costh);
            //std::cout<<"costh = "<<costh<<", val = "<<val<<std::endl;
        }
    }


    cov_mat.ResizeTo(ndof, ndof);
    cov_mat.Zero();

    for (int i=0;i<ndof;i++) // only consider data points with valid uncertainty
    {
            for (int j=0;j<ndof;j++)
            {
                cov_mat(i,j) = (*res_cov_matrix)[index_array[i]][index_array[j]];
            }
    }
    // for (int i=0;i<cov_mat.GetNrows();i++) for (int j=0;j<cov_mat.GetNrows();j++) {
    //     std::cout<<"cov_mat("<<i<<","<<j<<")="<<cov_mat(i, j)<<std::endl;
    // }

    double det = 0;
    double total_add = 0;
    TDecompLU inv_test;
    TMatrixD inv_matrix(cov_mat);
    while (!inv_test.InvertLU(inv_matrix, 1E-48, &det))
    {
        std::cerr << "Covariance matrix is non invertable. Determinant is " << det
                  << std::endl;
        for(int i = 0; i < inv_matrix.GetNrows(); ++i)
            inv_matrix(i,i) += 1e-24;
        total_add += 1e-24;
    }
    std::cout << "Covariance matrix inverted successfully." << std::endl;
    cov_mat.SetMatrixArray(inv_matrix.GetMatrixArray());
    // for (int i=0;i<cov_mat.GetNrows();i++) for (int j=0;j<cov_mat.GetNrows();j++) {
    //     std::cout<<"cov_mat("<<i<<","<<j<<")="<<cov_mat(i, j)<<std::endl;
    // }
    std::cout << "Added " << total_add << " to force positive-definite."
              << std::endl;

    pol_orders.assign(orders, orders+sizeof(orders)/sizeof(int));
    pol_range.assign(ranges, ranges+sizeof(ranges)/sizeof(double));
    int m_npar = pol_orders[0]+1;
    for (int i=1;i<pol_orders.size();i++)
        m_npar += pol_orders[i]+1-2;
    
    char* minName  = (char*)"Minuit2";
    char* algoName = (char*)"Migrad";

    // Setup the minimizer
    ROOT::Math::Minimizer* m_fitter = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    ROOT::Math::Functor m_fcn(&CalcLikelihood, m_npar);

    std::cout<<"Number of free parameters = "<<m_fcn.NDim()<<std::endl;

    m_fitter->SetFunction(m_fcn);
    m_fitter->SetStrategy(1);
    m_fitter->SetPrintLevel(2);
    m_fitter->SetTolerance(1.e-4);
    m_fitter->SetMaxIterations(1.e6);
    m_fitter->SetMaxFunctionCalls(1.e9);

    // A reasonable value for the fit to start
    m_fitter->SetVariable(0, "pol0_p0", 0.25, 0.1);
    int par_idx = 1;
    for (int i=0;i<pol_orders.size();i++)
    {
        if (i==0)
        {
            for (int j=1;j<=pol_orders[i];j++)
            {
                m_fitter->SetVariable(par_idx,Form("pol%i_p%i",i,j),0,0.1);
                par_idx++;
            }
        }
        else
        {
            for (int j=2;j<=pol_orders[i];j++)
            {
                m_fitter->SetVariable(par_idx,Form("pol%i_p%i",i,j),0,0.1);
                par_idx++;
            }
        }
    }
    
    bool did_converge = false;
    std::cout <<"Fit prepared." << std::endl;
    std::cout <<"Calling Minimize, running " << minName << ", "<< algoName << std::endl;
    did_converge = m_fitter->Minimize();

    if(!did_converge)
    {
        std::cout << "Fit did not converge."<< std::endl;
        std::cout << "Failed with status code: " << m_fitter->Status() << std::endl;
    }
    else
    {
        std::cout  << "Fit converged." << std::endl
                   << "Status code: " << m_fitter->Status() << std::endl;

        std::cout << "Calling HESSE." << std::endl;
        did_converge = m_fitter->Hesse();
    }

    if(!did_converge)
    {
        std::cout << "Hesse did not converge." << std::endl;
        std::cout << "Failed with status code: " << m_fitter->Status() << std::endl;
    }
    else
    {
        std::cout  << "Hesse converged." << std::endl
                   << "Status code: " << m_fitter->Status() << std::endl;
    }

    const double* par_val = m_fitter->X();
    const double* par_err = m_fitter->Errors();

    // for (int i=0;i<m_npar;i++) {
    //     std::cout<<m_fitter->VariableName(i)<<": "<<par_val[i]<<" +/- "<<par_err[i]<<std::endl;
    // }

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    hist_postfit->GetXaxis()->SetTitle("cos#theta");
    hist_postfit->GetYaxis()->SetTitle("A(cos#theta)");
    hist_postfit->SetLineWidth(2);
    hist_postfit->SetMinimum(0);
    hist_postfit->Draw("E0");
    // Also draw the fitted polynomials as histogram
    TH1D* h_fit_pol = new TH1D("","",1000,costh_min,costh_max);
    std::vector<double> pltpts;
    for (int i=1;i<=h_fit_pol->GetNbinsX();i++)
    {
        pltpts.push_back(h_fit_pol->GetBinCenter(i));
    }
    std::vector<double> valpts = CalcPol(par_val,pltpts);
    for (int i=1;i<=h_fit_pol->GetNbinsX();i++)
    {
        h_fit_pol->SetBinContent(i,valpts[i-1]);
    }
    h_fit_pol->SetLineColor(kRed);
    h_fit_pol->SetLineWidth(2);
    h_fit_pol->Draw("same");
    if (plot_extra_pol)
    {
        TH1D* h_fit_pol_fit = new TH1D("","",1000,costh_min,costh_max);
        valpts = CalcPol(par_pol,pltpts);
        for (int i=1;i<=h_fit_pol_fit->GetNbinsX();i++)
        {
            h_fit_pol_fit->SetBinContent(i,valpts[i-1]);
        }
        h_fit_pol_fit->SetLineColor(kGreen);
        h_fit_pol_fit->SetLineStyle(2);
        h_fit_pol_fit->SetLineWidth(2);
        h_fit_pol_fit->Draw("same");

        // h_fit_pol_fit->Divide(h_fit_pol);
        // h_fit_pol_fit->Draw();
    }
    for (int i=1;i<pol_range.size();i++)
    {
        TLine* l = new TLine(pol_range[i],hist_postfit->GetMinimum(),pol_range[i],hist_postfit->GetMaximum());
        l->SetLineStyle(2);
        l->Draw("same");
    }
    TLatex latex;
    latex.SetTextSize(0.06);
    for (int i=0;i<pol_orders.size();i++)
    {
        latex.DrawLatex(pol_range[i]+0.02,hist_postfit->GetMaximum(),Form("pol%i",pol_orders[i]));
    }
    c1->SaveAs("test.pdf");
}