double truth_alpha(double wavelength, double ABWFF=1.30, double RAYFF=0.75) {
    const int NUMENTRIES_water=60;
    const double GeV=1.e9;
    const double cm=1;

    double ENERGY_water[NUMENTRIES_water] =
     { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
       1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
       1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV, 1.82352e-09*GeV, 
       1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
       1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
       2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
       2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
       2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
       2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
       2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
       3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
       3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
       3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
       4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
       5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };

    double ABSORPTION_water[NUMENTRIES_water] =
    {
        16.1419*cm*ABWFF,  18.278*cm*ABWFF, 21.0657*cm*ABWFF, 24.8568*cm*ABWFF, 30.3117*cm*ABWFF, 
        38.8341*cm*ABWFF, 54.0231*cm*ABWFF, 81.2306*cm*ABWFF, 120.909*cm*ABWFF, 160.238*cm*ABWFF, 
        193.771*cm*ABWFF, 215.017*cm*ABWFF, 227.747*cm*ABWFF,  243.85*cm*ABWFF, 294.036*cm*ABWFF, 
        321.647*cm*ABWFF,  342.81*cm*ABWFF, 362.827*cm*ABWFF, 378.041*cm*ABWFF, 449.378*cm*ABWFF,
        739.434*cm*ABWFF, 1114.23*cm*ABWFF, 1435.56*cm*ABWFF, 1611.06*cm*ABWFF, 1764.18*cm*ABWFF, 
        2100.95*cm*ABWFF,  2292.9*cm*ABWFF, 2431.33*cm*ABWFF,  3053.6*cm*ABWFF, 4838.23*cm*ABWFF, 
        6539.65*cm*ABWFF, 7682.63*cm*ABWFF, 9137.28*cm*ABWFF, 12220.9*cm*ABWFF, 15270.7*cm*ABWFF, 
        19051.5*cm*ABWFF, 23671.3*cm*ABWFF, 29191.1*cm*ABWFF, 35567.9*cm*ABWFF,   42583*cm*ABWFF,
        49779.6*cm*ABWFF, 56465.3*cm*ABWFF,   61830*cm*ABWFF, 65174.6*cm*ABWFF, 66143.7*cm*ABWFF,   
        64820*cm*ABWFF,   61635*cm*ABWFF, 57176.2*cm*ABWFF, 52012.1*cm*ABWFF, 46595.7*cm*ABWFF, 
        41242.1*cm*ABWFF, 36146.3*cm*ABWFF, 31415.4*cm*ABWFF, 27097.8*cm*ABWFF, 23205.7*cm*ABWFF, 
        19730.3*cm*ABWFF, 16651.6*cm*ABWFF, 13943.6*cm*ABWFF, 11578.1*cm*ABWFF, 9526.13*cm*ABWFF
    };

    double RAYLEIGH_water[NUMENTRIES_water] = {
            386929*cm*RAYFF,  366249*cm*RAYFF,  346398*cm*RAYFF,  327355*cm*RAYFF,  309097*cm*RAYFF,  
            291603*cm*RAYFF,  274853*cm*RAYFF,  258825*cm*RAYFF,  243500*cm*RAYFF,  228856*cm*RAYFF,  
            214873*cm*RAYFF,  201533*cm*RAYFF,  188816*cm*RAYFF,  176702*cm*RAYFF,  165173*cm*RAYFF,
            154210*cm*RAYFF,  143795*cm*RAYFF,  133910*cm*RAYFF,  124537*cm*RAYFF,  115659*cm*RAYFF,  
            107258*cm*RAYFF, 99318.2*cm*RAYFF, 91822.2*cm*RAYFF,   84754*cm*RAYFF, 78097.3*cm*RAYFF, 
            71836.5*cm*RAYFF,   65956*cm*RAYFF, 60440.6*cm*RAYFF, 55275.4*cm*RAYFF, 50445.6*cm*RAYFF,
            45937*cm*RAYFF, 41735.2*cm*RAYFF, 37826.6*cm*RAYFF, 34197.6*cm*RAYFF, 30834.9*cm*RAYFF, 
            27725.4*cm*RAYFF, 24856.6*cm*RAYFF, 22215.9*cm*RAYFF, 19791.3*cm*RAYFF, 17570.9*cm*RAYFF,   
            15543*cm*RAYFF, 13696.6*cm*RAYFF, 12020.5*cm*RAYFF, 10504.1*cm*RAYFF, 9137.15*cm*RAYFF,
            7909.45*cm*RAYFF,  6811.3*cm*RAYFF, 5833.25*cm*RAYFF,  4966.2*cm*RAYFF, 4201.36*cm*RAYFF, 
            3530.28*cm*RAYFF, 2944.84*cm*RAYFF, 2437.28*cm*RAYFF, 2000.18*cm*RAYFF,  1626.5*cm*RAYFF, 
            1309.55*cm*RAYFF, 1043.03*cm*RAYFF, 821.016*cm*RAYFF,  637.97*cm*RAYFF, 488.754*cm*RAYFF
    };

    TGraph* gr_abs = new TGraph(NUMENTRIES_water,ENERGY_water,ABSORPTION_water);
    TGraph* gr_ray = new TGraph(NUMENTRIES_water,ENERGY_water,RAYLEIGH_water);
    double photoEnergy = 1239.84193/wavelength;
    double abslength = gr_abs->Eval(photoEnergy,0,"S");
    double raylength = gr_ray->Eval(photoEnergy,0,"S");
    double alpha = 1./(1./abslength+1./raylength);
    std::cout<<"Absorption length = "<<abslength<<" cm"<<std::endl;
    std::cout<<"Rayleigh length = "<<raylength<<" cm"<<std::endl;
    std::cout<<"Attenutation length = "<<alpha<<" cm"<<std::endl;

    return alpha;
}

double PoissonLLH (double mc, double w2, double data)
{

        // Standard Poisson LLH.
        double chi2 = 0.0;
        if(mc > 0.0)
        {
            chi2 = 2 * (mc - data);
            if(data > 0.0)
                chi2 += 2 * data * std::log(data / mc);
        }

        return (chi2 >= 0.0) ? chi2 : 0.0;

}

TH1D* hRate1; // number of PE per PMT
TH2D* hBinnedRate; // number of PE binned in R and costh
int m_calls;
std::vector<double> mPMT_R;
std::vector<double> mPMT_costh;
std::vector<double> costh_array;
std::vector<bool> mPMT_use;
double CalcLikelihood(const double* par)
{
    m_calls++;

    bool output_chi2 = false;
    if((m_calls < 1001 && (m_calls % 100 == 0 || m_calls < 20))
       || (m_calls > 1001 && m_calls % 1000 == 0))
        output_chi2 = true;

    double chi2_stat = 0;

    double* data_w  = hRate1->GetArray();
    double* data_w2 = hRate1->GetSumw2()->GetArray();

    TH2D* hBinnedRate_fit = (TH2D*)hBinnedRate->Clone();
    hBinnedRate_fit->Reset();

    for (int i=0;i<hRate1->GetNbinsX();i++) {
        if (!mPMT_use[i]) continue;
        double value = TMath::Exp(-mPMT_R[i]/par[0])/mPMT_R[i]/mPMT_R[i]*9000*9000; //an arbitrary normalization
        double costh_pmt = mPMT_costh[i];
        int costh_idx = hBinnedRate_fit->GetXaxis()->FindBin(costh_pmt);
        if (costh_idx>=1&&costh_idx<=hBinnedRate_fit->GetNbinsX()) {
            value*=par[costh_idx];
            hBinnedRate_fit->Fill(costh_pmt,mPMT_R[i],value);
        }
    }

    for (int i=1;i<=hBinnedRate_fit->GetNbinsX();i++) 
        for (int j=1;j<=hBinnedRate_fit->GetNbinsY();j++) {
            chi2_stat += PoissonLLH(hBinnedRate_fit->GetBinContent(i,j), 0, hBinnedRate->GetBinContent(i,j));
        }

    if(output_chi2)
    {
        std::cout << "Func Calls: " << m_calls << std::endl;
        std::cout << "Chi2 stat : " << chi2_stat << std::endl;
    }

    return chi2_stat;
}


void run_fit(const char* minName = "Minuit2", const char* algoName="Migrad"){
    int m_npar = hBinnedRate->GetNbinsX()+1; // number of costh bins + one alpha parameter
    m_calls = 0;
    ROOT::Math::Minimizer* m_fitter = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    ROOT::Math::Functor m_fcn(&CalcLikelihood, m_npar);

    std::cout<<"Number of free parameters = "<<m_fcn.NDim()<<std::endl;

    m_fitter->SetFunction(m_fcn);
    m_fitter->SetStrategy(1);
    m_fitter->SetPrintLevel(2);
    m_fitter->SetTolerance(1.e-4);
    m_fitter->SetMaxIterations(1.e6);
    m_fitter->SetMaxFunctionCalls(1.e9);

    m_fitter->SetVariable(0, "alpha", 10000, 100);
    for (int i=1;i<m_npar;i++)
        m_fitter->SetVariable(i, Form("norm%i",i), hRate1->GetMaximum()/100., hRate1->GetMaximum()/100.);
    
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

    for (int i=0;i<m_npar;i++) {
        std::cout<<m_fitter->VariableName(i)<<": "<<par_val[i]<<" +/- "<<par_err[i]<<std::endl;
    }

}


void fit_all(   std::string filename, int nmPMT_on=0, // number of mPMT modules used fit, 0 = using all
                double timetof_min = -952, double timetof_max = -945, // hit time window
                int nbins_costh = 10, double costh_min = 0.5, double costh_max = 1., // binning in costh
                int nbins_dist=100, double dist_min = 1000, double dist_max=9000, // binning R
                double cosths_min = 0.766 // limit due to source opening angle 
            ) 
{
    gROOT->Reset();
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    TChain* hitRate_pmtType1 = new TChain("hitRate_pmtType1");
    hitRate_pmtType1->Add(filename.c_str());

    //Only the first file is used to extract the PMT geometry
    TFile* f = hitRate_pmtType1->GetFile();

    double nHits, nPE, dist, costh, cosths, timetof;
    int mPMT_id;

    TH2D* hPMT1 = new TH2D("","",nbins_costh,costh_min,costh_max,nbins_dist,dist_min,dist_max);
    TTree* pmt_type1 = (TTree*)f->Get("pmt_type1");
    pmt_type1->SetBranchAddress("dist",&dist);
    pmt_type1->SetBranchAddress("costh",&costh);
    pmt_type1->SetBranchAddress("cosths",&cosths);
    pmt_type1->SetBranchAddress("mPMT_id",&mPMT_id);

    // uniformly masking mPMT modules when requested
    int nPMTpermPMT = 19;
    const int nmPMT_sim = pmt_type1->GetEntries();
    int mPMT_mask[nmPMT_sim];
    if (nmPMT_on>0){
        double mPMT_frac = (nmPMT_on+0.)/(nmPMT_sim/nPMTpermPMT);
        int mPMT_count = 0;
        for (int i=0;i<nmPMT_sim/nPMTpermPMT;i++){
            if ((mPMT_count+0.)/(i+1.)<mPMT_frac && mPMT_count<nmPMT_on) {
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    mPMT_mask[j]=0;
                }
                mPMT_count++;
            } else {
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    mPMT_mask[j]=1;
                }
            }
        }
    }

    int nmPMT_used=0;
    mPMT_use.clear();mPMT_R.clear();mPMT_costh.clear();
    mPMT_use.resize(nmPMT_sim,false);
    for (int i=0;i<pmt_type1->GetEntries();i++) {
        pmt_type1->GetEntry(i);
        mPMT_R.push_back(dist);
        mPMT_costh.push_back(-costh);
        if (nmPMT_on>0)
            if (mPMT_mask[mPMT_id]==1) continue; // ignore masked PMT
        if (cosths>cosths_min) // only include PMT within the source opening angle 
        {
            hPMT1->Fill(-costh,dist); 
            mPMT_use[i]=true;
            nmPMT_used++;
        }
    }
    
    TH1::SetDefaultSumw2(true);

    hRate1 = new TH1D("","",nmPMT_sim,0,nmPMT_sim);
    hBinnedRate = new TH2D("","",nbins_costh,costh_min,costh_max,nbins_dist,dist_min,dist_max);

    hitRate_pmtType1->SetBranchAddress("dist",&dist);
    hitRate_pmtType1->SetBranchAddress("costh",&costh);
    hitRate_pmtType1->SetBranchAddress("cosths",&cosths);
    hitRate_pmtType1->SetBranchAddress("nHits",&nHits);
    hitRate_pmtType1->SetBranchAddress("nPE",&nPE);
    hitRate_pmtType1->SetBranchAddress("timetof",&timetof);
    hitRate_pmtType1->SetBranchAddress("mPMT_id",&mPMT_id);

    for (ULong64_t i=0;i<hitRate_pmtType1->GetEntries();i++) {
        hitRate_pmtType1->GetEntry(i);
        if (nmPMT_on>0)
            if (mPMT_mask[mPMT_id]==1) continue;
        if (cosths>cosths_min&&timetof>timetof_min&&timetof<timetof_max) // only hits within the decided time window and the source opening angle  
        {
            double weight = nPE;
            hRate1->Fill(mPMT_id+0.5,weight);
            hBinnedRate->Fill(-costh,dist,weight);
        }
    }

    TCanvas* c1 = new TCanvas();
    hBinnedRate->GetXaxis()->SetTitle("cos(#theta_{PMT})");
    hBinnedRate->GetYaxis()->SetTitle("R (cm)");
    hBinnedRate->Draw("colz");
    c1->SaveAs(Form("R_theta_BinnedRate_%i.pdf",nmPMT_on));
    hPMT1->GetXaxis()->SetTitle("cos(#theta_{PMT})");
    hPMT1->GetYaxis()->SetTitle("R (cm)");
    c1->SaveAs(Form("PMT_R_theta_map_%i.pdf",nmPMT_on));

    int nonzerobins=0;
    for (int i=1;i<=hBinnedRate->GetNbinsX();i++)
        for (int j=1;j<=hBinnedRate->GetNbinsY();j++)
            if (hBinnedRate->GetBinContent(i,j)>0) nonzerobins++;

    run_fit();
    std::cout<<"Number of mPMT_used = "<<nmPMT_used<<std::endl;
    std::cout<<"Number of non-zero bins = "<<nonzerobins<<std::endl;
}

void fit_water_attenuation(){

    // TChain is used to load a number of files at the same time
    std::string filename = "/bundle/data/T2K/users/kmtsui/LI/diffuser32_350nm_x2scattering/absorption_diffuser*.root";
    fit_all(filename,0,-952,-940);
    double alpha = truth_alpha(350,1.3,1.5);

}