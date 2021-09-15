void makeLEDProfileInput_1D()
{
    const int nBins1=40;
    double bins1[nBins1+1]=
    {
      0.76655, 0.77700, 0.78761, 0.79875, 0.80872, 0.81934, 0.82883, 0.83808,
      0.84711, 0.85590, 0.86523, 0.87426, 0.88229, 0.89008, 0.89831, 0.90623,
      0.91259, 0.92055, 0.92700, 0.93319, 0.93911, 0.94477, 0.95063, 0.95618,
      0.96098, 0.96551, 0.96976, 0.97408, 0.97808, 0.98116, 0.98454, 0.98736,
      0.98989, 0.99233, 0.99444, 0.99607, 0.99742, 0.99857, 0.99933, 0.99981,
      1.00000
    };

    double Power_cosths[nBins1+1] =
    {
      0.50152, 0.54801, 0.59419, 0.63473, 0.67604, 0.71697, 0.75712, 0.80037,
      0.83589, 0.86910, 0.90307, 0.93628, 0.95584, 0.97296, 0.98300, 0.98918,
      0.99844, 1.00063, 1.00848, 1.01041, 1.01466, 1.01608, 1.01813, 1.02019,
      1.02200, 1.02122, 1.02277, 1.02084, 1.02007, 1.01813, 1.01698, 1.01402,
      1.01427, 1.01273, 1.01080, 1.00938, 1.00810, 1.00694, 1.00578, 1.00501,
      1.00308
    };

    TFile* fout = new TFile("led_prior_1D.root","RECREATE");
    TH1D* h_prior = new TH1D("h_prior","h_prior",nBins1,0,nBins1);
    for (int i=0;i<nBins1;i++)
        {
            double power = (Power_cosths[i]+Power_cosths[i+1])/2.;
            h_prior->SetBinContent(i+1,power);
        }
    h_prior->Write();

    TMatrixDSym cov_matrix(nBins1);
    for(int r = 0; r < nBins1; ++r)
    {
        for(int c = 0; c < nBins1; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.01,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_1percent");

    for(int r = 0; r < nBins1; ++r)
    {
        for(int c = 0; c < nBins1; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.05,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_5percent");

    for(int r = 0; r < nBins1; ++r)
    {
        for(int c = 0; c < nBins1; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.1,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_10percent");

    fout->Close();
}

void makeLEDProfileInput_2D()
{
    const int nBins1=40;
    double bins1[nBins1+1]=
    {
      0.76655, 0.77700, 0.78761, 0.79875, 0.80872, 0.81934, 0.82883, 0.83808,
      0.84711, 0.85590, 0.86523, 0.87426, 0.88229, 0.89008, 0.89831, 0.90623,
      0.91259, 0.92055, 0.92700, 0.93319, 0.93911, 0.94477, 0.95063, 0.95618,
      0.96098, 0.96551, 0.96976, 0.97408, 0.97808, 0.98116, 0.98454, 0.98736,
      0.98989, 0.99233, 0.99444, 0.99607, 0.99742, 0.99857, 0.99933, 0.99981,
      1.00000
    };

    double Power_cosths[nBins1+1] =
    {
      0.50152, 0.54801, 0.59419, 0.63473, 0.67604, 0.71697, 0.75712, 0.80037,
      0.83589, 0.86910, 0.90307, 0.93628, 0.95584, 0.97296, 0.98300, 0.98918,
      0.99844, 1.00063, 1.00848, 1.01041, 1.01466, 1.01608, 1.01813, 1.02019,
      1.02200, 1.02122, 1.02277, 1.02084, 1.02007, 1.01813, 1.01698, 1.01402,
      1.01427, 1.01273, 1.01080, 1.00938, 1.00810, 1.00694, 1.00578, 1.00501,
      1.00308
    };

    const int nBins2=12;
    double bins2[nBins2];
    for (int i=0;i<=nBins2;i++) bins2[i]=-TMath::Pi()+i*2.*TMath::Pi()/nBins2;

    double Variation_cosths[nBins1+1] =
    {
      1.52148, 1.32716, 1.60492, 1.51191, 1.32437, 1.07736, 0.80404, 0.84281,
      0.89956, 0.97436, 0.87347, 0.95811, 1.03983, 1.16010, 1.18226, 1.17236,
      1.14083, 1.09142, 1.10767, 1.02450, 0.96195, 0.91277, 0.88456, 0.85350,
      0.82181, 0.78467, 0.74606, 0.67583, 0.60713, 0.54087, 0.50040, 0.42365,
      0.37995, 0.32521, 0.27719, 0.23896, 0.19595, 0.17128, 0.13077, 0.08523,
      0.06337
    };

    TFile* fout = new TFile("led_prior_2D.root","RECREATE");
    TH1D* h_prior = new TH1D("h_prior","h_prior",nBins1*nBins2,0,nBins1*nBins2);
    for (int i=0;i<nBins1;i++)
        for (int j=0;j<nBins2;j++)
        {
            double phis = (bins2[j]+bins2[j+1])/2.;
            double power = (Power_cosths[i]+Power_cosths[i+1])/2.;
            double variation = (Variation_cosths[i]+Variation_cosths[i+1])/2.;
            double val = power*(1.+variation/100.*cos(phis));
            h_prior->SetBinContent(i*nBins2+j+1,val);
        }
    h_prior->Write();

    int ndims = nBins1*nBins2;
    TMatrixDSym cov_matrix(ndims);
    for(int r = 0; r < ndims; ++r)
    {
        for(int c = 0; c < ndims; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.01,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_1percent");

    for(int r = 0; r < ndims; ++r)
    {
        for(int c = 0; c < ndims; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.05,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_5percent");

    for(int r = 0; r < ndims; ++r)
    {
        for(int c = 0; c < ndims; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_prior->GetBinContent(r+1)*0.1,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_10percent");

    fout->Close();
}

void makeLEDProfileInput()
{
    makeLEDProfileInput_1D();
}