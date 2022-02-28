#include "../src/utils/CalcGroupVelocity.hh"
#include "truth_alpha.hh"

// Input: Processed files from WCSIM_TreeConvert  with -d option
// Output: template of indirect light for each PMT
Use with WCSIM_TreeConvert output with -d options to template a maping of indirect light
void build_template_and_spline(
    int pmttype = 0, double wavelength = 400,
    std::string nominalMC = "/bundle/data/T2K/users/kmtsui/LI/collimator32_400nm_nominal/out_collimator32_nominal_*.root", 
    double ABWFF=1.30, double RAYFF=0.75,
    std::vector<std::string> sctMC = std::vector<std::string>{"/bundle/data/T2K/users/kmtsui/LI/collimator32_400nm_x2scattering/out_collimator32_x2scattering_*.root"},
    std::vector<double> RAYFF_sct = std::vector<double>{1.5}, std::vector<double> ABWFF_sct = std::vector<double>{1.3}
    )
{
    // Define output according to PMT type
    std::string outname = pmttype == 0 ? "template_spline_BnL_collimator32_400nm.root" : "template2d_mPMT_diffuser4_400nm.root";
    std::string chainname = pmttype == 0 ? "hitRate_pmtType0" : "hitRate_pmtType1";
    std::string treename = pmttype == 0 ? "pmt_type0" : "pmt_type1";
    double timetof_indirect = pmttype == 0 ? 0.0 : 0.6 ;
    const int nbins = 5;
    double timetof_lo = pmttype == 0 ? 20 : 3 ;
    double timetof_hi = 200;

    double speed = CalcGroupVelocity(wavelength)*100./1.e9;
    std::cout<<"speed = "<<speed<<" cm/ns"<<std::endl;

    std::vector<double> lengths = truth_alpha(wavelength,ABWFF,RAYFF);
    double sct_length = lengths[2];
    double att_length = lengths[0];

    TChain* fChain_nominal = new TChain(chainname.c_str());
    fChain_nominal->Add(nominalMC.c_str());

    TFile* f = new TFile(fChain_nominal->GetFile()->GetName());

    TTree* t = (TTree*)f->Get(treename.c_str());
    double R, costh, cosths, omega, phim;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    const int nPMTs = t->GetEntries();
    std::cout<<"nPMTs = "<<nPMTs<<std::endl;
    std::vector<double> cosths_array, costh_array, R_array;

    for (int i=0;i<nPMTs;i++){
         t->GetEntry(i);
         cosths_array.push_back(cosths);
         costh_array.push_back(costh);
         R_array.push_back(R);
    }

    double nPE, timetof, nPE_digi, timetof_digi;
    int PMT_id, nRaySct, nReflec;
    fChain_nominal->SetBranchAddress("nPE",&nPE);
    fChain_nominal->SetBranchAddress("timetof",&timetof);
    fChain_nominal->SetBranchAddress("PMT_id",&PMT_id);
    fChain_nominal->SetBranchAddress("nRaySct",&nRaySct);
    fChain_nominal->SetBranchAddress("nReflec",&nReflec);
    fChain_nominal->SetBranchAddress("nPE_digi",&nPE_digi);
    fChain_nominal->SetBranchAddress("timetof_digi",&timetof_digi);

    const int nPoints_att = 11; const int nomIdx_att = 5; double step_att = 0.1;
    TH2D* hist_att_template[nPoints_att];
    std::cout<<"Initializing Latt template..."<<std::endl;
    for (int i=0;i<nPoints_att;i++)
    {
        hist_att_template[i] = new TH2D("","",nPMTs,0,nPMTs,nbins,timetof_lo,timetof_hi);

        hist_att_template[i]->Sumw2();
    }
    // for attenuation length, a simple reweight is good enough for spline generation
    for (unsigned long int i=0;i<fChain_nominal->GetEntries();i++){
        fChain_nominal->GetEntry(i);
        if (i%1000000==0) std::cout<<"Processing "<<i<<" out of "<<fChain_nominal->GetEntries()<<" events"<<std::endl;
        if (nPE_digi<=0) continue;
        if (timetof>timetof_indirect) // use timetof to get indirect photon
        {
            double dist = R_array[PMT_id] + timetof*speed;
            for (int j=0;j<nPoints_att;j++)
            {
                double att_length_new = att_length*((j-nomIdx_att)*step_att+1.);
                double reweight = exp(-dist/att_length_new)/exp(-dist/att_length);
                hist_att_template[j]->Fill(PMT_id+0.5,timetof_digi,nPE_digi*reweight);
            }
        }
    }

    // for scattering length, need separate MC to produce spline 
    const int nPoints_sct = sctMC.size();
    TH2D* hist_sct_template[nPoints_sct];
    for (int i=0;i<nPoints_sct;i++)
    {
        hist_sct_template[i] = new TH2D("","",nPMTs,0,nPMTs,nbins,timetof_lo,timetof_hi);

        hist_sct_template[i]->Sumw2();
    }
    for (int i=0;i<nPoints_sct;i++)
    {
        TChain* tIn = new TChain(chainname.c_str());
        tIn->Add(sctMC[i].c_str());
        tIn->SetBranchAddress("nPE",&nPE);
        tIn->SetBranchAddress("timetof",&timetof);
        tIn->SetBranchAddress("PMT_id",&PMT_id);
        tIn->SetBranchAddress("nRaySct",&nRaySct);
        tIn->SetBranchAddress("nReflec",&nReflec);
        tIn->SetBranchAddress("nPE_digi",&nPE_digi);
        tIn->SetBranchAddress("timetof_digi",&timetof_digi);

        for (unsigned long int j=0;j<tIn->GetEntries();j++)
        {
            tIn->GetEntry(j);
            if (j%1000000==0) std::cout<<"Processing "<<j<<" out of "<<tIn->GetEntries()<<" events"<<std::endl;
            if (nPE_digi<=0) continue;
            if (timetof>timetof_indirect) // use timetof to get indirect photon
            {
                double dist = R_array[PMT_id] + (timetof)*speed;
                double att_length_new = truth_alpha(wavelength,ABWFF_sct[i],RAYFF_sct[i])[0];
                double reweight = exp(-dist/att_length_new)/exp(-dist/att_length);
                hist_sct_template[i]->Fill(PMT_id+0.5,timetof_digi,nPE_digi/reweight);
            }
        }
    }

    TFile* fout = new TFile(outname.c_str(),"RECREATE");
    hist_att_template[nomIdx_att]->Write("timetof_template");
    TDirectory *spline_att_dir = fout->mkdir("Spline_Att");
    spline_att_dir->cd();    // make the "tof" directory the current directory
    for (int p=0;p<nPMTs;p++)
    {
        TDirectory *pmt_dir = spline_att_dir->mkdir(Form("PMT%i",p));
        pmt_dir->cd();
        for (int i=1;i<=nbins;i++)
        {
            char nameHisto[256];
            //sprintf(nameHisto,"Spline_Att_PMT%i_Bin_%i",p,i);
            sprintf(nameHisto,"Bin_%i",i);
            TGraph* g = new TGraph(nPoints_att);
            g->SetName(nameHisto);
            g->SetTitle(nameHisto);
            g->SetMarkerStyle(20);
            g->SetMarkerColor(2);
            for (int j=0;j<nPoints_att;j++)
            {
                double x = att_length*((j-nomIdx_att)*step_att+1.);
                double y = hist_att_template[nomIdx_att]->GetBinContent(p+1,i)>0 ? hist_att_template[j]->GetBinContent(p+1,i)/hist_att_template[nomIdx_att]->GetBinContent(p+1,i) : 1;
                g->SetPoint(j,x,y);
            }
            g->Write();
            delete g;
        }

    }

    std::vector<double> sortSct (RAYFF_sct.begin(),RAYFF_sct.end());
    sortSct.push_back(RAYFF);
    std::sort(sortSct.begin(),sortSct.end());
    TDirectory *spline_sct_dir = fout->mkdir("Spline_Sct");
    spline_sct_dir->cd();    // make the "tof" directory the current directory
    for (int p=0;p<nPMTs;p++)
    {
        TDirectory *pmt_dir = spline_sct_dir->mkdir(Form("PMT%i",p));
        pmt_dir->cd();
        for (int i=1;i<=nbins;i++)
        {
            char nameHisto[256];
            //sprintf(nameHisto,"Spline_Att_PMT%i_Bin_%i",p,i);
            sprintf(nameHisto,"Bin_%i",i);
            TGraph* g = new TGraph(nPoints_sct+1);
            g->SetName(nameHisto);
            g->SetTitle(nameHisto);
            g->SetMarkerStyle(20);
            g->SetMarkerColor(2);
            for (int j=0;j<nPoints_sct+1;j++)
            {
                double x = sct_length*sortSct[j]/RAYFF;
                double y = 1;
                for (int k=0;k<nPoints_sct;k++)
                {
                    if (fabs(RAYFF_sct[k]-sortSct[j])<1.e-6)
                    {
                        y = hist_att_template[nomIdx_att]->GetBinContent(p+1,i)>0 ? hist_sct_template[k]->GetBinContent(p+1,i)/hist_att_template[nomIdx_att]->GetBinContent(p+1,i) : 1;
                        break;
                    }
                }
                g->SetPoint(j,x,y);
            }
            g->Write();
            delete g;
        }

    }

    fout->Close();
}