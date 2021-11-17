// Use with WCSIM_TreeConvert output with -d options to produce a maping of indirect light background
void make_scatter_map(std::string filename = "/bundle/data/T2K/users/kmtsui/LI/production/out_diffuser4_400nm_nominal_*.root", int pmttype = 0){

    std::string chainname = pmttype == 0 ? "hitRate_pmtType0" : "hitRate_pmtType1";
    std::string treename = pmttype == 0 ? "pmt_type0" : "pmt_type1";
    double timetof_indirect = pmttype == 0 ? 1.0 : 0.6 ;
    double timetof_cut1 = pmttype == 0 ? 15 : 2 ;
    double timetof_cut2 = 200;

    std::string outname = pmttype == 0 ? "scattering_map_BnL_diffuser4_400nm.root" : "scattering_map_mPMT_diffuser4_400nm.root";

    TChain* fChain_digitized = new TChain(chainname.c_str());
    fChain_digitized->Add(filename.c_str());

    TFile* f = new TFile(fChain_digitized->GetFile()->GetName());

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
    fChain_digitized->SetBranchAddress("nPE",&nPE);
    fChain_digitized->SetBranchAddress("timetof",&timetof);
    fChain_digitized->SetBranchAddress("PMT_id",&PMT_id);
    fChain_digitized->SetBranchAddress("nRaySct",&nRaySct);
    fChain_digitized->SetBranchAddress("nReflec",&nReflec);
    fChain_digitized->SetBranchAddress("nPE_digi",&nPE_digi);
    fChain_digitized->SetBranchAddress("timetof_digi",&timetof_digi);
    TH1D* h_signal_region = new TH1D("","",nPMTs,0,nPMTs);
    TH1D* h_control_region = new TH1D("","",nPMTs,0,nPMTs);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        fChain_digitized->GetEntry(i);
        if (nPE_digi<=0) continue;
        if (timetof>timetof_indirect) // use timetof to get indirect photon
        {
            if (timetof_digi<timetof_cut1) h_signal_region->Fill(PMT_id+0.5,nPE_digi);
            else if (timetof_digi<timetof_cut2) h_control_region->Fill(PMT_id+0.5,nPE_digi);
        }
    }

    TFile* fout = new TFile(outname.c_str(),"RECREATE");
    for (int i=1;i<=h_signal_region->GetNbinsX();i++)
    {
        double x = h_control_region->GetBinContent(i);
        double y = h_signal_region->GetBinContent(i);
        if (x>0&&y>0)
        {
            double val = y/x;
            double err = sqrt(1./x+1./y);
            h_signal_region->SetBinContent(i,val);
            h_signal_region->SetBinError(i,err);
        }
        else
        {
            h_signal_region->SetBinContent(i,0);
            h_signal_region->SetBinError(i,0);
        }
    }
    h_signal_region->Write("scattering_map");
    fout->Close();
}