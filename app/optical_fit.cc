#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "OPTICALFIT/AnaSample.hh"
#include "OPTICALFIT/AnaTree.hh"
#include "OPTICALFIT/BinManager.hh"
#include "OPTICALFIT/Fitter.hh"

int main(int argc, char** argv)
{
    std::string fname_input;
    std::string fname_output;

    char option;
    while((option = getopt(argc, argv, "j:f:o:s:t:nh")) != -1)
    {
        switch(option)
        {
            case 'f':
                fname_input = optarg;
                break;
            case 'o':
                fname_output = optarg;
                break;
            case 'h':
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
                          << "-f : Input file\n"
                          << "-o : Output file\n";
            default:
                return 0;
        }
    }

    // Add analysis samples:
    std::vector<AnaSample*> samples;
    BinManager* bm = new BinManager(10,0.5,1,10,1000,9000);
    auto s0 = new AnaSample(0, "B&L_PMT", bm, 0);
    s0->SetCut("timetof",-950,-940);
    s0->SetCut("cosths",0.766,1);
    s0->SetCut("costh",0.5,1);
    s0->LoadEventsFromFile(fname_input,"hitRate_pmtType0", "pmt_type0");
    samples.push_back(s0);
    auto s1 = new AnaSample(1, "mPMT", bm, 1);
    s1->SetCut("timetof",-950,-940);
    s1->SetCut("cosths",0.766,1);
    s1->SetCut("costh",0.5,1);
    s1->LoadEventsFromFile(fname_input,"hitRate_pmtType1", "pmt_type1");
    samples.push_back(s1);

    samples.at(0)->InitEventMap();
    samples.at(1)->InitEventMap();

    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");

    Fitter fitter(fout,0,8);

    std::vector<FitParameter> fitpara;
    const int npar = 11;
    std::string parname[npar]={"alpha","norm0_1","norm0_2","norm0_3","norm0_4","norm0_5","norm1_1","norm1_2","norm1_3","norm1_4","norm1_5"};
    int par_type[npar]={0,1,1,1,1,1,1,1,1,1,1}; // 0: attenuation length, 1: angular response
    int pmt_type[npar]={-1,0,0,0,0,0,1,1,1,1,1}; // 0: apply to B&L PMTs, 1: apply to mPMTs, -1: apply to both
    std::string var[npar]={"R","costh","costh","costh","costh","costh","costh","costh","costh","costh","costh"};
    double var_low[npar]={0,0.5,0.6,0.7,0.8,0.9,0.5,0.6,0.7,0.8,0.9};
    double var_high[npar]={1.e4,0.6,0.7,0.8,0.9,1.0,0.6,0.7,0.8,0.9,1.0};
    double prior[npar]={1.e4,100,100,100,100,100,2.,2.,2.,2.,2};
    double step[npar]={1.e4,100,100,100,100,100,1.,1.,1.,1.,1};
    double low[npar]={0,0,0,0,0,0,0,0,0,0,0};
    double high[npar]={1.e5,1.e3,1.e3,1.e3,1.e3,1.e3,1.e3,1.e3,1.e3,1.e3,1.e3};
    int random[npar]={0,0,0,0,0,0,0,0,0,0,0};
    bool fixed[npar]={false,false,false,false,false,false,false,false,false,false,false};
    for (int i=0;i<npar;i++) {
        FitParameter fp;
        fp.name = parname[i];
        fp.par_type=par_type[i];
        fp.pmt_type=pmt_type[i];
        fp.var=var[i];
        fp.var_low=var_low[i];
        fp.var_high=var_high[i];
        fp.prior=prior[i];
        fp.step=step[i];
        fp.low=low[i];
        fp.high=high[i];
        fp.random=random[i];
        fp.fixed=fixed[i];
        fitpara.push_back(fp);
    }
    fitter.InitFitter(fitpara);

    bool did_converge = false;

    did_converge = fitter.Fit(samples, false);
    if(!did_converge)
        std::cout << "Fit did not coverge." << std::endl;


    fout->Close();

    return 0;
}    