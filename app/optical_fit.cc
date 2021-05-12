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
    auto s = new AnaSample(0, "Sample mPMT", bm, 1);
    samples.push_back(s);

    //read events
    AnaTree selTree(fname_input.c_str(), "hitRate_pmtType1", "pmt_type1");
    selTree.SetTimetofCut(-950,-940);
    selTree.SetCosthsCut(0.766,1);
    selTree.GetEvents(samples.at(0));

    samples.at(0)->InitEventMap();

    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");

    Fitter fitter(fout,0,8);

    std::vector<FitParameter> fitpara;
    std::string parname[6]={"alpha","norm1","norm2","norm3","norm4","norm5"};
    int par_type[6]={0,1,1,1,1,1}; // 0: attenuation length, 1: angular response
    int pmt_type[6]={-1,1,1,1,1,1}; // 0: apply to B&L PMTs, 1: apply to mPMTs, -1: apply to both
    std::string var[6]={"R","costh","costh","costh","costh","costh"};
    double var_low[6]={0,0.5,0.6,0.7,0.8,0.9};
    double var_high[6]={1.e4,0.6,0.7,0.8,0.9,1.0};
    double prior[6]={1.e4,2.,2.,2.,2.,2};
    double step[6]={1.e4,1.,1.,1.,1.,1};
    double low[6]={0,0,0,0,0,0,};
    double high[6]={1.e5,1.e2,1.e2,1.e2,1.e2,1.e2};
    int random[6]={0,0,0,0,0,0};
    bool fixed[6]={false,false,false,false,false,false};
    for (int i=0;i<6;i++) {
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