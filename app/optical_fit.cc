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

#include "toml/toml_helper.h"

int main(int argc, char** argv)
{
    std::string fname_input;
    std::string fname_output;
    std::string config_file;

    char option;
    while((option = getopt(argc, argv, "j:f:o:c:s:t:nh")) != -1)
    {
        switch(option)
        {
            case 'f':
                fname_input = optarg;
                break;
            case 'o':
                fname_output = optarg;
                break;
            case 'c':
                config_file = optarg;
                break;
            case 'h':
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
                          << "-f : Input file\n"
                          << "-o : Output file\n"
                          << "-c : Config file\n";
            default:
                return 0;
        }
    }

    auto const &card_toml = toml_h::parse_card(config_file);
    auto const &samples_config = toml_h::find(card_toml, "samples");

    // Add analysis samples:
    std::vector<AnaSample*> samples;
    BinManager* bm = new BinManager(10,0.5,1,10,1000,9000);

    for (auto const &name : toml_h::find<std::vector<std::string>>(samples_config, "names"))
    {
        std::cout << "Sample name: " << name << std::endl;
        auto const &ele = toml_h::find<toml::array>(samples_config, name);
        auto pmttype = toml_h::find<int>(ele,0);
        auto filename = toml_h::find<std::string>(ele,1);
        auto ev_tree_name = toml_h::find<std::string>(ele,2);
        auto pmt_treename = toml_h::find<std::string>(ele,3);
        std::cout << " pmttype =  "<<pmttype<<", filename =  "<<filename<<", ev_tree_name =  "<<ev_tree_name<<", pmt_treename =  "<<pmt_treename<<std::endl;

        auto binning = toml_h::find<toml::array>(ele,5);
        auto binning_file = toml_h::find<std::string>(binning,0);
        if (binning_file.find_first_of("/")!=0) {
            std::size_t found = config_file.find_last_of("/");
            if (found!= std::string::npos) {
                std::string prefix = config_file.substr(0,found+1);
                binning_file = prefix+binning_file;
            }
        }
        auto binning_var = toml_h::find<std::vector<std::string>>(binning,1);
        std::cout<< "binning_file = "<<binning_file<<", binning_var = [ ";
        for (auto t : binning_var) std::cout<< t <<", ";
        std::cout<<"]"<<std::endl;

        auto s = new AnaSample(samples.size(), name , binning_file, pmttype);
        s->SetBinVar(binning_var);

        for (auto const &cut : toml_h::find<toml::array>(ele,4)){
            auto cutvar = toml_h::find<std::string>(cut,0);
            auto cutlow = toml_h::find<double>(cut,1);
            auto cuthigh = toml_h::find<double>(cut,2);
            std::cout << " cutvar =  "<<cutvar<<", cutlow =  "<<cutlow<<", cuthigh =  "<<cuthigh<<std::endl;
            s->SetCut(cutvar,cutlow,cuthigh);
        }

        s->LoadEventsFromFile(filename,ev_tree_name,pmt_treename);
        s->InitEventMap();
        samples.push_back(s);
    }

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