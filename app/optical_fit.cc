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
#include "OPTICALFIT/AnaFitParameters.hh"

#include "toml/toml_helper.h"

int main(int argc, char** argv)
{
    std::string fname_input;
    std::string fname_output;
    std::string config_file;
    int num_threads = 1;

    char option;
    while((option = getopt(argc, argv, "j:f:o:c:n:s:t:h")) != -1)
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
            case 'n':
                num_threads = std::stoi(optarg);
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
    auto const &fitparameters_config = toml_h::find(card_toml, "fitparameters");

    // Add analysis samples:
    std::vector<AnaSample*> samples;
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

    //Add fit parameters
    std::vector<AnaFitParameters*> fitparas;
    for (auto const &name : toml_h::find<std::vector<std::string>>(fitparameters_config, "names"))
    {
        std::cout << "Parameter name: " << name << std::endl;
        auto const &ele = toml_h::find<toml::array>(fitparameters_config, name);
        auto pmttype = toml_h::find<int>(ele,0);
        auto functype = toml_h::find<std::string>(ele,1);
        auto npar = toml_h::find<int>(ele,2);
        
        auto par_setup = toml_h::find<toml::array>(ele,3);
        std::vector<std::string> parnames;
        std::vector<double> priors;
        std::vector<double> steps;
        std::vector<double> lows;
        std::vector<double> highs;
        std::vector<bool> fixeds;
        for (auto const &par : par_setup){
            auto parname = toml_h::find<std::string>(par,0);
            auto prior = toml_h::find<double>(par,1);
            auto step = toml_h::find<double>(par,2);
            auto low = toml_h::find<double>(par,3);
            auto high = toml_h::find<double>(par,4);
            auto fixed = toml_h::find<bool>(par,5);
            parnames.push_back(parname);
            priors.push_back(prior);
            steps.push_back(step);
            lows.push_back(low);
            highs.push_back(high);
            fixeds.push_back(fixed);
        }
        if (parnames.size()<npar)
        {   int lastidx = parnames.size()-1;
            for (int i=parnames.size();i<npar;i++)
            {
                parnames.push_back(Form("%s_%i",parnames[lastidx].c_str(),i));
                priors.push_back(priors[lastidx]);
                steps.push_back(steps[lastidx]);
                lows.push_back(lows[lastidx]);
                highs.push_back(highs[lastidx]);
                fixeds.push_back(fixeds[lastidx]);
            }
        }
        for (int i=0;i<npar;i++)
            std::cout<<parnames[i]<<" "<<priors[i]<<" "<<steps[i]<<" "<<lows[i]<<" "<<highs[i]<<" "<<fixeds[i]<<" "<<std::endl;

        auto binning = toml_h::find<toml::array>(ele,4);
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

        auto fitpara = new AnaFitParameters(name,pmttype);
        fitpara->InitParameters(parnames,priors,steps,lows,highs,fixeds);
        fitpara->SetParameterFunction(functype);
        fitpara->SetBinVar(binning_var);
        fitpara->SetBinning(binning_file);
        fitpara->InitEventMap(samples);

        fitpara->SetNumThreads(num_threads);

        fitparas.push_back(fitpara);
    }

    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");

    Fitter fitter(fout,0,num_threads);

    fitter.InitFitter(fitparas);

    bool did_converge = false;

    did_converge = fitter.Fit(samples, false);
    if(!did_converge)
        std::cout << "Fit did not coverge." << std::endl;


    fout->Close();

    return 0;
}    