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
    std::string fname_output;
    std::string config_file;
    int num_threads = 1;
    int seed = 0;

    char option;
    while((option = getopt(argc, argv, "j:f:o:c:n:s:t:h")) != -1)
    {
        switch(option)
        {
            case 'o':
                fname_output = optarg;
                break;
            case 'c':
                config_file = optarg;
                break;
            case 'n':
                num_threads = std::stoi(optarg);
                if (num_threads<1) num_threads = 1;
                break;
            case 's':
                seed = std::stoi(optarg);
                if (seed<0) seed = 0;
                break;
            case 'h':
                std::cout << "USAGE: "
                          << argv[0] << "\nOPTIONS:\n"
                          << "-o : Output file\n"
                          << "-c : Config file\n"
                          << "-s : RNG seed \n"
                          << "-n : Number of threads\n";
            default:
                return 0;
        }
    }

    auto const &card_toml = toml_h::parse_card(config_file);
    auto const &samples_config = toml_h::find(card_toml, "samples");
    auto const &fitparameters_config = toml_h::find(card_toml, "fitparameters");
    auto const &minimizer_config = toml_h::find(card_toml, "Minimizer");

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
        if (binning_file.find_first_of("/")!=0) 
        {  //assume the binning file lives in the same directory of config_file
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
        if (binning_file.find_first_of("/")!=0) 
        {   //assume the binning file lives in the same directory of config_file
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

        fitparas.push_back(fitpara);
    }

    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");

    // Load minimizer config
    MinSettings min_settings;
    min_settings.minimizer = toml_h::find<std::string>(minimizer_config, "minimizer");
    min_settings.algorithm = toml_h::find<std::string>(minimizer_config, "algorithm");
    min_settings.likelihood = toml_h::find<std::string>(minimizer_config, "likelihood");
    min_settings.print_level = toml_h::find<int>(minimizer_config, "print_level");
    min_settings.strategy = toml_h::find<int>(minimizer_config, "strategy");
    min_settings.tolerance = toml_h::find<double>(minimizer_config, "tolerance");
    min_settings.max_iter = toml_h::find<double>(minimizer_config, "max_iter");
    min_settings.max_fcn = toml_h::find<double>(minimizer_config, "max_fcn");

    Fitter fitter(fout,seed,num_threads);
    fitter.SetMinSettings(min_settings);
    fitter.InitFitter(fitparas);

    bool did_converge = false;

    did_converge = fitter.Fit(samples);
    if(!did_converge)
        std::cout << "Fit did not coverge." << std::endl;


    fout->Close();

    return 0;
}    