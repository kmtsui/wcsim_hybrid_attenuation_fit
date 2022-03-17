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
#include "OPTICALFIT/ColorOutput.hh"

#include "toml/toml_helper.h"

const std::string TAG = color::GREEN_STR + "[optical_fit]: " + color::RESET_STR;
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::RED_STR + "[WARNING]: " + color::RESET_STR;

void HelpMessage()
{
    std::cout   << TAG << "USAGE: "
                << "optical_fit" << "\nOPTIONS:\n"
                << "-o : Output file\n"
                << "-c : Config file\n"
                << "-s : RNG seed \n"
                << "-n : Number of threads\n"
                << "-t : Number of toy fits \n";
}

int main(int argc, char** argv)
{
    std::string fname_output = "fitoutput.root";
    std::string config_file;
    int num_threads = 1;
    int seed = 0;
    int toys = 0;

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
                std::cout << TAG<<"Set seed = "<<seed<<std::endl;
                break;
            case 't':
                toys = std::stoi(optarg);
                if (toys<0) toys = 0;
                std::cout << TAG<<"Number of toys = "<<toys<<std::endl;
                break;
            case 'h':
                HelpMessage();
            default:
                return 0;
        }
    }

    if (config_file.size()==0)
    {
        std::cout<< ERR << "No config file!" << std::endl;
        HelpMessage();
        return -1;
    }

    TFile* fout = TFile::Open(fname_output.c_str(), "RECREATE");
    Fitter fitter(fout,seed,num_threads);

    auto const &card_toml = toml_h::parse_card(config_file);
    auto const &samples_config = toml_h::find(card_toml, "samples");
    auto const &fitparameters_config = toml_h::find(card_toml, "fitparameters");
    auto const &minimizer_config = toml_h::find(card_toml, "Minimizer");

    // Add analysis samples:
    std::vector<AnaSample*> samples;
    for (auto const &name : toml_h::find<std::vector<std::string>>(samples_config, "names"))
    {
        std::cout << TAG << "Sample name: " << name << std::endl;
        auto const &ele = toml_h::find<toml::array>(samples_config, name);
        auto pmttype = toml_h::find<int>(ele,0);
        auto filename = toml_h::find<std::string>(ele,1);
        auto ev_tree_name = toml_h::find<std::string>(ele,2);
        auto pmt_treename = toml_h::find<std::string>(ele,3);
        std::cout << TAG << " pmttype =  "<<pmttype<<", filename =  "<<filename<<", ev_tree_name =  "<<ev_tree_name<<", pmt_treename =  "<<pmt_treename<<std::endl;

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
        std::cout << TAG<< "binning_file = "<<binning_file<<", binning_var = [ ";
        for (auto t : binning_var) std::cout << t <<", ";
        std::cout <<"]"<<std::endl;

        auto s = new AnaSample(samples.size(), name , binning_file, pmttype);
        s->SetBinVar(binning_var);

        for (auto const &cut : toml_h::find<toml::array>(ele,4))
        {
            auto cutvar = toml_h::find<std::string>(cut,0);
            auto cutlow = toml_h::find<double>(cut,1);
            auto cuthigh = toml_h::find<double>(cut,2);
            std::cout << TAG << " cutvar =  "<<cutvar<<", cutlow =  "<<cutlow<<", cuthigh =  "<<cuthigh<<std::endl;
            s->SetCut(cutvar,cutlow,cuthigh);
        }

        if (ele.size()>6)
        {
            for (int i=6; i<ele.size();i++ )
            {
                auto opt = toml_h::find<toml::array>(ele,i);
                auto optname = toml_h::find<std::string>(opt,0);
                if (optname=="norm")
                {
                    auto norm = toml_h::find<double>(opt,1);
                    std::cout << TAG<<"Setting norm = "<<norm<<std::endl;
                    s->SetNorm(norm);
                } 
                else if (optname=="mask")
                {
                    auto mask = toml_h::find<int>(opt,1);
                    std::cout << TAG<<"Masking to use only "<< mask << " PMTs" <<std::endl;
                    s->MaskPMT(mask);
                    if (pmttype==1)
                    {
                        auto nPMTpermPMT = toml_h::find<int>(samples_config, "nPMTpermPMT");
                        std::cout << TAG << "Set nPMTpermPMT =  "<< nPMTpermPMT <<std::endl;
                        s->SetnPMTpermPMT(nPMTpermPMT);
                    }
                }
                else if (optname=="mask_mPMT")
                {
                    if (pmttype==1)
                    {
                        auto mask = toml_h::find<std::vector<int>>(opt,1);
                        std::cout << TAG << "Masking the small PMTs inside mPMT: ";
                        for (auto m : mask) std::cout << m <<" ";
                        std::cout << std::endl;
                        auto nPMTpermPMT = toml_h::find<int>(samples_config, "nPMTpermPMT");
                        std::cout << TAG << "Set nPMTpermPMT =  "<< nPMTpermPMT <<std::endl;
                        s->MaskmPMT(mask);
                        s->SetnPMTpermPMT(nPMTpermPMT);
                    }
                }
                else if (optname=="scatter_control")
                {
                    auto time1 = toml_h::find<double>(opt,1);
                    auto time2 = toml_h::find<double>(opt,2);
                    auto time3 = toml_h::find<double>(opt,3);
                    auto factor = toml_h::find<double>(opt,4);
                    std::cout << TAG<<"Fitting with scattering control region. timetof region = ["<<time1<<","<<time2<<","<<time3<<"], scaling factor = "<< factor <<std::endl;
                    s->SetScatter(time1,time2,time3,factor);
                }
                else if (optname=="scatter_map")
                {
                    auto fname = toml_h::find<std::string>(opt,1);
                    auto hname = toml_h::find<std::string>(opt,2);
                    auto time1 = toml_h::find<double>(opt,3);
                    auto time2 = toml_h::find<double>(opt,4);
                    auto time3 = toml_h::find<double>(opt,5);
                    std::cout << TAG<<"Applying scattering map "<<hname<<" from "<<fname<<", for timetof region = ["<<time1<<","<<time2<<","<<time3<<"]"<<std::endl;
                    TFile fs(fname.c_str());
                    TH1D* hs = (TH1D*)fs.Get(hname.c_str());
                    s->SetScatterMap(time1,time2,time3,*hs);
                }
                else if (optname=="time_offset")
                {
                    auto width = toml_h::find<double>(opt,1);
                    std::cout << TAG<<"Random timetof offset per PMT, sampled from a Gaussian of width = "<< width <<std::endl;
                    s->SetTimeOffset(true, width);
                }
                else if (optname=="time_smear")
                {
                    auto mean = toml_h::find<double>(opt,1);
                    auto width = toml_h::find<double>(opt,2);
                    std::cout << TAG<<"Random timetof smearing per PMT, sampled from a Gaussian of mean = "<< mean << ", width = "<< width <<std::endl;
                    s->SetTimeSmear(true, mean, width);
                }
                else if (optname=="z0")
                {
                    auto z0 = toml_h::find<double>(opt,1);
                    std::cout << TAG<<"Setting diffuser Z-position to "<< z0 <<std::endl;
                    s->SetZ0(z0);
                }
                else if (optname=="template")
                {
                    // auto fname = toml_h::find<std::string>(opt,1);
                    // auto hname = toml_h::find<std::string>(opt,2);
                    // auto offset = toml_h::find<double>(opt,3);
                    // auto lo = toml_h::find<double>(opt,4);
                    // auto hi = toml_h::find<double>(opt,5);
                    // std::cout << TAG<<"Use timetof template "<< hname <<" from "<< fname <<std::endl;
                    // TFile fs(fname.c_str());
                    // TH1D* hs = (TH1D*)fs.Get(hname.c_str());
                    // s->SetTemplate(*hs, offset, lo, hi);

                    auto fname = toml_h::find<std::string>(opt,1);
                    auto hname = toml_h::find<std::string>(opt,2);
                    auto offset = toml_h::find<double>(opt,3);
                    auto combine = toml_h::find<bool>(opt,4);
                    auto template_only = toml_h::find<bool>(opt,5);
                    std::cout << TAG<<"Use timetof template "<< hname << " from "<< fname <<std::endl;
                    TFile fs(fname.c_str());
                    TH2D* hist = (TH2D*)fs.Get(hname.c_str());
                    s->SetTemplate(*hist, offset, combine, template_only);
                }
                else if (optname=="pmt_eff")
                {
                    auto fname = toml_h::find<std::string>(opt,1);
                    auto hname = toml_h::find<std::string>(opt,2);
                    std::cout << TAG << "Use PMT efficiency histogram "<< hname << " from " << fname << std::endl;
                    TFile fs(fname.c_str());
                    TH1D* hist = (TH1D*)fs.Get(hname.c_str());
                    s->SetPMTEff(*hist);
                }
                else if (optname=="pmt_eff_var")
                {
                    auto sigma = toml_h::find<double>(opt,1);
                    std::cout << TAG << "Apply PMT efficiency variation with 1-sigma =  "<< sigma << std::endl;
                    s->SetPMTEffVar(sigma);
                }
            }
        }

        s->LoadEventsFromFile(filename,ev_tree_name,pmt_treename);
        s->InitEventMap();
        samples.push_back(s);
    }

    //Add fit parameters
    std::vector<AnaFitParameters*> fitparas;
    for (auto const &name : toml_h::find<std::vector<std::string>>(fitparameters_config, "names"))
    {
        std::cout << TAG << "Parameter name: " << name << std::endl;
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
        // Read parameter setup vector
        for (auto const &par : par_setup)
        {
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
        // If length of parameter setup vector < npar, copy setup from last entry
        if (parnames.size()<npar)
        {   
            int lastidx = parnames.size()-1;
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
            std::cout << TAG<<parnames[i]<<" "<<priors[i]<<" "<<steps[i]<<" "<<lows[i]<<" "<<highs[i]<<" "<<fixeds[i]<<" "<<std::endl;

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
        std::cout << TAG<< "binning_file = "<<binning_file<<", binning_var = [ ";
        for (auto t : binning_var) std::cout << t <<", ";
        std::cout <<"]"<<std::endl;

        auto fitpara = new AnaFitParameters(name,pmttype);
        // Optional config
        if (ele.size()>5)
        {
            for (int i=5; i<ele.size();i++ )
            {
                auto opt = toml_h::find<toml::array>(ele,i);
                auto optname = toml_h::find<std::string>(opt,0);
                if (optname=="covariance") // load prior covariance matrix
                {
                    auto fname = toml_h::find<std::string>(opt,1);
                    auto matname = toml_h::find<std::string>(opt,2);
                    TFile f(fname.c_str());
                    std::cout << TAG<<"Using covariance matrix "<<matname<<" from "<<fname<<std::endl;
                    TMatrixDSym* cov_mat = (TMatrixDSym*)f.Get(matname.c_str());
                    fitpara->SetCovarianceMatrix(*cov_mat);
                } 
                else if (optname=="prior") // set prior central values
                {
                    auto fname = toml_h::find<std::string>(opt,1);
                    auto hname = toml_h::find<std::string>(opt,2);
                    TFile f(fname.c_str());
                    std::cout << TAG<<"Set prior central values to  "<<hname<<" from "<<fname<<std::endl;
                    TH1D* hist = (TH1D*)f.Get(hname.c_str());
                    for (int j=0;j<npar;j++)
                        priors[j] = hist->GetBinContent(j+1);
                } 
                else if (optname=="spline") // set spline reweight
                {
                    auto fname = toml_h::find<std::vector<std::string>>(opt,1);
                    auto sname = toml_h::find<std::vector<std::string>>(opt,2);
                    //auto num = toml_h::find<int>(opt,3);
                    fitpara->SetSpline(fname,sname);
                } 
            }
        }
        fitpara->SetParameterFunction(functype);
        fitpara->InitParameters(parnames,priors,steps,lows,highs,fixeds);
        fitpara->SetBinVar(binning_var);
        fitpara->SetBinning(binning_file);
        fitpara->InitEventMap(samples);

        fitparas.push_back(fitpara);
    }

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

    fitter.SetMinSettings(min_settings);
    fitter.InitFitter(fitparas);

    bool stat_fluc = toml_h::find<bool>(minimizer_config, "stat_fluc");
    bool did_converge = false;

    did_converge = fitter.Fit(samples, stat_fluc);
    if(!did_converge)
        std::cout << TAG << "Fit did not coverge." << std::endl;

    // optonal MCMC for error estimation
    int MCMCSteps = toml_h::find<int>(minimizer_config, "MCMCSteps");
    if (MCMCSteps>0)
    {
        double MCMCStepSize = toml_h::find<double>(minimizer_config, "MCMCStepSize");
        std::cout << TAG << "Running MCMC for " << MCMCSteps << " steps, with step size = " << MCMCStepSize << std::endl;
        fitter.RunMCMCScan(MCMCSteps,MCMCStepSize);
    }

    // optional 1D parameter scan
    int ScanSteps = toml_h::find<int>(minimizer_config, "ScanSteps");
    if (ScanSteps>0)
    {
        std::vector<int> ParameterScans = toml_h::find<std::vector<int>>(minimizer_config, "ParameterScans");
        std::cout << TAG << "Running 1D scan for parameters ";
        for (int i=0;i<ParameterScans.size();i++)
            std::cout << ParameterScans[i]<<", ";
        std::cout << "for " << ScanSteps <<" steps" << std::endl;
        fitter.ParameterScans(ParameterScans,ScanSteps);
    }

    fout->Close();

    std::cout << TAG << "Output saved to " << fname_output << std::endl;

    if (toys>0)
    {
        std::cout << TAG << "Running " << toys << " toy fits..." << std::endl;

        for (int i=0;i<toys;i++)
        {
            std::cout << TAG << "Processing Toy Fit " << i << " :" << std::endl;

            std::string toy_output = fname_output;
            toy_output.insert(toy_output.size()-5,Form("_toy%i",i));
            TFile* fout_toy = TFile::Open(toy_output.c_str(), "RECREATE");
            fitter.SetDirectory(fout_toy);

            for (int j=0;j<samples.size();j++)
            {
                samples[j]->InitToy();
            }

            bool toy_converge = fitter.Fit(samples, stat_fluc);
            if(!toy_converge)
                std::cout << TAG << "Toy Fit " << i <<" did not coverge." << std::endl;

            fout_toy->Close();
        }
    }
    return 0;
}    