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
    BinManager* bm = new BinManager(10,0.5,1,100,1000,9000);
    auto s = new AnaSample(0, "Sample 0", bm);
    samples.push_back(s);

    //read events
    AnaTree selTree(fname_input.c_str(), "hitRate_pmtType1", "pmt_type1");
    selTree.SetTimetofCut(-950,-940);
    selTree.SetCosthsCut(0.766,1);
    selTree.GetEvents(samples.at(0));

    return 0;
}    