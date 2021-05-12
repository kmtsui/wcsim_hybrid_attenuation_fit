#ifndef BINMANAGER_H
#define BINMANAGER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class BinManager
{
    public:
        BinManager();
        BinManager(const std::string& filename);
        BinManager(const int nx, const double xmin, const double xmax, const int ny, const double ymin, const double ymax);

        int GetNbins() const;
        int SetBinning(const std::string& filename);
        int GetBinIndex(const std::vector<double>& val) const;
        double GetBinWidth(const unsigned int i) const;
        double GetBinWidth(const int i, const int d) const;
        std::vector<double> GetBinVector(const double d) const;
        std::vector<std::vector<std::pair<double, double>>> GetEdgeVector() const { return bin_edges; };
        std::vector<std::pair<double, double>> GetEdgeVector(const int d) const { return bin_edges.at(d); };
        void Print() const;

    private:
        bool CheckBinIndex(const int i, const int d, const double val) const;

        unsigned int dimension;
        unsigned int nbins;
        std::string fname_binning;
        std::vector<std::vector<std::pair<double, double>>> bin_edges;
};

#endif
