#ifndef PARAMETERFUNCTION_HH
#define PARAMETERFUNCTION_HH

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TMath.h>

#include "AnaEvent.hh"

class ParameterFunction
{
public:
    virtual ~ParameterFunction() {};
    virtual double operator()(double par, AnaEvent ev)
    {
        return 0.0;
    }
};

class Identity : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        return par;
    }
};

class Attenuation : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        double R = ev.GetR();
        double omega = ev.GetOmega();
        double val = TMath::Exp(-R/par)*omega;

        return val;
    }
};

#endif
