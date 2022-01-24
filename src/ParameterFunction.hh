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

enum FunctionType
{
    kIdentity    = 0,
    kAttenuation = 1,
    kScatter = 2,
    kPolynomialCosth = 3,
    kAttenuationZ = 4,
    kSourcePhiVar = 5,
    kSpline = 6, 
};


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

class AttenuationZ : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        double R = ev.GetR();
        double omega = ev.GetOmega();
        double dz = ev.GetDz();
        double z0 = ev.GetZ0();
        double alpha_z0 = alpha0 + slopeA*z0;
        double val;
        double da = slopeA*dz;
        if (fabs(da)>1.e-9)
        {
            val = TMath::Power(1+da/alpha_z0,-R/da)*omega;
        }
        else
        {
            val = TMath::Exp(-R/alpha_z0)*omega;
        }

        return val;
    }
    double alpha0;
    double slopeA;
};

class Scatter : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        return 1.;
    }
};

class SourcePhiVar : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        return 1+par*cos(ev.GetPhis());
    }
};

class PolynomialCosth : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        double costh = ev.GetCosth();
        double val = 0;
        for (int i=0;i<pol_orders.size();i++)
        {
            if (costh>=pol_range[i] && costh<pol_range[i+1])
            {
                //std::cout<<"Using pol"<<i<<std::endl;
                // for (int j=0;j<=pol_orders[i];j++)
                // {
                //     val += pol_coeff[i][j]*TMath::Power(costh-pol_range[i],j);
                // }
                double x = costh-pol_range[i];
                for (int j=pol_orders[i];j>=0;j--)
                {
                    val = val*x + pol_coeff[i][j];
                }
                break;
            }
        }
        //std::cout<<"CalcPol:: costh = "<<costh<<", val = "<<val<<std::endl;
        return val;
    }
    std::vector<int> pol_orders; // order of polynomial in each piece 
    std::vector<double> pol_range; // applicable range for each polynomial
    // Coefficients and boundary conditions for each polynomial
    std::vector<std::vector<double>> pol_coeff;
    std::vector<double> pol_p0, pol_p1;
    void SetPolOrders(std::vector<int>& vec) { pol_orders = vec; }
    void SetPolRange(std::vector<double>& vec) { pol_range = vec; }
    void Print()
    {
        std::cout<<"Polynomial orders: ";
        for (int i=0;i<pol_orders.size();i++)
            std::cout<<pol_orders[i]<<" ";
        std::cout<<"\nPolynomial ranges: ";
        for (int i=0;i<pol_range.size();i++)
            std::cout<<pol_range[i]<<" ";
        std::cout<<std::endl;
    }
    void SetPolynomial(std::vector<double>& params)
    {
        pol_coeff.clear(); pol_p0.clear(); pol_p1.clear(); 
        int par_index = 0;
        for (int i=0;i<pol_orders.size();i++)
        {
            std::vector<double> coeff;
            if (i==0) // for the first polynomial, the coefficients are unconstrained
            {
                for (int j=0;j<=pol_orders[i];j++)
                {
                    coeff.push_back(params[j]);
                    par_index++;
                } 
            }
            else // for others, we need to match the 0-th and 1-st order derivatives 
            {
                coeff.push_back(pol_p0[i-1]);
                coeff.push_back(pol_p1[i-1]);
                for (int j=2;j<=pol_orders[i];j++)
                {
                    coeff.push_back(params[par_index]);
                    par_index++;
                } 
            }
            pol_coeff.emplace_back(coeff);
            double p0 = 0;
            double p1 = 0;
            for (int j=0;j<=pol_orders[i];j++) // store the 0-th and 1-st order derivatives at end-point as boundary conditions
            {
                p0 += coeff[j]*TMath::Power(pol_range[i+1]-pol_range[i],j);
                p1 += coeff[j]*j*TMath::Power(pol_range[i+1]-pol_range[i],j-1);
            } 
            pol_p0.push_back(p0);
            pol_p1.push_back(p1);
        }
    }
};

class Spline : public ParameterFunction
{
public:
    double operator()(double par, AnaEvent ev)
    {
        return 1.;
    }
};

#endif
