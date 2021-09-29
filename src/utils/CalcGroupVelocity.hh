#ifndef CALCGROUPVELOCITY_HH
#define CALCGROUPVELOCITY_HH

#include <TGraph.h>

double CalcGroupVelocity(double wavelength) {
    const int NUMENTRIES_water=60;
    const double GeV=1.e9;

    const double c_light   = 2.99792458e+8;

    // arrray of refractive index copied from WCSIM
    double ENERGY_water[NUMENTRIES_water] =
     { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
       1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
       1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV, 1.82352e-09*GeV, 
       1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
       1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
       2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
       2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
       2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
       2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
       2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
       3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
       3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
       3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
       4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
       5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };

    double RINDEX1[NUMENTRIES_water] = 
     {1.32885, 1.32906, 1.32927, 1.32948, 1.3297, 1.32992, 1.33014, 
      1.33037, 1.3306, 1.33084, 1.33109, 1.33134, 1.3316, 1.33186, 1.33213,
      1.33241, 1.3327, 1.33299, 1.33329, 1.33361, 1.33393, 1.33427, 1.33462,
      1.33498, 1.33536, 1.33576, 1.33617, 1.3366, 1.33705, 1.33753, 1.33803,
      1.33855, 1.33911, 1.3397, 1.34033, 1.341, 1.34172, 1.34248, 1.34331,
      1.34419, 1.34515, 1.3462, 1.34733, 1.34858, 1.34994, 1.35145, 1.35312,
      1.35498, 1.35707, 1.35943, 1.36211, 1.36518, 1.36872, 1.37287, 1.37776,
      1.38362, 1.39074, 1.39956, 1.41075, 1.42535};

    // Calculate group velocity
    // Copy from G4MaterialPropertiesTable.cc
    double energy_vg[NUMENTRIES_water];
    double groupvel[NUMENTRIES_water];

    double E0 = ENERGY_water[0];
    double n0 = RINDEX1[0];

    double E1 = ENERGY_water[1];
    double n1 = RINDEX1[1];

    // add entry at first photon energy
    double vg = c_light/(n0+(n1-n0)/std::log(E1/E0));
    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    if((vg<0) || (vg>c_light/n0))  { vg = c_light/n0; }
    energy_vg[0] = E0; 
    groupvel[0] = vg;

    // add entries at midpoints between remaining photon energies
    for (int i=2;i<NUMENTRIES_water;i++)
    {
      vg = c_light/( 0.5*(n0+n1)+(n1-n0)/std::log(E1/E0));
      if((vg<0) || (vg>c_light/(0.5*(n0+n1))))  { vg = c_light/(0.5*(n0+n1)); }
      energy_vg[i-1] =  0.5*(E0+E1);
      groupvel[i-1] = vg;

      E0 = E1;
      n0 = n1;
      E1 = ENERGY_water[i];
      n1 = RINDEX1[i];
    }

    vg = c_light/(n1+(n1-n0)/std::log(E1/E0));
    if((vg<0) || (vg>c_light/n1))  { vg = c_light/n1; }
    energy_vg[NUMENTRIES_water-1] = E1;
    groupvel[NUMENTRIES_water-1] = vg;

    TGraph* gr_groupvel = new TGraph(NUMENTRIES_water,energy_vg,groupvel);
    double photoEnergy = 1239.84193/wavelength;

    return gr_groupvel->Eval(photoEnergy,0,"S");
}

#endif // CALCGROUPVELOCITY_HH