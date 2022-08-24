#include<vector>
double truth_alpha(double wavelength=401.9, double ABWFF=1.3000000000, double RAYFF=1000000000) {
    const int NUMENTRIES_water=60;
    const double GeV=1.e9;
    const double cm=1;
    //const double MIEFF=0.0;

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

    double ABSORPTION_water[NUMENTRIES_water] =
    {
        16.1419*cm*ABWFF,  18.278*cm*ABWFF, 21.0657*cm*ABWFF, 24.8568*cm*ABWFF, 30.3117*cm*ABWFF, 
        38.8341*cm*ABWFF, 54.0231*cm*ABWFF, 81.2306*cm*ABWFF, 120.909*cm*ABWFF, 160.238*cm*ABWFF, 
        193.771*cm*ABWFF, 215.017*cm*ABWFF, 227.747*cm*ABWFF,  243.85*cm*ABWFF, 294.036*cm*ABWFF, 
        321.647*cm*ABWFF,  342.81*cm*ABWFF, 362.827*cm*ABWFF, 378.041*cm*ABWFF, 449.378*cm*ABWFF,
        739.434*cm*ABWFF, 1114.23*cm*ABWFF, 1435.56*cm*ABWFF, 1611.06*cm*ABWFF, 1764.18*cm*ABWFF, 
        2100.95*cm*ABWFF,  2292.9*cm*ABWFF, 2431.33*cm*ABWFF,  3053.6*cm*ABWFF, 4838.23*cm*ABWFF, 
        6539.65*cm*ABWFF, 7682.63*cm*ABWFF, 9137.28*cm*ABWFF, 12220.9*cm*ABWFF, 15270.7*cm*ABWFF, 
        19051.5*cm*ABWFF, 23671.3*cm*ABWFF, 29191.1*cm*ABWFF, 35567.9*cm*ABWFF,   42583*cm*ABWFF,
        49779.6*cm*ABWFF, 56465.3*cm*ABWFF,   61830*cm*ABWFF, 65174.6*cm*ABWFF, 66143.7*cm*ABWFF,   
        64820*cm*ABWFF,   61635*cm*ABWFF, 57176.2*cm*ABWFF, 52012.1*cm*ABWFF, 46595.7*cm*ABWFF, 
        41242.1*cm*ABWFF, 36146.3*cm*ABWFF, 31415.4*cm*ABWFF, 27097.8*cm*ABWFF, 23205.7*cm*ABWFF, 
        19730.3*cm*ABWFF, 16651.6*cm*ABWFF, 13943.6*cm*ABWFF, 11578.1*cm*ABWFF, 9526.13*cm*ABWFF
    };

    double RAYLEIGH_water[NUMENTRIES_water] = {
            386929*cm*RAYFF,  366249*cm*RAYFF,  346398*cm*RAYFF,  327355*cm*RAYFF,  309097*cm*RAYFF,  
            291603*cm*RAYFF,  274853*cm*RAYFF,  258825*cm*RAYFF,  243500*cm*RAYFF,  228856*cm*RAYFF,  
            214873*cm*RAYFF,  201533*cm*RAYFF,  188816*cm*RAYFF,  176702*cm*RAYFF,  165173*cm*RAYFF,
            154210*cm*RAYFF,  143795*cm*RAYFF,  133910*cm*RAYFF,  124537*cm*RAYFF,  115659*cm*RAYFF,  
            107258*cm*RAYFF, 99318.2*cm*RAYFF, 91822.2*cm*RAYFF,   84754*cm*RAYFF, 78097.3*cm*RAYFF, 
            71836.5*cm*RAYFF,   65956*cm*RAYFF, 60440.6*cm*RAYFF, 55275.4*cm*RAYFF, 50445.6*cm*RAYFF,
            45937*cm*RAYFF, 41735.2*cm*RAYFF, 37826.6*cm*RAYFF, 34197.6*cm*RAYFF, 30834.9*cm*RAYFF, 
            27725.4*cm*RAYFF, 24856.6*cm*RAYFF, 22215.9*cm*RAYFF, 19791.3*cm*RAYFF, 17570.9*cm*RAYFF,   
            15543*cm*RAYFF, 13696.6*cm*RAYFF, 12020.5*cm*RAYFF, 10504.1*cm*RAYFF, 9137.15*cm*RAYFF,
            7909.45*cm*RAYFF,  6811.3*cm*RAYFF, 5833.25*cm*RAYFF,  4966.2*cm*RAYFF, 4201.36*cm*RAYFF, 
            3530.28*cm*RAYFF, 2944.84*cm*RAYFF, 2437.28*cm*RAYFF, 2000.18*cm*RAYFF,  1626.5*cm*RAYFF, 
            1309.55*cm*RAYFF, 1043.03*cm*RAYFF, 821.016*cm*RAYFF,  637.97*cm*RAYFF, 488.754*cm*RAYFF
    };

//    double MIE_water[NUMENTRIES_water] = {
//     7790020*cm*MIEFF, 7403010*cm*MIEFF, 7030610*cm*MIEFF, 6672440*cm*MIEFF, 6328120*cm*MIEFF,
//     5997320*cm*MIEFF, 5679650*cm*MIEFF, 5374770*cm*MIEFF, 5082340*cm*MIEFF, 4802000*cm*MIEFF,
//     4533420*cm*MIEFF, 4276280*cm*MIEFF, 4030220*cm*MIEFF, 3794950*cm*MIEFF, 3570120*cm*MIEFF,
//     3355440*cm*MIEFF, 3150590*cm*MIEFF, 2955270*cm*MIEFF, 2769170*cm*MIEFF, 2592000*cm*MIEFF,
//     2423470*cm*MIEFF, 2263300*cm*MIEFF, 2111200*cm*MIEFF, 1966900*cm*MIEFF, 1830120*cm*MIEFF,
//     1700610*cm*MIEFF, 1578100*cm*MIEFF, 1462320*cm*MIEFF, 1353040*cm*MIEFF, 1250000*cm*MIEFF,
//     1152960*cm*MIEFF, 1061680*cm*MIEFF,  975936*cm*MIEFF,  895491*cm*MIEFF,  820125*cm*MIEFF,
//      749619*cm*MIEFF,  683760*cm*MIEFF,  622339*cm*MIEFF,  565152*cm*MIEFF,  512000*cm*MIEFF,
//      462688*cm*MIEFF,  417027*cm*MIEFF,  374832*cm*MIEFF,  335923*cm*MIEFF,  300125*cm*MIEFF,
//      267267*cm*MIEFF,  237184*cm*MIEFF,  209715*cm*MIEFF,  184704*cm*MIEFF,  162000*cm*MIEFF,
//      141456*cm*MIEFF,  122931*cm*MIEFF,  106288*cm*MIEFF, 91395.2*cm*MIEFF,   78125*cm*MIEFF,
//     66355.2*cm*MIEFF, 55968.2*cm*MIEFF, 46851.2*cm*MIEFF, 38896.2*cm*MIEFF,   32000*cm*MIEFF
//   };


    TGraph* gr_abs = new TGraph(NUMENTRIES_water,ENERGY_water,ABSORPTION_water);
    TGraph* gr_ray = new TGraph(NUMENTRIES_water,ENERGY_water,RAYLEIGH_water);
//    TGraph* gr_mie = new TGraph(NUMENTRIES_water,ENERGY_water,MIE_water);
    double photoEnergy = 1239.84193/wavelength;
    double alpha = 1./(1./gr_abs->Eval(photoEnergy)+1./gr_ray->Eval(photoEnergy));// + 1./gr_mie->Eval(photoEnergy));
    std::cout<<"Absorption length = "<<gr_abs->Eval(photoEnergy)<<" cm"<<std::endl;
    std::cout<<"Scattering length = "<<gr_ray->Eval(photoEnergy)<<" cm"<<std::endl;
    std::cout<<"*Attenutation length* = "<<alpha<<" cm"<<std::endl;
    std::cout<< "photon energy = " << photoEnergy << "eV"<<std::endl;
    std::cout<<"wavelength = " << wavelength << "nm" << std::endl;
    return alpha;//,gr_ray->Eval(photoEnergy)};//, gr_mie->Eval(photoEnergy)};
}
