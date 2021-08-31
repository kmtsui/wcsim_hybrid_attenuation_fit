#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TRandom3.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

// Mock up digitization code copied from WCSIM
// Currently contains BoxandLine20inchHQE and PMT3inchR14374 options
#include "OPTICALFIT/utils/WCSIMDigitization.hh"

using namespace std;

TRandom3* rng;

WCSimRootGeom *geo = 0; 

// Simple example of reading a generated Root file
const int nPMTtypes = 2;
double PMTradius[nPMTtypes];

double CalcSolidAngle(double r, double R, double costh)
{
  double weight = 2*TMath::Pi()*(1-R/sqrt(R*R+r*r));
  //weight *= 1.-0.5*sqrt(1-costh*costh);
  //weight *= 0.5+0.5*costh;

  return weight;
}

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

TGraph *gr_power_cosths, *gr_variation_cosths;
void SetLEDProfile()
{
  std::cout<<"Set up LED profile for event reweight"<<std::endl;
  const int nCosths = 41;
  double Cosths[nCosths] = 
    {
      0.76655, 0.77700, 0.78761, 0.79875, 0.80872, 0.81934, 0.82883, 0.83808,
      0.84711, 0.85590, 0.86523, 0.87426, 0.88229, 0.89008, 0.89831, 0.90623,
      0.91259, 0.92055, 0.92700, 0.93319, 0.93911, 0.94477, 0.95063, 0.95618,
      0.96098, 0.96551, 0.96976, 0.97408, 0.97808, 0.98116, 0.98454, 0.98736,
      0.98989, 0.99233, 0.99444, 0.99607, 0.99742, 0.99857, 0.99933, 0.99981,
      1.00000
    };
  double Power_cosths[nCosths] =
    {
      0.50152, 0.54801, 0.59419, 0.63473, 0.67604, 0.71697, 0.75712, 0.80037,
      0.83589, 0.86910, 0.90307, 0.93628, 0.95584, 0.97296, 0.98300, 0.98918,
      0.99844, 1.00063, 1.00848, 1.01041, 1.01466, 1.01608, 1.01813, 1.02019,
      1.02200, 1.02122, 1.02277, 1.02084, 1.02007, 1.01813, 1.01698, 1.01402,
      1.01427, 1.01273, 1.01080, 1.00938, 1.00810, 1.00694, 1.00578, 1.00501,
      1.00308
    };
  gr_power_cosths = new TGraph(nCosths,Cosths,Power_cosths);

  double Variation_cosths[nCosths] =
    {
      1.52148, 1.32716, 1.60492, 1.51191, 1.32437, 1.07736, 0.80404, 0.84281,
      0.89956, 0.97436, 0.87347, 0.95811, 1.03983, 1.16010, 1.18226, 1.17236,
      1.14083, 1.09142, 1.10767, 1.02450, 0.96195, 0.91277, 0.88456, 0.85350,
      0.82181, 0.78467, 0.74606, 0.67583, 0.60713, 0.54087, 0.50040, 0.42365,
      0.37995, 0.32521, 0.27719, 0.23896, 0.19595, 0.17128, 0.13077, 0.08523,
      0.06337
    };
  gr_variation_cosths = new TGraph(nCosths,Cosths,Variation_cosths);
}

double GetLEDWeight(double cosths, double phis)
{
  if (cosths<0.76655) return 0.;
  return gr_power_cosths->Eval(cosths)*(1.+gr_variation_cosths->Eval(cosths)/100.*cos(phis));
}

void HelpMessage()
{
  std::cout << "USAGE: "
            << "WCSIM_TreeConvert" << "\nOPTIONS:\n"
            << "-f : Intput file\n"
            << "-o : Output file\n"
            << "-l : Laser wavelength\n"
            << "-w : Apply diffuser profile reweight\n"
            << "-b : Use only B&L PMTs\n"
            << "-d : Run with raw Cherenkov hits and perform ad-hoc digitization\n"
            << "-t : Use separated triggers\n"
            << "-v : Verbose\n"
            << "-s : Start event\n"
            << "-e : End event\n"
            << "-r : RNG seed\n";
}

int main(int argc, char **argv){
  
  char * filename=NULL;
  char * outfilename=NULL;
  bool verbose=false;
  bool hybrid = true;
  double cvacuum = 3e8 / 1e9;//speed of light, in meter per ns.
  double nindex = 1.373;//refraction index of water
  bool plotDigitized = true; //using digitized hits
  bool separatedTriggers=false; //Assume two independent triggers, one for mPMT, one for B&L
  bool diffuserProfile = false; //Reweigh PMT hits by source angle

  double wavelength = 400; //wavelength in nm

  int nPMTpermPMT=19;

  int startEvent=0;
  int endEvent=0;
  int seed = 0;
  char c;
  while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:hdtvw")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
      case 'f':
        filename = optarg;
        break;
      case 'd':
        plotDigitized = false; //using raw hits
        break;
      case 'b':
        hybrid = false; // no mPMT
        break;
      case 't':
        separatedTriggers = true;
        break;
      case 'v':
        verbose = true;
        break;
      case 'o':
	      outfilename = optarg;
	      break;
      case 's':
	      startEvent = std::stoi(optarg);
        if (startEvent<0) startEvent = 0; 
	      break;
      case 'e':
	      endEvent = std::stoi(optarg);
	      break;
      case 'r':
	      seed = std::stoi(optarg);
        if (seed<0) seed=0;
        std::cout<<"Set RNG seed = "<<seed<<std::endl;
	      break;
      case 'l':
        wavelength = std::stod(optarg);
        if (wavelength<0) {
          std::cout<<"Wavelength < 0, using default = 400 nm"<<std::endl;
          wavelength = 400;
        }
	      break;
      case 'w':
        diffuserProfile = true;
        break;
      case 'h':
        HelpMessage();
      default:
        return 0;
    }
  }
  
  double vg = CalcGroupVelocity(wavelength);
  std::cout<<"Using wavelength = "<<wavelength<<" nm, group velocity = "<<vg<<" m/s, n = "<<cvacuum/vg<<std::endl;
  vg /= 1.e9; // convert to m/ns

  rng = new TRandom3(seed);
  gRandom = rng;
  
  TFile *file = TFile::Open(filename);
  // Open the file
  if (filename==NULL){
    std::cout << "Error, no input file: " << std::endl;
    HelpMessage();
    return -1;
  }
  if (!file->IsOpen()){
    std::cout << "Error, could not open input file: " << filename << std::endl;
    return -1;
  }
  
  SetLEDProfile();

  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");

  // Get the number of events
  int nevent = ((int)tree->GetEntries());//std::min(((int)tree->GetEntries()),100000);
  if(endEvent!=0 && endEvent<=nevent) nevent = endEvent;
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  TBranch *branch2;
  if(hybrid){
    branch2 = tree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
  // Force deletion to prevent memory leak 
    tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
      exit(9);
  }
  geotree->GetEntry(0);
  PMTradius[0]=geo->GetWCPMTRadius();
  PMTradius[1]=geo->GetWCPMTRadius(true);

  // Options tree - only need 1 "event"
  TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
  WCSimRootOptions *opt = 0; 
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  if(outfilename==NULL) outfilename = (char*)"out.root";
  
  TFile * outfile = new TFile(outfilename,"RECREATE");
  std::cout<<"File "<<outfilename<<" is open for writing"<<std::endl;

  double nHits, nPE, timetof;
  double weight;
  int nReflec, nRaySct, nMieSct;
  double photonStartTime;
  double nPE_digi, timetof_digi;
  int PMT_id;
  // TTree for storing the hit information. One for B&L PMT, one for mPMT
  TTree* hitRate_pmtType0 = new TTree("hitRate_pmtType0","hitRate_pmtType0");
  hitRate_pmtType0->Branch("nHits",&nHits); // dummy variable, always equal to 1
  hitRate_pmtType0->Branch("nPE",&nPE); // number of PE
  hitRate_pmtType0->Branch("timetof",&timetof); // hittime-tof
  hitRate_pmtType0->Branch("PMT_id",&PMT_id);
  hitRate_pmtType0->Branch("weight",&weight);
  // Branches below only filled for raw hits
  if (!plotDigitized)
  {
    hitRate_pmtType0->Branch("nReflec",&nReflec); // Number of reflection experienced by a photon before reaching the sensitive detector
    hitRate_pmtType0->Branch("nRaySct",&nRaySct); // Number of Rayleigh scattering
    hitRate_pmtType0->Branch("nMieSct",&nMieSct); // Number of Mie scattering
    hitRate_pmtType0->Branch("photonStartTime",&photonStartTime); // True photon start time
    hitRate_pmtType0->Branch("nPE_digi",&nPE_digi); // nPE after ad-hoc digitization
    hitRate_pmtType0->Branch("timetof_digi",&timetof_digi); // hittime-tof after ad-hoc digitization
  }
  TTree* hitRate_pmtType1 = new TTree("hitRate_pmtType1","hitRate_pmtType1");
  hitRate_pmtType1->Branch("nHits",&nHits);
  hitRate_pmtType1->Branch("nPE",&nPE);
  hitRate_pmtType1->Branch("timetof",&timetof);
  hitRate_pmtType1->Branch("PMT_id",&PMT_id);
  hitRate_pmtType1->Branch("weight",&weight);
  if (!plotDigitized)
  {
    hitRate_pmtType1->Branch("nReflec",&nReflec);
    hitRate_pmtType1->Branch("nRaySct",&nRaySct);
    hitRate_pmtType1->Branch("nMieSct",&nMieSct);
    hitRate_pmtType1->Branch("photonStartTime",&photonStartTime);
    hitRate_pmtType1->Branch("nPE_digi",&nPE_digi);
    hitRate_pmtType1->Branch("timetof_digi",&timetof_digi);
  }

  double vtxpos[3];
  tree->GetEntry(0);
  wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
  for (int i=0;i<3;i++) vtxpos[i]=wcsimrootevent->GetVtx(i);

  // Save the PMT geometry information relative to source
  double dist, costh, cosths, phis, omega, phim, costhm;
  int mPMT_id;
  TTree* pmt_type0 = new TTree("pmt_type0","pmt_type0");
  pmt_type0->Branch("R",&dist);          // distance to source
  pmt_type0->Branch("costh",&costh);     // photon incident angle relative to PMT
  pmt_type0->Branch("cosths",&cosths);   // PMT costheta angle relative to source
  pmt_type0->Branch("phis",&phis);       // PMT phi angle relative to source
  pmt_type0->Branch("costhm",&costhm);   // costhm = costh
  pmt_type0->Branch("phim",&phim);       // dummy
  pmt_type0->Branch("omega",&omega);     // solid angle subtended by PMT
  pmt_type0->Branch("PMT_id",&PMT_id);   // unique PMT id
  pmt_type0->Branch("mPMT_id",&mPMT_id); // dummy 
  pmt_type0->Branch("weight",&weight); 
  TTree* pmt_type1 = new TTree("pmt_type1","pmt_type1");
  pmt_type1->Branch("R",&dist);
  pmt_type1->Branch("costh",&costh);
  pmt_type1->Branch("cosths",&cosths);
  pmt_type1->Branch("phis",&phis);
  pmt_type1->Branch("costhm",&costhm);   // photon incident theta angle relative to central PMT 
  pmt_type1->Branch("phim",&phim);       // photon incident phi angle relative to central PMT 
  pmt_type1->Branch("omega",&omega);
  pmt_type1->Branch("PMT_id",&PMT_id);
  pmt_type1->Branch("mPMT_id",&mPMT_id); //sub-ID of PMT inside a mPMT module
                                         // 0 -11 : outermost ring
                                         // 12 - 17: middle ring
                                         // 18: central PMT
  pmt_type1->Branch("weight",&weight); 

  // LI source direction is always perpendicular to the wall
  // Separate treatment for barrel and endcap 
  double vDirSource[3];
  double vSource_localXaxis[3];
  double vSource_localYaxis[3];
  double endcapZ=3000;
  // Barrel injector
  if (abs(vtxpos[2])<endcapZ) {
    vDirSource[0]=vtxpos[0];
    vDirSource[1]=vtxpos[1];
    double norm = sqrt(vtxpos[0]*vtxpos[0]+vtxpos[1]*vtxpos[1]);
    vDirSource[0]/=-norm;
    vDirSource[1]/=-norm;
    vDirSource[2]=0;
    vSource_localXaxis[0]=0;vSource_localXaxis[1]=0;vSource_localXaxis[2]=1;
    vSource_localYaxis[0]=-vtxpos[1]/norm;vSource_localYaxis[1]=vtxpos[0]/norm;vSource_localYaxis[2]=0;
  } else {
    vDirSource[0]=0;
    vDirSource[1]=0;
    if (vtxpos[2]>endcapZ) 
    {
      vDirSource[2]=-1;
      vSource_localXaxis[0]=1;vSource_localXaxis[1]=0;vSource_localXaxis[2]=0;
      vSource_localYaxis[0]=0;vSource_localYaxis[1]=-1;vSource_localYaxis[2]=0;
    }
    else 
    {
      vDirSource[2]=1;
      vSource_localXaxis[0]=1;vSource_localXaxis[1]=0;vSource_localXaxis[2]=0;
      vSource_localYaxis[0]=0;vSource_localYaxis[1]=1;vSource_localYaxis[2]=0;
    }
  }

  int nPMTs_type0=geo->GetWCNumPMT();
  int nPMTs_type1=0; if (hybrid) nPMTs_type1=geo->GetWCNumPMT(true);
  std::vector<double> ledweight_type0(nPMTs_type0,-1);
  std::vector<double> ledweight_type1(nPMTs_type1,-1);
  for (int pmtType=0;pmtType<nPMTtypes;pmtType++) 
  {
    int nPMTs_type = pmtType==0 ? nPMTs_type0 : nPMTs_type1;
    for (int i=0;i<nPMTs_type;i++) 
    {
      WCSimRootPMT pmt;
      if (pmtType==0) pmt = geo->GetPMT(i,false);
      else pmt = geo->GetPMT(i,true);
      PMT_id = i;
      if (pmtType == 0) mPMT_id = 0;
      else mPMT_id = PMT_id%nPMTpermPMT;
      double PMTpos[3];
      double PMTdir[3];                   
      for(int j=0;j<3;j++){
        PMTpos[j] = pmt.GetPosition(j);
        PMTdir[j] = pmt.GetOrientation(j);
      }
      double particleRelativePMTpos[3];
      for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vtxpos[j];
      double vDir[3];double vOrientation[3];
      for(int j=0;j<3;j++){
        vDir[j] = particleRelativePMTpos[j];
        vOrientation[j] = PMTdir[j];
      }
      double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
      double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
      for(int j=0;j<3;j++){
        vDir[j] /= Norm;
        vOrientation[j] /= NormOrientation;
      }
      dist = Norm;
      costh = -(vDir[0]*vOrientation[0]+vDir[1]*vOrientation[1]+vDir[2]*vOrientation[2]);
      cosths = vDir[0]*vDirSource[0]+vDir[1]*vDirSource[1]+vDir[2]*vDirSource[2];
      double localx = vDir[0]*vSource_localXaxis[0]+vDir[1]*vSource_localXaxis[1]+vDir[2]*vSource_localXaxis[2];
      double localy = vDir[0]*vSource_localYaxis[0]+vDir[1]*vSource_localYaxis[1]+vDir[2]*vSource_localYaxis[2];
      phis = atan2(localy,localx);
      weight = GetLEDWeight(cosths,phis);
      if (pmtType==0) ledweight_type0[i]=weight;
      if (pmtType==1) ledweight_type1[i]=weight;
      double pmtradius = pmtType==0 ? PMTradius[0] : PMTradius[1]; 
      omega = CalcSolidAngle(pmtradius,dist,costh);
      costhm = costh;
      phim = 0;
      if (pmtType==1 && mPMT_id!=nPMTpermPMT-1)
      { 
        // Calculate photon incident cos(phi) angle relative to central PMT 
        int idx_centralpmt = int(PMT_id/nPMTpermPMT)*nPMTpermPMT + nPMTpermPMT-1; // locate the central PMT
        double PMTpos_central[3], PMTdir_central[3];
        WCSimRootPMT pmt_central = geo->GetPMT(idx_centralpmt,true);
        // central PMT position relative to current PMT
        for(int j=0;j<3;j++){
          PMTpos_central[j] = pmt_central.GetPosition(j)-PMTpos[j];
          PMTdir_central[j] = pmt_central.GetOrientation(j);
        }
        costhm = -(vDir[0]*PMTdir_central[0]+vDir[1]*PMTdir_central[1]+vDir[2]*PMTdir_central[2]);
        TVector3 v_central(PMTpos_central);
        TVector3 v_dir(vDir);
        TVector3 v_orientation(vOrientation);
        // Use cross product to extract the perpendicular component
        TVector3 v_dir1 = v_orientation.Cross(v_central);
        TVector3 v_dir2 = v_orientation.Cross(-v_dir);
        // Use dot cross to calculate the phi angle
        phim = v_dir1.Angle(v_dir2);
      }
      if (pmtType==0) pmt_type0->Fill();
      if (pmtType==1) pmt_type1->Fill();
    }
  }
  outfile->cd();
  pmt_type0->Write();
  pmt_type1->Write();

  // Ad-hoc digitizer, PMT specific 
  BoxandLine20inchHQE_Digitizer* BnLDigitizer = new BoxandLine20inchHQE_Digitizer();
  PMT3inchR14374_Digitizer* mPMTDigitizer = new PMT3inchR14374_Digitizer();

  // Now loop over events
  for (int ev=startEvent; ev<nevent; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);

    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    //wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    }

    if(verbose){
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
    }

    std::vector<double> triggerInfo;
    triggerInfo.clear();
    triggerInfo = wcsimrootevent->GetTriggerInfo();

    std::vector<double> triggerInfo2;
    triggerInfo2.clear();
    if(hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();


    if(verbose){
      for(unsigned int v=0;v<triggerInfo.size();v++){
	      cout << "Trigger entry #" << v << ", info = " << triggerInfo[v] << endl;
      }
      if(hybrid){
        for(unsigned int v=0;v<triggerInfo2.size();v++){
          cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << endl;
        }
      }
    }

    double triggerShift[nPMTtypes];
    double triggerTime[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      triggerShift[pmtType]=0;
      triggerTime[pmtType]=0;
      if(triggerInfo.size()>=3){
        if(pmtType==0){
          triggerShift[pmtType] = triggerInfo[1];
          triggerTime[pmtType] = triggerInfo[2];
        }
      }
      if(triggerInfo2.size()>=3){
        if(pmtType==1 && hybrid){
          triggerShift[pmtType] = triggerInfo2[1];
          triggerTime[pmtType] = triggerInfo2[2];
        }
      }
    }    
  

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
    int ncherenkovhits2 = 0; if(hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
    int ncherenkovdigihits2 = 0;if(hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits(); 
    
    if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      printf("Ncherenkovhits2 %d\n",     ncherenkovhits2);
      printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
      cout << "RAW HITS:" << endl;
    }

    if(!plotDigitized) for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(separatedTriggers){
        if(triggerInfo2.size()!=0 && pmtType==0) continue;
        if(triggerInfo.size()!=0 && pmtType==1) continue;
      }
      if(verbose) cout << "PMT Type = " << pmtType << endl;
 
      // Grab the big arrays of times and parent IDs
      
      TClonesArray *timeArray;//An array of pointers on CherenkovHitsTimes.
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      
      double particleRelativePMTpos[3];
      double totalPe = 0;
      //int totalHit = 0;

      int nhits;
      if(pmtType == 0) nhits = ncherenkovhits;
      else nhits = ncherenkovhits2;

      // Get the number of Cherenkov hits
      // Loop over sub events
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (int i=0; i< nhits ; i++)
      {
        //if(verbose) cout << "Hit #" << i << endl;

        WCSimRootCherenkovHit *wcsimrootcherenkovhit;
        if(pmtType==0) wcsimrootcherenkovhit = (WCSimRootCherenkovHit*) (wcsimrootevent->GetCherenkovHits())->At(i);
        else wcsimrootcherenkovhit = (WCSimRootCherenkovHit*) (wcsimrootevent2->GetCherenkovHits())->At(i);
        
        int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
        int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
        int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);

        WCSimRootPMT pmt;
        if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
        else pmt  = geo->GetPMT(tubeNumber-1,true);

        PMT_id = tubeNumber-1;

        double PMTpos[3];
        double PMTdir[3];
        for(int j=0;j<3;j++){
          PMTpos[j] = pmt.GetPosition(j);
          PMTdir[j] = pmt.GetOrientation(j);
        }
        
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        //double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
        //double phiAngle;
        //if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
        //else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
        
        // Use vertex positions for LI events
        for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vtxpos[j];
          
        double vDir[3];double vOrientation[3];
        for(int j=0;j<3;j++){
          vDir[j] = particleRelativePMTpos[j];
          vOrientation[j] = PMTdir[j];
        }
        double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
        double tof = Norm*1e-2/(vg);
        //if(verbose) cout << "Time of flight = " << tof << endl;
        double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
        for(int j=0;j<3;j++){
          vDir[j] /= Norm;
          vOrientation[j] /= NormOrientation;
        }
        
        WCSimRootCherenkovHitTime * HitTime = (WCSimRootCherenkovHitTime*) timeArray->At(timeArrayIndex);//Takes the first hit of the array as the timing, It should be the earliest hit
        //WCSimRootCherenkovHitTime HitTime = (WCSimRootCherenkovHitTime) timeArray->At(j);		  
        double time = HitTime->GetTruetime();

        photonStartTime = HitTime->GetPhotonStartTime();
        nReflec = 0;
        nRaySct = 0;
        nMieSct = 0;
        for (int idx = timeArrayIndex; idx<timeArrayIndex+peForTube; idx++) {
          WCSimRootCherenkovHitTime * cht = (WCSimRootCherenkovHitTime*) timeArray->At(idx);

          // only works for peForTube = 1
          // if peForTube > 1, you don't know whether reflection and scattering happens at the same time for a single photon
          if (cht->GetReflection()>0) nReflec++;
          if (cht->GetRayScattering()>0) nRaySct++;
          if (cht->GetMieScattering()>0) nMieSct++;
        }
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        timetof = time-tof;
        nHits = 1; nPE = peForTube; 

        // A simple way to mimic the digitization procedure
        nPE_digi = nPE;
        timetof_digi = timetof;
        if (pmtType==0) BnLDigitizer->Digitize(nPE_digi,timetof_digi);
        if (pmtType==1) mPMTDigitizer->Digitize(nPE_digi,timetof_digi);

        if (pmtType==0) hitRate_pmtType0->Fill();
        if (pmtType==1) hitRate_pmtType1->Fill();

      } // End of loop over Cherenkov hits
      if(verbose) cout << "Total Pe : " << totalPe << endl;
    }


    // Get the number of digitized hits
    // Loop over sub events
    if(verbose) cout << "DIGITIZED HITS:" << endl;

    if(plotDigitized) for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(separatedTriggers){
        if(triggerInfo2.size()!=0 && pmtType==0) continue;
        if(triggerInfo.size()!=0 && pmtType==1) continue;
      }
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      // Grab the big arrays of times and parent IDs
      /*
      // Not used
      TClonesArray *timeArray;
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      */
      double particleRelativePMTpos[3];
      double totalPe = 0;
      //int totalHit = 0;

      int nhits;
      if(pmtType == 0) nhits = ncherenkovdigihits;
      else nhits = ncherenkovdigihits2;
      
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (int i=0; i< nhits ; i++)
      { 
        TObject *Hit;
        if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
        else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

        WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
          dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
        
        int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId();
        double peForTube      = wcsimrootcherenkovdigihit->GetQ();

        //int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
        //int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
        //int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
        WCSimRootPMT pmt;
        if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
        else pmt  = geo->GetPMT(tubeNumber-1,true); 

        PMT_id = (tubeNumber-1.);
        
        double PMTpos[3];
        double PMTdir[3];                   
        for(int j=0;j<3;j++){
          PMTpos[j] = pmt.GetPosition(j);
          PMTdir[j] = pmt.GetOrientation(j);
        }
        
        
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        //double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
        //double phiAngle;
        //if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
        //else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));

        // Use vertex positions for LI events
        for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vtxpos[j];
        
        double vDir[3];double vOrientation[3];
        for(int j=0;j<3;j++){
          vDir[j] = particleRelativePMTpos[j];
          vOrientation[j] = PMTdir[j];	  
        }
        double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
        double tof = Norm*1e-2/(vg);
        double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
        for(int j=0;j<3;j++){
          vDir[j] /= Norm;
          vOrientation[j] /= NormOrientation;
        }

        double time = wcsimrootcherenkovdigihit->GetT();
      
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

        timetof = time-tof+triggerTime[pmtType]-triggerShift[pmtType];

        nHits = 1; nPE = peForTube; 
        weight = 0;
        if (diffuserProfile) 
        {
          weight = pmtType==0 ? ledweight_type0[PMT_id] : ledweight_type1[PMT_id];   
          nPE *= pmtType==0 ? ledweight_type0[PMT_id] : ledweight_type1[PMT_id];  
        }
        if (pmtType==0) hitRate_pmtType0->Fill();
        if (pmtType==1) hitRate_pmtType1->Fill();


      } // End of loop over Cherenkov hits
      if(verbose) cout << "Total Pe : " << totalPe << std::endl; //", total hit : " << totalHit << endl;
    }

    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    if(hybrid) wcsimrootsuperevent2->ReInitialize();
    
  } // End of loop over events

  outfile->cd();
  hitRate_pmtType0->Write();
  hitRate_pmtType1->Write();

  outfile->Close();
  
  return 0;
 }
