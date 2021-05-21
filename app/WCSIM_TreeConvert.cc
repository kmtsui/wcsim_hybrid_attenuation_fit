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

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

#include "utils/heuman_lambda.hh"

using namespace std;

WCSimRootGeom *geo = 0; 

// Simple example of reading a generated Root file
const int nPMTtypes = 2;
double PMTradius[nPMTtypes];

double CalcSolidAngle(double r, double R, double costh)
{
  //double weight = 2*TMath::Pi()*(1-R/sqrt(R*R+r*r));
  //weight *= 1.-0.5*sqrt(1-costh*costh);
  //weight *= 0.5+0.5*costh;

  //return weight;

  if (costh<0) return 0.;
  // analytic way to calculate solid angle subtended by PMT
  // r = radius of PMT, R = distance to source, costh = cos angle of source relative to PMT
  // assume the PMT is a perfect hemisphere
  // calculate the slope of line that connects the source and tangent to the PMT surface
  // then calculate the angle subtended
  double sinth = sqrt(1-costh*costh);
  double cosc;
  double costh1, costh2;
  if (R*sinth-r>0) // one of the line limited by the edge of hemisphere
  {
    // tangent line to the hemisphere
    double m1 = (R*R*sinth*costh-r*sqrt(R*R-r*r))/(R*R*sinth*sinth-r*r);
    double c1 = -m1*R*sinth+R*costh;
    double a1 = R;
    double b1 = sqrt(R*sinth*R*sinth+(R*costh-c1)*(R*costh-c1));
    costh1 = std::min((a1*a1+b1*b1-c1*c1)/2./a1/b1,1.);

    // half of the solid angle 
    double omega = TMath::Pi()*(1-costh1);

    // another limited by the edge of hemisphere
    double a2 = R;
    double b2 = sqrt(R*costh*R*costh+(R*sinth-r)*(R*sinth-r));
    double c2 = r;
    costh2 = std::min((a2*a2+b2*b2-c2*c2)/2./a2/b2,1.);

    // the source completely perpendicular to the PMT direction
    // only get the half solid angle
    if (costh2>=1)
      return omega;

    // solid angle subtended by an ellipse
    // formula copied from this paper https://www.sciencedirect.com/science/article/abs/pii/S0168900209022918
    double phi = asin(costh2);
    double k2 = 1-costh1*costh1/costh2*costh2;

    double heumanlambda = heuman_lambda(phi,k2);

    omega += TMath::Pi()*(1-heumanlambda);

    return omega;
  }
  else 
  {
    double m1 = (R*R*sinth*costh+r*sqrt(R*R-r*r))/(R*R*sinth*sinth-r*r);
    double m2 = (R*R*sinth*costh-r*sqrt(R*R-r*r))/(R*R*sinth*sinth-r*r);
    double c1 = -m1*R*sinth+R*costh;
    double c2 = R*sinth-R*costh/m2;
    double c = sqrt(c1*c1+c2*c2);
    double a = sqrt(R*sinth*R*sinth+(R*costh-c1)*(R*costh-c1));
    double b = sqrt(R*costh*R*costh+(R*sinth-c2)*(R*sinth-c2));
    cosc = (a*a+b*b-c*c)/2/a/b;

    double coshalfc = sqrt((cosc+1)/2);
    return (1-coshalfc)*2*TMath::Pi();
  }

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

int main(int argc, char **argv){
  
  char * filename=NULL;
  char * outfilename=NULL;
  bool verbose=false;
  bool hybrid = true;
  double cvacuum = 3e8 / 1e9;//speed of light, in meter per ns.
  double nindex = 1.373;//refraction index of water
  bool plotDigitized = true; //using digitized hits
  bool separatedTriggers=false;//Assume two independent triggers, one for mPMT, one for B&L

  double wavelength = 400; //wavelength in nm

  int startEvent=0;
  int endEvent=0;
  char c;
  while( (c = getopt(argc,argv,"f:o:b:s:e:l:hrtv")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
      case 'f':
        filename = optarg;
        break;
      case 'r':
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
	      break;
      case 'e':
	      endEvent = std::stoi(optarg);
	      break;
      case 'l':
        wavelength = std::stod(optarg);
        if (wavelength<0) {
          std::cout<<"Wavelength < 0, using default = 400 nm"<<std::endl;
          wavelength = 400;
        }
	      break;
      case 'h':
        std::cout << "USAGE: "
                  << argv[0] << "\nOPTIONS:\n"
                  << "-f : Intput file\n"
                  << "-o : Output file\n"
                  << "-l : Laser wavelength\n"
                  << "-b : Use only B&L PMTs\n"
                  << "-r : Use raw hits\n"
                  << "-t : Use separated triggers\n"
                  << "-v : Verbose\n"
                  << "-s : Start event\n"
                  << "-e : End event\n";
      default:
        return 0;
    }
  }
  
  double vg = CalcGroupVelocity(wavelength);
  std::cout<<"Using wavelength = "<<wavelength<<" nm, group velocity = "<<vg<<" m/s, n = "<<cvacuum/vg<<std::endl;
  vg /= 1.e9; // convert to m/ns

  
  TFile *file = TFile::Open(filename);
  // Open the file
  if (filename==NULL){
    std::cout << "Error, no input file: " << std::endl;
    return -1;
  }
  if (!file->IsOpen()){
    std::cout << "Error, could not open input file: " << filename << std::endl;
    return -1;
  }
  
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

  double nHits, nPE, dist, costh, timetof, cosths, omega;
  int PMT_id, mPMT_PMTNo; //mPMT_id
  // TTree for storing the hit information. One for B&L PMT<, one for mPMT
  TTree* hitRate_pmtType0 = new TTree("hitRate_pmtType0","hitRate_pmtType0");
  hitRate_pmtType0->Branch("nHits",&nHits); // dummy variable, always equal to 1
  hitRate_pmtType0->Branch("nPE",&nPE); // number of PE
  hitRate_pmtType0->Branch("R",&dist); // distance to source
  hitRate_pmtType0->Branch("costh",&costh); // photon incident angle relative to PMT
  hitRate_pmtType0->Branch("cosths",&cosths); // PMT angle relative to source
  hitRate_pmtType0->Branch("omega",&omega); // solid angle subtended by PMT
  hitRate_pmtType0->Branch("timetof",&timetof); // hittime-tof
  hitRate_pmtType0->Branch("PMT_id",&PMT_id);
  TTree* hitRate_pmtType1 = new TTree("hitRate_pmtType1","hitRate_pmtType1");
  hitRate_pmtType1->Branch("nHits",&nHits);
  hitRate_pmtType1->Branch("nPE",&nPE);
  hitRate_pmtType1->Branch("R",&dist);
  hitRate_pmtType1->Branch("costh",&costh);
  hitRate_pmtType1->Branch("cosths",&cosths);
  hitRate_pmtType1->Branch("omega",&omega);
  hitRate_pmtType1->Branch("timetof",&timetof);
  hitRate_pmtType1->Branch("PMT_id",&PMT_id);
  hitRate_pmtType1->Branch("mPMT_PMTNo",&mPMT_PMTNo); //sub-ID of PMT inside a mPMT module

  double vtxpos[3];
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

    for (int i=0;i<3;i++) vtxpos[i]=wcsimrootevent->GetVtx(i);

    // LI source direction is always perpendicular to the wall
    // Separate treatment for barrel and endcap 
    double vDirSource[3];
    double endcapZ=3000;
    if (abs(vtxpos[2])<endcapZ) {
      vDirSource[0]=vtxpos[0];
      vDirSource[1]=vtxpos[1];
      double norm = sqrt(vtxpos[0]*vtxpos[0]+vtxpos[1]*vtxpos[1]);
      vDirSource[0]/=-norm;
      vDirSource[1]/=-norm;
      vDirSource[2]=0;
    } else {
      vDirSource[0]=0;
      vDirSource[1]=0;
      if (vtxpos[2]>endcapZ) vDirSource[2]=-1;
      else vDirSource[2]=1;
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
        //int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
        int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);

        WCSimRootPMT pmt;
        if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
        else pmt  = geo->GetPMT(tubeNumber-1,true);

        //if (pmtType == 0) PMT_id = tubeNumber-1;
        //else mPMT_id = (tubeNumber-1.);
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

        if(pmtType == 1) mPMT_PMTNo = pmt.GetmPMT_PMTNo();
        
        WCSimRootCherenkovHitTime * HitTime = (WCSimRootCherenkovHitTime*) timeArray->At(i);//Takes the first hit of the array as the timing, It should be the earliest hit
        //WCSimRootCherenkovHitTime HitTime = (WCSimRootCherenkovHitTime) timeArray->At(j);		  
        double time = HitTime->GetTruetime();
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        timetof = time-tof;
        nHits = 1; nPE = peForTube; dist = Norm; costh = -(vDir[0]*vOrientation[0]+vDir[1]*vOrientation[1]+vDir[2]*vOrientation[2]);
        cosths = vDir[0]*vDirSource[0]+vDir[1]*vDirSource[1]+vDir[2]*vDirSource[2];
        double pmtradius = pmtType==0 ? PMTradius[0] : PMTradius[1]; 
        omega = 0; //CalcSolidAngle(pmtradius,dist,costh);
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

        //if(pmtType == 0) PMT_id = (tubeNumber-1.);
        //else mPMT_id = (tubeNumber-1.);
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

        if(pmtType == 1) mPMT_PMTNo = pmt.GetmPMT_PMTNo();

        timetof = time-tof+triggerTime[pmtType]-triggerShift[pmtType];

        nHits = 1; nPE = peForTube; dist = Norm; costh = -(vDir[0]*vOrientation[0]+vDir[1]*vOrientation[1]+vDir[2]*vOrientation[2]);
        cosths = vDir[0]*vDirSource[0]+vDir[1]*vDirSource[1]+vDir[2]*vDirSource[2];
        double pmtradius = pmtType==0 ? PMTradius[0] : PMTradius[1]; 
        omega = 0; //CalcSolidAngle(pmtradius,dist,costh);
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
  // Save also PMT geometry information
  TTree* pmt_type0 = new TTree("pmt_type0","pmt_type0");
  pmt_type0->Branch("R",&dist);
  pmt_type0->Branch("costh",&costh);
  pmt_type0->Branch("cosths",&cosths);
  pmt_type0->Branch("omega",&omega);
  pmt_type0->Branch("PMT_id",&PMT_id);
  TTree* pmt_type1 = new TTree("pmt_type1","pmt_type1");
  pmt_type1->Branch("R",&dist);
  pmt_type1->Branch("costh",&costh);
  pmt_type1->Branch("cosths",&cosths);
  pmt_type1->Branch("omega",&omega);
  pmt_type1->Branch("PMT_id",&PMT_id);
  pmt_type1->Branch("mPMT_PMTNo",&mPMT_PMTNo);
  double vDirSource[3];
  double endcapZ=3000;
  if (abs(vtxpos[2])<endcapZ) {
    vDirSource[0]=vtxpos[0];
    vDirSource[1]=vtxpos[1];
    double norm = sqrt(vtxpos[0]*vtxpos[0]+vtxpos[1]*vtxpos[1]);
    vDirSource[0]/=-norm;
    vDirSource[1]/=-norm;
    vDirSource[2]=0;
  } else {
    vDirSource[0]=0;
    vDirSource[1]=0;
    if (vtxpos[2]>endcapZ) vDirSource[2]=-1;
    else vDirSource[2]=1;
  }

  int nPMTs_type0=geo->GetWCNumPMT();
  int nPMTs_type1=0; if (hybrid) nPMTs_type1=geo->GetWCNumPMT(true);
  for (int pmtType=0;pmtType<nPMTtypes;pmtType++) 
  {
    int nPMTs_type = pmtType==0 ? nPMTs_type0 : nPMTs_type1;
    for (int i=0;i<nPMTs_type;i++) 
    {
      WCSimRootPMT pmt;
      if (pmtType==0) pmt = geo->GetPMT(i,false);
      else pmt = geo->GetPMT(i,true);
      if (pmtType == 0) PMT_id = i;
      else {
          PMT_id = i;
          mPMT_PMTNo = pmt.GetmPMT_PMTNo();
      }
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
      double pmtradius = pmtType==0 ? PMTradius[0] : PMTradius[1]; 
      omega = CalcSolidAngle(pmtradius,dist,costh);
      if (pmtType==0) pmt_type0->Fill();
      if (pmtType==1) pmt_type1->Fill();
    }
  }
  pmt_type0->Write();
  pmt_type1->Write();
  outfile->Close();
  
  return 0;
 }
