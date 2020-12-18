#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <math.h>       
#include <cmath>        
#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "TVector3.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH3D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile3D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>
#include <fftw3.h>
#include "TVirtualFFT.h"
using namespace std;
#define NProj 1

struct Particle{
  Float_t x0,y0,z0,px0,py0,pz0;
  Float_t x1,y1,z1,px1,py1,pz1;
  Float_t Einit,Estop;
  Int_t   Id;

};

vector<double> Energy;
vector<double> dEdXBins;

double findWET(double, double);
void ComputeSpline(Particle *, TProfile3D*);
void findBin(double,double,double,int*,int*);
std::string* part_name;
std::string* proc_name;

double Xmin, Xmax;
double Ymin, Ymax;
double Zmin, Zmax;
TH3D* CT;
int NbinsX  = 300;
int NbinsY  = 300;
int NbinsZ  = 300; 
double dPhi, dTheta;
int main(int argc, char** argv){
  Particle Point;
  char* phaseFileName = Form("%s",argv[1]);
  /*///--------------------------------------
  char* phantomFileName = Form("%s",argv[2]);    
  // Load Phantom data
  TFile* phantomFile = new TFile(phantomFileName,"update");
  CT = (TH3D*)phantomFile->Get("rsp");
  ///--------------------------------------
  */
  
  // Load Particle data
  TFile* phaseFile = new TFile(phaseFileName,"update");
  TTree* t = (TTree*)phaseFile->Get("phase");
  t->SetBranchAddress("x0",&Point.x0);
  t->SetBranchAddress("y0",&Point.y0);
  t->SetBranchAddress("z0",&Point.z0);

  t->SetBranchAddress("px0",&Point.px0);
  t->SetBranchAddress("py0",&Point.py0);
  t->SetBranchAddress("pz0",&Point.pz0);

  t->SetBranchAddress("x1",&Point.x1);
  t->SetBranchAddress("y1",&Point.y1);
  t->SetBranchAddress("z1",&Point.z1);

  t->SetBranchAddress("px1",&Point.px1);
  t->SetBranchAddress("py1",&Point.py1);
  t->SetBranchAddress("pz1",&Point.pz1);

  t->SetBranchAddress("Einit",&Point.Einit);
  t->SetBranchAddress("Estop",&Point.Estop);
  t->SetBranchAddress("proc_name",&proc_name);
  t->SetBranchAddress("part_name",&part_name);
  t->GetEntry(0);

  Xmin = -150; Xmax = 150;//Xmin =t->GetMinimum("x0") ;   Xmax = t->GetMaximum("x1") ; 
  Ymin = -150; Ymax = 150;//Ymin = t->GetMinimum("y1") ;   Ymax = t->GetMaximum("y1") ; 
  Zmin = -150; Zmax = 150;//Zmin = t->GetMinimum("z1") ;   Zmax = t->GetMaximum("z1") ; 
  
  std::string line;
  std::ifstream SPWater ("dEdX/Water_Geant4_P.dat");
  double data[3];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<3;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }

  int NEntries = t->GetEntries();
  //Fill the Distance Driven Binning
  TProfile3D* ddb= new TProfile3D(Form("ddb"),Form("ddb"),NbinsX, Xmin, Xmax, NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  for(int i=0;i<100;i++){//Loop over all protons
    t->GetEntry(i);
    if(i%20000 == 0) cout<<i<<endl;
    if(part_name->compare("proton")==0){
      ComputeSpline(&Point,ddb);
    }
  }  
  //prepare the plan
  int N = max(64, int( pow(2,ceil(log2(2*NbinsY))))) ; // for padding
  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  ddb->Write("",TObject::kOverwrite);
  for(int ix=0; ix<NbinsX;ix++){
    for(int iz=0; iz<NbinsZ;iz++){
      

    }
  }

  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
  phaseFile->Close();
  return 0;
}

////////////////////////////////////////////
// Extract WET
////////////////////////////////////////////
double findWET(double Einit,double Estop){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
  int it_Estop = lower_bound(Energy.begin(), Energy.end(), Estop)-Energy.begin();
  double WET   = 0 ;
  for(int i = it_Estop; i<it_Einit; i++){
    //WET += 0.01/dEdXBins[i];
    WET += 0.1/dEdXBins[i];
  }
  return WET;
}
////////////////////////////////////////////
// Compute Spline
////////////////////////////////////////////
void ComputeSpline(Particle *Point, TProfile3D* ddb){
  TVector3 p0(Point->px0,Point->py0,Point->pz0);
  TVector3 p1(Point->px1,Point->py1,Point->pz1);
  TVector3 m0(Point->x0,Point->y0,Point->z0);
  TVector3 m1(Point->x1,Point->y1,Point->z1);
  int NStep = 512;
  TVector3 m, mi, mf;
  mi = m0;
  mf = m1;
  double WEPL   = findWET(Point->Einit,0.1);
  double WET    = findWET(Point->Einit,Point->Estop);
  double alpha1 = 1.01+0.43*pow(WET/WEPL,2);
  double alpha2 = 0.99-0.46*pow(WET/WEPL,2);
  double Length = TVector3(mf-mi).Mag();

  p0.SetMag(Length*alpha1);
  p1.SetMag(Length*alpha2);

  TVector3 A,B,C,D;
  A     =    2*mi - 2*mf + p0+p1;
  B     =   -3*mi + 3*mf - 2*p0-p1;
  C     =    p0;
  D     =    mi;
  double TotL    = 0;
  TVector3 m_old = mi;
  std::map<pair<int,int>,double> Lengthmap;
  std::pair<std::map<std::pair<int,int>,double>::iterator,bool> ret;
  int binx,biny,binz;

  for(int i=0;i<NStep;i++){
    double t=double(i)/NStep;
    m = D+t*(C+t*(B+t*A));
    ddb->Fill(m.x(), m.y(), m.z(), WET);

  }
}

