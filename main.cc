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
#include "TObjString.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>
#include <fftw3.h>
#include "TVirtualFFT.h"
using namespace std;

struct Particle{
  Float_t x0,y0,z0,px0,py0,pz0;
  Float_t x1,y1,z1,px1,py1,pz1;
  Float_t Einit,wepl;
  Float_t WET_prob, Y_prob, Z_prob;
  Int_t   Id, dEE, MaxEnergyTransFilter, ThresholdFilter, TS, include;

};
vector<double> Energy;
vector<double> dEdXBins;
double findWET(double, double);
void ComputeSplinePrior(Particle *, TProfile3D*, TH3D*, double);
void ComputeSpline(Particle *, TProfile3D*);
void findBin(double,double,double,int*,int*);
std::string* part_name;
std::string* proc_name;
double Xmin, Xmax;
double Ymin, Ymax;
double Zmin, Zmax;
TH3D* CT;
TH3D* Hull;
int NbinsX  = 320;
int NbinsY  = 320;
int NbinsZ  = 20; 
double dPhi, dTheta;
int main(int argc, char** argv){
  if(argc <4){
    cout<<"Please enter: "<<endl;
    cout<<"1. the projection file"<<endl;
    cout<<"2. Filter type - (1) none, (2) ThreeSigma, (3) Prior Filter"<<endl;
    cout<<"3. whether an image (0) or a noise (1) reconstruction is required"<<endl;
    cout<<"4. Which tomographic reconstruction Fourier Filter to use: Cosine Shepp-Logan Hamming Hann"<<endl;
    return 0;
  }
  int prior_use =0; 
  float P_t = 0; //Threshold for prior filter
  int prior_arg = 5; //The argument in which the prior file lies
  Particle Point;
  TString* phaseFileName = new TString(argv[1]);
  int filter = atoi(argv[2]);
  int noise = atoi(argv[3]);
  string CT_filter = argv[4];
  if(filter == 3 && argc >6){ //Prior filter + hull detection
    prior_use = 1;
    P_t = atof(argv[5]);
    prior_arg = 6;
    cout<<"Prior Filter threshold: "<<P_t<<endl;    
  }
  if(filter != 3 && argc >5){ //Not prior filter + hull detection
    prior_use = 1;
  }
  if(filter == 3 && argc <7){ //Prior Filter without hull detection
    P_t=atof(argv[5]);
    cout<<"Prior Filter threshold: "<<P_t<<endl;
  }
  
  if(prior_use == 1){
  cout<<"using prior"<<endl;
  ///--------------------------------------
  char* phantomFileName = Form("%s",argv[prior_arg]);    
  cout<<phantomFileName<<endl;
  // Load Phantom data
  TFile* phantomFile = new TFile(phantomFileName,"update");
  CT = (TH3D*)phantomFile->Get("RSP");
  Hull = (TH3D*)CT->Clone("Hull");
  Hull->Reset();
  ///--------------------------------------
  }
  TObjArray *x = phaseFileName->Tokenize("_");
  TObjString* ang_s = (TObjString*)x->At(x->GetLast());
  double ang = atof(ang_s->String());
  cout<<argv[1]<<" "<<ang<<endl;
  // Load Particle data
  TFile* phaseFile = new TFile(phaseFileName->Data(),"update");
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
  t->SetBranchAddress("dEEFilter", &Point.dEE);
  t->SetBranchAddress("ThresholdFilter", &Point.ThresholdFilter);
  t->SetBranchAddress("MaxEnergyTransFilter", &Point.MaxEnergyTransFilter);
  t->SetBranchAddress("ThreeSigma_30", &Point.TS);
  t->SetBranchAddress("WET_prob", &Point.WET_prob);
  t->SetBranchAddress("Y_prob", &Point.Y_prob);
  t->SetBranchAddress("Z_prob", &Point.Z_prob);
  //t->SetBranchAddress("Einit",&Point.Einit);
  t->SetBranchAddress("wepl",&Point.wepl);
  t->SetBranchAddress("Include_Particle",&Point.include);
  //t->SetBranchAddress("proc_name",&proc_name);
  //t->SetBranchAddress("part_name",&part_name);
  t->GetEntry(0);

  Xmin = -100; Xmax = 100;//Xmin =t->GetMinimum("x0") ;   Xmax = t->GetMaximum("x1") ; 
  Ymin = -100; Ymax = 100;//Ymin = t->GetMinimum("y1") ;   Ymax = t->GetMaximum("y1") ; 
  Zmin = -10; Zmax = 10;//Zmin = t->GetMinimum("z1") ;   Zmax = t->GetMaximum("z1") ; 
  
  std::string line;
  std::ifstream SPWater ("dEdX/Water_Geant4_He.dat");
  double data[3];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<3;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }

  int NEntries = t->GetEntries();
  //Fill the Distance Driven Binning
  TProfile3D* ddb= new TProfile3D(Form("ddb"),Form("ddb"),NbinsX, Xmin, Xmax, NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax,"s");
  TProfile3D* ddb_filtered= new TProfile3D(Form("ddb_filtered"),Form("ddb_filtered"),NbinsX, Xmin, Xmax, NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  srand(20); // seed random number generator
  for(int i=0;i<NEntries;i++){//Loop over all protons
  t->GetEntry(i);
  if(i%20000 == 0) cout<<i<<endl;
  //if(Point.include == 1){ 

      if(prior_use == 1){
        if(filter == 0){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSplinePrior(&Point, ddb, CT, ang);
            }
         }
        

        if(filter == 1){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter){
            if(Point.wepl >0 && Point.wepl <260){
	      ComputeSplinePrior(&Point, ddb, CT, ang);
            }
          }
        }
        if(filter == 2){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter && Point.TS ==1){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSplinePrior(&Point, ddb, CT, ang);
            }
          }
        }

        if(filter == 3){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter && Point.WET_prob >P_t && Point.Y_prob >P_t && Point.Z_prob >P_t){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSplinePrior(&Point, ddb, CT, ang);
            }
          }
        }
      }

      else{
        if(filter == 0){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSpline(&Point, ddb);
            }
         }
        
        if(filter == 1){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSpline(&Point, ddb);
            }
          }
        }

        if(filter == 2){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter &&  Point.TS ==1){
            if(Point.wepl >0 && Point.wepl <260){
              ComputeSpline(&Point, ddb);
            }
          }
        }

        if(filter == 3){
          if (Point.dEE && Point.ThresholdFilter && Point.MaxEnergyTransFilter && Point.WET_prob >P_t && Point.Y_prob >P_t && Point.Z_prob >P_t){
          if(Point.wepl >0 && Point.wepl <260){
              ComputeSpline(&Point, ddb);
            }
          }
        }
      }

  }
  //}
  ddb->Write("ddb",TObject::kOverwrite);
  ddb_filtered->Write("",TObject::kOverwrite);
  //Prepare the plan for the data
  cout<<"Prepare the FFT Plan"<<endl;
  int N = max(64, int( pow(2,ceil(log2(2*NbinsY+1))))) ; // for padding
  fftw_complex *out_fft, *in_ifft, *out_filter;
  double       *in_fft , *out_ifft, *in_filter;  
  fftw_plan fft, ifft, fft_filter;
  in_fft           = (double*) fftw_malloc(sizeof(double) * N);
  in_filter        = (double*) fftw_malloc(sizeof(double) * N);  
  out_ifft         = (double*) fftw_malloc(sizeof(double) * N);    
  out_fft          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1)); // only half the frequencies because it's real-only input
  in_ifft          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1)); // Negative frequencies are symmetrical and redundant
  out_filter       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

  fft              = fftw_plan_dft_r2c_1d(N, in_fft, out_fft, FFTW_ESTIMATE);
  fft_filter       = fftw_plan_dft_r2c_1d(N, in_filter, out_filter, FFTW_ESTIMATE);  
  ifft             = fftw_plan_dft_c2r_1d(N, in_ifft, out_ifft, FFTW_ESTIMATE);

  // Filter in Fourier Frequency (from real data to minimize artefact)
  //------------------------------------------------------------------------------------------------------------------
  cout<<"Define the Filter:"<<CT_filter<<endl;

  TH1D* Filter_Real    = new TH1D("Filter_Real", "Filter_Real", N, 0,N); // To save
  TH1D* Filter_Fourier = new TH1D("Filter_Fourier", "Filter_Fourier", (N/2+1), 0, (N/2+1)); // To save  

  // Ramp Filter
    in_filter[0] = 0.25;
    for(int i =1; i<N; i++){
    if(i%2==0) in_filter[i] = 0;
      else{
        if(i<(N/2+1)) in_filter[i]= -1 / pow( M_PI*i,2);
        else in_filter[i] = -1 / pow( M_PI*(N-i),2);
      }
    } //Checked

    fftw_execute(fft_filter); // Transform to Fourier domain
    for(int i=0; i<(N/2+1); i++){
       out_filter[i][0] = 2*out_filter[i][0];
      
    }    
  
  if(CT_filter == "Cosine"){
    for(int i=1; i<(N/2)+1; i++){
     if( i <= N/2+1) out_filter[i][0] = out_filter[i][0]*cos(((double)i*M_PI)/(N));
     else{
       out_filter[i][0] = 0;
       //out_filter[i][1] = 0;
      } 
    }    
  }
  if(CT_filter =="Shepp-Logan"){
    for(int i=1; i<(N/2)+1; i++){ 
      out_filter[i][0]= out_filter[i][0]*((N)/M_PI)*sin(((double)i*M_PI)/N);
    }
  }

  if(CT_filter == "Hamming"){
    for(int i=1; i<(N/2)+1; i++){ 
      out_filter[i][0]= out_filter[i][0]*(0.54+(0.46*cos(2*((double)i*M_PI)/N)));
    }
  }
  if(CT_filter == "Hann"){
    for(int i=1; i<(N/2)+1; i++){ 
      out_filter[i][0]= out_filter[i][0]*(0.5+(0.5*cos(2*((double)i*M_PI)/N)));
    }
  }
  // Plot the filter funcions
  for(int i =0; i<=N; i++) Filter_Real->SetBinContent(i+1, in_filter[i]);

  for(int i=0; i<(N/2+1); i++) Filter_Fourier->SetBinContent(i+1,out_filter[i][0]);
  

  // Fourier Transform the DDB, filter and perform the inverse transform
  //------------------------------------------------------------------------------------------------------------------
  TVector3 Source_Position(0,0,0); // At infinitiy
  TVector3 Source_Ref_Position(0,0,0); // Position of the voxel in the source referential
  TVector3 Lab_Ref_Position(0,0,0); // Position of the voxel in the lab referential
  TVector3 uX(1,0,0);     
  TVector3 uY(0,1,0);     
  TVector3 uZ(0,0,1);
  double rad = ang*TMath::Pi()/180.;
  Source_Position.RotateZ(rad);
  uX.RotateZ(rad);
  uY.RotateZ(rad);  
 // TH3D* ddb_fourier = new TH3D(Form("ddb_fourier"),Form("ddb_fourier"),(N/2+1),0,(N/2+1),(N/2+1),0,(N/2+1),(N/2+1),0,N/2+1);  
  cout<<rad<<" Filter Slice by Slice"<<endl;
  double Norm = 1./N;  // Inverse transform scale the value by 1/N 
  for(int ix=0; ix<NbinsX;ix++){ // Cross through X  
    for(int iz=0; iz<NbinsZ;iz++){ // Rotate through Z
      
      // Do the Fourier Filtering
      if(noise == 0){
        for(int iy=0; iy<NbinsY-1; iy++) in_fft[iy] = ddb->GetBinContent(ix,iy,iz);     // Fill the fft array
        fftw_execute(fft); // Transform to Fourier Domain
        for(int j=0; j<(N/2+1); j++){
        
          in_ifft[j][0] = out_fft[j][0]*out_filter[j][0]; // Real 
	  in_ifft[j][1] = out_fft[j][1]*out_filter[j][0];  // Imaginary
         
     
        }
      }
      if(noise == 1){
        for(int iy=0; iy<NbinsY-1; iy++) in_fft[iy] = pow(ddb->GetBinError(ix,iy,iz),2);     // Fill the fft array
        fftw_execute(fft); // Transform to Fourier Domain
        for(int j=0; j<(N/2+1); j++){
          in_ifft[j][0] = out_fft[j][0]*pow(out_filter[j][0],2); // Real 
          in_ifft[j][1] = out_fft[j][1]*pow(out_filter[j][0],2); // Imaginary
        }
      }
      fftw_execute(ifft); // Inverse transform
      for(int iy=0; iy<NbinsY; iy++){ // Rotate and Fill
	//ddb_filtered->SetBinContent(ix,iy,iz, Norm*out_ifft[iy]); //Normalize -- Checked	  
	Source_Ref_Position.SetX(ddb_filtered->GetXaxis()->GetBinCenter(ix));
	Source_Ref_Position.SetY(ddb_filtered->GetYaxis()->GetBinCenter(iy));      
	Source_Ref_Position.SetZ(ddb_filtered->GetZaxis()->GetBinCenter(iz));
	ddb_filtered->Fill(Source_Ref_Position.x(), Source_Ref_Position.y(), Source_Ref_Position.z(),(Norm*out_ifft[iy]));

      }
    }
  }
  /*

  ****This section was part of an experiment in interpolation that may be revisited later*****

  TH3F* ddb_rotated = new TH3F("ddb_rotated","ddb_rotated",430,-215,215,430,-215,215,430,-215,215);
  cout<<"Starting Rotation and Interpolation"<<endl;

  
  TVector3 ddb_rot_center(0,0,0);
  TVector3 ddb_filtered_rot(0,0,0);
  TVector3 ZeroBin(0,0,0);
  TVector3 LastBin(0,0,0);
  ZeroBin.SetX(ddb_filtered->GetXaxis()->GetBinCenter(1));
  ZeroBin.SetY(ddb_filtered->GetYaxis()->GetBinCenter(1));
  ZeroBin.SetZ(ddb_filtered->GetZaxis()->GetBinCenter(1));
  LastBin.SetX(ddb_filtered->GetXaxis()->GetBinCenter(NbinsX));
  LastBin.SetY(ddb_filtered->GetYaxis()->GetBinCenter(NbinsY));
  LastBin.SetZ(ddb_filtered->GetZaxis()->GetBinCenter(NbinsZ));

  for(int ix=0; ix<NbinsX;ix++){ // Cross through X  
    for(int iz=0; iz<NbinsZ;iz++){ // Rotate through Z
      for(int iy=0; iy<NbinsY; iy++){

        Source_Ref_Position.SetX(ddb_filtered->GetXaxis()->GetBinCenter(ix));
        Source_Ref_Position.SetY(ddb_filtered->GetYaxis()->GetBinCenter(iy));
        Source_Ref_Position.SetZ(ddb_filtered->GetZaxis()->GetBinCenter(iz));
        Source_Ref_Position.RotateZ(-1*rad);

        int binx=ddb_rotated->GetXaxis()->FindBin(Source_Ref_Position.x());
        int biny=ddb_rotated->GetYaxis()->FindBin(Source_Ref_Position.y());
        int binz=ddb_rotated->GetZaxis()->FindBin(Source_Ref_Position.z());

        ddb_rot_center.SetX(ddb_rotated->GetXaxis()->GetBinCenter(binx));
        ddb_rot_center.SetY(ddb_rotated->GetYaxis()->GetBinCenter(biny));
        ddb_rot_center.SetZ(ddb_rotated->GetZaxis()->GetBinCenter(binz));


        ddb_filtered_rot.SetX(ddb_rotated->GetXaxis()->GetBinCenter(binx));
        ddb_filtered_rot.SetY(ddb_rotated->GetXaxis()->GetBinCenter(biny));
        ddb_filtered_rot.SetZ(ddb_rotated->GetZaxis()->GetBinCenter(binz));
        ddb_filtered_rot.RotateZ(rad);

        if(ddb_filtered_rot.x() > ZeroBin.x() && ddb_filtered_rot.x() < LastBin.x() &&
           ddb_filtered_rot.y() > ZeroBin.y() && ddb_filtered_rot.y() < LastBin.y() &&
           ddb_filtered_rot.z() > ZeroBin.z() && ddb_filtered_rot.z() < LastBin.z()){


          ddb_rotated->Fill(ddb_rot_center.x(), ddb_rot_center.y(), ddb_rot_center.z(),ddb_filtered->Interpolate(ddb_filtered_rot.x(),ddb_filtered_rot.y(),ddb_filtered_rot.z()));
        }

      }
    }
  }*/

                                                                     
 
  Filter_Real->Write("",TObject::kOverwrite);
  Filter_Fourier->Write("",TObject::kOverwrite);
  ddb_filtered->Write("",TObject::kOverwrite);
  Hull->Write("",TObject::kOverwrite);
  //ddb_rotated->Write("",TObject::kOverwrite);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);  
  fftw_free(in_fft); fftw_free(out_fft);
  fftw_free(in_ifft); fftw_free(out_ifft);  
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
  for(int i = it_Estop; i<it_Einit; i++) WET += 0.1/dEdXBins[i];
  return WET;
}
////////////////////////////////////////////
// Compute Spline Prior
////////////////////////////////////////////
void ComputeSplinePrior(Particle *Point, TProfile3D* ddb,TH3D* RSPMap, double ang){
  TVector3 p0(Point->px0,Point->py0,Point->pz0);
  TVector3 p1(Point->px1,Point->py1,Point->pz1);
  TVector3 m0(Point->x0,Point->y0,Point->z0);
  TVector3 m1(Point->x1,Point->y1,Point->z1);
  Point->Einit=512;
  double TrackLength = TVector3(m1-m0).Mag();
  int NStep = TrackLength*4;//512;
  TVector3 m_entry, m_exit; // position of the path at the entry and exit
  double t, t_entry = 0., t_exit = 1.; // fraction of the path at which point the spline enter the Hull
  TVector3 m, mi, mf;
  mi = m0;
  mf = m1;

  double Xmax = ddb->GetXaxis()->GetXmax();
  double Xmin = ddb->GetXaxis()->GetXmin();
  double Ymax = ddb->GetYaxis()->GetXmax();
  double Ymin = ddb->GetYaxis()->GetXmin();
  double Zmax = ddb->GetZaxis()->GetXmax();
  double Zmin = ddb->GetZaxis()->GetXmin();

  TVector3 origin(0, 0, 0); // cm
  m0 -= origin; m1 -= origin;  
  // Negative rotation because we change the coordinates and not the phantom
  m0.RotateZ(-1*(ang)*M_PI/180.);
  m1.RotateZ(-1*(ang)*M_PI/180.);
  p0.RotateZ(-1*(ang)*M_PI/180.);
  p1.RotateZ(-1*(ang)*M_PI/180.);
  m0 += origin; m1 += origin;
  p0.SetMag(TVector3(m1-m0).Mag()); 
  p1.SetMag(TVector3(m1-m0).Mag());
  int recon_radius = 80;
  double RSP = 0;
  int outside = 1;
  //Propagate from the entrance to the Hull
  for(int k =1; k<NStep-1; k++){
    t_entry = double(k)/NStep;
    m_entry = m0 + t_entry*p0;
    if (sqrt(pow(m_entry.x(),2)+pow(m_entry.y(),2))<recon_radius){
      int global = RSPMap->FindBin(m_entry.x(), m_entry.y(),m_entry.z());
      RSP = RSPMap->GetBinContent(global);
      outside = 0;
    }
    else RSP = 0;
    if(outside == 0.0){
      t_entry = (double(k))/NStep; //step back 10mm from hull to avoid issues at the edges
      m_entry = m0 + t_entry*p0;
      break;
    }
  }
  outside =1; 
  // Retro Propagate from the exit to the Hull
  for(int k =NStep-1; k>=1; k--){
    t_exit = double(k)/NStep;
    m_exit = m1 - (1-t_exit)*p1;
    if (sqrt(pow(m_exit.x(),2)+pow(m_exit.y(),2))<recon_radius){
      int global = RSPMap->FindBin(m_exit.x(),m_exit.y(),m_exit.z());
      RSP = RSPMap->GetBinContent(global);
      outside = 0;
    }
    else RSP = 0;
    if(outside == 0){
      t_exit = (double(k))/NStep; //step back 10mm from hull to avoid issues at the edges
      m_exit = m1 - (1-t_exit)*p1;
      break;
    }
  }

  if(outside == 1) { // No Hull we are in air
    m_exit = m1;
    m_entry = m0;
    t_entry = 0.;
    t_exit = 1.;
  }

  m_entry -= origin; m_exit -= origin;  
  // Negative rotation because we change the coordinates and not the phantom

  m_entry.RotateZ(1*(ang)*M_PI/180.);
  m_exit.RotateZ(1*(ang)*M_PI/180.);
  m1.RotateZ(1*(ang)*M_PI/180.);
  m0.RotateZ(1*(ang)*M_PI/180.);
  p0.RotateZ(1*(ang)*M_PI/180.);
  p1.RotateZ(1*(ang)*M_PI/180.);
  m0 += origin; m1 += origin;

  /// For sanity
  Hull->Fill(m_entry.x(), m_entry.y(), m_entry.z());
  Hull->Fill(m_exit.x(), m_exit.y(), m_exit.z());  
  double WEPL   = 260;//findWET(Point->Einit,0.1);
  double WET    = (double) Point->wepl;
  double alpha1 = 1.01+0.43*pow(WET/WEPL,2);
  double alpha2 = 0.99-0.46*pow(WET/WEPL,2);
  double HullLength = TVector3(m_exit-m_entry).Mag(); 

  p0.SetMag(HullLength*alpha1);
  p1.SetMag(HullLength*alpha2);

  TVector3 A,B,C,D;
  A       =    2*m_entry - 2*m_exit + p0+p1;
  B       =   -3*m_entry + 3*m_exit - 2*p0-p1;
  C       =    p0;
  D       =    m_entry;
  
  std::map<pair<int,int>,double> Lengthmap;
  std::pair<std::map<std::pair<int,int>,double>::iterator,bool> ret;
  int binx,biny,binz;
  int step_start = (int) (t_entry*NStep);
  int step_end = (int) (t_exit*NStep); 
  int NStep_hull = step_end - step_start;
  double wepl_step;
  double Length = TVector3(m0-m1).Mag();
  if( NStep != NStep_hull){
    p0.SetMag(TVector3(m_entry-m0).Mag()); 
    p1.SetMag(TVector3(m1-m_exit).Mag());
     
  //Before the phantom - straight line in air
    for(int i=0;i<step_start;i++){
      t = double(i)/step_start;
      m = m0 + (t*p0);
      if(m.x() >= Xmin && m.x() <= Xmax && m.y() >= Ymin && m.y() <= Ymax && m.z() >= Zmin && m.z() <= Zmax){
        ddb->Fill(m.x(), m.y(), m.z(), WET);
      }
    }

    // In the phantom - MLP in water
    for(int i=0;i<=NStep_hull;i++){
      t = double(i)/NStep_hull;
      m = (D+t*(C+t*(B+t*A)));
      if(m.x() >= Xmin && m.x() <= Xmax && m.y() >= Ymin && m.y() <= Ymax && m.z() >= Zmin && m.z() <= Zmax){
        ddb->Fill(m.x(), m.y(), m.z(), WET);
      }
    }  

    // After the phantom - Straight line in air
    for(int i=0;i<(NStep - step_end);i++){
      t = double(i)/(NStep - step_end);  
      m = m1 - (t*p1);
      if(m.x() >= Xmin && m.x() <= Xmax && m.y() >= Ymin && m.y() <= Ymax && m.z() >= Zmin && m.z() <= Zmax){
        ddb->Fill(m.x(), m.y(), m.z(), WET);
      }
    }
  }
  else{//doesn't enter hull - straight line
    p0.SetMag(TVector3(m1-m0).Mag()); 
    for(int i=0;i<=NStep;i++){
      t = double(i)/NStep;
      m = m0 + t*p0;    
      if(m.x() >= Xmin && m.x() <= Xmax && m.y() >= Ymin && m.y() <= Ymax && m.z() >= Zmin && m.z() <= Zmax){
        ddb->Fill(m.x(), m.y(), m.z(), WET);
      }
    }
  }
}

////////////////////////////////////////////
//// Compute Spline
//////////////////////////////////////////////
void ComputeSpline(Particle *Point, TProfile3D* ddb){
  TVector3 p0(Point->px0,Point->py0,Point->pz0);
  TVector3 p1(Point->px1,Point->py1,Point->pz1);
  TVector3 m0(Point->x0,Point->y0,Point->z0);
  TVector3 m1(Point->x1,Point->y1,Point->z1);
  int NStep = 512;
  TVector3 m, mi, mf;
  mi = m0;
  mf = m1;
  double Xmax = ddb->GetXaxis()->GetXmax();
  double Xmin = ddb->GetXaxis()->GetXmin();
  double Ymax = ddb->GetYaxis()->GetXmax();
  double Ymin = ddb->GetYaxis()->GetXmin();
  double Zmax = ddb->GetZaxis()->GetXmax();
  double Zmin = ddb->GetZaxis()->GetXmin();
  double WEPL   = 260;
  double WET    = Point->wepl;
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
  double wepl_step;
  for(int i=0;i<NStep;i++){
    double t=double(i)/NStep;
    m = D+t*(C+t*(B+t*A));
    if(m.x() >= Xmin && m.x() <= Xmax && m.y() >= Ymin && m.y() <= Ymax && m.z() >= Zmin && m.z() <= Zmax){
      ddb->Fill(m.x(), m.y(), m.z(), WET);
    }
  }
}

