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
#include "TMatrixD.h"
#include "TVector3.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH3D.h"
#include "TF2.h"
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

struct Particle{
  Float_t x0,y0,z0,px0,py0,pz0;
  Float_t x1,y1,z1,px1,py1,pz1;
  Float_t Einit,Estop;
  Float_t WET, WEPL;
  Int_t   Id;
  Float_t angle;
};

vector<double> Energy;
vector<double> dEdXBins;

// DDB filling functions
void ComputeSplineNoPrior(Particle *, TProfile3D*);
void ComputeSplinePrior(Particle *, TH3D*, TProfile3D*);
void ComputeLikelihood(Particle *, TMatrixD&, TMatrixD&, TMatrixD&, TMatrixD&, TProfile3D*);
std::string* part_name;
std::string* proc_name;

//Math operations
TMatrixD inverse(TMatrix );
TMatrixD reverse(TMatrix );
TMatrixD Div(TMatrixD , TMatrixD);
TMatrixD Mult(TMatrixD , TMatrixD);
TMatrixD Pow(TMatrixD , float);
TMatrixD linspace(float, float, int );
double Gauss(double, double, double);
double Gauss2D(TMatrixD, TMatrixD, TMatrixD);

//MLP functions
void GatherEnergyRadLength(TMatrixD &, TMatrixD&, TMatrixD&, TMatrixD&, double, TH3D*);
void CalculateMLPSigma(Particle*, TMatrixD&, TMatrixD&, TMatrixD&, TMatrixD&, TMatrixD&, TMatrixD&);
double findWET(double, double);
double E2beta(double);
double pv(double);
double gauss(double , double, double, double, double, double , double);
double simps(TMatrixD, TMatrixD );
void simps_mult(TMatrixD , TMatrixD , double , double &, double &, double &);

double Xmin, Xmax;
double Ymin, Ymax;
double Zmin, Zmax;
TH3D* Prior;
int NbinsX  = 300;
int NbinsY  = 300;
int NbinsZ  = 300; 
double dPhi, dTheta;
int main(int argc, char** argv){
  Particle Point;
  char* phaseFileName = Form("%s",argv[1]);

  // Load Phantom data
  /*char* phantomFileName = Form("%s",argv[2]);    
  TFile* phantomFile = new TFile(phantomFileName,"update");
  Prior = (TH3D*)phantomFile->Get("rsp");*/
  
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

  // For the MLP
  TMatrixD tracks_X  = linspace(Xmin,Xmax, NbinsX+1); // mm
  TMatrixD tracks_Y(NbinsX+1, 1); // mm
  TMatrixD tracks_Z(NbinsX+1, 1); // mm
  TMatrixD tracks_L_inv(NbinsX+1, 1); // inverse radiation length -> mm
  TMatrixD tracks_E(NbinsX+1, 1); // energy -> MeV
  TMatrixD tracks_pv_inv(NbinsX+1, 1); // momentum-velocity
  TMatrixD tracks_sig(NbinsX+1, 3); // lateral position sigma at each point
  
  int NEntries = 100000;//t->GetEntries();

  //Fill the Distance Driven Binning matrix 
  TProfile3D* ddb= new TProfile3D(Form("ddb"),Form("ddb"),NbinsX, Xmin, Xmax, NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  TH3D* ddb_filtered= new TH3D(Form("ddb_filtered"),Form("ddb_filtered"),NbinsX, Xmin, Xmax, NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  for(int i=0;i<NEntries;i++){//Loop over all protons
    t->GetEntry(i);
    Point.WET    = findWET(Point.Einit,Point.Estop);
    Point.WEPL   = findWET(Point.Einit,0.1);


    if(i%20000 == 0) cout<<i<<endl;
    if(part_name->compare("proton")==0){
      
      //Single Event -- CSP
      //ComputeSplineNoPrior(&Point,ddb);
      //ComputeSplinePrior(&Point,Prior,ddb);      
      
      //Single Event -- MLP
      GatherEnergyRadLength(tracks_X, tracks_L_inv, tracks_E, tracks_pv_inv , Point.Einit, Prior);
      CalculateMLPSigma(&Point, tracks_X, tracks_L_inv, tracks_pv_inv, tracks_sig, tracks_Y, tracks_Z);
      ComputeLikelihood(&Point, tracks_X, tracks_Y, tracks_Z, tracks_sig, ddb);
      
    }
  }

  ddb->Write("ddb",TObject::kOverwrite);
  //Prepare the plan for the data
  cout<<"Prepare the FFT Plan"<<endl;
  int N = max(64, int( pow(2,ceil(log2(2*NbinsY))))) ; // for padding
  fftw_complex *out_fft, *in_ifft, *out_filter;
  double       *in_fft , *out_ifft, *in_filter;  
  fftw_plan fft, ifft, fft_filter;;
  in_fft           = (double*) fftw_malloc(sizeof(double) * N);
  in_filter        = (double*) fftw_malloc(sizeof(double) * N);  
  out_ifft         = (double*) fftw_malloc(sizeof(double) * N);    
  out_fft          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1)); // only half the frequencies because it's real-only input
  in_ifft          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1)); // Negative frequencies are symmetrical and redundant
  out_filter       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft         = fftw_plan_dft_r2c_1d(N, in_fft, out_fft, FFTW_ESTIMATE);
  fft_filter  = fftw_plan_dft_r2c_1d(N, in_filter, out_filter, FFTW_ESTIMATE);  
  ifft        = fftw_plan_dft_c2r_1d(N, in_ifft, out_ifft, FFTW_ESTIMATE);

  // Filter in Fourier Frequency (from real data to minimize artefact)
  //------------------------------------------------------------------------------------------------------------------
  cout<<"Define the Filter (Ramp so far)"<<endl;
  in_filter[0] = 0.25;
  for(int i =1; i<N; i++){
    if(i%2==0) in_filter[i] = 0;
    else{
      if(i<(N/2+1)) in_filter[i]= -1 / pow( M_PI*i,2);
      else in_filter[i] = -1 / pow( M_PI*(N-i),2);
    }
  } //Checked

  fftw_execute(fft_filter); // Transform to Fourier domain
  for(int i=0; i<(N/2+1); i++) out_filter[i][0] = 2*out_filter[i][0];

  // Fourier Transform the DDB, filter and perform the inverse transform
  //------------------------------------------------------------------------------------------------------------------

  TVector3 Source_Position(0,0,0); // At infinitiy
  TVector3 Source_Ref_Position(0,0,0); // Position of the voxel in the source referential
  TVector3 Lab_Ref_Position(0,0,0); // Position of the voxel in the lab referential

  TVector3 uX(1,0,0);     
  TVector3 uY(0,1,0);     
  TVector3 uZ(0,0,1);

  // Here we should do the loop -- change it for parallelization
  int NProj = 30;
  int Step  = (int) 360./NProj;
  for(int deg=0; deg<360; deg+=Step){
    double rad = float(deg)*TMath::Pi()/180.;
    Source_Position.RotateZ(rad);
    uX.RotateZ(rad);
    uY.RotateZ(rad);  

    cout<<rad<<" Filter Slice by Slice"<<endl;
    double Norm = 1./N;  // Inverse transform scale the value by 1/N 
    for(int ix=0; ix<NbinsX;ix++){ // Cross through X
      for(int iz=0; iz<NbinsZ;iz++){ // Rotate through Z

	// Do the Fourier Filtering
	for(int iy=0; iy<NbinsY; iy++) in_fft[iy] = ddb->GetBinContent(ix,iy,iz);     // Fill the fft array
	fftw_execute(fft); // Transform to Fourier Domain
	for(int j=0; j<(N/2+1); j++){
	  in_ifft[j][0] = out_fft[j][0]*out_filter[j][0]; // Real 
	  in_ifft[j][1] = out_fft[j][1]*out_filter[j][0]; // Imaginary
	}
	fftw_execute(ifft); // Inverse transform
	for(int iy=0; iy<NbinsY; iy++){ // Rotate and Fill
	  //ddb_filtered->SetBinContent(ix,iy,iz, Norm*out_ifft[iy]); //Normalize -- Checked	  

	  Source_Ref_Position.SetX(ddb_filtered->GetXaxis()->GetBinCenter(ix)); // This could be adjusted before for improved speed
	  Source_Ref_Position.SetY(ddb_filtered->GetYaxis()->GetBinCenter(iy));      
	  Source_Ref_Position.SetZ(ddb_filtered->GetZaxis()->GetBinCenter(iz));
	  Source_Ref_Position.RotateZ(-1*rad); // In the inertial referential need to rotate backward

	  // For Cone Beam -- NOT WORKING!
	  //Lab_Ref_Position.SetZ( (Source_Ref_Position - Source_Position).Dot(uZ) );
	  //Lab_Ref_Position.SetX( (Source_Ref_Position - Source_Position).Dot(uX) );
	  //Lab_Ref_Position.SetY( (Source_Ref_Position - Source_Position).Dot(uY) );
	  //ddb_filtered->Fill(Lab_Ref_Position.x(), Lab_Ref_Position.y(), Lab_Ref_Position.z(), Norm*out_ifft[iy]); //Normalized -- Checked
	  
	  ddb_filtered->Fill(Source_Ref_Position.x(), Source_Ref_Position.y(), Source_Ref_Position.z(), Norm*out_ifft[iy]); //Normalized -- Checked
	}
      }
    }
  }
  ddb_filtered->Scale(TMath::Pi()/(2.*NProj));
  ddb_filtered->Write("",TObject::kOverwrite);
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
// Compute Spline without a prior (initial)
////////////////////////////////////////////
void ComputeSplineNoPrior(Particle *Point, TProfile3D* ddb){
  TVector3 p0(Point->px0,Point->py0,Point->pz0);
  TVector3 p1(Point->px1,Point->py1,Point->pz1);
  TVector3 m0(Point->x0,Point->y0,Point->z0);
  TVector3 m1(Point->x1,Point->y1,Point->z1);
  int NStep = 512;
  TVector3 m, mi, mf;
  mi = m0;
  mf = m1;
  double alpha1 = 1.01+0.43*pow(Point->WET/Point->WEPL,2);
  double alpha2 = 0.99-0.46*pow(Point->WET/Point->WEPL,2);
  double Length = TVector3(mf-mi).Mag();

  p0.SetMag(Length*alpha1);
  p1.SetMag(Length*alpha2);

  TVector3 A,B,C,D;
  A     =    2*mi - 2*mf + p0+p1;
  B     =   -3*mi + 3*mf - 2*p0-p1;
  C     =    p0;
  D     =    mi;
  for(int i=0;i<NStep;i++){
    double t=double(i)/NStep;
    m = D+t*(C+t*(B+t*A));
    ddb->Fill(m.x(), m.y(), m.z(), Point->WET);
  }
}

////////////////////////////////////////////
// Compute Spline with a prior (Hull) 
////////////////////////////////////////////
void ComputeSplinePrior(Particle *Point, TH3D* Prior, TProfile3D* ddb){
  TVector3 p0(Point->px0,   Point->py0,   Point->pz0);
  TVector3 p1(Point->px1,   Point->py1,   Point->pz1);
  TVector3 m0(Point->x0/10, Point->y0/10, Point->z0/10); // mm -> cm
  TVector3 m1(Point->x1/10, Point->y1/10, Point->z1/10); // mm -> cm

  int NStep = 512;
  TVector3 m, m_entry, m_exit, m_old; // position of the path at the entry and exit
  double t, t_entry = 0., t_exit = 1.; // fraction of the path at which point the spline enter the Hull
  TVector3 origin(0, 0, 0); // cm
  Point->angle = 0;
  m0 -= origin; m1 -= origin;  
  // Negative rotation because we change the coordinates and not the phantom
  m0.RotateZ(-1*(Point->angle)*M_PI/180.);
  m1.RotateZ(-1*(Point->angle)*M_PI/180.);
  p0.RotateZ(-1*(Point->angle)*M_PI/180.);
  p1.RotateZ(-1*(Point->angle)*M_PI/180.);
  m0 += origin; m1 += origin;
  
  m_entry = m0;
  m_exit  = m1;
  
  double RSP    = 0.0;  
  double recon_radius = 10.0;
  // Propagate from the entrance to the Hull
  for(int k =1; k<NStep-1; k++){
    t_entry = double(k)/NStep;
    m_entry = m0 + t_entry*p0;
    int global = Prior->FindBin(m_entry.x(), m_entry.y(),m_entry.z());
    double RSP = Prior->GetBinContent(global);
    if(RSP>0.4) break;
  }

  // Retro Propagate from the exit to the Hull
  for(int k =NStep-1; k>=1; k--){
    t_exit = double(k)/NStep;
    m_exit = m1 - t_exit*p1;
    int global = Prior->FindBin(m_exit.x(),m_exit.y(),m_exit.z());
    double RSP = Prior->GetBinContent(global);
    if(RSP>0.4) break;
  }

  if(t_entry > t_exit) { // No Hull we are in air
    m_exit = m1;
    m_entry = m0;
    t_entry = 0.;
    t_exit = 1.;
  }

  double alpha1 = 1.01+0.43*pow(Point->WET/Point->WEPL,2);
  double alpha2 = 0.99-0.46*pow(Point->WET/Point->WEPL,2);
  double TotLength  = TVector3(m1-m0).Mag();
  double HullLength = TVector3(m_exit-m_entry).Mag(); 

  double Einit  = 200;  
  double WET = 0;
  p0.SetMag(alpha1*HullLength);
  p1.SetMag(alpha2*HullLength);
  TVector3 A,B,C,D;
  A       =    2*m_entry - 2*m_exit + p0+p1;
  B       =   -3*m_entry + 3*m_exit - 2*p0-p1;
  C       =    p0;
  D       =    m_entry;
  m_old   =    m0;
  for(int i=0;i<NStep;i++){
    
    t = double(i)/NStep;
    // before the Hull -- straight line
    if(t < t_entry) m = m0 + t*(m_entry-m0);
    
    // after the Hull -- straight line
    else if (t > t_exit) m = m_exit+ t*(m1-m_exit);

    // in the Hull -- cubic spline
    else m = D+t*(C+t*(B+t*A));

    ddb->Fill(m.x(), m.y(), m.z(), Point->WET);

  }
}



////////////////////////////////////////////
// Compute Likelihood 
////////////////////////////////////////////
void ComputeLikelihood(Particle *Point, TMatrixD& tracks_X, TMatrixD& tracks_Y, TMatrixD& tracks_Z, TMatrixD& tracks_sig, TProfile3D* ddb){
  Int_t biny,binz;
  Float_t posY, posZ;
  for(int ix=0;ix<NbinsX;ix++){
    TF2 *f2 = new TF2("f2","xygaus", Xmin, Xmax, Ymin, Ymax);
    f2->SetParameters(1.,tracks_Y(ix,0),tracks_sig(ix,0),tracks_Z(ix,0),tracks_sig(ix,0));
    biny = ddb->GetYaxis()->FindBin(tracks_Y(ix,0));
    binz = ddb->GetZaxis()->FindBin(tracks_Z(ix,0));    
    for(int iy = biny-5; iy<biny+5; biny++){
      for(int iz = binz-5; iy<binz+5; binz++){
	posY = ddb->GetYaxis()->GetBinCenter(iy);
	posZ = ddb->GetZaxis()->GetBinCenter(iz);	
	ddb->Fill(tracks_X(ix,0), posY, posZ, f2->Eval(posY,posZ)*Point->WET);
      }
    }
  }
}

////////////////////////////////////////////
// Compute MLP Uncertainty -- Working but not checked
////////////////////////////////////////////
void CalculateMLPSigma(Particle *Point, TMatrixD& tracks_X, TMatrixD& tracks_L_inv, TMatrixD& tracks_pv_inv, TMatrixD& tracks_sig, TMatrixD& tracks_Y, TMatrixD& tracks_Z){

  TMatrixD Y0(2,1), Y1(2,1), Z0(2,1), Z1(2,1);
  TMatrixD R0(2,2), R1(2,2);   // Transvection matrix
  TMatrixD SI0(2,2), SI1(2,2); // Covariance matrix
  
  int N           = tracks_X.GetNoElements();
  double x0       = tracks_X(0,0);
  double x1       = tracks_X(N-1,0);
  double Epsilon  = 13.6;
  TVector3 p_old  = TVector3(x0, 0, 0);
  double sigma_t0, sigma_theta0, sigma_tt0;
  double sigma_t1, sigma_theta1, sigma_tt1;
  Y0(0,0) = Point->y0; Y0(1,0) = Point->py0;
  Y1(0,0) = Point->y1; Y1(1,0) = Point->py1;
  Z0(0,0) = Point->z0; Z0(1,0) = Point->pz0;
  Z1(0,0) = Point->z1; Z1(1,0) = Point->pz1;
  for (int idx = 1; idx <= N - 2; idx++){ // supposed to be from 1 to -1 but tracks_idx is defined from 0 to N-1 so we from 1 to N-2    
    double t     = tracks_X(idx,0);

    // SIGMA 0
    TMatrixD L_0 = tracks_L_inv.GetSub(0,idx+1,0,0); // matrix of the inverse
    TMatrixD pv0 = tracks_pv_inv.GetSub(0,idx+1,0,0); // matrix of 1/pv squared
    TMatrixD us0 = tracks_X.GetSub(0,idx+1,0,0);
    TMatrixD y0  = Mult( pv0, L_0);
    R0(0,0)  = 1        ; R0(0,1)  = t-Point->x0;
    R0(1,0)  = 0        ; R0(1,1)  = 1;

    double E0s0  = pow(1+0.038* TMath::Log(simps(us0, (t-x0)*L_0)),2);
    simps_mult(us0, y0, t, sigma_t0, sigma_tt0, sigma_theta0);
    SI0(0,0) = sigma_t0 ; SI0(0,1) = sigma_tt0   ;
    SI0(1,0) = sigma_tt0; SI0(1,1) = sigma_theta0;

    // SIGMA 1
    TMatrixD L_1 = tracks_L_inv.GetSub(idx,tracks_X.GetNoElements()-1,0,0);
    TMatrixD pv1 = tracks_pv_inv.GetSub(idx,tracks_X.GetNoElements()-1,0,0);
    TMatrixD us1 = tracks_X.GetSub(idx,tracks_X.GetNoElements()-1,0,0);
    TMatrixD y1  = Mult( pv1, L_1);
    R1(0,0)  = 1        ; R1(0,1)  = t-Point->x1;
    R1(1,0)  = 0        ; R1(1,1)  = 1;
    
    double E0s1  = pow(1+0.038* TMath::Log(simps(us1, (x1-t)*L_1)),2);
    simps_mult(us1,y1,t,sigma_t1,sigma_tt1,sigma_theta1);
    SI1(0,0) = sigma_t1 ; SI1(0,1) = sigma_tt1   ;
    SI1(1,0) = sigma_tt1; SI1(1,1) = sigma_theta1;

    //Mean and error MLP
    SI0 = inverse(E0s0*SI0);
    SI1 = inverse(E0s1*SI1);

    TMatrixD y_part1      = SI0*(R0*Y0) + SI1*(R1*Y1);
    TMatrixD z_part1      = SI0*(R0*Z0) + SI1*(R1*Z1);
    TMatrixD part2        = inverse(SI0+SI1); // Error
    tracks_Y(idx,0)       = (part2*y_part1)(0,0);
    tracks_Z(idx,0)       = (part2*z_part1)(0,0);
    tracks_sig(idx,2)     = Epsilon*TMath::Sqrt(part2(0,0)); // Position Variance
    tracks_sig(idx,1)     = Epsilon*TMath::Sqrt(part2(0,1)); // Covariance 
    tracks_sig(idx,0)     = Epsilon*TMath::Sqrt(part2(1,1)); // Direction Variance
    // Uncertainty with a reconstruction with both detector
    
    // Uncertainty with a reconstruction with uncertainty from the entrance
    /*tracks_sig(idx,0) = TMath::Sqrt(Epsilon*E0s0*sigma_theta0/2);
    tracks_sig(idx,1) = TMath::Sqrt(Epsilon*E0s0*sigma_tt0/2);
    tracks_sig(idx,2) = TMath::Sqrt(Epsilon*E0s0*sigma_t0/2);*/
    
    // Uncertainty with a reconstruction at the exit -- rear
    /*tracks_sig(idx,0) = TMath::Sqrt(Epsilon*E0s1*sigma_theta1*2);
    tracks_sig(idx,1) = TMath::Sqrt(Epsilon*E0s1*sigma_tt1*2);
    tracks_sig(idx,2) = TMath::Sqrt(Epsilon*E0s1*sigma_t1*2);*/
  }                                                                                                                                                                           

  // Fill the missing value caused by the asymptotic behaviour of the log law

  tracks_Y(0,0)     = tracks_Y(1,0);
  tracks_Z(0,0)     = tracks_Z(1,0);  
  tracks_sig(0,0)   = tracks_sig(1,0);
  tracks_sig(0,1)   = tracks_sig(1,1);
  tracks_sig(0,2)   = tracks_sig(1,2);

  tracks_Y(N-1,0)   = tracks_Y(N-2,0);
  tracks_Z(N-1,0)   = tracks_Z(N-2,0);    
  tracks_sig(N-1,0) = tracks_sig(N-2,0);
  tracks_sig(N-1,1) = tracks_sig(N-2,1);
  tracks_sig(N-1,2) = tracks_sig(N-2,2);
}                                                                                                                                                                                                 

////////////////////////////////////////////
// Gather EnergyRadLength  -- NOT WORKING
//////////////////////////////////////////// 
void GatherEnergyRadLength(TMatrixD &tracks_X, TMatrixD& tracks_L_inv, TMatrixD& tracks_E, TMatrixD& tracks_pv_inv , double Einit, TH3D* Prior){

  double dx           = Prior->GetXaxis()->GetBinWidth(0);
  tracks_E(0,0)       = Einit;
  tracks_pv_inv(0,0)  = pow( 1./pv(Einit),2);
  tracks_L_inv(0,0)   = 1./(3.5E5); // 3.5E5, 261
  double Eout = Einit;
  for(int idx = 1; idx <=  NbinsX ; idx++){

    double RSP = 1.0;//Prior->GetBinContent(idx,idy,idz);

    //Gather the rad-length information - only relevant between air/water (CACF 2017)
    if(RSP>0.5) tracks_L_inv(idx,0)      = 1./361;   // in water in mm
    else if(RSP<0.5) tracks_L_inv(idx,0) = 1./3.5E5;  // in air in mm
    
    //Gather the energy-momentum information
    int idE          = lower_bound(Energy.begin(), Energy.end(), Eout) -Energy.begin();
    Eout            -= (dx/10)*RSP*dEdXBins[idE];
    if(Eout<=0) Eout = 0.;
    tracks_E(idx,0)  = Eout; // in MeV
    tracks_pv_inv(idx,0) = pow( 1./pv(Eout), 2); // 1/pv^2
  }
}
////////////////////////////////////////////
// linspace
////////////////////////////////////////////
TMatrixD linspace(float x0, float x1, int NStep){
    int nrows = NStep;
    int ncols = 1;
    TMatrixD result(nrows,ncols);
    for(int i=0;i<nrows ;i++){
        result(i,0) = x0 + double(i)*double(x1 - x0)/(NStep-1);
    }
    return result;
}
////////////////////////////////////////////                                                                                                                                                                      
// Inverse
////////////////////////////////////////////
TMatrixD inverse(TMatrix A){

  int nrows = A.GetNrows();
  int ncols = A.GetNcols();
  TMatrixD result(nrows,ncols);
  result(0,0) = A(1,1); result(1,1) = A(0,0);
  result(1,0) = -1*A(0,1); result(0,1) = -1*A(1,0);
  double det  =  A(0,0)*A(1,1) - A(1,0)*A(0,1);
  result *= 1./det;
  return result;
}

////////////////////////////////////////////
// Inverse
////////////////////////////////////////////
TMatrixD reverse(TMatrix A){
  int nrows = A.GetNrows();
  int ncols = A.GetNcols();
  TMatrixD result(nrows,ncols);
  for(int i = 0; i<nrows; i++) result(i,0) = A(nrows-1-i,0);
  return result;
}

////////////////////////////////////////////
// element-wise division
////////////////////////////////////////////
TMatrixD Div(TMatrixD A, TMatrixD B){
    int nrows = A.GetNrows();
    int ncols = A.GetNcols();
    TMatrixD result(nrows,ncols);
    for(int i=0;i<nrows ;i++){
        for(int j=0;j<ncols;j++){
            result(i,j) = A(i,j)/B(i,j);
        }
    }
    return result;
}

////////////////////////////////////////////
// element-wise power
////////////////////////////////////////////
TMatrixD Pow(TMatrixD A, float p){
    int nrows = A.GetNrows();
    int ncols = A.GetNcols();
    TMatrixD result(nrows,ncols);
    for(int i=0;i<nrows ;i++){
        for(int j=0;j<ncols;j++){
            result(i,j) = pow(A(i,j),p);
        }
    }
    return result;
}

////////////////////////////////////////////
// element-wise multiplication
////////////////////////////////////////////
TMatrixD Mult(TMatrixD A, TMatrixD B){
    int nrows  = A.GetNrows();
    int ncols = A.GetNcols();
    TMatrixD result(nrows,ncols);
    for(int i=0;i<nrows ;i++){
        for(int j=0;j<ncols;j++){
            result(i,j) = A(i,j)*B(i,j);
        }
    }
    return result;
}

////////////////////////////////////////////
// Simpson rule of integration for the three sigma
////////////////////////////////////////////
void simps_mult(TMatrixD x, TMatrixD y, double t, double &sigma_t0, double &sigma_tt0, double &sigma_theta0){
  double h = 0;
  sigma_tt0 = 0; sigma_t0 = 0; sigma_theta0 = 0;
  for (int i=0; i < x.GetNoElements(); i++){

    if ( i == 0 || i == x.GetNoElements()-1 ){ // for the first and last elements
      h = (x(1,0)-x(0,0));
      sigma_theta0 += y(i,0)*h/3.;
      sigma_tt0    += (t-x(i,0))*y(i,0)*h/3.;
      sigma_t0     += (t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
    }
    else
      {
        h = (x(i+1,0)-x(i,0));
        if (i%2==0)
          {
            sigma_theta0 += 2*y(i,0)*h/3.;
            sigma_tt0    += 2*(t-x(i,0))*y(i,0)*h/3.;
            sigma_t0     += 2*(t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
          }
        else
          {
            h = (x(i+1,0)-x(i,0));
            sigma_theta0 += 4*y(i,0)*h/3.;
            sigma_tt0    += 4*(t-x(i,0))*y(i,0)*h/3.;
            sigma_t0     += 4*(t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
          }
      }
  }
}

////////////////////////////////////////////
// Simpson rule of integration
////////////////////////////////////////////
double simps(TMatrixD x, TMatrixD y){
  double sum=0.0;
  double h =0.0;
  for (int i=0; i < x.GetNoElements(); i+=1){
      if ( i == 0 || i == x.GetNoElements()-1 ) // for the first and last elements
        {
          h = x(1,0)- x(0,0);
          sum += y(i,0)*h/3.;
        }
      else
        {
          h = x(i+1,0)- x(i,0);
          if (i%2==0) sum += 2*y(i,0)*h/3.;
          else        sum += 4*y(i,0)*h/3.; // the rest of data
        }
  }
  return sum;
}

////////////////////////////////////////////
// Velocity divided by light velocity
////////////////////////////////////////////
double E2beta(double E){
  double mc2   = 938.27; // relativistic mass for protons
  double tau   = E/mc2;
  return (tau+2)*tau/pow(tau+1,2);
}

////////////////////////////////////////////
// momentum-velocity
////////////////////////////////////////////
double pv(double E){
  double mc2   = 938.27; // relativistic mass for protons
  double tau   = E/mc2;
  return E*(tau+2)/(tau+1);
}
////////////////////////////////////////////
// Single-variate normal distribution
////////////////////////////////////////////
double Gauss(double x, double x0, double sigma){
  double num = TMath::Exp(-pow(x-x0,2)/(2*pow(sigma,2)));
  double den = TMath::Sqrt(2*TMath::Pi()*pow(sigma,2));
  return  num/den;
}
////////////////////////////////////////////
// Bi-variate normal distribution
////////////////////////////////////////////
double Gauss2D(TMatrixD Y0, TMatrixD Y1, TMatrixD Sigma){
  TMatrixD Sigma_I = Sigma.Invert();
  TMatrixD Diff  = Y1 - Y0;
  TMatrixD Diff_t(TMatrixD::kTransposed,Diff);
  TMatrixD Part1 = Sigma_I*Diff;
  TMatrixD num   = (Diff_t*Part1);
  double den     = (2*TMath::Pi())*TMath::Sqrt(Sigma.Determinant());
  return TMath::Exp(-0.5*num(0,0))/den;
}


