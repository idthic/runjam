#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include "HydroSpectrum.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

#ifdef ALPHA
# define abs  fabs
inline double max(double a, double b) { return a>=b ? a:b;}
inline double min(double a, double b) { return a<=b ? a:b;}
# include <strstream>
#else
# include <sstream>
#endif

#include <iostream>

using namespace std;
const double HydroSpectrum::mPion = 139.0;
//const double HydroSpectrum::mPion  = 1020.0;//For phi meson
//const double HydroSpectrum::mPion  = 3096.9;//For J/psi meson
const double HydroSpectrum::mKaon    = 493.6;
const double HydroSpectrum::mProton  = 939.0;

//...scale transformation 1 [fm^-1]=197.32 [MeV]
const double sctr = 197.32;

//HydroSpectrum::HydroSpectrum(int kint)
HydroSpectrum::HydroSpectrum(int kint, int eos_pce)
{
  kineticTemp = kint;

  pi   = 3.1415926;
  fin  = -1;


  hbarc = 0.19732;

  //...degree of freedom  [pi+ or pi- or pi0]
  degpi = 1;
  //      degpi  = 3;//For  phi meson
  degk  = 1;
  degp  = 2;

  //...pion mass [fm^-1]
	mpi   = 0.0;
	mk    = 0.0;
	mpro  = 0.0;

  if(eos_pce !=10){
    mpi   = mPion/sctr;
    mk    = mKaon/sctr;
    mpro  = mProton/sctr;
  }
  //
	mupi = 0.0;
	muk = 0.0;
	mup = 0.0;
  //                                       chemical potential at f.o.(pion)
  //                                        from pi to delta(1232)
  //                                       chemical potential at f.o.(kaon)
  //                                       from pi to delta(1232)
  //                                       chemical potential at f.o.
  //                                       from pi to delta(1232)

  if(eos_pce == 1){
    switch(kineticTemp) {

      //   (tf=80mev)
    case 1:
      mupi = 0.951925e+02/sctr;
      muk = 0.233670e+03/sctr;
      mup = 0.456089e+03/sctr;
      cout << "HydroSpectrum Tf=80MeV" << endl;
      break;
	
      //     (tf=100mev)
    case 2:
      mupi = 0.833141e+02/sctr;
      //	  mupi = 0.356941e+03/sctr;//For phi-meson
      //	  mupi = 0.0/sctr;//For J/psi
      muk = 0.180805e+03/sctr;
      mup = 0.348810e+03/sctr;
      cout << "HydroSpectrum Tf=100MeV" << endl;
      break;
	
	
      //     (tf=120mev)
    case 3:
      mupi = 0.646814e+02/sctr;
      muk = 0.128598e+03/sctr;
      mup = 0.245865e+03/sctr;
      cout << "HydroSpectrum Tf=120MeV" << endl;
      break;
	
      //     (tf=140mev)
    case 4:
      mupi = 0.406648e+02/sctr;
      muk = 0.633578e+02/sctr;
      mup = 0.145518e+03/sctr;
      cout << "HydroSpectrum Tf=140MeV" << endl;
      break;
	
      //     (tf=160mev)
    case 5:
      mupi = 0.137867e+02/sctr;
      muk = 0.249125e+02/sctr;
      mup = 0.476744e+02/sctr;
      cout << "HydroSpectrum Tf=160MeV" << endl;
      break;
	
    default:
      cerr << "HydroSpectrum::  Not avaiable sorry"<<endl;
      cerr << " kinetic temperature = " << kineticTemp << endl;
      exit(1);
    }
  }



  //      co = 1;
  //     ^^^^^^^---> V1
  co = 2;
  //     ^^^^^^^---> V2
  //      co = 3;
  //     ^^^^^^^---> V3
  //      co = 4;
  //     ^^^^^^^---> V3


  // KM, 2013/04/20, initializes the flag to rotate the data in freezeout.dat
  {
    const char* env=std::getenv("HydroSpectrum_RotateFreezeoutData");
    this->fRotateFreezeoutData=env&&std::atoi(env);
    if(this->fRotateFreezeoutData)
      std::cout<<"HydroSpectrum: RotateFreezeoutData mode enabled!"<<std::endl;
  }
}

//**********************************************************************
//...sgn=1: fermion  -1:boson
/*
double HydroSpectrum::thermaldist(double mt,double yy,double mu,int iw,int sgn)
{
	double pds1,pds4,pds6,pds7;
	double coshhhp = cosh(yy+hh);
	double coshhhm = cosh(yy-hh);
	double sinhhhp = sinh(hh+yy);
	double sinhhhm = sinh(hh-yy);
	double coshyvp = cosh(yy+yv);
	double coshyvm = cosh(yy-yv);
	double invcoshyv  = 1.0/cosh(yv);
	double beta = 1.0/tf;
	
    if((iw == 1) || (iw == 5)) {
        if(bulk == 1){
	pds1 = mt*dss*coshhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds4 = mt*dss*coshhhm + pti*cosp*dsx - pti*sinp*dsy;
	pds6 = mt*dss*coshhhp - pti*cosp*dsx + pti*sinp*dsy;
	pds7 = mt*dss*coshhhp - pti*cosp*dsx - pti*sinp*dsy;
	}else{
	pds1 = mt*dss*sinhhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds4 = mt*dss*sinhhhm + pti*cosp*dsx - pti*sinp*dsy;
	pds6 = mt*dss*sinhhhp - pti*cosp*dsx + pti*sinp*dsy;
	pds7 = mt*dss*sinhhhp - pti*cosp*dsx - pti*sinp*dsy;
	}
	double pu1 = gam*(mt*coshyvm*invcoshyv - pti*cosp*vx - pti*sinp*vy);
	double pu4 = gam*(mt*coshyvm*invcoshyv - pti*cosp*vx + pti*sinp*vy);
	double pu6 = gam*(mt*coshyvp*invcoshyv + pti*cosp*vx - pti*sinp*vy);
	double pu7 = gam*(mt*coshyvp*invcoshyv + pti*cosp*vx + pti*sinp*vy);

#ifdef ALPHA
	double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
	double bose4 = exp(max(-50.0,min(50.0,(pu4-mu)/tf)));
	double bose6 = exp(max(-50.0,min(50.0,(pu6-mu)/tf)));
	double bose7 = exp(max(-50.0,min(50.0,(pu7-mu)/tf)));
	bose1 = pds1/(bose1+sgn);
	bose4 = pds4/(bose4+sgn);
	bose6 = pds6/(bose6+sgn);
	bose7 = pds7/(bose7+sgn);
#else
	double bose1=0.0;
	double bose4=0.0;
	double bose6=0.0;
	double bose7=0.0;
	double bse1=beta*(pu1-mu);
	double bse4=beta*(pu4-mu);
	double bse6=beta*(pu6-mu);
	double bse7=beta*(pu7-mu);
	if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
	if(bse4<30.0) bose4 = pds4/(exp(bse4)+sgn);
	if(bse6<30.0) bose6 = pds6/(exp(bse6)+sgn);
	if(bse7<30.0) bose7 = pds7/(exp(bse7)+sgn);
#endif
return bose1+bose4+bose6+bose7;


    } else if((iw == 3) || (iw == 7)) {
      if(bulk == 1){
	pds1 = mt*dss*coshhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds6 = mt*dss*coshhhp - pti*cosp*dsx + pti*sinp*dsy;
	}else	{
	pds1 = mt*dss*sinhhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds6 = mt*dss*sinhhhp - pti*cosp*dsx + pti*sinp*dsy;
	}
      //	double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
      //	double pds6 = ee*ds0 - pti*cosp*dsx + pti*sinp*dsy - pz*dsz;
	double pu1 = gam*(mt*coshyvm*invcoshyv - pti*cosp*vx - pti*sinp*vy);
	double pu6 = gam*(mt*coshyvp*invcoshyv + pti*cosp*vx - pti*sinp*vy);
#ifdef ALPHA
	double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
	double bose6 = exp(max(-50.0,min(50.0,(pu6-mu)/tf)));
	bose1 = pds1/(bose1+sgn);
	bose6 = pds6/(bose6+sgn);
#else
	double bose1=0.0;
	double bose6=0.0;
	double bse1=beta*(pu1-mu);
	double bse6=beta*(pu6-mu);

	if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
	if(bse6<30.0) bose6 = pds6/(exp(bse6)+sgn);
#endif
	return bose1+bose6;

    } else if((iw == 2) || (iw == 6)) {
      if(bulk == 1){
	pds1 = mt*dss*coshhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds4 = mt*dss*coshhhm + pti*cosp*dsx - pti*sinp*dsy;
	}else	{
	pds1 = mt*dss*sinhhhm + pti*cosp*dsx + pti*sinp*dsy;
	pds4 = mt*dss*sinhhhm + pti*cosp*dsx - pti*sinp*dsy;
	}
      //	double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
      //	double pds4 = ee*ds0 + pti*cosp*dsx - pti*sinp*dsy + pz*dsz;
	double pu1 = gam*(mt*coshyvm*invcoshyv - pti*cosp*vx - pti*sinp*vy);
	double pu4 = gam*(mt*coshyvp*invcoshyv - pti*cosp*vx + pti*sinp*vy);
#ifdef ALPHA
	double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
	double bose4 = exp(max(-50.0,min(50.0,(pu4-mu)/tf)));
	bose1 = pds1/(bose1+sgn);
	bose4 = pds4/(bose4+sgn);
#else
	double bose1=0.0;
	double bose4=0.0;
	double bse1=beta*(pu1-mu);
	double bse4=beta*(pu4-mu);
	
	if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
	if(bse4<30.0) bose4 = pds4/(exp(bse4)+sgn);
#endif
	return bose1+bose4;


    } else if((iw == 4) || (iw == 8)) {
        if(bulk == 1){
	pds1 = mt*dss*coshhhm + pti*cosp*dsx + pti*sinp*dsy;
	}else	{
	pds1 = mt*dss*sinhhhm + pti*cosp*dsx + pti*sinp*dsy;
	}
      //	double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
	double pu1 = gam*(mt*coshyvm*invcoshyv - pti*cosp*vx - pti*sinp*vy);
#ifdef ALPHA
	double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
	bose1 = pds1/(bose1+sgn);
#else
	double bose1=0.0;
	double bse1=beta*(pu1-mu);
	if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
#endif
	return bose1;

    } else {
	cerr << "HydroSpec::thermaldist funny iw " << iw << endl;
	exit(1);
    }

  }
*/


double HydroSpectrum::thermaldist(double ee,double pz,double mu,int iw,int sgn){
  if((iw == 1) || (iw == 5)) {
    double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
    double pds4 = ee*ds0 + pti*cosp*dsx - pti*sinp*dsy + pz*dsz;
    double pds6 = ee*ds0 - pti*cosp*dsx + pti*sinp*dsy - pz*dsz;
    double pds7 = ee*ds0 - pti*cosp*dsx - pti*sinp*dsy - pz*dsz;
    double pu1 = gam*(ee - pti*cosp*vx - pti*sinp*vy - pz*vz);
    double pu4 = gam*(ee - pti*cosp*vx + pti*sinp*vy - pz*vz);
    double pu6 = gam*(ee + pti*cosp*vx - pti*sinp*vy + pz*vz);
    double pu7 = gam*(ee + pti*cosp*vx + pti*sinp*vy + pz*vz);
#ifdef ALPHA
    double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
    double bose4 = exp(max(-50.0,min(50.0,(pu4-mu)/tf)));
    double bose6 = exp(max(-50.0,min(50.0,(pu6-mu)/tf)));
    double bose7 = exp(max(-50.0,min(50.0,(pu7-mu)/tf)));
    bose1 = pds1/(bose1+sgn);
    bose4 = pds4/(bose4+sgn);
    bose6 = pds6/(bose6+sgn);
    bose7 = pds7/(bose7+sgn);
#else
    double bose1=0.0;
    double bose4=0.0;
    double bose6=0.0;
    double bose7=0.0;
    double bse1=beta*(pu1-mu);
    double bse4=beta*(pu4-mu);
    double bse6=beta*(pu6-mu);
    double bse7=beta*(pu7-mu);
    if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
    if(bse4<30.0) bose4 = pds4/(exp(bse4)+sgn);
    if(bse6<30.0) bose6 = pds6/(exp(bse6)+sgn);
    if(bse7<30.0) bose7 = pds7/(exp(bse7)+sgn);
#endif
    return bose1+bose4+bose6+bose7;
    //	return (bose1+abs(bose1)+bose4+abs(bose4)
    //		+bose6+abs(bose6)+bose7+abs(bose7))/2.0;

  } else if((iw == 3) || (iw == 7)) {
    double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
    double pds6 = ee*ds0 - pti*cosp*dsx + pti*sinp*dsy - pz*dsz;
    double pu1 = gam*(ee - pti*cosp*vx - pti*sinp*vy - pz*vz);
    double pu6 = gam*(ee + pti*cosp*vx - pti*sinp*vy + pz*vz);
#ifdef ALPHA
    double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
    double bose6 = exp(max(-50.0,min(50.0,(pu6-mu)/tf)));
    bose1 = pds1/(bose1+sgn);
    bose6 = pds6/(bose6+sgn);
#else
    double bose1=0.0;
    double bose6=0.0;
    double bse1=beta*(pu1-mu);
    double bse6=beta*(pu6-mu);
    if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
    if(bse6<30.0) bose6 = pds6/(exp(bse6)+sgn);
#endif
    return bose1+bose6;
    //	return (bose1+abs(bose1)+bose6+abs(bose6))/2.0;

  } else if((iw == 2) || (iw == 6)) {
    double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
    double pds4 = ee*ds0 + pti*cosp*dsx - pti*sinp*dsy + pz*dsz;
    double pu1 = gam*(ee - pti*cosp*vx - pti*sinp*vy - pz*vz);
    double pu4 = gam*(ee - pti*cosp*vx + pti*sinp*vy - pz*vz);
#ifdef ALPHA
    double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
    double bose4 = exp(max(-50.0,min(50.0,(pu4-mu)/tf)));
    bose1 = pds1/(bose1+sgn);
    bose4 = pds4/(bose4+sgn);
#else
    double bose1=0.0;
    double bose4=0.0;
    double bse1=beta*(pu1-mu);
    double bse4=beta*(pu4-mu);
    if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
    if(bse4<30.0) bose4 = pds4/(exp(bse4)+sgn);
#endif
    return bose1+bose4;
    //	return (bose1+abs(bose1)+bose4+abs(bose4))/2.0;

  } else if((iw == 4) || (iw == 8)) {
    double pds1 = ee*ds0 + pti*cosp*dsx + pti*sinp*dsy + pz*dsz;
    double pu1 = gam*(ee - pti*cosp*vx - pti*sinp*vy - pz*vz);
#ifdef ALPHA
    double bose1 = exp(max(-50.0,min(50.0,(pu1-mu)/tf)));
    bose1 = pds1/(bose1+sgn);
#else
    double bose1=0.0;
    double bse1=beta*(pu1-mu);
    if(bse1<30.0) bose1 = pds1/(exp(bse1)+sgn);
#endif
    return bose1;
    //	return (bose1+abs(bose1))/2.0;
  } else {
    cerr << "HydroSpec::thermaldist funny iw " << iw << endl;
    exit(1);
  }
}

void HydroSpectrum::openDataFile(string fname)
{
  //fdata.open(unit=20,file='../src/25/freezeout.dat',status='OLD')
  fdata.open(fname.c_str(), ios::in);
  if(!fdata){
    cerr << "HydroSpectrum::openDataFile(fn)! unable to open file " << fname << endl;
    exit(1);
  }
}
void HydroSpectrum::openDataFile(string fname,string fname_ecc)
{

  if(fdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! file already opened "
         << fname << endl;
    exit(1);
  }
  if(eccdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! file already opened "
         << fname_ecc << endl;
    exit(1);
  }

  //fdata.open(unit=20,file='../src/25/freezeout.dat',status='OLD')
  fdata.open(fname.c_str(), ios::in);
  if (!fdata)  {
    cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! unable to open file " << fname << endl;
    exit(1);
  }

  //eccdata.open(fname_ecc.c_str(), ios::in);
  eccdata.open(fname_ecc.c_str(), ifstream::in);
  if (!eccdata)  {
    cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! unable to open file <" << fname_ecc
	       << ">" << endl;
    exit(1);
  }

}

void HydroSpectrum::openDataFile2(string fname1,string fname2)
{
  if(fdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile2! file already opened"
         << fname1 << endl;
    exit(1);
  }
  if(pdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile2! file already opened"
         << fname2 << endl;
    exit(1);
  }

  fdata.open(fname1.c_str(), ios::in);
  if (!fdata)  {
    cerr << "HydroSpectrum::openDataFile2! unable to open file " << fname1 << endl;
    exit(1);
  }
  pdata.open(fname2.c_str(), ios::in);
  if (!pdata)  {
    cerr << "HydroSpectrum::openDataFile2! unable to open file " << fname2 << endl;
    exit(1);
  }

}
void HydroSpectrum::openDataFile(string fname1,string fname2,string fname_ecc)
{

  if(fdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname1 << endl;
    exit(1);
  }
  if(pdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname2 << endl;
    exit(1);
  }
  if(eccdata.is_open()) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname_ecc << endl;
    exit(1);
  }

  fdata.open(fname1.c_str(), ios::in);
  if (!fdata) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file " << fname1 << endl;
    exit(1);
  }
  pdata.open(fname2.c_str(), ios::in);
  if (!pdata) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file " << fname2 << endl;
    exit(1);
  }

  eccdata.open(fname_ecc.c_str(), ifstream::in);
  if (!eccdata) {
    cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file <" << fname_ecc
	       << ">" << endl;
    exit(1);
  }

}

// input data of hydrodynamic flow
int HydroSpectrum::readData()
{
  if(!fdata){
    std::cerr<<"HydroSpectrum::readData: unexpected end of freezeout.dat file"<<std::endl;
    std::exit(1);
  }

  fdata >> bulk;
  if(bulk == fin) return 1;

  fdata >> dss; // $2 time component of the surface
  fdata >> dsx; // $3 x component of the surface
  fdata >> dsy; // $4 y component of the surface
  // fdata >> dsz;
  fdata >> hh;  // $5 eta of the cell?
  fdata >> tf;  // $6 temperature of the fluid
  fdata >> nbf; // $7 baryon number density of the fluid
  fdata >> vx;  // $8 velocity along the x axis
  fdata >> vy;  // $9 velocity along the y axis

  fdata >> yv;  // $10 rapidity
  // vz = tanh(yv);

  fdata >> iw;

  // 2013/04/20, KM, reverse data
  if(this->fRotateFreezeoutData){
    // hh=-hh;
    // if(bulk!=1)
    //   dss=-dss;

    if(bulk==1)
      hh=-hh;
    yv=-yv;
    if(iw!=8){
      std::cerr<<"HydroSpectrum_RotateFreezeoutData: not supported for the case that iw=="<<iw<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // if(tf*197.32 > 100.0) cout << tf*197.32 << endl;

  beta = 1.0/tf;

  vz = tanh(yv);

  if(bulk == 1){
    ds0 = dss*cosh(hh);
    dsz = -dss*sinh(hh);
  }else{
    ds0 = dss*sinh(hh);
    dsz = -dss*cosh(hh);
  }

  return 0;
}
int HydroSpectrum::readPData()
{

	pdata >> tau;
	pdata >> xx;
	pdata >> yy;
	pdata >> eta;

  // 2013/04/20, KM, reverse data
  if(this->fRotateFreezeoutData){
    eta=-eta;
  }

	return 0;
}
int HydroSpectrum::readEData()
{

  int nlines = 0;

  char buffer[1280];
  while (!eccdata.eof()) {
    eccdata.getline(buffer,sizeof(buffer),'\n');
    //	if(eccdata.eof() || strlen(buffer) == 0 ) break;
    if(nlines>=100) {
	    cerr << " too many data " << endl;
	    exit(1);
    }
    //sscanf(buffer,"%lf %lf %lf %lf %lf",&hh_ecc[nlines],&coe_ecc[nlines],
    //       &ecc_ecc[nlines],&area_ecc[nlines],&aveene_ecc[nlines]);
    sscanf(buffer,"%lf %lf %lf %lf %lf %lf",&hh_ecc[nlines],
    &coe_ecc[nlines],&ecc_ecc[nlines],&area_ecc[nlines],
    &aveene_ecc[nlines],&eccp_ecc[nlines]);

    nlines++;
  }

/* //
  string templine;
  while(eccdata) {
    getline(eccdata,templine);
    if(templine.substr(0,1) == "#") continue;
    if(eccdata.eof() || templine.size() == 0 ) break;

    if(nlines>=50) {
      cerr << " too many data " << endl;
      exit(1);
    }
    istringstream ist(templine);
    ist >> hh_ecc[nlines] >> coe_ecc[nlines] >> ecc_ecc[nlines];
    if(ist.fail()) {
      cerr << "funny " << endl;
      exit(1);
    }
    nlines++;
  }
// */

  return nlines;
}
/**********************************************************
C
C     SUBROUTINE GAUS12
C       SET WEIGHT VSLUE FOR INTEGRATION
C
C**********************************************************/
void HydroSpectrum::Gauss12(double xini,double xfin,double* xn,double* wn)
{
  //double xn[12],wn[12];

  double x[6]={0.98156063424671925e00, 0.90411725637047486e00,
               0.76990267419430469e00, 0.58731795428661745e00,
               0.36783149899818019e00, 0.12523340851146892e00};

  double w[6]={0.47175336386511827e-01, 0.10693932599531843e00,
               0.16007832854334623e00, 0.20316742672306592e00,
               0.23349253653835481e00, 0.24914704581340279e00};

  double a1=0.5*(xfin-xini);
  double a2=0.5*(xfin+xini);

  for(int i=0;i<6;i++) {
    xn[i  ]= x[i  ]*a1+a2;
    xn[i+6]=-x[5-i]*a1+a2;
    wn[i  ]= w[i  ]*a1;
    wn[i+6]= w[5-i]*a1;
  }

}

//******************************************************************
//
//        SUBROUTINE GAUSS38-LAGUERRE20
//              DATA IS GIVEN BY H.NAKAMURA
//              THIS SUBROUTINE IS WRITTEN BY T.HIRANO
//                                           '96 11 29
//******************************************************************
void HydroSpectrum::GauLag(double xini,double xmid,double* xn,double* wn)
{
  //double xn(58),wn(58),xini,xmid

  static double x[38]={- 9.980499305357e-01, - 9.897394542664e-01,
                       - 9.748463285902e-01, - 9.534663309335e-01,
                       - 9.257413320486e-01, - 8.918557390046e-01,
                       - 8.520350219324e-01, - 8.065441676053e-01,
                       - 7.556859037540e-01, - 6.997986803792e-01,
                       - 6.392544158297e-01, - 5.744560210478e-01,
                       - 5.058347179279e-01, - 4.338471694324e-01,
                       - 3.589724404794e-01, - 2.817088097902e-01,
                       - 2.025704538921e-01, - 1.220840253379e-01,
                       - 4.078514790458e-02,   4.078514790458e-02,
                       1.220840253379e-01,   2.025704538921e-01,
                       2.817088097902e-01,   3.589724404794e-01,
                       4.338471694324e-01,   5.058347179279e-01,
                       5.744560210478e-01,   6.392544158297e-01,
                       6.997986803792e-01,   7.556859037540e-01,
                       8.065441676053e-01,   8.520350219324e-01,
                       8.918557390046e-01,   9.257413320486e-01,
                       9.534663309335e-01,   9.748463285902e-01,
                       9.897394542664e-01,   9.980499305357e-01};

  static double y[20]={7.053988969198e-02,   3.721268180016e-01,
                       9.165821024832e-01,   1.707306531028e00,
                       2.749199255310e00,    4.048925313850e00,
                       5.615174970848e00,    7.459017453700e00,
                       9.594392870116e00,    1.203880254193e+01,
                       1.481429346483e+01,   1.794889545992e+01,
                       2.147878834852e+01,   2.545170266868e+01,
                       2.993255471277e+01,   3.501343423054e+01,
                       4.083305702714e+01,   4.761999407237e+01,
                       5.581079574152e+01,   6.652441652668e+01};

  static double w[38]={5.002880749632e-3,    1.161344471647e-2,
                       1.815657770961e-2,    2.457973973823e-2,
                       3.083950054518e-2,    3.689408159400e-2,
                       4.270315850467e-2,    4.822806186076e-2,
                       5.343201991033e-2,    5.828039914700e-2,
                       6.274093339213e-2,    6.678393797914e-2,
                       7.038250706690e-2,    7.351269258474e-2,
                       7.615366354845e-2,    7.828784465821e-2,
                       7.990103324353e-2,    8.098249377060e-2,
                       8.152502928039e-2,    8.152502928039e-2,
                       8.098249377060e-2,    7.990103324353e-2,
                       7.828784465821e-2,    7.615366354845e-2,
                       7.351269258474e-2,    7.038250706690e-2,
                       6.678393797914e-2,    6.274093339213e-2,
                       5.828039914700e-2,    5.343201991033e-2,
                       4.822806186076e-2,    4.270315850467e-2,
                       3.689408159400e-2,    3.083950054518e-2,
                       2.457973973823e-2,    1.815657770961e-2,
                       1.161344471647e-2,    5.002880749632e-3};

  static double u[20]={1.687468018511e-1,    2.912543620060e-1,
                       2.666861028666e-1,    1.660024532694e-1,
                       7.482606467538e-2,    2.496441732002e-2,
                       6.202550828442e-3,    1.144962383071e-3,
                       1.557417806303e-4,    1.540143980270e-5,
                       1.086486325120e-6,    5.330127202618e-8,
                       1.757985532862e-9,    3.725505729478e-11,
                       4.767526861010e-13,   3.372848343678e-15,
                       1.155014009180e-17,   1.539522135110e-20,
                       5.286443248824e-24,   1.656456602226e-28};


  for(int i=0;i<38;i++) {
    xn[i] = (xmid-xini)*x[i]/2. +(xini+xmid)/2.;
    wn[i] = (xmid-xini)*w[i]/2.;
  }

  //       do 998 j = 1,20
  for(int j=0;j<20;j++) {
    xn[j+38] = y[j] + xmid;
    wn[j+38] = exp(y[j])*u[j];
  }

}

//**********************************************************************
void HydroSpectrum::Gauss38(double xini,double xfin,double* xn,double* wn)
{
  //double xn[38],wn[38];
  double  x[38], w[38];

  x[37]=9.980499305357e-1;
  x[36]=9.897394542664e-1;
  x[35]=9.748463285902e-1;
  x[34]=9.534663309335e-1;
  x[33]=9.257413320486e-1;
  x[32]=8.918557390046e-1;
  x[31]=8.520350219324e-1;
  x[30]=8.065441676053e-1;
  x[29]=7.556859037540e-1;
  x[28]=6.997986803792e-1;
  x[27]=6.392544158297e-1;
  x[26]=5.744560210478e-1;
  x[25]=5.058347179279e-1;
  x[24]=4.338471694324e-1;
  x[23]=3.589724404794e-1;
  x[22]=2.817088097902e-1;
  x[21]=2.025704538921e-1;
  x[20]=1.220840253379e-1;
  x[19]=4.078514790458e-2;

  //    .....   WEIGHT       ...........
  w[37]=5.002880749632e-3;
  w[36]=1.161344471647e-2;
  w[35]=1.815657770961e-2;
  w[34]=2.457973973823e-2;
  w[33]=3.083950054518e-2;
  w[32]=3.689408159400e-2;
  w[31]=4.270315850467e-2;
  w[30]=4.822806186076e-2;
  w[29]=5.343201991033e-2;
  w[28]=5.828039914700e-2;
  w[27]=6.274093339213e-2;
  w[26]=6.678393797914e-2;
  w[25]=7.038250706690e-2;
  w[24]=7.351269258474e-2;
  w[23]=7.615366354845e-2;
  w[22]=7.828784465821e-2;
  w[21]=7.990103324353e-2;
  w[20]=8.098249377060e-2;
  w[19]=8.152502928039e-2;

  for(int i=0;i<19;i++) {
    x[i] = -x[37-i];
    w[i] =  w[37-i];
  }
  for(int i=0;i<38;i++) {
    xn[i] =(xfin-xini)*x[i]/2.0+(xini+xfin)/2.0;
    wn[i] =(xfin-xini)*w[i]/2.0;
  }

}
