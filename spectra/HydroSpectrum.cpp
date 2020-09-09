#include "HydroSpectrum.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <util/Constants.hpp>

#ifdef ALPHA
# include <strstream>
#else
# include <sstream>
#endif

#include <iostream>

#ifdef ALPHA
inline double clamp(double value, double lower, double upper) {
  return value < lower ? lower : value > uppper ? upper : value;
}

// 1 / (exp(x) + sgn)
inline double thermal_distribution(double x, double sgn) {
  double const value = 1.0 / (std::exp(clamp(x, -50.0, 50.0)) + sgn);
  return value < 0.0 ? 0.0 : value;
}
#else
// 1 / (exp(x) + sgn)
inline double thermal_distribution(double x, double sgn) {
  double const value = x < 30.0 ? 1.0 / (std::exp(x) + sgn) : 0.0;
  return value < 0.0 ? 0.0 : value;
}
#endif


namespace idt {
namespace runjam {

using namespace std;

const double HydroSpectrum::mPion = 139.0;
//const double HydroSpectrum::mPion  = 1020.0;//For phi meson
//const double HydroSpectrum::mPion  = 3096.9;//For J/psi meson
const double HydroSpectrum::mKaon    = 493.6;
const double HydroSpectrum::mProton  = 939.0;

//HydroSpectrum::HydroSpectrum(int kint)
HydroSpectrum::HydroSpectrum(int kint, int eos_pce) {
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
  if (eos_pce != 10) {
    mpi   = mPion / hbarc_MeVfm;
    mk    = mKaon / hbarc_MeVfm;
    mpro  = mProton / hbarc_MeVfm;
  }

	mupi = 0.0; //!< chemical potential at f.o. (pion)   from pi to delta(1232)
	muk = 0.0;  //!< chemical potential at f.o. (kaon)   from pi to delta(1232)
	mup = 0.0;  //!< chemical potential at f.o. (proton) from pi to delta(1232)

  if (eos_pce == 1) {
    switch (kineticTemp) {
      //   (tf=80mev)
    case 1:
      mupi = 0.951925e+02 / hbarc_MeVfm;
      muk = 0.233670e+03 / hbarc_MeVfm;
      mup = 0.456089e+03 / hbarc_MeVfm;
      cout << "HydroSpectrum Tf=80MeV" << std::endl;
      break;

      //     (tf=100mev)
    case 2:
      mupi = 0.833141e+02 / hbarc_MeVfm;
      //	  mupi = 0.356941e+03 / hbarc_MeVfm;//For phi-meson
      //	  mupi = 0.0 / hbarc_MeVfm;//For J/psi
      muk = 0.180805e+03 / hbarc_MeVfm;
      mup = 0.348810e+03 / hbarc_MeVfm;
      cout << "HydroSpectrum Tf=100MeV" << std::endl;
      break;

      //     (tf=120mev)
    case 3:
      mupi = 0.646814e+02 / hbarc_MeVfm;
      muk = 0.128598e+03 / hbarc_MeVfm;
      mup = 0.245865e+03 / hbarc_MeVfm;
      cout << "HydroSpectrum Tf=120MeV" << std::endl;
      break;

      //     (tf=140mev)
    case 4:
      mupi = 0.406648e+02 / hbarc_MeVfm;
      muk = 0.633578e+02 / hbarc_MeVfm;
      mup = 0.145518e+03 / hbarc_MeVfm;
      cout << "HydroSpectrum Tf=140MeV" << std::endl;
      break;

      //     (tf=160mev)
    case 5:
      mupi = 0.137867e+02 / hbarc_MeVfm;
      muk = 0.249125e+02 / hbarc_MeVfm;
      mup = 0.476744e+02 / hbarc_MeVfm;
      cout << "HydroSpectrum Tf=160MeV" << std::endl;
      break;

    default:
      std::cerr << "HydroSpectrum::  Not avaiable sorry" << std::endl;
      std::cerr << " kinetic temperature = " << kineticTemp << std::endl;
      std::exit(1);
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
    const char* env = std::getenv("HydroSpectrum_RotateFreezeoutData");
    this->fRotateFreezeoutData = env && std::atoi(env);
    if (this->fRotateFreezeoutData)
      std::cout << "HydroSpectrum: RotateFreezeoutData mode enabled!" << std::endl;
  }
}

//**********************************************************************
double HydroSpectrum::thermaldist(double ee, double pz, double mu, int iw, int sgn) {
  if (iw == 1 || iw == 5) {
    double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
    double pds4 = ee * ds0 + pti * cosp * dsx - pti * sinp * dsy + pz * dsz;
    double pds6 = ee * ds0 - pti * cosp * dsx + pti * sinp * dsy - pz * dsz;
    double pds7 = ee * ds0 - pti * cosp * dsx - pti * sinp * dsy - pz * dsz;
    double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
    double pu4 = gam * (ee - pti * cosp * vx + pti * sinp * vy - pz * vz);
    double pu6 = gam * (ee + pti * cosp * vx - pti * sinp * vy + pz * vz);
    double pu7 = gam * (ee + pti * cosp * vx + pti * sinp * vy + pz * vz);
    double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
    double bose4 = pds4 * thermal_distribution(beta * (pu4 - mu), sgn);
    double bose6 = pds6 * thermal_distribution(beta * (pu6 - mu), sgn);
    double bose7 = pds7 * thermal_distribution(beta * (pu7 - mu), sgn);
    return bose1 + bose4 + bose6 + bose7;
  } else if (iw == 3 || iw == 7) {
    double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
    double pds6 = ee * ds0 - pti * cosp * dsx + pti * sinp * dsy - pz * dsz;
    double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
    double pu6 = gam * (ee + pti * cosp * vx - pti * sinp * vy + pz * vz);
    double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
    double bose6 = pds6 * thermal_distribution(beta * (pu6 - mu), sgn);
    return bose1 + bose6;
  } else if (iw == 2 || iw == 6) {
    double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
    double pds4 = ee * ds0 + pti * cosp * dsx - pti * sinp * dsy + pz * dsz;
    double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
    double pu4 = gam * (ee - pti * cosp * vx + pti * sinp * vy - pz * vz);
    double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
    double bose4 = pds4 * thermal_distribution(beta * (pu4 - mu), sgn);
    return bose1 + bose4;
  } else if (iw == 4 || iw == 8) {
    double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
    double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
    double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
    return bose1;
  } else {
    std::cerr << "HydroSpec::thermaldist funny iw " << iw << std::endl;
    std::exit(1);
  }
}

void HydroSpectrum::openFDataFile(string fname) {
  if (fdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile! file already opened '" << fname << "'" << std::endl;
    std::exit(1);
  }
  fdata_fname = fname;
  fdata.open(fname.c_str());
  if (!fdata) {
    std::cerr << "HydroSpectrum::openDataFile! unable to open file '" << fname << "'" << std::endl;
    std::exit(1);
  }
}
void HydroSpectrum::openPDataFile(string fname2) {
  if (pdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile2! file already opened '" << fname2 << "'" << std::endl;
    std::exit(1);
  }
  pdata_fname = fname2;
  pdata.open(fname2.c_str(), ios::in);
  if (!pdata)  {
    std::cerr << "HydroSpectrum::openDataFile2! unable to open file '" << fname2 << "'" << std::endl;
    std::exit(1);
  }
}

void HydroSpectrum::openEDataFile(string fname_ecc) {
  if (eccdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! file already opened '" << fname_ecc << "'" << std::endl;
    std::exit(1);
  }
  eccdata_fname = fname_ecc;
  eccdata.open(fname_ecc.c_str());
  if (!eccdata)  {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! unable to open file '" << fname_ecc << "'" << std::endl;
    std::exit(1);
  }
}

// input data of hydrodynamic flow
int HydroSpectrum::readFData() {
  if (!fdata) {
    std::cerr << "HydroSpectrum::readFData: unexpected end of freezeout.dat file" << std::endl;
    std::exit(1);
  }

  fdata >> bulk;
  if (bulk == fin) return 1;

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
  if (this->fRotateFreezeoutData) {
    // hh = -hh;
    // if (bulk != 1)
    //   dss = -dss;

    if (bulk == 1)
      hh = -hh;
    yv = -yv;
    if (iw != 8) {
      std::cerr << "HydroSpectrum_RotateFreezeoutData: not supported for the case that iw==" << iw << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // if (tf * 197.32 > 100.0) cout << tf * 197.32 << std::endl;
  beta = 1.0 / tf;

  vz = tanh(yv);

  if (bulk == 1) {
    ds0 = dss * cosh(hh);
    dsz = -dss * sinh(hh);
  } else {
    ds0 = dss * sinh(hh);
    dsz = -dss * cosh(hh);
  }

  return 0;
}
int HydroSpectrum::readPData() {
	pdata >> tau;
	pdata >> xx;
	pdata >> yy;
	pdata >> eta;

  // 2013/04/20, KM, reverse data
  if (this->fRotateFreezeoutData) {
    eta=-eta;
  }

	return 0;
}
int HydroSpectrum::readEData() {
  char buffer[1280];
  std::string line;
  std::size_t index = 0, iline = 0;
  while (std::getline(eccdata, line)) {
    iline++;
    if (index >= 100) {
	    std::cerr << eccdata_fname << ":" << iline << ":"
                << " too many data" << std::endl;
	    std::exit(1);
    }

    int r = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf",
      &hh_ecc[index],
      &coe_ecc[index],
      &ecc_ecc[index],
      &area_ecc[index],
      &aveene_ecc[index],
      &eccp_ecc[index]);
    if (r != 6) {
	    std::cerr << eccdata_fname << ":" << iline << ":"
                << " invalid format" << std::endl;
      std::exit(1);
    }
    index++;
  }

  return index;
}

}
}
