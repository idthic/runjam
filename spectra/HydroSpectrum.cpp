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
namespace hydro2jam {

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
/// @fn double HydroSpectrum::thermaldist(double ee, double pz, double mu, int iw, int sgn);
/// @param[in] ee
/// @param[in] pz
/// @param[in] mu
/// @param[in] iw
/// @param[in] sgn 1: fermion  -1:boson
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

void HydroSpectrum::openDataFile(string fname) {
  fdata_fname = fname;
  fdata.open(fname.c_str());
  if (!fdata) {
    std::cerr << "HydroSpectrum::openDataFile(fn)! unable to open file " << fname << std::endl;
    std::exit(1);
  }
}
void HydroSpectrum::openDataFile(string fname, string fname_ecc) {
  if (fdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! file already opened "
         << fname << std::endl;
    std::exit(1);
  }
  if (eccdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! file already opened "
         << fname_ecc << std::endl;
    std::exit(1);
  }

  fdata_fname = fname;
  fdata.open(fname.c_str());
  if (!fdata)  {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! unable to open file " << fname << std::endl;
    std::exit(1);
  }

  eccdata_fname = fname_ecc;
  eccdata.open(fname_ecc.c_str());
  if (!eccdata)  {
    std::cerr << "HydroSpectrum::openDataFile(fn,fnEcc)! unable to open file <" << fname_ecc
	       << ">" << std::endl;
    std::exit(1);
  }
}

void HydroSpectrum::openDataFile2(string fname1, string fname2) {
  if (fdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile2! file already opened" << fname1 << std::endl;
    std::exit(1);
  }
  if (pdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile2! file already opened" << fname2 << std::endl;
    std::exit(1);
  }

  fdata_fname = fname1;
  fdata.open(fname1.c_str(), ios::in);
  if (!fdata)  {
    std::cerr << "HydroSpectrum::openDataFile2! unable to open file " << fname1 << std::endl;
    std::exit(1);
  }
  pdata_fname = fname2;
  pdata.open(fname2.c_str(), ios::in);
  if (!pdata)  {
    std::cerr << "HydroSpectrum::openDataFile2! unable to open file " << fname2 << std::endl;
    std::exit(1);
  }
}

void HydroSpectrum::openDataFile(string fname1, string fname2, string fname_ecc) {
  if (fdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname1 << std::endl;
    std::exit(1);
  }
  if (pdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname2 << std::endl;
    std::exit(1);
  }
  if (eccdata.is_open()) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! file already opened"
         << fname_ecc << std::endl;
    std::exit(1);
  }

  fdata_fname = fname1;
  fdata.open(fname1.c_str(), ios::in);
  if (!fdata) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file " << fname1 << std::endl;
    std::exit(1);
  }
  pdata_fname = fname2;
  pdata.open(fname2.c_str(), ios::in);
  if (!pdata) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file " << fname2 << std::endl;
    std::exit(1);
  }
  eccdata_fname = fname_ecc;
  eccdata.open(fname_ecc.c_str(), ifstream::in);
  if (!eccdata) {
    std::cerr << "HydroSpectrum::openDataFile(fn1,fn2,fnEcc)! unable to open file <" << fname_ecc
	       << ">" << std::endl;
    std::exit(1);
  }
}

// input data of hydrodynamic flow
int HydroSpectrum::readData() {
  if (!fdata) {
    std::cerr << "HydroSpectrum::readData: unexpected end of freezeout.dat file" << std::endl;
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
/**********************************************************
C
C     SUBROUTINE GAUS12
C       SET WEIGHT VSLUE FOR INTEGRATION
C
C**********************************************************/
void HydroSpectrum::Gauss12(double xini, double xfin, double* xn, double* wn)
{
  //double xn[12], wn[12];

  double x[6] = {0.98156063424671925e00, 0.90411725637047486e00,
               0.76990267419430469e00, 0.58731795428661745e00,
               0.36783149899818019e00, 0.12523340851146892e00};

  double w[6] = {0.47175336386511827e-01, 0.10693932599531843e00,
               0.16007832854334623e00, 0.20316742672306592e00,
               0.23349253653835481e00, 0.24914704581340279e00};

  double a1 = 0.5 * (xfin - xini);
  double a2 = 0.5 * (xfin + xini);

  for (int i = 0; i < 6; i++) {
    xn[i] = x[i] * a1 + a2;
    xn[i + 6] = -x[5 - i] * a1 + a2;
    wn[i] = w[i] * a1;
    wn[i + 6] = w[5 - i] * a1;
  }
}

//******************************************************************
//
//        SUBROUTINE GAUSS38-LAGUERRE20
//              DATA IS GIVEN BY H.NAKAMURA
//              THIS SUBROUTINE IS WRITTEN BY T.HIRANO
//                                           '96 11 29
//******************************************************************
void HydroSpectrum::GauLag(double xini, double xmid, double* xn, double* wn)
{
  //double xn(58),wn(58),xini,xmid

  static double x[38] = {
    -9.980499305357e-01, -9.897394542664e-01,
    -9.748463285902e-01, -9.534663309335e-01,
    -9.257413320486e-01, -8.918557390046e-01,
    -8.520350219324e-01, -8.065441676053e-01,
    -7.556859037540e-01, -6.997986803792e-01,
    -6.392544158297e-01, -5.744560210478e-01,
    -5.058347179279e-01, -4.338471694324e-01,
    -3.589724404794e-01, -2.817088097902e-01,
    -2.025704538921e-01, -1.220840253379e-01,
    -4.078514790458e-02,  4.078514790458e-02,
    1.220840253379e-01,   2.025704538921e-01,
    2.817088097902e-01,   3.589724404794e-01,
    4.338471694324e-01,   5.058347179279e-01,
    5.744560210478e-01,   6.392544158297e-01,
    6.997986803792e-01,   7.556859037540e-01,
    8.065441676053e-01,   8.520350219324e-01,
    8.918557390046e-01,   9.257413320486e-01,
    9.534663309335e-01,   9.748463285902e-01,
    9.897394542664e-01,   9.980499305357e-01};

  static double y[20] = {
    7.053988969198e-02,   3.721268180016e-01,
    9.165821024832e-01,   1.707306531028e00,
    2.749199255310e00,    4.048925313850e00,
    5.615174970848e00,    7.459017453700e00,
    9.594392870116e00,    1.203880254193e+01,
    1.481429346483e+01,   1.794889545992e+01,
    2.147878834852e+01,   2.545170266868e+01,
    2.993255471277e+01,   3.501343423054e+01,
    4.083305702714e+01,   4.761999407237e+01,
    5.581079574152e+01,   6.652441652668e+01};

  static double w[38] = {
    5.002880749632e-3,    1.161344471647e-2,
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

  static double u[20] = {
    1.687468018511e-1,    2.912543620060e-1,
    2.666861028666e-1,    1.660024532694e-1,
    7.482606467538e-2,    2.496441732002e-2,
    6.202550828442e-3,    1.144962383071e-3,
    1.557417806303e-4,    1.540143980270e-5,
    1.086486325120e-6,    5.330127202618e-8,
    1.757985532862e-9,    3.725505729478e-11,
    4.767526861010e-13,   3.372848343678e-15,
    1.155014009180e-17,   1.539522135110e-20,
    5.286443248824e-24,   1.656456602226e-28};


  for (int i = 0; i < 38; i++) {
    xn[i] = (xmid - xini) * x[i] / 2. +(xini + xmid) / 2.;
    wn[i] = (xmid - xini) * w[i] / 2.;
  }

  //       do 998 j = 1,20
  for (int j = 0; j < 20; j++) {
    xn[j + 38] = y[j] + xmid;
    wn[j + 38] = exp(y[j]) * u[j];
  }

}

//**********************************************************************
void HydroSpectrum::Gauss38(double xini, double xfin, double* xn, double* wn)
{
  //double xn[38],wn[38];
  double  x[38], w[38];

  x[37] = 9.980499305357e-1;
  x[36] = 9.897394542664e-1;
  x[35] = 9.748463285902e-1;
  x[34] = 9.534663309335e-1;
  x[33] = 9.257413320486e-1;
  x[32] = 8.918557390046e-1;
  x[31] = 8.520350219324e-1;
  x[30] = 8.065441676053e-1;
  x[29] = 7.556859037540e-1;
  x[28] = 6.997986803792e-1;
  x[27] = 6.392544158297e-1;
  x[26] = 5.744560210478e-1;
  x[25] = 5.058347179279e-1;
  x[24] = 4.338471694324e-1;
  x[23] = 3.589724404794e-1;
  x[22] = 2.817088097902e-1;
  x[21] = 2.025704538921e-1;
  x[20] = 1.220840253379e-1;
  x[19] = 4.078514790458e-2;

  //    .....   WEIGHT       ...........
  w[37] = 5.002880749632e-3;
  w[36] = 1.161344471647e-2;
  w[35] = 1.815657770961e-2;
  w[34] = 2.457973973823e-2;
  w[33] = 3.083950054518e-2;
  w[32] = 3.689408159400e-2;
  w[31] = 4.270315850467e-2;
  w[30] = 4.822806186076e-2;
  w[29] = 5.343201991033e-2;
  w[28] = 5.828039914700e-2;
  w[27] = 6.274093339213e-2;
  w[26] = 6.678393797914e-2;
  w[25] = 7.038250706690e-2;
  w[24] = 7.351269258474e-2;
  w[23] = 7.615366354845e-2;
  w[22] = 7.828784465821e-2;
  w[21] = 7.990103324353e-2;
  w[20] = 8.098249377060e-2;
  w[19] = 8.152502928039e-2;

  for (int i = 0; i < 19; i++) {
    x[i] = -x[37 - i];
    w[i] =  w[37 - i];
  }
  for (int i = 0; i <38; i++) {
    xn[i] = (xfin - xini) * x[i] / 2.0 + (xini + xfin) / 2.0;
    wn[i] = (xfin - xini) * w[i] / 2.0;
  }

}

}
}
