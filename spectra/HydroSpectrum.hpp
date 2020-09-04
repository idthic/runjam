// -*- mode: c++ -*-
#ifndef hydro2jam_spectra_HydroSpectrum_hpp
#define hydro2jam_spectra_HydroSpectrum_hpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

namespace idt {
namespace hydro2jam {

class HydroSpectrum {
protected:
  double hbarc, pi;
  int degpi, degk, degp; // degree of freedom.
  double  mpi, mk, mpro; // masses of pions, kaons and protons.
  double  mupi, muk, mup; // chemical potentials.
  int co; //=1: v1  =2:v2  =3: v3
  int baryonfree;
  int kineticTemp;

  int fin;

  // FData
  int bulk, iw;
  double dss, ds0, dsx, dsy, dsz, tf, nbf, vx, vy, vz, yv, hh, beta;

  // PData
  double tau, xx, yy, eta;

  // EData
  double hh_ecc[100], coe_ecc[100], ecc_ecc[100];
  double area_ecc[100], aveene_ecc[100], eccp_ecc[100];

  std::string fdata_fname;
  std::string pdata_fname;
  std::string eccdata_fname;
  std::ifstream fdata;
  std::ifstream pdata;
  std::ifstream eccdata;

private:
  bool fRotateFreezeoutData;
public:
  static const double mPion, mKaon, mProton;
  HydroSpectrum(int kint, int eos_pce);

  void openFDataFile(std::string fn_freezeout_dat);
  void openPDataFile(std::string fname2);
  void openEDataFile(std::string fname_ecc);
  void closeFDataFile() { fdata.clear(); fdata.close(); }
  void closePDataFile() { pdata.clear(); pdata.close(); }
  void closeEDataFile() { eccdata.clear();eccdata.close(); }
  int  readFData();
  int  readPData();
  int  readEData();

  void setCoef(int c) {co = c;}
  static void Gauss12(double xini, double xfin, double* xn, double* wn);
  static void Gauss38(double xini, double xfin, double* xn, double* wn);
  static void GauLag(double xini, double xmid, double* xn, double* wn);

protected:
  double pti, cosp, sinp, gam;

  /// @fn double thermaldist(double ee, double pz, double mu, int iw, int sgn);
  /// @param[in] ee
  /// @param[in] pz
  /// @param[in] mu
  /// @param[in] iw
  /// @param[in] sgn 1: fermion  -1:boson
  /// @remarks pti, cosp, sinp, and gam need to be set before call of
  ///   this function.
  double thermaldist(double ee, double pz, double mu, int iw, int sgn);
  //double thermaldist(double mt, double yy, double mu, int iw, int sgn);
};

}
}

#endif
