// -*- mode:c++;indent-tabs-mode:nil -*-
#pragma once
#ifndef ElementReso_h
#define ElementReso_h
#include <vector>
#include "HydroSpectrum.hpp"
#include "IResonanceList.hpp"

namespace idt {
namespace hydro2jam {

class ElementReso : public HydroSpectrum{
protected:
  ResonanceListPCE rlist;
  //double  *mass, *deg, *degeff, *mu,  *anti;
  //int bf[151];

  double  phi[12], phiw[12], y[38], yw[38];
  double  pt[58], ptw[58];
  std::string elemFile[151];
  std::ifstream fdata;
  std::ifstream edata;
  std::ofstream outdat[151];
  std::ofstream outdatPos[151];
  std::ofstream outdatNeg[151];

  int    baryonfree;
  double tmpf;
  double mubf;
  double meanf;
  int    nreso_loop;

public:
  double ymin, ymax;
  ElementReso(std::string dir, std::string* fname, int kint, int eos_pce, std::string fname2);
  //    ElementReso(std::string dir, std::string* fname, int kint);
  ~ElementReso();
  void setBaryonFree(int i) {baryonfree=i;}
  //    void setNresoLoop(int i) {nreso_loop=i;}
  void setTMPF(double t) {tmpf=t / hbarc_MeVfm * 1000.0;}
  void setMUBF(double m) {mubf=m / hbarc_MeVfm * 1000.0;}
  void initialize();
  void analyze(std::string fn);
private:
  void integrateForResonance(std::string const& fname_freezeout_dat, int ireso);
};

}
}

#endif
