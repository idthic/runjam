// -*- mode: c++ -*-
#pragma once
#ifndef runjam_spectra_ElementReso_hpp
#define runjam_spectra_ElementReso_hpp
#include <vector>
#include "HydroSpectrum.hpp"
#include "IResonanceList.hpp"

namespace idt {
namespace runjam {

class ElementReso: public HydroSpectrum {
protected:
  ResonanceListPCE rlist;

  std::ifstream fdata;
  std::ifstream edata;
  std::vector<std::string>   elemFile;
  std::vector<std::ofstream> outdat;
  std::vector<std::ofstream> outdatPos;
  std::vector<std::ofstream> outdatNeg;

  int    baryonfree;
  double tmpf;
  double mubf;
  double meanf;

public:
  double ymin, ymax;
  ElementReso(std::string dir, std::string* fname, int kint, int eos_pce, std::string fname2);
  ~ElementReso();
  void setBaryonFree(int i) { baryonfree = i; }
  void setTMPF(double t) { tmpf = t / hbarc_MeVfm * 1000.0; }
  void setMUBF(double m) { mubf = m / hbarc_MeVfm * 1000.0; }
  void initialize();
  void analyze(std::string fn);
private:
  void integrateForResonance(std::string const& fname_freezeout_dat, int ireso);
};

}
}

#endif
