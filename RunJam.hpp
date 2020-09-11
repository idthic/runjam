// -*- mode: c++ -*-
#ifndef runjam_RunJam_hpp
#define runjam_RunJam_hpp

#include <string>
#include <fstream>
#include "spectra/ParticleSample.hpp"
#include "args.hpp"

namespace idt {
namespace runjam {

class RunJam {
private:
  double aveNumberPart1, aveNumberPart2;
  int nevent;
  int numberTestParticle;
  int nv;
  int nbary;
  int nmeson;

  std::ofstream ofs;   // file for phase space data output.
  std::ofstream ofs0;   // file for phase space data output before rescattering.

  std::ofstream ofs_bin0;
  std::ofstream ofs_bin;

private:
  void initialize(runjam_context const& iparam);

public:
  RunJam(runjam_context const& iparam);
  ~RunJam();

  void   setNumberOfTestParticle(int i) { numberTestParticle = i; }
  double getIniAverageParticleNumber1() { return aveNumberPart1; }
  double getIniAverageParticleNumber2() { return aveNumberPart2; }
  void   setWeakDecay();
  void   unsetWeakDecay();
  void   setMSTC(int i, int j);
  void   generateEvent(ParticleSampleBase* psamp, std::string const& cascadeMode);
  void   initJam(ParticleSampleBase* psamp);
  void   cmCorrection();
  void   printPhaseSpaceData(std::ofstream& output);

  static int sampleJamID(int ir);
};

}
}
#endif
