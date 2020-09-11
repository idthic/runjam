// -*- mode: c++ -*-
#ifndef runjam_RunJam_hpp
#define runjam_RunJam_hpp

#include <string>
#include <fstream>
#include "spectra/ParticleSample.hpp"
#include "jam/Jam1.hpp"
#include "args.hpp"

namespace idt {
namespace runjam {

class RunJam {
private:
  Jam1* jam;
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

  void   setNumberOfTestParticle(int i) {numberTestParticle=i;}
  double getIniAverageParticleNumber1() {return aveNumberPart1;}
  double getIniAverageParticleNumber2() {return aveNumberPart2;}
  void   setWeakDecay() {jam->setMSTC(42,0);}  //=0: allow weak decays
  void   unsetWeakDecay() {jam->setMSTC(42,1);}  //=1:no weak decays
  void   setMSTC(int i,int j) {jam->setMSTC(i,j);}
  void   generateEvent(IParticleSample* psamp, std::string const& cascadeMode);
  void   initJam(IParticleSample* psamp);
  void   cmCorrection();
  void   printPhaseSpaceData(std::ofstream& output);

  static int sampleJamID(int ir);
};

}
}
#endif
