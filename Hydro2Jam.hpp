// -*- mode:c++ -*-
#ifndef Hydro2Jam_h
#define Hydro2Jam_h

#include <string>
#include <fstream>
#include "spectra/IParticleSample.hpp"
#include "jam/Jam1.hpp"
#include "args.hpp"

namespace idt {
namespace hydro2jam {

class Hydro2Jam {
private:
  Jam1* jam;
  double   aveNumberPart1, aveNumberPart2;
  std::string dirReso;
  int kinTmp;
  int eosPCE;
  std::string resodata;
  int nevent;
  int numberTestParticle;
  int nv;
  int nbary;
  int nmeson;

  std::ofstream ofs;   // file for phase space data output.
  std::ofstream ofs0;   // file for phase space data output before rescattering.
  int dumpPhaseSpaceData;

  bool flag_decayOnly;

private:
  void initialize(hydro2jam_context const& iparam);

public:
  Hydro2Jam(hydro2jam_context const& iparam);
  ~Hydro2Jam();

  void   setResoData(std::string d) {resodata=d;}
  void   setNumberOfTestParticle(int i) {numberTestParticle=i;}
  double getIniAverageParticleNumber1() {return aveNumberPart1;}
  double getIniAverageParticleNumber2() {return aveNumberPart2;}
  void   setWeakDecay() {jam->setMSTC(42,0);}  //=0: allow weak decays
  void   unsetWeakDecay() {jam->setMSTC(42,1);}  //=1:no weak decays
  void   setMSTC(int i,int j) {jam->setMSTC(i,j);}
  void   generateEventFromHypersurfaceFiles(
    std::string const& fn_freezeout_dat,
    std::string const& fn_position_dat,
    int baryonfree, double deltat, double deltax, double deltay, double deltah);
  void   generateEvent(IParticleSample* psamp);
  void   initJam(IParticleSample* psamp);
  void   cmCorrection();
  void   printPhaseSpaceData(std::ofstream& output);

  static int sampleJamID(int ir);
};

}
}
#endif
