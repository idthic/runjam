// -*- mode:c++;indent-tabs-mode:nil -*-
#ifndef ParticleSample_h
#define ParticleSample_h
#include <vector>
#include <string>
#include "ElementReso.h"
#include "HydroParticleCF.h"
#include "IParticleSample.h"

class ParticleSample:public ElementReso,public IParticleSample
{
private:
  std::ifstream *resDataPos;
  std::string   outfilepos;
  std::ofstream outdatpos;
  std::ifstream *resDataNeg;
  std::string   outfileneg;
  std::ofstream outdatneg;
  int di; // iteration number in bisection method
  int isOutput;
  std::vector<HydroParticleCF*> plist;
  int    baryonfree;
  double tmpf;
  double mubf;
  double meanf;

  bool mode_delayed_cooperfrye;
  bool fReverseParticleList;
  bool fShuffleParticleList;
public:
  ParticleSample(std::string dir, std::string* outf, int kin, int eos_pce,
  std::string fname);
  ~ParticleSample();
  void setBaryonFree(int i) {baryonfree=i;}
  void setTMPF(double t) {tmpf=t/sctr*1000.0;}
  void setMUBF(double m) {mubf=m/sctr*1000.0;}
  void initialize(std::string const& fn,std::string const& fn_p);
  void analyze(std::string fn, std::string fn_p);
  void finish();
  double dx,dy,dh,dtau;
  void setDx(double d) {dx=d;}
  void setDy(double d) {dy=d;}
  void setDh(double d) {dh=d;}
  void setDtau(double d) {dtau=d;}
  double getDtau() {return dtau;}
  double getDx() {return dx;}
  double getDy() {return dy;}
  double getDh() {return dh;}
  void setIsOutput(int i){this->isOutput=i;}
  int getIsOutput() const{return this->isOutput;}
  std::string const& getFileNamePos() const{return outfilepos;}
  std::vector<HydroParticleCF*> const& getParticleList() const{return plist;}

private:
  std::string fn_freezeout_dat;
  std::string fn_position_dat;
public:
  void setHypersurfaceFilenames(std::string const& fn_freezeout_dat,std::string const& fn_position_dat){
    this->fn_freezeout_dat=fn_freezeout_dat;
    this->fn_position_dat=fn_position_dat;
  }
  void update(){
    this->analyze(this->fn_freezeout_dat,this->fn_position_dat);
  }

private:
  void getSample(
    double vx,double vy,double vz,
    double ds0,double dsx,double dsy,
    double dsz,int ir, int ipos,
    double tau,double xx,double yy,double eta);
  void putParticle(
    double px,double py,double pz,
    double e,double m, int ir, double tau,double x,
    double y, double eta,int ipos);
  void outputData(
    double prx,double pry,double prz,
    double er,double mres, int ir, double tau,double xx,
    double yy, double eta, int ipos);
};

#endif
