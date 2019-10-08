#include <cmath>
#include <vector>
#include "IParticleSample.h"

//-----------------------------------------------------------------------------
// implementation of class ParticleSampleBase

void ParticleSampleBase::addParticleMinkowski(int iReso,double px,double py,double pz,double m,double x,double y,double z,double t){
  //const double hbarC=0.197327053;
  HydroParticleCF* particle=new HydroParticleCF(iReso);
  double const e=m<0.0?-1.0:std::sqrt(px*px+py*py+pz*pz+m*m);
  particle->setPx(px);
  particle->setPy(py);
  particle->setPz(pz);
  particle->setPe(e );
  particle->setX(x);
  particle->setY(y);
  particle->setT(t);
  particle->setZ(z);
  this->plist.push_back(particle);
}

void ParticleSampleBase::addParticleTauEta(int iReso,double px,double py,double pz,double m,double x,double y,double eta,double tau){
  double const t=tau*std::cosh(eta);
  double const z=tau*std::sinh(eta);
  this->addParticleMinkowski(iReso,px,py,pz,m,x,y,z,t);
}

void ParticleSampleBase::clearParticleList(){
  for(std::vector<HydroParticleCF*>::const_iterator i=plist.begin();i!=plist.end();++i)
    delete (*i);
  this->plist.clear();
}
