#include <cmath>
#include <vector>
#include "IParticleSample.hpp"

//-----------------------------------------------------------------------------
// implementation of class ParticleSampleBase

namespace idt {
namespace runjam {

void ParticleSampleBase::addParticleMinkowski(int iReso, double px, double py, double pz, double m, double x, double y, double z, double t) {
  double const e = m < 0.0 ? -1.0 : std::sqrt(px * px + py * py + pz * pz + m * m);
  Particle* particle = new Particle(iReso);
  particle->px = px;
  particle->py = py;
  particle->pz = pz;
  particle->e = e;
  particle->x = x;
  particle->y = y;
  particle->z = z;
  particle->t = t;
  this->plist.push_back(particle);
}

void ParticleSampleBase::addParticleTauEta(int iReso, double px, double py, double pz, double m, double x, double y, double eta, double tau) {
  double const t = tau * std::cosh(eta);
  double const z = tau * std::sinh(eta);
  this->addParticleMinkowski(iReso, px, py, pz, m, x, y, z, t);
}

void ParticleSampleBase::clearParticleList() {
  for(std::vector<Particle*>::const_iterator i = plist.begin(); i != plist.end(); ++i)
    delete *i;
  this->plist.clear();
}

  static std::vector<IParticleSampleFactory*> particleSampleFactories;

  IParticleSample* CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
    std::size_t n = particleSampleFactories.size();
    for (std::size_t i = 0; i < n; i++)
      if (IParticleSample* ret = particleSampleFactories[i]->CreateInstance(ctx, type, inputfile))
        return ret;
    return 0;
  }

  void IParticleSampleFactory::Register(IParticleSampleFactory* factory) {
    particleSampleFactories.push_back(factory);
  }


}
}
