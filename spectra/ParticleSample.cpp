#include <cmath>
#include <vector>
#include <algorithm>
#include "ParticleSample.hpp"
#include "util.hpp"

//-----------------------------------------------------------------------------
// implementation of class ParticleSampleBase

namespace idt {
namespace runjam {

  void ParticleSampleBase::shuffleParticleList() {
#if __cplusplus >= 201703L
    std::shuffle(this->plist.begin(), this->plist.end(), idt::util::random_engine());
#else
    //std::random_shuffle(this->plist.begin(), this->plist.end());
    std::size_t i  = plist.size();
    while (i-- > 1) std::swap(plist[i], plist[idt::util::irand(i)]);
#endif
  }

  void ParticleSampleBase::addParticleMinkowski(int pdg, double px, double py, double pz, double m, double x, double y, double z, double t) {
    double const e = m < 0.0 ? -1.0 : std::sqrt(px * px + py * py + pz * pz + m * m);
    this->plist.emplace_back();
    Particle& particle = this->plist.back();
    particle.pdg = pdg;
    particle.mass = m;
    particle.px = px;
    particle.py = py;
    particle.pz = pz;
    particle.e = e;
    particle.x = x;
    particle.y = y;
    particle.z = z;
    particle.t = t;
  }

  void ParticleSampleBase::addParticleTauEta(int pdg, double px, double py, double pz, double m, double x, double y, double eta, double tau) {
    double const t = tau * std::cosh(eta);
    double const z = tau * std::sinh(eta);
    this->addParticleMinkowski(pdg, px, py, pz, m, x, y, z, t);
  }

  void ParticleSampleBase::clearParticleList() {
    this->plist.clear();
  }

  static std::vector<ParticleSampleFactoryBase*> particleSampleFactories;

  ParticleSampleBase* CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
    std::size_t n = particleSampleFactories.size();
    for (std::size_t i = 0; i < n; i++)
      if (ParticleSampleBase* ret = particleSampleFactories[i]->CreateInstance(ctx, type, inputfile))
        return ret;
    return 0;
  }

  void ParticleSampleFactoryBase::Register(ParticleSampleFactoryBase* factory) {
    particleSampleFactories.push_back(factory);
  }


}
}
