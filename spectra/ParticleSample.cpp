#include <cmath>
#include <vector>
#include <algorithm>
#include "ParticleSample.hpp"
#include "util.hpp"

//-----------------------------------------------------------------------------
// implementation of class ParticleSampleBase

namespace idt {
namespace runjam {

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

  void ParticleSampleBase::shuffleParticleList() {
    for (std::size_t i  = plist.size(); i >= 2; i--)
      std::swap(plist[i - 1], plist[idt::util::irand(i)]);

// C++ version で振る舞いが変わると再現性に欠くので自分で書く。
// #if __cplusplus >= 201703L
//     std::shuffle(this->plist.begin(), this->plist.end(), idt::util::random_engine());
// #else
//     std::random_shuffle(this->plist.begin(), this->plist.end());
// #endif
  }

  void ParticleSampleBase::adjustCenterOfMassByGalileiBoost() {
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    double s = 0.0;
    for (Particle& part: *this) {
      px += part.px;
      py += part.py;
      pz += part.pz;
      cx += part.x * part.mass;
      cy += part.y * part.mass;
      cz += part.z * part.mass;
      s += part.mass;
    }

    std::size_t c = this->size();
    cx = -cx / s;
    cy = -cy / s;
    cz = -cz / s;
    px = -px / c;
    py = -py / c;
    pz = -pz / c;

    for (Particle& part: *this) {
      part.x += cx;
      part.y += cy;
      part.z += cz;
      part.px += px;
      part.py += py;
      part.pz += pz;
      double const m = part.mass;
      part.e = std::sqrt(m * m + px * px + py * py + pz * pz);
    }
  }

  void ParticleSampleBase::adjustCenterOfMassByLorentzBoost() {
    double u0 = 0.0, u1 = 0.0, u2 = 0.0, u3 = 0.0;

    vector4 u = {};
    for (Particle& part: *this) u += part.mom;
    double norm = u * u;
    if (norm <= 0.0) return;
    u /= std::sqrt(norm);

    double gw = 0.0, g1 = 0.0, g2 = 0.0, g3 = 0.0;
    for (Particle& part: *this) {
      part.mom.boost(~u);
      part.pos.boost(~u);
      gw += part.e;
      g1 += part.e * part.x;
      g2 += part.e * part.y;
      g3 += part.e * part.z;
    }
    g1 /= gw;
    g2 /= gw;
    g3 /= gw;
    for (Particle& part: *this) {
      part.x -= g1;
      part.y -= g2;
      part.z -= g3;
    }
  }

  void OversampledParticleSampleBase::update() {
    // 一括生成済の時
    if (this->pcache.size() > 1) {
      if (++this->indexOfCachedEvents < this->pcache.size() - 1)
        return;

      this->pcache.clear();
      this->indexOfCachedEvents = -1;
    }

    // 一括生成要求がある時
    if (this->numberOfExpectedEvents > 0) {
      int ncache = this->numberOfExpectedEvents;
      this->updateWithOverSampling(this->m_overSamplingFactor * ncache);

      // 二項分布で各イベントの粒子数を決定する
      this->pcache.clear();
      this->pcache.reserve(ncache + 1);
      this->pcache.emplace_back(0);
      std::size_t pos = 0, size = base::plist.size();
      for (; ncache > 1; ncache--) {
        pos += idt::util::irand_binomial(size - pos, 1.0 / ncache);
        this->pcache.emplace_back(pos);
      }
      this->pcache.emplace_back(size);

      this->shuffleParticleList();
      this->numberOfExpectedEvents = 0;
      this->indexOfCachedEvents = 0;
      return;
    }

    this->updateWithOverSampling(this->m_overSamplingFactor);
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
