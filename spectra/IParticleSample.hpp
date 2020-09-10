// -*- mode: c++ -*-
#ifndef runjam_spectra_IParticleSample_hpp
#define runjam_spectra_IParticleSample_hpp
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../args.hpp"

namespace idt {
namespace runjam {

struct Particle {
public:
  int pdg;
  double mass;

  double x;
  double y;
  double z;
  double t;

  double e;
  double px;
  double py;
  double pz;
public:
  Particle(int pdg): pdg(pdg) {
    this->mass = -1.0;
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->t = 0.0;
    this->px = 0.0;
    this->py = 0.0;
    this->pz = 0.0;
    this->e  = 0.0;
  }
  Particle() {}

};

class IParticleSample {
public:
  virtual ~IParticleSample() {}

  virtual void setAdviceNumberOfExpectedEvents(int nEvents) {}

  /// @fn void update();
  /// \~en generates a resonance distribution
  /// \~ja 粒子分布の生成を実行します。
  virtual void update() = 0;

  /// @fn std::vector<Particle*> const& getParticleList() const;
  /// \~en retrieves the generated resonance distribution.
  /// \~ja 生成した粒子分布を取得します。
  virtual std::vector<Particle*> const& getParticleList() const = 0;
};

class ParticleSampleBase: public IParticleSample {
protected:
  std::vector<Particle*> plist;

protected:
  void clearParticleList();

  /// @param[in] px [GeV/c] momentum in the x-direction
  /// @param[in] py [GeV/c] momentum in the y-direction
  /// @param[in] pz [GeV/c] momentum in the z-direction
  /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
  void addParticleMinkowski(int iReso, double px, double py, double pz, double m, double x, double y, double z, double t);

  /// @param[in] px [GeV/c] momentum in the x-direction
  /// @param[in] py [GeV/c] momentum in the y-direction
  /// @param[in] pz [GeV/c] momentum in the z-direction
  /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
  void addParticleTauEta(int iReso, double px, double py, double pz, double m, double x, double y, double tau, double eta);

public:
  ParticleSampleBase() {}
  virtual ~ParticleSampleBase() { this->clearParticleList(); }
  virtual std::vector<Particle*> const& getParticleList() const { return this->plist; }

private:
  // コピー禁止
  ParticleSampleBase(ParticleSampleBase const&) {
    std::cerr << "ParticleSampleBase(copy ctor): copy not supported!" << std::endl;
    std::exit(1);
  }
  ParticleSampleBase& operator=(ParticleSampleBase const&) {
    std::cerr << "ParticleSampleBase(copy assign): copy not supported!" << std::endl;
    std::exit(1);
  }
};

  class IParticleSampleFactory {
  public:
    virtual IParticleSample* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) = 0;
    virtual ~IParticleSampleFactory() {}
  protected:
    void Register(IParticleSampleFactory* factory);
  };

  class ParticleSampleFactoryRegistered: public IParticleSampleFactory {
  protected:
    ParticleSampleFactoryRegistered() { Register(this); }
  };

  IParticleSample* CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile);

}
}

#endif
