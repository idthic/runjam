// -*- mode:c++;coding:utf-8;indent-tabs-mode:nil -*-
//
// 2014-05-01 KM,
//   * Created
//
#pragma once
#ifndef spectra_ParticleSampleViscous_h
#define spectra_ParticleSampleViscous_h
#include <vector>
#include <string>
#include "IParticleSample.hpp"
#include "IResonanceList.hpp"

//-----------------------------------------------------------------------------
// prototype declarations

struct HypersurfaceElementC0Lrf;

void SampleParticlesC0lrf(
  std::vector<Particle*>& plist,
  HypersurfaceElementC0Lrf const& surface,
  IResonanceList const* rlist,
  double overSamplingFactor = 1.0,
  bool turnsOffViscousEffect = false
);

class OversampledParticleSampleBase;
class ParticleSampleViscous;
class ParticleSampleFromHydrojet;

int checkViscousCooperFryeInterpolated(bool debug);

//-----------------------------------------------------------------------------
// class definitions

struct HypersurfaceElementC0Lrf {
  double m_ds[4];           //!< 曲線座標での面素ベクトル
  double m_pos[4];          //!< 曲線座標での位置 (面素の中心位置)
  double m_dx[4];           //!< 曲線座標での幅
  double m_velocity[4];     //!< チルダ座標での流速
  double m_stress[6];       //!< 局所静止系での応力
  double m_stressMax;       //!< 第一主応力

  double m_temperature;     //!< 温度
  double m_energy;          //!< エネルギー密度
  double m_pressure;        //!< 圧力

public:
  double const& surfaceElement(int covariantIndex)     const { return this->m_ds[covariantIndex]; }
  double const& position      (int contravariantIndex) const { return this->m_pos[contravariantIndex]; }
  //double const& flowVelocity  (int contravariantIndex) const { return m_velocity[contravariantIndex]; }
  double temperature() const { return this->m_temperature; }
  //double  chemicalPotential(int ireso) const { return 0.0; } // 未実装
};

class OversampledParticleSampleBase: public ParticleSampleBase {
  typedef ParticleSampleBase base;

private:
  // default value = 1.0
  double m_overSamplingFactor;
public:
  void setOverSamplingFactor(double value) {
    this->m_overSamplingFactor = value;
  }
  double getOverSamplingFactor() const {
    return this->m_overSamplingFactor;
  }

  // 実装: 複数事象一括生成について。
  //
  // 1 numberOfExpectedEvents が有限の値に設定されている時、一括生成が要求されている事を意味する。
  //   numberOfExpectedEvents は setAdviceNumberOfExpectedEvents を通して設定できる。
  //   一括生成が要求されている時に update が呼ばれると一括生成が実行され、
  //   numberOfExpectedEvents は 0 にクリアされる。
  // 2 一括生成された粒子は base::plist に保持され寿命が管理される。
  //   同時に、事象 #ievent の粒子一覧は pcache[ievent] に記録される。
  //   pcache.size()>0 の時、一括生成された事象が未だ残っている事を表す。
  //   (pcache.size()>0 の間 base::plist には一括生成された全粒子が格納されている事になる。)
  // 3 pcache の事象を使い切ると pcache はクリアされる。
  //   この場合は通常の 1 事象の生成が行われる。
  //   その過程で、今迄一括生成の全粒子 plist も解放・クリアされる。
  //

private:
  int numberOfExpectedEvents;
  int indexOfCachedEvents;
  std::vector<std::vector<Particle*> > pcache;

public:
  virtual std::vector<Particle*> const& getParticleList() const {
    if (this->indexOfCachedEvents >= 0)
      return this->pcache[indexOfCachedEvents];
    else
      return this->plist;
  }

public:
  virtual void setAdviceNumberOfExpectedEvents(int nEvents) {
    this->numberOfExpectedEvents = nEvents;
    this->indexOfCachedEvents = -1;
  }

public:
  OversampledParticleSampleBase(): m_overSamplingFactor(1.0) {
    this->setAdviceNumberOfExpectedEvents(0);
  }

  virtual void updateWithOverSampling(double overSamplingFactor) = 0;

  virtual void update();
};

class ParticleSampleViscous: public OversampledParticleSampleBase {
  typedef OversampledParticleSampleBase base;

  IResonanceList* rlist;
  std::string fname_hypersurface;

private:
  bool m_turnsOffViscousEffect;
public:
  bool getTurnsOffViscousEffect() const {
    return this->m_turnsOffViscousEffect;
  }
  void setTurnsOffViscousEffect(bool value) {
    this->m_turnsOffViscousEffect = value;
  }

private:
  double m_switchingTemperature;
public:
  double getSwitchingTemperature() const {
    return this->m_switchingTemperature;
  }
  void setSwitchingTemperature(double value) {
    this->m_switchingTemperature = value;
  }

public:
  ParticleSampleViscous(IResonanceList* rlist, std::string const& fname_hypersurface);

  virtual void updateWithOverSampling(double overSamplingFactor);
};

class ParticleSampleFromHydrojet: public OversampledParticleSampleBase {
  typedef OversampledParticleSampleBase base;
  IResonanceList* rlist;
  std::string fname_freezeout;
  std::string fname_position;

private:
  double m_switchingTemperature;
public:
  double getSwitchingTemperature() const {
    return this->m_switchingTemperature;
  }
  void setSwitchingTemperature(double value) {
    this->m_switchingTemperature = value;
  }

private:
  double dx, dy, dh, dtau;
public:
  void setDx(double d) { dx = d; }
  void setDy(double d) { dy = d; }
  void setDh(double d) { dh = d; }
  void setDtau(double d) { dtau = d; }
  double getDtau() { return dtau; }
  double getDx() { return dx; }
  double getDy() { return dy; }
  double getDh() { return dh; }

public:
  ParticleSampleFromHydrojet(
    IResonanceList* rlist,
    std::string const& fname_freezeout,
    std::string const& fname_position);
  ParticleSampleFromHydrojet(
    IResonanceList* rlist,
    std::string const& dname_hydro);

private:
  virtual void updateWithOverSampling(double overSamplingFactor);
private:
  bool readHypersurfaceElement(HypersurfaceElementC0Lrf& surface, std::ifstream& ifsf, std::ifstream& ifsp) const;
};

#endif
