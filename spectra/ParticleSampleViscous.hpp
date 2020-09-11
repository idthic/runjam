// -*- mode: c++ -*-
//
// 2014-05-01 KM,
//   * Created
//
#ifndef runjam_spectra_ParticleSampleViscous_hpp
#define runjam_spectra_ParticleSampleViscous_hpp
#include <vector>
#include <string>
#include "ParticleSample.hpp"
#include "ResonanceList.hpp"

namespace idt {
namespace runjam {

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

  void SampleParticlesC0lrf(
    std::vector<Particle*>& plist,
    HypersurfaceElementC0Lrf const& surface,
    IResonanceList const* rlist,
    double overSamplingFactor = 1.0,
    bool turnsOffViscousEffect = false
  );

  int checkViscousCooperFryeInterpolated(bool debug);

}
}

#endif
