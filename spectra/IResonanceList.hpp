// -*- mode: c++ -*-
//
// 2015-03-25 KM
//   * class IResonanceList: moved from ParticleSampleViscous.h
//
#ifndef runjam_spectra_IResonanceData_hpp
#define runjam_spectra_IResonanceData_hpp
#include <vector>
#include <string>
#include <args.hpp>

namespace idt {
namespace runjam {

class IResonanceList {
public:
  //! 共鳴の数を返します。
  virtual int numberOfResonances() const = 0;

  //! 指定した共鳴の質量を返します。単位は [fm^{-1}] です。
  //! @param[in] ireso
  virtual double mass(int ireso) const = 0;
  //! 指定した共鳴の統計の符号を返します。
  //! boson に対して +1 を、fermion に対して -1 を返します。
  //! @param[in] ireso
  virtual double statisticsSign(int ireso) const = 0;
  //! 指定した共鳴の化学ポテンシャルを返します。
  //! @param[in] ireso
  virtual double chemicalPotential(int ireso) const = 0;
  //! 指定した共鳴の自由度 (isospin, spin, etc.) の数 g を返します。
  //! @param[in] ireso
  virtual int numberOfDegrees(int ireso) const = 0;
  //! 指定した共鳴に対応する PDG Monte-Carlo code を生成します。
  //! @param[in] ireso
  virtual int generatePDGCode(int ireso) const = 0;

  virtual ~IResonanceList() {}
};

class ResonanceListPCE: public IResonanceList {
public:
  struct resonance {
    double mass;
    double deg;
    double degeff;
    double mu;
    int    bf;
    int    anti;
    std::vector<int> pdg_codes;
  };

private:
  static resonance resT[5][21];

  std::vector<resonance> data;

  void initialize(int kineticTemp, int eos_pce,std::string const& fname_rlist);
public:
  ResonanceListPCE(runjam_context const& ctx);
  ResonanceListPCE(int kineticTemp, int eos_pce,std::string const& fname_rlist);

  resonance const& operator[](int ireso) const {
    return this->data[ireso];
  }
  resonance& operator[](int ireso) {
    return this->data[ireso];
  }

public:
  virtual double mass(int ireso) const override { return this->data[ireso].mass; }
  virtual double statisticsSign(int ireso) const override { return -this->data[ireso].bf; } // +1 for boson, -1 for fermion
  virtual double chemicalPotential(int ireso) const override { return this->data[ireso].mu; }
  virtual int numberOfDegrees(int ireso) const override { return this->data[ireso].deg; }
  virtual int generatePDGCode(int ireso) const override;
  virtual int numberOfResonances() const override { return this->data.size(); }
};

}
}

#endif
