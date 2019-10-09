// -*- mode:c++;coding:utf-8;indent-tabs-mode:nil -*-
//
// 2015-03-25 KM
//   * class IResonanceList: moved from ParticleSampleViscous.h
//
#ifndef spectra_IResonanceData_h
#define spectra_IResonanceData_h
#include <vector>
#include <string>
#include <util/Constants.hpp>
#include "../args.hpp"

namespace idt {
namespace hydro2jam {

class IResonanceList{
public:
  //! 指定した共鳴の質量を返します。単位は [fm^{-1}] です。
  //! @param[in] ireso
  virtual double mass             (int ireso) const=0;
  //! 指定した共鳴の統計の符号を返します。
  //! boson に対して +1 を、fermion に対して -1 を返します。
  //! @param[in] ireso
  virtual double statisticsSign   (int ireso) const=0;
  //! 指定した共鳴の化学ポテンシャルを返します。
  //! @param[in] ireso
  virtual double chemicalPotential(int ireso) const=0;
  //! 指定した共鳴の自由度 (isospin, spin, etc.) の数 g を返します。
  //! @param[in] ireso
  virtual int numberOfDegrees  (int ireso) const=0;
  //! 共鳴の数を返します。
  virtual int numberOfResonances() const=0;

  virtual ~IResonanceList(){}
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
  };

  int m_numberOfResonances;

  static int nreso;

private:
  static resonance resT[5][21];

  std::vector<resonance> data;

  void initialize(int kineticTemp, int eos_pce,std::string const& fname_rlist);
public:
  ResonanceListPCE(hydro2jam_context const& ctx);
  ResonanceListPCE(int kineticTemp, int eos_pce,std::string const& fname_rlist);

  resonance const& operator[](int ireso) const{
    return this->data[ireso];
  }
  resonance& operator[](int ireso){
    return this->data[ireso];
  }

public:
  double mass             (int ireso) const{return this->data[ireso].mass;}
  double statisticsSign   (int ireso) const{return -this->data[ireso].bf;} // +1 for boson, -1 for fermion
  double chemicalPotential(int ireso) const{return this->data[ireso].mu;}
  int numberOfDegrees  (int ireso) const{return this->data[ireso].deg;}
  int numberOfResonances() const{return this->m_numberOfResonances;}
};

}
}

#endif
