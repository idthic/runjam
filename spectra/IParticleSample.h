// -*- mode:c++;coding:utf-8;indent-tabs-mode:nil -*-
//
// 2014-05-09 KM,
//   * Created
//   * class IParticleSample was moved from ParticleSample.h.
//   * Added comments.
//
#ifndef spectra_IParticleSample_h
#define spectra_IParticleSample_h
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "HydroParticleCF.h"

struct ParticleIDType{
  enum value_type{
    HydroParticleID, // usually denoted as id
    PDGCode,         // usually denoted as kf
    JamInternalCode, // usually denoted as kc
  };
};

class IParticleSample{
public:
  virtual ~IParticleSample(){}

  /// @fn int getIsOutput() const
  /// \~en returns if the resonance distributions are output or not.
  /// - the return value 0 means that the distributions are passed to initJam with IParticleSample::getParticleList().
  /// - the return value 1 means that the distributions are passed to initJam with a file.
  ///   The filename to whom the distributions is written is obtained by IParticleSample::getFileNamePos().
  /// \~ja 粒子分布を出力するかどうかを返します。
  /// - 0 を返した場合は、 IParticleSample::getParticleList() を使用して粒子分布が initJam に渡されます。
  /// - 1 を返した場合は、一旦粒子分布をファイルに書き込んだ後 initJam にそのファイル名が渡されます。
  ///   書き込み先ファイル名は IParticleSample::getFileNamePos() で取得します。
  virtual int getIsOutput() const=0;

  /// @fn std::string const& getFileNamePos() const;
  /// \~en returns the filename containing the output resonance distributions when getIsOutput()==1.
  /// \~ja getIsOutput()==true の時に、粒子分布を格納したファイル名を返します。
  virtual std::string const& getFileNamePos() const=0;

  virtual void setAdviceNumberOfExpectedEvents(int nEvents){}

  /// @fn void update();
  /// \~en generates a resonance distribution
  /// \~ja 粒子分布の生成を実行します。
  virtual void update()=0;

  /// @fn std::vector<HydroParticleCF*> const& getParticleList() const;
  /// \~en retrieves the generated resonance distribution.
  /// \~ja 生成した粒子分布を取得します。
  virtual std::vector<HydroParticleCF*> const& getParticleList() const=0;

  virtual ParticleIDType::value_type getParticleIdType() const{return ParticleIDType::HydroParticleID;}
};

class ParticleSampleBase:public IParticleSample{
protected:
  std::vector<HydroParticleCF*> plist;

protected:
  void clearParticleList();

  /// @param[in] px [GeV/c] momentum in the x-direction
  /// @param[in] py [GeV/c] momentum in the y-direction
  /// @param[in] pz [GeV/c] momentum in the z-direction
  /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
  void addParticleMinkowski(int iReso,double px,double py,double pz,double m,double x,double y,double z,double t);

  /// @param[in] px [GeV/c] momentum in the x-direction
  /// @param[in] py [GeV/c] momentum in the y-direction
  /// @param[in] pz [GeV/c] momentum in the z-direction
  /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
  void addParticleTauEta(int iReso,double px,double py,double pz,double m,double x,double y,double tau,double eta);

public:
  ParticleSampleBase(){}
  virtual ~ParticleSampleBase(){this->clearParticleList();}

  virtual int getIsOutput() const{return false;}
  virtual std::string const& getFileNamePos() const{static std::string dummy;return dummy;}
  virtual std::vector<HydroParticleCF*> const& getParticleList() const{return this->plist;}

private:
  // コピー禁止
  ParticleSampleBase(ParticleSampleBase const&){
    std::cerr<<"ParticleSampleBase(copy ctor): copy not supported!"<<std::endl;
    std::exit(1);
  }
  ParticleSampleBase& operator=(ParticleSampleBase const&){
    std::cerr<<"ParticleSampleBase(copy assign): copy not supported!"<<std::endl;
    std::exit(1);
  }
};

#endif
