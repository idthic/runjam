// -*- mode:c++;indent-tabs-mode:nil -*-
#pragma once
#ifndef spectra_PARTICLE_SAMPLE_FROM_OVERSAMPLED_PHASESPACE_H
#define spectra_PARTICLE_SAMPLE_FROM_OVERSAMPLED_PHASESPACE_H
#include <vector>
#include <string>
#include "IParticleSample.hpp"

class ParticleSampleFromOversampledPhasespace:public ParticleSampleBase{
  std::string fname_phasespace_dat;
public:
  ParticleSampleFromOversampledPhasespace(std::string const& fname_phasespace_dat)
    :fname_phasespace_dat(fname_phasespace_dat)
  {
    this->m_currentSampleIndex=-1;
    this->m_overSamplingFactor=-1;
  }

private:
  int m_overSamplingFactor;
  int m_currentSampleIndex;
  std::vector<std::vector<Particle*> > pcache;
public:
  void setOverSamplingFactor(int value) {
    this->m_overSamplingFactor = value;
  }
  int getOverSamplingFactor() const {
    return m_overSamplingFactor;
  }

public:
  virtual ParticleIDType::value_type getParticleIdType() const{return ParticleIDType::PDGCode;}
  virtual void update();

  virtual std::vector<Particle*> const& getParticleList() const{
    if(m_currentSampleIndex<0){
      std::cerr<<"ParticleSampleFromOversampledPhasespace: A phasespace file has not been read."<<std::endl;
      std::exit(EXIT_FAILURE);
    }else if(m_currentSampleIndex>=this->pcache.size()){
      std::cerr<<"ParticleSampleFromOversampledPhasespace: No more samples."<<std::endl;
      std::exit(EXIT_FAILURE);
    }

    return this->pcache[m_currentSampleIndex];
  }

private:
  void readPhasespaceDat();
};

#endif
