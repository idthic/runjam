// -*- mode:c++;indent-tabs-mode:nil -*-
#pragma once
#ifndef spectra_PHASESPACE_DATA_SAMPLE_H
#define spectra_PHASESPACE_DATA_SAMPLE_H
#include <vector>
#include <string>
#include "IParticleSample.h"

class ParticleSampleFromPhasespaceDat:public ParticleSampleBase{
  std::string fname_phasespace_dat;
public:
  ParticleSampleFromPhasespaceDat(std::string const& fname_phasespace_dat)
    :fname_phasespace_dat(fname_phasespace_dat)
  {}

  virtual ParticleIDType::value_type getParticleIdType() const{return ParticleIDType::PDGCode;}
  virtual void update();
};

#endif
