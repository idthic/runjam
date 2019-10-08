// -*- mode:c++;indent-tabs-mode:nil -*-
#pragma once
#ifndef DEBUG_SAMPLE_H
#define DEBUG_SAMPLE_H
#include <vector>
#include <string>
#include "HydroParticleCF.h"
#include "ParticleSample.h"

class ParticleSampleDebugJamAsymmetry:public ParticleSampleBase{
public:
  virtual void update();
};

#endif
