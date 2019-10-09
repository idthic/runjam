// -*- mode:c++;indent-tabs-mode:nil -*-
#ifndef spectra_ParticleSampleRead_hpp
#define spectra_ParticleSampleRead_hpp
#include <vector>
#include <string>
#include "IParticleSample.hpp"

namespace idt {
namespace hydro2jam {

  class ParticleSampleRead: public ParticleSampleBase {
    std::string fname_particlesample_dat;
  public:
    ParticleSampleRead(std::string const& fname_particlesample_dat):
      fname_particlesample_dat(fname_particlesample_dat)
    {}

  public:
    virtual void update();

  private:
    void readFile();
  };

}
}

#endif
