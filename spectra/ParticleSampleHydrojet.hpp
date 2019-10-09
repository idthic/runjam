// -*- mode:c++ -*-
#ifndef spectra_ParticleSampleHydrojet_hpp
#define spectra_ParticleSampleHydrojet_hpp
#include <vector>
#include <string>
#include "IParticleSample.hpp"
#include "../args.hpp"

namespace idt {
namespace hydro2jam {

IParticleSample* CreateParticleSampleHydrojet(hydro2jam_context const& ctx);

}
}

#endif
