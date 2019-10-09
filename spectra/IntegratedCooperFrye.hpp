// -*- C++ -*-
#ifndef hydro2jam_spectra_IntegratedCooperFrye_hpp
#define hydro2jam_spectra_IntegratedCooperFrye_hpp
#include <ksh/phys/Minkowski.hpp>
namespace idt {
namespace hydro2jam {
  void IntegrateBosonCooperFrye(
    double& dNPos,double& dNNeg,
    const kashiwa::phys::vector4& u,
    const kashiwa::phys::vector4& ds,
    double beta,double mass,double mu
    );
  void IntegrateFermionCooperFrye(
    double& dNPos,double& dNNeg,
    const kashiwa::phys::vector4& u,
    const kashiwa::phys::vector4& ds,
    double beta,double mass,double mu
    );
}
}

#endif
