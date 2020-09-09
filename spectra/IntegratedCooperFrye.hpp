// -*- C++ -*-
#ifndef runjam_spectra_IntegratedCooperFrye_hpp
#define runjam_spectra_IntegratedCooperFrye_hpp
#include <ksh/phys/Minkowski.hpp>
namespace idt {
namespace runjam {
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
