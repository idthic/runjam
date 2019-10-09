// -*- C++ -*-
#ifndef HYDROJET_SPETRA_IntegratedCooperFrye_H
#define HYDROJET_SPETRA_IntegratedCooperFrye_H
#include <ksh/phys/Minkowski.hpp>
namespace hydrojet{
namespace spectra{
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
