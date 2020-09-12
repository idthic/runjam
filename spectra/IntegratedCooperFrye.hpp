// -*- C++ -*-
#ifndef runjam_spectra_IntegratedCooperFrye_hpp
#define runjam_spectra_IntegratedCooperFrye_hpp
#include "ParticleSample.hpp"
namespace idt {
namespace runjam {
  void IntegrateBosonCooperFrye(
    double& dNPos,double& dNNeg,
    const vector4& u, const vector4& ds,
    double beta,double mass,double mu
    );
  void IntegrateFermionCooperFrye(
    double& dNPos,double& dNNeg,
    const vector4& u, const vector4& ds,
    double beta,double mass,double mu
    );
}
}

#endif
