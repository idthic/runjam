#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstdio>
#include <ksh/integrator.hpp>
#include "IntegratedCooperFrye.hpp"

namespace idt {
namespace runjam {
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

  static double const sqrtTangentLowerBound = 0;
  static double const sqrtTangentUpperBound = std::sqrt(M_PI/2);

  template<int BF>
  struct BFTraits {
    static double Integral2(double xsig, double bmu) {
      return kashiwa::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound, sqrtTangentUpperBound, [=](double t) -> double {
        double tantt = std::tan(t * t);
        double jacob = 2*t*(tantt * tantt + 1);
        double x = tantt + xsig;
        double fE2 = x * x / (std::exp(x - bmu) - BF);
        return jacob * fE2;
      });
    }
    static double IntegralP(double xsig, double bmu, double bmass) {
      return kashiwa::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound, sqrtTangentUpperBound, [=](double t) -> double {
        double tantt = std::tan(t * t);
        double jacob = 2 * t * (tantt * tantt + 1);
        double x = tantt + xsig;
        double fEP = x * std::sqrt(x * x - bmass * bmass) / (std::exp(x - bmu) - BF);
        return jacob * fEP;
      });
    }
    static double Integral0(double xsig, double bmu) {
      return -BF * std::log(1 - BF * std::exp(bmu - xsig));
    }
  };

  template<typename BFTr>
  void IntegrateCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
    dNPos = 0;
    dNNeg = 0;

    double const pi_beta3 = M_PI / (beta * beta * beta) / (8.0 * M_PI * M_PI * M_PI);
    double const bmu = beta * mu;
    double const bmass = beta * mass;
    double const dS0_ = u * ds;
    double const dS0 = std::abs(dS0_);
    (dS0_ >= 0 ? dNPos : dNNeg) = 4 * pi_beta3 * dS0 * BFTr::IntegralP(bmass, bmu, bmass);

    double const dsds = ds * ds;
    if (dsds >= 0) return; // if(timelike)return;
    double const dSz = std::sqrt(dS0 * dS0 - dsds);
    double const vsig = dS0/dSz;
    double const xsig = bmass / std::sqrt(1 - vsig * vsig);
    double const dNSec = pi_beta3 * (
      dSz * ((vsig * vsig + 1) * BFTr::Integral2(xsig, bmu) - bmass * bmass * BFTr::Integral0(xsig, bmu))
      - 2 * dS0 * BFTr::IntegralP(xsig, bmu, bmass));

    dNPos += dNSec;
    dNNeg += dNSec;
  }

  void IntegrateBosonCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
    //mwg_assert(mu<mass);
    return IntegrateCooperFrye<BFTraits<1> >(dNPos, dNNeg, u, ds, beta, mass, mu);
  }
  void IntegrateFermionCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
    return IntegrateCooperFrye<BFTraits<-1> >(dNPos, dNNeg, u, ds, beta, mass, mu);
  }

  double IntegrateCFMomentum(int sign, double xsig, double bmu, double bmass) {
    if (sign == -1)
      return BFTraits<-1>::IntegralP(xsig, bmu, bmass);
    else
      return BFTraits<+1>::IntegralP(xsig, bmu, bmass);
  }

//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
}
}
