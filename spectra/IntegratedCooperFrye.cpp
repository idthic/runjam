#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstdio>
#include <ksh/integrator.hpp>
#include <ksh/phys/Minkowski.hpp>
#include "IntegratedCooperFrye.hpp"

//#define MWG_STD_LAMBDA

namespace ksh=kashiwa;

namespace hydrojet{
namespace spectra{
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

  typedef kashiwa::phys::vector4 vector4;
  static double const sqrtTangentLowerBound=0;
  static double const sqrtTangentUpperBound=std::sqrt(M_PI/2);

  template<int BF>
  struct BFTraits{
#ifdef MWG_STD_LAMBDA
#ifdef DBG20110803
    static double Integral0(double xsig,double bmu){
      return ksh::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound,sqrtTangentUpperBound,[=](double t) -> double{
        double tantt=std::tan(t*t);
        double jacob=2*t*(tantt*tantt+1);
        double x=tantt+xsig;
        double fE0=1/(std::exp(x-bmu)-BF);
        return jacob*fE0;
      });
    }
#endif
    static double Integral2(double xsig,double bmu){
      return ksh::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound,sqrtTangentUpperBound,[=](double t) -> double{
        double tantt=std::tan(t*t);
        double jacob=2*t*(tantt*tantt+1);
        double x=tantt+xsig;
        double fE2=x*x/(std::exp(x-bmu)-BF);
        return jacob*fE2;
      });
    }
    static double IntegralP(double xsig,double bmu,double bmass){
      return ksh::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound,sqrtTangentUpperBound,[=](double t) -> double{
        double tantt=std::tan(t*t);
        double jacob=2*t*(tantt*tantt+1);
        double x=tantt+xsig;
        double fEP=x*std::sqrt(x*x-bmass*bmass)/(std::exp(x-bmu)-BF);
        return jacob*fEP;
      });
    }
#else
    struct LambdaIntegral2{
      double const& xsig;
      double const& bmu;
      LambdaIntegral2(const double& xsig,const double& bmu):xsig(xsig),bmu(bmu){}
      double operator()(double t) const{
        double tantt=std::tan(t*t);
        double jacob=2*t*(tantt*tantt+1);
        double x=tantt+xsig;
        double fE2=x*x/(std::exp(x-bmu)-BF);
        return jacob*fE2;
      }
    };
    static double Integral2(double xsig,double bmu){
      return ksh::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound,sqrtTangentUpperBound,LambdaIntegral2(xsig,bmu));
    }
    struct LambdaIntegralP{
      double const& xsig;
      double const& bmu;
      double const& bmass;
      LambdaIntegralP(const double& xsig,const double& bmu,const double& bmass):xsig(xsig),bmu(bmu),bmass(bmass){}
      double operator()(double t) const{
        double tantt=std::tan(t*t);
        double jacob=2*t*(tantt*tantt+1);
        double x=tantt+xsig;
        double fEP=x*std::sqrt(x*x-bmass*bmass)/(std::exp(x-bmu)-BF);
        return jacob*fEP;
      }
    };
    static double IntegralP(double xsig,double bmu,double bmass){
      return ksh::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound,sqrtTangentUpperBound,LambdaIntegralP(xsig,bmu,bmass));
    }
#endif
    static double Integral0(double xsig,double bmu){
      return -BF*std::log(1-BF*std::exp(bmu-xsig));
    }
  };

  template<typename BFTr>
  void IntegrateCooperFrye(double& dNPos,double& dNNeg,const vector4& u,const vector4& ds,double beta,double mass,double mu){
    dNPos=0;
    dNNeg=0;

    double const pi_beta3=M_PI/(beta*beta*beta)/(8.0*M_PI*M_PI*M_PI);
    double const bmu=beta*mu;
    double const bmass=beta*mass;
    double const dS0_=u*ds;
    double const dS0=std::abs(dS0_);
    (dS0_>=0?dNPos:dNNeg)=4*pi_beta3*dS0*BFTr::IntegralP(bmass,bmu,bmass);

    double const dsds=ds*ds;
    if(dsds>=0)return; // if(timelike)return;
    double const dSz=std::sqrt(dS0*dS0-dsds);
    double const vsig=dS0/dSz;
    double const xsig=bmass/std::sqrt(1-vsig*vsig);
    double const dNSec=pi_beta3*(
      dSz*((vsig*vsig+1)*BFTr::Integral2(xsig,bmu)-bmass*bmass*BFTr::Integral0(xsig,bmu))
      -2*dS0*BFTr::IntegralP(xsig,bmu,bmass)
      );

    dNPos+=dNSec;
    dNNeg+=dNSec;
    //std::fprintf(stderr,"dbg: pi_beta3=%lf bmu=%lf bmass=%lf dS0_=%lf dSz=%lf dN+=%lf dN-=%lf\n",pi_beta3,bmu,bmass,dS0_,dSz,dNPos,dNNeg);
    //std::fprintf(stderr,"dbg: dS0=%.15lf dS0*dS0=%.15lf ds*ds=%.15lf dS0*dS0-ds*ds=%.15lf dSz=%.15lf\n",dS0_,dS0*dS0,ds*ds,dS0*dS0-ds*ds,dSz);
  }

  void IntegrateBosonCooperFrye(double& dNPos,double& dNNeg,const vector4& u,const vector4& ds,double beta,double mass,double mu){
    //mwg_assert(mu<mass);
    return IntegrateCooperFrye<BFTraits<1> >(dNPos,dNNeg,u,ds,beta,mass,mu);
  }
  void IntegrateFermionCooperFrye(double& dNPos,double& dNNeg,const vector4& u,const vector4& ds,double beta,double mass,double mu){
    return IntegrateCooperFrye<BFTraits<-1> >(dNPos,dNNeg,u,ds,beta,mass,mu);
  }

  double IntegrateCFMomentum(int sign,double xsig,double bmu,double bmass){
    if(sign==-1)
      return BFTraits<-1>::IntegralP(xsig,bmu,bmass);
    else
      return BFTraits<+1>::IntegralP(xsig,bmu,bmass);
  }

//-----------------------------------------------------------------------------

  // template<int BF>
  // struct BFTraitsViscous{
  //   struct LambdaForIntegralBulk{
  //     double const& xsig;
  //     double const& bmu;
  //     double const& bmass;
  //     double const& stressMax;
  //     LambdaForIntegralBulk(const double& xsig,const double& bmu,double const& bmass,const double& stressMax)
  //       :xsig(xsig),bmu(bmu),bmass(bmass),stressMax(stressMax)
  //     {}
  //     double operator()(double t) const{
  //       double const tantt=std::tan(t*t);
  //       double const jacob=2*t*(tantt*tantt+1);
  //       double const x=tantt+xsig;

  //       double const p2=x*x-bmass*bmass;
  //       double const f0=1.0/(std::exp(x-bmu)-BF);
  //       double const fu=f0*(1.0+stressmax*p2*(1.0+BF*f0));
  //       return jacob*x*std::sqrt(p2)*fu;
  //     }
  //   };
  //   static double IntegralBulk(double xsig,double bmu,double bmass,double stressMax){
  //     return ksh::IntegrateByGaussLegendre<100>(
  //       sqrtTangentLowerBound,sqrtTangentUpperBound,
  //       LambdaForIntegralBulk(xsig,bmu,bmass,stressMax));
  //   }

  //   struct LambdaForIntegralSurf{
  //     double const& xsig;
  //     double const& vsig;
  //     double const& bmu;
  //     double const& bmass;
  //     double const& stressMax;
  //     LambdaForIntegralSurf(double const& xsig,double const& vsig,const double& bmu,double const& bmass,const double& stressMax)
  //       :xsig(xsig),vsig(vsig),bmu(bmu),bmass(bmass),stressMax(stressMax)
  //     {}
  //     double operator()(double t) const{
  //       double const tantt=std::tan(t*t);
  //       double const jacob=2*t*(tantt*tantt+1);
  //       double const x=tantt+xsig;

  //       double const p2=x*x-bmass*bmass;
  //       double const f0=1.0/(std::exp(x-bmu)-BF);
  //       double const fu=f0*(1.0+stressmax*p2*(1.0+BF*f0));
  //       double const ac=x*vsig-std::sqrt(p2);
  //       return jacob*ac*ac*fu;
  //     }
  //   };
  //   static double IntegralSurf(double xsig,double vsig,double bmu,double bmass,double stressMax){
  //     return ksh::IntegrateByGaussLegendre<100>(
  //       sqrtTangentLowerBound,sqrtTangentUpperBound,
  //       LambdaForIntegralSurf(xsig,vsig,bmu,bmass,stressMax));
  //   }
  // };

  // template<typename BFTr>
  // void IntegrateViscousCF(double& dNPos,double& dNNeg,const vector4& u,const vector4& ds,double beta,double mass,double mu,double stressMax){
  //   dNPos=0;
  //   dNNeg=0;

  //   double const head=(1.0/(8.0*M_PI*M_PI))/(beta*beta*beta); // T^3/4pi^2
  //   double const bmu=beta*mu;
  //   double const bmass=beta*mass;

  //   double const dS0_=u*ds;
  //   double const dS0=std::abs(dS0_);
  //   (dS0_>=0?dNPos:dNNeg)=head*4*dS0*BFTr::IntegralBulk(bmass,bmu,bmass,stressMax);

  //   double const dsds=ds*ds;
  //   if(dsds>=0)return; // if(timelike)return;
  //   double const dSz=std::sqrt(dS0*dS0-dsds);
  //   double const vsig=dS0/dSz;
  //   double const xsig=bmass/std::sqrt(1-vsig*vsig);
  //   double const dNSurf=head*dSz*BFTr::IntegralSurf(xsig,vsig,bmu,bmass,stressMax);
  //   dNPos+=dNSurf;
  //   dNNeg+=dNSurf;
  // }



//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
}
}
