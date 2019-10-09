#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <ksh/integrator.hpp>

#include "uty/Random.hpp"
#include "ParticleSampleViscous.hpp"

namespace idt {
namespace hydro2jam {
namespace {

  // sqrt(pi/2)
  static const double SQRT_TANGENT_ASYMPTOTE = M_SQRT2 / M_2_SQRTPI;

  //---------------------------------------------------------------------------
  // settings

  //...scale transformation 1 [fm^-1]=197.32 [MeV]
  //sctr = 197.32;
  // static const double CONST_INVERSE_MEVFM=197.32;
  // static const double CONST_INVERSE_MEVFM=197.327053 ; // 197.327053 hydrojet/src/physicsbase/Const.h
  static const double CONST_INVERSE_MEVFM = 197.3269718; // 197.3269718(44) idt/rfh/i2/common/def.h
  static const double CONST_MEVFM = 1.0 / CONST_INVERSE_MEVFM;

  // // double facranmax = 1.6;
  // // double facranmax = 1.7;//06/28/2010, lower switching T, larger radial flow
  // // double facranmax = 1.8;//08/27/2019, LHC, larger radial flow
  // static const double momentumRandomMaxFactor = 1.8;
  // static const double maximalMomentum = 6.0e3 * CONST_MEVFM;

	//di = 10;
  static const int CONST_NumberOfBisectionIteration = 20; // di = 20;

  // freezeout を実行する温度の最低を指定します。
  static const double CONST_FreezeoutSkipTemperature = 0.01; // unit: [/fm]

  static const double CONST_ENERGY_MAX_FACTOR = 100.0;

  //! @def CFInterpolationType
  //! - 1 Lagrange polynomials と局所 cubic spline の組合せ
  //! - 2 局所 7次 Langrange polynomials
#define CFInterpolationType 1

  //! @def ParticleSampleViscous20150330_InterpolateUCDF
  //! UCDF を計算する際に補間テーブルを用いる事を指定します。
  //!
  //! (これを指定しても、そう速くならない。
  //! テーブル生成も、補間も時間はかかっていない様だ。
  //! 別の所に律速段階があるという事か。)
#define ParticleSampleViscous20150330_InterpolateUCDF

  //---------------------------------------------------------------------------
  // utilities

  // 反変ベクトルを boost します。
  // 共変ベクトルを逆 boost します。
  void LorentzBoost(double* p, double const* u) {
    double const p0_ = u[0] * p[0] + u[1] * p[1] + u[2] * p[2] + u[3] * p[3];
    double const fac = (p0_ + p[0]) / (1.0 + u[0]);
    p[0] = p0_;
    p[1] += fac*u[1];
    p[2] += fac*u[2];
    p[3] += fac*u[3];
  }

//-----------------------------------------------------------------------------
// test implementation for viscous case

  namespace interpolation {

    // cooper-frye integration without chemical potential
    class CFIntegrator {
    public:
      double bmass;
      double sign;
      double xsig;
      double vsig;
      void calculateIntegrandBulk(double t, double& idealPart, double& viscousPart) const {
        double const tantt = std::tan(t * t);
        double const x = tantt + this->xsig;
        if (x >= CONST_ENERGY_MAX_FACTOR) {
          idealPart = 0.0, viscousPart = 0.0;
          return;
        }
        double const jacob = 2 * t * (tantt * tantt + 1);

        // f_0(x), \delta f(x)
        double const bp2 = x * x - bmass * bmass;
        double const exp_ = std::exp(x);
        double const f0 = 1.0 / (exp_ - sign);
        double const delf = f0 * f0 * exp_ * bp2;

        // bE*bp
        double const pE = x * std::sqrt(bp2);
        idealPart   = jacob * pE * f0;
        viscousPart = jacob * pE * delf;
      }

      void integrateBulk(double& ideal, double& viscous) const {
        typedef kashiwa::integrator_detail::GaussLegendre<100> gl_nodes;
        double const lower  = 0.0;
        double const upper  = SQRT_TANGENT_ASYMPTOTE;
        double const center = 0.5 * (upper + lower);
        double const dxds   = 0.5 * (upper - lower);
        ideal = 0.0, viscous = 0.0;
        for (int i = 0; i < gl_nodes::data_size; i++) {
          double const s = gl_nodes::data[i].t;
          double const w = gl_nodes::data[i].w;

          double i1, v1;
          this->calculateIntegrandBulk(center + dxds * s, i1, v1);
          double i2, v2;
          this->calculateIntegrandBulk(center - dxds * s, i2, v2);

          ideal   += w * (i1 + i2);
          viscous += w * (v1 + v2);
        }
        ideal   *= dxds;
        viscous *= dxds;
      }

      void calculateIntegrandSurface(double t, double& idealPart, double& viscousPart) const {
        double const tantt = std::tan(t * t);
        double const x = tantt + this->xsig;
        if (x >= CONST_ENERGY_MAX_FACTOR) {
          idealPart = 0.0, viscousPart = 0.0;
          return;
        }
        double const jacob = 2 * t * (tantt * tantt + 1);

        // f_0(x), \delta f(x)
        double const bp2 = x * x - bmass * bmass;
        double const exp_ = std::exp(x);
        double const f0 = 1.0 / (exp_ - sign);
        double const delf = f0 * f0 * exp_ * bp2;

        // (x*vsig-bp)^2
        double const bp = std::sqrt(bp2);
        double const prod = x * vsig - bp;
        double const prod2 = prod * prod;
        idealPart   = jacob * prod2 * f0  ;
        viscousPart = jacob * prod2 * delf;
      }

      void integrateSurface(double& ideal, double& viscous) const {
        typedef kashiwa::integrator_detail::GaussLegendre<100> gl_nodes;
        double const lower  = 0.0;
        double const upper  = SQRT_TANGENT_ASYMPTOTE;
        double const center = 0.5 * (upper + lower);
        double const dxds   = 0.5 * (upper - lower);
        ideal = 0.0, viscous = 0.0;
        for (int i = 0; i < gl_nodes::data_size; i++) {
          double const s = gl_nodes::data[i].t;
          double const w = gl_nodes::data[i].w;

          double i1, v1;
          this->calculateIntegrandSurface(center + dxds * s, i1, v1);
          double i2, v2;
          this->calculateIntegrandSurface(center - dxds * s, i2, v2);

          ideal   += w * (i1 + i2);
          viscous += w * (v1 + v2);
        }
        ideal   *= dxds;
        viscous *= dxds;
      }

      // xsig = bmass + tan(t*t)
      void integrateSurfaceByParamTsig(double t, double& ipart, double& vpart) {
        if (t > SQRT_TANGENT_ASYMPTOTE * (1.0 - DBL_EPSILON)) {
          // avoid tan(t^2) jump error
          ipart = 0.0;
          vpart = 0.0;
        } else {
          double const tantt = std::tan(t * t);
          this->xsig = tantt + this->bmass;
          this->vsig = std::sqrt(1.0 - (bmass / xsig) * (bmass / xsig));
          this->integrateSurface(ipart, vpart);
        }
      }

      //-----------------------------------------------------------------------
      // E^2, pE, p^2

      void calculateIntegrandMoment(double t, double* buff) const {
        double const tantt = std::tan(t * t);
        double const x = tantt + this->xsig;
        if (x >= CONST_ENERGY_MAX_FACTOR) {
          std::fill(buff, buff + 6, 0.0);
          return;
        }
        double const jacob = 2 * t * (tantt * tantt + 1);

        // f_0(x), \delta f(x)
        double const bp2 = x * x - bmass * bmass;
        double const exp_ = std::exp(x);
        double const f0 = 1.0 / (exp_ - sign);
        double const delf = f0 * f0 * exp_ * bp2;

        double const _xbp = x * std::sqrt(bp2);
        double const _x2  = x * x;
        buff[0] = jacob * f0   * _xbp;  //  f0   E p
        buff[1] = jacob * delf * _xbp; // delf E p
        buff[2] = jacob * f0   * _x2 ; // f0   E^2
        buff[3] = jacob * delf * _x2 ; // delf E^2
        buff[4] = jacob * f0   * bp2 ; // f0   p^2
        buff[5] = jacob * delf * bp2 ; // delf p^2
      }

      void integrateMoment(double* buff) const {
        typedef kashiwa::integrator_detail::GaussLegendre<100> gl_nodes;
        double const lower  = 0.0;
        double const upper  = SQRT_TANGENT_ASYMPTOTE;
        double const center = 0.5 * (upper + lower);
        double const dxds   = 0.5 * (upper - lower);

        std::fill(buff, buff + 6, 0.0);

        for (int i = 0; i < gl_nodes::data_size; i++) {
          double const s = gl_nodes::data[i].t;
          double const w = gl_nodes::data[i].w;

          double buffL[6], buffR[6];
          this->calculateIntegrandMoment(center + dxds * s, buffL);
          this->calculateIntegrandMoment(center - dxds * s, buffR);
          for (int k = 0; k < 6; k++)
            buff[k] += w * (buffL[k] + buffR[k]);
        }

        for (int k = 0; k < 6; k++) buff[k] *= dxds;
      }

      void integrateMomentByParamTsig(double t, double* buff) {
        if (t > SQRT_TANGENT_ASYMPTOTE * (1.0 - DBL_EPSILON)) {
          // avoid tan(t^2) jump error
          std::fill(buff, buff + 6, 0.0);
        } else {
          double const tantt = std::tan(t * t);
          this->xsig = tantt + this->bmass;
          this->integrateMoment(buff);
        }
      }
    };

    template<int ORDER = 30>
    class LagrangeInterpolation {
    public:
      typedef double interpolation_data[ORDER];

    private:
      double domain_begin;
      double domain_end  ;

    private:
      //! sampling points of Lagrange interpolation
      std::vector<double> nodes;
      void initializeNodes() {
        this->nodes.resize(ORDER, 0.0);

        // Chebychev nodes
        {
          for (int i = 0; i < ORDER; i++) nodes[i] = std::cos(M_PI - (0.5 * M_PI / ORDER) * (2 * i + 1));
          double const scale = 1.0 / nodes[ORDER - 1];
          for (int i = 0; i < ORDER; i++) nodes[i] *= scale;
        }

        // Gauss-Legendre points (Legendre nodes)
        // {
        //   typedef kashiwa::integrator_detail::GaussLegendre<ORDER> gl_interp;
        //   double const scale = 1.0 / gl_interp::data[gl_interp::data_size - 1].t;
        //   for (int i = 0; i < gl_interp::data_size; i++) {
        //     double const s = gl_interp::data[i].t * scale;
        //     nodes[gl_interp::order / 2     + i] = +s;
        //     nodes[gl_interp::data_size - 1 - i] = -s;
        //   }
        // }

        // rescale
        double const center = 0.5 * (domain_end + domain_begin);
        double const half   = 0.5 * (domain_end - domain_begin);
        for (int i = 0; i < ORDER; i++) nodes[i] = center + half * nodes[i];
      }
    public:
      void initialize(double domainBegin = 0.0, double domainEnd = SQRT_TANGENT_ASYMPTOTE) {
        this->domain_begin = domainBegin;
        this->domain_end   = domainEnd;
        this->initializeNodes();
      }

    private:
      void calculateNormalization(double (&normalization)[ORDER]) const {
        std::fill(normalization, normalization + ORDER, 1.0);
        for (int i = 0; i < ORDER; i++) {
          for (int j = i + 1; j < ORDER; j++) {
            double const primary = 1.0 / (nodes[i] - nodes[j]);
            normalization[i] *= +primary;
            normalization[j] *= -primary;
          }
        }
      }
    public:
      template<typename F>
      void initializeData(int count, interpolation_data* data, F const& functor) const {
        double normalization[ORDER];
        this->calculateNormalization(normalization);

        std::vector<double> buff;
        buff.resize(count, 0.0);
        for (int is = 0; is < ORDER; is++) {
          functor(&buff[0], nodes[is]);
          for (int k = 0; k < count; k++)
            data[k][is] = buff[k] * normalization[is];
        }
      }

    public:
      LagrangeInterpolation(double domainBegin = 0.0, double domainEnd = SQRT_TANGENT_ASYMPTOTE) {
        this->initialize(domainBegin, domainEnd);
      }

    private:
      void checkDomain(double tsig) const {
        const double margin = 1e-10;
        if (tsig < domain_begin - margin || domain_end + margin < tsig) {
          std::fprintf(
            stderr, "ParticleSampleViscous.cxx(LagrangeInterpolation): out of domain (tsig=%g not in [%g, %g])\n",
            tsig, domain_begin, domain_end);
          std::exit(EXIT_FAILURE);
        }
      }
    public:
      // 複数の関数の同一点 tsig の上での補間値を得る。
      // 複数の interpolation_data を受け取り複数の double を返す。
      void interpolate(int count, double* result, interpolation_data const* data, double tsig) const {
        this->checkDomain(tsig);

        std::fill(result, result + count, 0.0);
        double accumulated = 1.0;
        for (int i = 0; i < ORDER; i++) {
          double const primary = tsig - nodes[i];
          for (int k = 0; k < count; k++)
            result[k] = result[k] * primary + accumulated * data[k][i];
          accumulated *= primary;
        }
      }
    };

    // Cubic spline interpolation with equi-partitioned intervals
    template<int NPART>
    class CubicSplineInterpolation {
    public:
      struct interpolation_data {
        double y[NPART + 1];
        double h[NPART + 1];
      };

    public:
      template<typename F>
      void initializeData(int count, interpolation_data* data, F const& functor) const {
        double const lower  = 0.0;
        double const upper  = SQRT_TANGENT_ASYMPTOTE;
        double const dxds   = (upper - lower) / double(NPART);

        std::vector<double> buff;
        buff.resize(count, 0.0);
        for (int s = 0; s <= NPART; s++) {
          double const x = lower + dxds * s;
          functor(&buff[0], x);
          for (int k = 0; k < count; k++)
            data[k].y[s] = buff[k];
        }

        for (int k = 0; k < count; k++)
          spline3_determine_natural_hvalues(NPART, data[k].y, data[k].h);
      }
    private:
#if 0
      // global に決めようとすると数値誤差が蓄積して係数が発散する。
      static void spline3_determine_natural_hvalues__global(int npart, double const* _yval, double* _hval) {
        // (a-1) 順方向にする時
        double const*& yval = _yval;
        double*& hval = _hval;
        // // (a-2) 逆方向
        // std::reverse_iterator<double const*> yval(_yval + npart + 1);
        // std::reverse_iterator<double*> hval(_hval + npart + 1);

        // (a) 通常
        {
          hval[0] = 0;
          hval[1] = 0;
          for (int s = 2; s <= npart; s++) {
            double const rhs = yval[s] - 2.0 * yval[s - 1] + yval[s - 2];
            hval[s] = rhs - (4.0 * hval[s - 1] + hval[s - 2]);
          }

          double const sqrt3 = std::sqrt(3);
          double const lam1 = -2.0 + sqrt3;
          double const lam2 = -2.0 - sqrt3;
          hval[1] = -2 * sqrt3 * hval[npart] / (std::pow(lam1, npart) - std::pow(lam2, npart));

          // // (b-1)
          // double lam1s = lam1;
          // double lam2s = lam2;
          // for (int s = 2; s < npart; s++) {
          //   lam1s *= lam1;
          //   lam2s *= lam2;
          //   hval[s] += (hval[1] / (2.0 * sqrt3)) * (lam1s - lam2s);
          // }

          // // (b-2)
          // for (int s = 2; s < npart; s++)
          //   hval[s] += (hval[1] / (2.0 * sqrt3)) * (std::pow(lam1, s) - std::pow(lam2, s));

          // (b-3)
          for (int s = 2; s < npart; s++)
            hval[s] = (yval[s] - 2.0 * yval[s - 1] + yval[s - 2]) - (4.0 * hval[s - 1] + hval[s - 2]);

          hval[npart] = 0.0;
        }

        // // (c) 逐次的に改善 (発散は大して減らない。2, 3桁減る程度)
        // {
        //   double const sqrt3 = std::sqrt(3);
        //   double const lam1 = -2.0 + sqrt3;
        //   double const lam2 = -2.0 - sqrt3;
        //   hval[0] = 0;
        //   hval[1] = 0;
        //   for (int n = 3; n <= npart; n++) {
        //     for (int s = 2; s <= n; s++) {
        //       double const rhs = yval[s] - 2.0 * yval[s - 1] + yval[s - 2];
        //       hval[s] = rhs - (4.0 * hval[s - 1] + hval[s - 2]);
        //     }
        //     hval[1] -= 2 * sqrt3 * hval[n] / (std::pow(lam1, n) - std::pow(lam2, n));
        //   }

        //   for (int s = 2; s <= npart; s++) {
        //     double const rhs = yval[s] - 2.0 * yval[s - 1] + yval[s - 2];
        //     hval[s] = rhs - (4.0 * hval[s - 1] + hval[s - 2]);
        //   }
        //   hval[npart] = 0;
        // }
      }
#endif
      static void spline3_determine_natural_hvalues__local(int npart, double const* _yval, double* _hval) {
        // (a-1) 順方向
        double const*& yval = _yval;
        double*& hval = _hval;
        // // (a-2) 逆方向にしても大して変わらない
        // std::reverse_iterator<double const*> yval(_yval + npart + 1);
        // std::reverse_iterator<double*> hval(_hval + npart + 1);

        // 大局的に spline の方程式を解くと数値的に不安定なので局所的に解く。
        // 具体的には、現在地点 s0 での h 係数 (h[s0+1]) を、N個先で natural spline 条件を課して求める。
        // (当然ながらこれにより、境界上での1,2次微係数の連続性は失われるが、ずれは小さいと期待する。)
        // N は色々試した結果
        const int N = 10;
        // とする事にした。これより多くしても余り精度が改善しない様子だからである。

        hval[0] = 0.0;
        hval[1] = 0.0;
        double const sqrt3 = std::sqrt(3);
        double const lam2  = -2.0 - sqrt3;
        double const lam2N = std::pow(lam2, N);
        double const lam1N = 1.0 / lam2N;
        double const scale = -2.0 * sqrt3 / (lam1N - lam2N);
        for (int s0 = 0; s0 <= npart - N; s0++) {
          for (int s = s0 + 2, sM = s0 + N; s <= sM; s++)
            hval[s] = (yval[s] - 2.0 * yval[s - 1] + yval[s - 2]) - (4.0 * hval[s - 1] + hval[s - 2]);

          double const h1del = scale * hval[s0 + N];
          hval[s0 + 1] += h1del;
          hval[s0 + 2] -= h1del * 4.0;
        }

        // この時点での状態:
        //   s0(last) = npart - N
        //   h[s0(last) + N] = h[npart] を 0 にする様に調整済
        //   h[s0(last) + 2] = h[npart - N + 2] まで確定済

        for (int s = npart - N + 3; s <= npart; s++)
          hval[s] = (yval[s] - 2.0 * yval[s - 1] + yval[s - 2]) - (4.0 * hval[s - 1] + hval[s - 2]);

        hval[npart] = 0;
      }
      static void spline3_determine_natural_hvalues(int npart, double const* yval, double* hval) {
        return spline3_determine_natural_hvalues__local(npart, yval, hval);
      }
      static double spline3_interpolate(double const a0, double const a1, double const* yval, double const* hval) {
        return yval[0] * a0 + yval[1] * a1 - a0 * a1 * (hval[0] + hval[1] + hval[0] * a0 + hval[1] * a1);
      }

    public:
      // 複数の関数の同一点 tsig の上での補間値を得る。
      // 複数の interpolation_data を受け取り複数の double を返す。
      void interpolate(int count, double* result, interpolation_data const* data, double tsig) const {
        double const lower = 0.0;
        double const upper = SQRT_TANGENT_ASYMPTOTE;
        double const dxds  = (upper - lower) / NPART;

        double const s = (NPART / SQRT_TANGENT_ASYMPTOTE) * tsig;
        int    const is = int(s);
        if (is < 0) {
          for (int k = 0; k < count; k++)
            result[k] = data[k].y[0];
        } else if (is >= NPART) {
          for (int k = 0; k < count; k++)
            result[k] = data[k].y[NPART];
        } else {
          double const a1 = s - is;
          double const a0 = 1.0 - a1;
          for (int k = 0; k < count; k++)
            result[k] = spline3_interpolate(a0, a1, data[k].y + is, data[k].h + is);
        }
      }
    };

    class CFInterpolater1 {
      double bmass;
      double integ_bulk_ideal;
      double integ_bulk_viscous;

#ifdef ParticleSampleViscous20150330_InterpolateUCDF
      static const int number_of_interpolated_functions = 8;
#else
      static const int number_of_interpolated_functions = 2;
#endif

#if CFInterpolationType == 1
      typedef LagrangeInterpolation<10>     interpL_type;
      typedef CubicSplineInterpolation<200> interpS_type;
      interpL_type interpL;
      interpS_type interpS;
      interpL_type::interpolation_data dataL[number_of_interpolated_functions];
      interpS_type::interpolation_data dataS[number_of_interpolated_functions];
#elif CFInterpolationType == 2
      static const int NINTERP = 30;
      LagrangeInterpolation<7> interpL[NINTERP];
      LagrangeInterpolation<7>::interpolation_data dataL[NINTERP][number_of_interpolated_functions];
#else
#  error hydrojet: invalid macro value of CFInterpolationType
#endif

    private:
      struct initializeData_lambda1 {
        CFIntegrator& eval;
        initializeData_lambda1(CFIntegrator& eval): eval(eval) {}
        void operator()(double* buff, double tsig) const {
          eval.integrateSurfaceByParamTsig(tsig, buff[0], buff[1]);
#ifdef ParticleSampleViscous20150330_InterpolateUCDF
          eval.integrateMomentByParamTsig(tsig, buff + 2);
#endif
        }
      };
    public:
      void initializeData(IResonanceList* reso, int ir, double beta) {
        CFIntegrator eval;
        this->bmass = beta * reso->mass(ir);
        eval.bmass = this->bmass;
        eval.sign = reso->statisticsSign(ir);
        eval.xsig = this->bmass;
        eval.vsig = 0.0;
        eval.integrateBulk(
          this->integ_bulk_ideal,
          this->integ_bulk_viscous);

        initializeData_lambda1 lambda1(eval);
#if CFInterpolationType == 1
        interpL.initialize(0.0, 0.1);
        interpL.initializeData(number_of_interpolated_functions, dataL, lambda1);
        interpS.initializeData(number_of_interpolated_functions, dataS, lambda1);
#elif CFInterpolationType == 2
        const double width=SQRT_TANGENT_ASYMPTOTE/NINTERP;
        for (int ipart=0; ipart<NINTERP; ipart++) {
          interpL[ipart].initialize(width*ipart, width*(ipart+1));
          interpL[ipart].initializeData(number_of_interpolated_functions, dataL[ipart], lambda1);
        }
#endif
      }

    private:
      void interpolateAt(int offset, int count, double* buff, double tvalue) const {
#if CFInterpolationType == 1
        if (tvalue<0.1)
          interpL.interpolate(count, buff, this->dataL + offset, tvalue);
        else
          interpS.interpolate(count, buff, this->dataS + offset, tvalue);
#elif CFInterpolationType == 2
        const double width = SQRT_TANGENT_ASYMPTOTE/NINTERP;
        int const ipart = std::max(0, std::min(NINTERP - 1, int(tvalue / width)));
        interpL[ipart].interpolate(count, buff, this->dataL[ipart] + offset, tvalue);
#endif
      }

    public:
      double IsotropicPartTotalCDF(double xsig, double dsig0, double dsig_abs, double piMax_2enthalpy) const {
        double ret = 0.0;

        if (dsig0 > 0) {
          double const ideal    =  this->integ_bulk_ideal  ;
          double const viscous  =  this->integ_bulk_viscous;
          ret += 4.0 * dsig0 * (ideal + piMax_2enthalpy * viscous);
        }

        if (std::abs(dsig0) < dsig_abs) {
          double const tsig = std::sqrt(std::atan(xsig - this->bmass));

          double result[2];
          this->interpolateAt(0, 2, result, tsig);

          double const ideal = result[0];
          double const viscous = result[1];
          ret += dsig_abs * (ideal + piMax_2enthalpy * viscous);
        }

        return ret;
      }

    public:
#ifdef ParticleSampleViscous20150330_InterpolateUCDF
      double IsotropicPartUCDF(double xsig, double dsig0, double dsig_abs, double piMax_2enthalpy, double x) const {
        double ret = 0.0;

        double const t = std::sqrt(std::atan(x - this->bmass));

        double integ_pE[2];
        this->interpolateAt(2, 2, integ_pE, t);
        // integ_pE[0] : E p f0
        // integ_pE[1] : E p delf

        if (dsig0 > 0) {
          double const ideal   = integ_pE[0];
          double const viscous = integ_pE[1];
          ret += 4.0 * dsig0 * (ideal + piMax_2enthalpy * viscous);
        }

        if (std::abs(dsig0) < dsig_abs) {
          double ideal, viscous;
          if (x <= xsig) {
            // 積分範囲 [tsig, ∞)
            double const tsig = std::sqrt(std::atan(xsig - this->bmass));

            double integ_surf[2];
            this->interpolateAt(0, 2, integ_surf, tsig);
            ideal   = integ_surf[0];
            viscous = integ_surf[1];
          } else {
            // 積分範囲 [t, ∞)

            double integ_mom[4];
            this->interpolateAt(4, 4, integ_mom, t);
            // integ_sum[0] : E^2 f0
            // integ_sum[1] : E^2 delf
            // integ_sum[2] : p^2 f0
            // integ_sum[3] : p^2 delf

            // (x*vsig-bp)^2
            double const vsig = std::sqrt(1.0 - (bmass / xsig) * (bmass / xsig));
            ideal   = vsig * vsig * integ_mom[0] - 2.0 * vsig * integ_pE[0] + integ_mom[2];
            viscous = vsig * vsig * integ_mom[1] - 2.0 * vsig * integ_pE[1] + integ_mom[3];
          }
          ret += dsig_abs * (ideal + piMax_2enthalpy * viscous);
        }

        return ret;
      }
#endif
    };

    class TotalCDFInterpolater {
      double m_beta;
      std::vector<CFInterpolater1> m_data;
      void initialize(IResonanceList* reso, double beta) {
        this->m_beta = beta;
        this->m_data.resize(reso->numberOfResonances());
        for (int ir = 0, irN = reso->numberOfResonances(); ir < irN; ir++) {
          m_data[ir].initializeData(reso, ir, beta);
        }
      }
    public:
      TotalCDFInterpolater(IResonanceList* reso, double beta) {
        this->initialize(reso, beta);
      }

      double beta() const { return this->m_beta; }

      double IsotropicPartTotalCDF(int ir, double xsig, double dsig0, double dsig_abs, double piMax_2enthalpy) const {
        return this->m_data[ir].IsotropicPartTotalCDF(xsig, dsig0, dsig_abs, piMax_2enthalpy);
      }
#ifdef ParticleSampleViscous20150330_InterpolateUCDF
      double IsotropicPartUCDF(int ir, double xsig, double dsig0, double dsig_abs, double piMax_2enthalpy, double x) const {
        return this->m_data[ir].IsotropicPartUCDF(xsig, dsig0, dsig_abs, piMax_2enthalpy, x);
      }
#endif
    };
  } /* end of namespace interpolation */

  struct ViscousCooperFryeIntegrator {
    double dsig0;
    double dsig_abs;
    double vsig; // in [0.0, 1.0] if the surface is spacelike, -1.0 if the surface is timelike
    double piMax_2enthalpy;

    // for each resonance
    int ireso; // used by m_interp
    double bmass;
    double sign; // +1 for boson, -1 for fermion
    double xsig;

  private:
    double calculateIntegrandBulk(double x) const {
      if (x >= CONST_ENERGY_MAX_FACTOR)return 0;
      //if (!isfinite(exp_)) { std::cerr << "exp: x=" << x << std::endl; std::exit(2); }

      // f(x)
      double const bp2 = x * x - bmass * bmass;
      double const exp_ = std::exp(x);
      double const f0 = 1.0 / (exp_ - sign);
      double const f = f0 * (1 + piMax_2enthalpy * exp_ * f0 * bp2);

      // bE*bp
      double const pE = x * std::sqrt(bp2);
      //if (!isfinite(pE * f))std::cerr << "x=" << x << std::endl;
      return pE * f;
    }

    struct IntegrandBulk {
      ViscousCooperFryeIntegrator const& parent;
      double const& xsig;

      IntegrandBulk(ViscousCooperFryeIntegrator const& parent, double const& xsig): parent(parent), xsig(xsig) {}
      double operator()(double t) const {
        // x, dx/dt
        double const tantt = std::tan(t * t);
        double const jacob = 2 * t * (tantt * tantt + 1);
        double const x = tantt + xsig;
        return jacob * parent.calculateIntegrandBulk(x);
      }
    };

    double calculateIntegrandSurface(double x) const {
      if (x >= CONST_ENERGY_MAX_FACTOR)return 0;
      // f(x)
      double const bp2 = x * x - bmass * bmass;
      double const exp_ = std::exp(x);
      double const f0 = 1.0 / (exp_ - sign);
      double const f = f0 * (1 + piMax_2enthalpy * exp_ * f0 * bp2);

      // (x*vsig-bp)^2
      double const bp = std::sqrt(bp2);
      double const prod = x * vsig - bp;
      return f * prod * prod;
    }

    struct IntegrandSurface {
      ViscousCooperFryeIntegrator const& parent;
      double const& xsig;

      IntegrandSurface(ViscousCooperFryeIntegrator const& parent, double const& xsig): parent(parent), xsig(xsig) {}
      double operator()(double t) const {
        // x, dx/dt
        double const tantt = std::tan(t * t);
        double const jacob = 2 * t * (tantt * tantt + 1);
        double const x = tantt + xsig;
        return jacob * parent.calculateIntegrandSurface(x);
      }
    };

  public:
    double internalIsotropicPartUCDF(double xlow) const {
      double ret = 0.0;
      if (dsig0 > 0.0) {
        IntegrandBulk const integrand(*this, xlow);
        ret += 4.0 * dsig0 * kashiwa::IntegrateByGaussLegendre<100>(0.0, SQRT_TANGENT_ASYMPTOTE, integrand);
      }
      if (vsig >= 0.0) {
        IntegrandSurface const integrand(*this, std::max(xsig, xlow));
        ret += dsig_abs * kashiwa::IntegrateByGaussLegendre<100>(0.0, SQRT_TANGENT_ASYMPTOTE, integrand);
      }

      return ret;
    }
    double internalIsotropicPartTotalCDF() const {
      return internalIsotropicPartUCDF(this->bmass);
    }

  private:
    interpolation::TotalCDFInterpolater const* m_interp;
  public:
    void setCDFInterpolater(interpolation::TotalCDFInterpolater* interp) {
      this->m_interp = interp;
    }

    ViscousCooperFryeIntegrator() {
      this->m_interp = NULL;
    }

    double IsotropicPartTotalCDF() const {
      if (this->m_interp)
        return this->m_interp->IsotropicPartTotalCDF(ireso, xsig, dsig0, dsig_abs, piMax_2enthalpy);
      else
        return this->internalIsotropicPartTotalCDF();
    }

    /// F1 の上側累積関数
    /// @param xlow 条件 bmass<=xlow
    double IsotropicPartUCDF(double xlow) const {
#ifdef ParticleSampleViscous20150330_InterpolateUCDF
      // if (this->m_interp) {
      //   double const interpolated = this->m_interp->IsotropicPartUCDF(ireso, xsig, dsig0, dsig_abs, piMax_2enthalpy, xlow);
      //   double const integrated = this->internalIsotropicPartUCDF(xlow);
      //   //if (std::abs((interpolated - integrated) / (interpolated + integrated + 1e-8)) > 1e-5)
      //   if (std::abs(interpolated - integrated) > 1e-5)
      //     std::cerr << "IsotropicPartUCDF: " << interpolated << " (" << integrated << ")" << std::endl;
      // }

      if (this->m_interp)
        return this->m_interp->IsotropicPartUCDF(ireso, xsig, dsig0, dsig_abs, piMax_2enthalpy, xlow);
      else
#endif
        return this->internalIsotropicPartUCDF(xlow);
    }

    /// F1 の上側累積分布関数の逆関数
    /// @param value 累積確率を指定します。
    /// @return 対応する x = beta*energy を返します。
    double IsotropicPartInverseUCDF(double value) const {
      static const double asymptote = SQRT_TANGENT_ASYMPTOTE;

      // 二分法、目的関数が単調減少である事に注意。
      double alow = 0.0;
      double aupp = SQRT_TANGENT_ASYMPTOTE;
      for (int i = 0; i < CONST_NumberOfBisectionIteration && aupp - alow > 1e-5; i++) {
        double const asec = (alow + aupp) * 0.5;
        double const tsec = std::tan(asec * asec);
        double const xsec = bmass + tsec;
        // assert F(alow) >= value > F(aupp)
        if (IsotropicPartUCDF(xsec) >= value)
          alow = asec;
        else
          aupp = asec;
      }

      double const a = (alow + aupp) * 0.5;
      double const t = std::tan(a * a);
      double const x = bmass + t;
      return x;
    }
  };

  struct SurfaceParticleSampler: ViscousCooperFryeIntegrator {
    typedef ViscousCooperFryeIntegrator base;
    using base::dsig0;
    using base::dsig_abs;
    using base::vsig;
    using base::piMax_2enthalpy;

    // change for each resonance
    using base::bmass;
    using base::sign;
    using base::xsig;

    // surface information
    HypersurfaceElementC0Lrf const* surface;
    double temperature;
    double beta;
    double dsig[4];
    double _1_2enthalpy;
    double const* stress;

  private: // preferences
    double m_overSamplingFactor;
  public:
    void setOverSamplingFactor(double value) {
      this->m_overSamplingFactor = value;
    }
    // double getOverSamplingFactor() const {
    //   return this->m_overSamplingFactor;
    // }

  private:
    bool m_turnsOffViscousEffect;
  public:
    void setTurnsOffViscousEffect(bool value) {
      if ((this->m_turnsOffViscousEffect = value)) {
        this->_1_2enthalpy    = 0.0;
        this->piMax_2enthalpy = 0.0;
      } else {
        // piMax = \pi_max = max{固有値 of \pi}
        double const piMax = surface->m_stressMax;
        this->_1_2enthalpy = 0.5 / (surface->m_energy + surface->m_pressure);
        this->piMax_2enthalpy = piMax * _1_2enthalpy;
      }
    }

  public:
    SurfaceParticleSampler(HypersurfaceElementC0Lrf const* surface) {
      this->m_overSamplingFactor=1.0;
      this->initialize(surface);

      this->m_dominating = NULL;
      this->totalIntegral = 0.0;
    }

  public:
    void initialize(HypersurfaceElementC0Lrf const* surface) {
      this->surface = surface;
      this->temperature = surface->temperature();
      this->stress = surface->m_stress;

      double const tau = surface->position(0);
      //double const tau = 1.0;
      this->beta    = 1.0 / temperature;
      this->dsig[0] = surface->surfaceElement(0) * tau;
      this->dsig[1] = surface->surfaceElement(1);
      this->dsig[2] = surface->surfaceElement(2) * tau;
      this->dsig[3] = surface->surfaceElement(3) * tau;

      // 面素 dσ を局所静止系での表示に持ってくる。
      // 共変ベクトル dsig を逆ブーストする (= 反変ベクトルの順ブーストと同じ操作)。
      LorentzBoost(this->dsig, surface->m_velocity);

      this->dsig0 = this->dsig[0];
      this->dsig_abs = std::sqrt(dsig[1] * dsig[1] + dsig[2] * dsig[2] + dsig[3] * dsig[3]);
      double const dsig0_abs = std::abs(dsig0);
      if (dsig0_abs < dsig_abs) {
        // spacelike surface
        vsig = dsig0_abs / dsig_abs;
      } else {
        // timelike surface
        vsig = -1.0;
      }

      this->sign  =  0;
      this->bmass = -1;
      this->setTurnsOffViscousEffect(false);
    }


  private:
    // stress =
    //   [ stress[0] stress[1] stress[2] ]
    //   [ stress[1] stress[3] stress[4] ]
    //   [ stress[2] stress[4] stress[5] ]
    double UnisotropicPartProbability(double* bp) {
      double const pppi
        = stress[0] * bp[1] * bp[1]
        + stress[3] * bp[2] * bp[2]
        + stress[5] * bp[3] * bp[3]
        + 2.0*(
          stress[1] * bp[1] * bp[2]
          +stress[2] * bp[1] * bp[3]
          +stress[4] * bp[2] * bp[3]);

      double const x = bp[0];
      double const bp2 = x * x - bmass * bmass;

      // f0*exp(x) = exp(x)/(exp(x)-sgn) = 1/(1-sgn*exp(-x))
      double const f0exp = 1.0 / (1.0 - sign * std::exp(-x));

      double const prob1 = (1.0 + _1_2enthalpy * pppi * f0exp) / (1 + piMax_2enthalpy * bp2 * f0exp);
      return std::max(prob1, 0.0);
    }

  private:
    double totalIsotropicPartCDF;

    /// 局所静止系での beta*momentum を生成します。
    /// @return 棄却された場合に false を返します。
    bool generateMomentum(double* bpout) {
      // (1) 運動量分布(等方部分)

      // (0, 1] の範囲で一様乱数を生成する為に 1.0-getRand() を用いる。
      // 理由: ucdf==0 は E=p=∞ を意味するので除外したい。
      //   また、ucdf==1 は p=0, E=m を意味するので (位相体積0だが) 一応可能である。
      //   従って、ucdf in (0, totalIsotropicPartCDF] となる。従って (0, 1] の乱数を用いたい。
      //   一方で、getRand() は [0, 1) の一様乱数を生成するので 1.0-getRand() とすれば良い。
      // To generate uniform random numbers in the interval (0, 1], we use 1.0-getRand() instead of getRand().
      // Because:
      //   the case ucdf==0 means E, p is infinity, so it should be excluded.
      //   the case ucdf==1 means E=m, and p=0, and this is physically acceptible even if its phasespace volume is zero.
      //   Therefore, the ucdf should be generated in the interval (0, totalIsotropicPartCDF],
      //   and we want to use the uniform random numbers in (0, 1]. Since getrand() is the random number generator
      //   in the interval [0, 1), we can use 1.0-getRand() to generate the desired numbers.
      double const x = IsotropicPartInverseUCDF((1.0 - Random::getRand()) * totalIsotropicPartCDF);
      double const bp = std::sqrt(x * x - bmass * bmass);
      if (!std::isfinite(bp)) {
        std::cerr << "ParticleSampleViscous.cxx(generateMomentum): beta*p nan/inf" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // (2) [pdS]^+ の角度分布を生成(等方ではない)
      {
        // weight ∝ ∫[c = -min(1, A) ～ 1] dc (a+b*c); c = cost
        double cost;
        double const a = x * dsig[0];
        double const b = bp * dsig_abs;
        if (b < -a) {
          // assert(0);
          // a < -b の時は積分値が0になるのでこの函数は呼ばれない筈。
          return false;
        } else if (a <= b) {
          // chi = (a+b*c)^2
          double const chi = (a + b) * (a + b) * Random::getRand();
          cost = (std::sqrt(chi) - a) / b;
        } else {
          // chi = (a+b*c)^2
          // psi = (chi-a*a)/b
          double const psi = b + (4.0 * Random::getRand() - 2.0) * a;
          double const chi = a * a + b * psi;
          cost = psi / (a + std::sqrt(chi));
        }

        double const sint = std::sqrt(1.0 - cost * cost);
        double const phd = 2 * M_PI * Random::getRand();

        // dsig∥z な座標系での運動量
        bpout[0] = x;
        bpout[1] = bp * sint * std::cos(phd);
        bpout[2] = bp * sint * std::sin(phd);
        bpout[3] = bp * cost;
      }

      // dsig の方向に回転
      double const dsigt = std::sqrt(dsig[1] * dsig[1] + dsig[2] * dsig[2]);
      {
        double const cos_the_sig = dsig[3] / dsig_abs;
        double const sin_the_sig = dsigt / dsig_abs;
        double const vz = bpout[3];
        double const vx = bpout[1];
        bpout[3] = cos_the_sig * vz - sin_the_sig * vx;
        bpout[1] = sin_the_sig * vz + cos_the_sig * vx;
      }
      {
        double const cos_phi_sig = dsig[1] / dsigt;
        double const sin_phi_sig = dsig[2] / dsigt;
        double const vx = bpout[1];
        double const vy = bpout[2];
        bpout[1] = cos_phi_sig * vx - sin_phi_sig * vy;
        bpout[2] = sin_phi_sig * vx + cos_phi_sig * vy;
      }

      // (3) 非等方部分の確率的棄却
      //   棄却した時は、再生成は行わない。
      return Random::getRand() < UnisotropicPartProbability(bpout);
    }

    bool generateParticleSample(Particle* particle) {
      // 局所静止系での運動量
      double p[4];
      if (!generateMomentum(p)) return false;
      p[0] *= temperature;
      p[1] *= temperature;
      p[2] *= temperature;
      p[3] *= temperature;

      // チルダ座標での運動量
      LorentzBoost(p, surface->m_velocity);

      // 曲線座標での位置
      double x[4];
      for (int i = 0; i < 4; i++) {
        x[i] = surface->position(i);
        if (surface->m_dx[i] > 0.0)
          x[i] += (Random::getRand() - 0.5) * surface->m_dx[i];
      }

      double const cosh_ = std::cosh(x[1]);
      double const sinh_ = std::sinh(x[1]);

      // 実験座標(t, x, y, z)での運動量(tau-eta)
      double const laboratory_momentum[4] = {
        p[0] * cosh_ + p[1] * sinh_,
        p[2],
        p[3],
        p[1] * cosh_ + p[0] * sinh_
      };

      // 実験座標(t, x, y, z)での位置(tau-eta)
      double const laboratory_coord[4] = {
        x[0] * cosh_,
        x[2],
        x[3],
        x[0] * sinh_
      };

      // 粒子の追加
      const double hbarC = 0.197327053;
      particle->e  = laboratory_momentum[0] * hbarC;
      particle->px = laboratory_momentum[1] * hbarC;
      particle->py = laboratory_momentum[2] * hbarC;
      particle->pz = laboratory_momentum[3] * hbarC;
      particle->t  = laboratory_coord[0];
      particle->x  = laboratory_coord[1];
      particle->y  = laboratory_coord[2];
      particle->z  = laboratory_coord[3];
      return true;
    }

    void dbg_checkIsotropicPartUCDF() {
      std::FILE* f = std::fopen("dbg_IsotropicPartUCDF.txt", "w");
      const int N = 1000;
      std::fprintf(f, "# x=bete*E UCDF(in) UCDF(out) p\n");
      for (int i = 0; i < N; i++) {
        double const cdf = (1.0 / N) * (i + 0.5);
        double const x = IsotropicPartInverseUCDF(totalIsotropicPartCDF * cdf);
        std::fprintf(
          f, "%g %g %g %g\n", x, cdf, IsotropicPartUCDF(x) / totalIsotropicPartCDF,
          std::sqrt(x * x - bmass * bmass) * temperature
        );
      }
      std::fclose(f);
    }

  private:
    double totalIntegral;
    interpolation::TotalCDFInterpolater* m_dominating;
  public:
    void setDominatingCDFInterpolater(interpolation::TotalCDFInterpolater* dominatingCDF) {
      this->m_dominating = dominatingCDF;
    }

  public:
    void SampleResonance(std::vector<Particle*>& plist, IResonanceList const* reso, int ireso) {
      double const _sign  = reso->statisticsSign(ireso);
      double const _bmass = beta * reso->mass(ireso);
      //double const _bmu   = beta * surface->chemicalPotential(ireso); //□mu not supported

      // update ireso, sign, bmass, xsig, totalIntegral
      //   反粒子の場合など、前の共鳴と同じパラメータになっている場合がある。
      //   その場合には totalIsotropicPartCDF の計算を省略できる (151 → 90 に減る)。
      //   □化学ポテンシャルに対応する時は、bmu の一致も確かめる。
      //     但し化学ポテンシャルが入っている場合は基本的には bmu の符号が異なるので計算は省略できない。
      base::ireso = ireso;
      if (_sign != this->sign || _bmass != this->bmass) {
        // update parameters and totalIsotropicPartCDF

        this->sign  = _sign;
        this->bmass = _bmass;
        //this->bmu   = _bmu; //□mu not supported
        this->xsig = bmass / std::sqrt(1.0 - vsig * vsig);

        if (m_dominating) {
          // 計算コストの小さい優関数で totalIntegral を評価し、n >= 1 の時に本当の積分を実行する。
          double const betaRatio = beta / m_dominating->beta();
          this->totalIntegral
            = betaRatio * betaRatio * betaRatio
            * m_dominating->IsotropicPartTotalCDF(
              ireso,
              xsig / betaRatio,
              dsig0, dsig_abs,
              betaRatio * betaRatio * piMax_2enthalpy);
        } else
          this->totalIntegral = base::IsotropicPartTotalCDF();
      }

      if (std::isnan(this->totalIntegral)) {
        std::cerr << "ParticleSampleViscous.cxx: IsotropicPartUCDF([bmass,infty]) is nan" << std::endl;
        std::fprintf(stderr, "IsotropicPartCDF([bmass,infty])=%g bmass=%g xsig=%g vsig=%g\n", totalIsotropicPartCDF, bmass, xsig, vsig);
        std::fflush(stderr);
        std::exit(1);
      }

      // nlambda = <n> = Poisson 分布期待値
      double const nlambda1 = totalIntegral * (1.0 / (8 * M_PI * M_PI)) * (temperature * temperature * temperature);
      double const nlambda = nlambda1 * reso->numberOfDegrees(ireso) * this->m_overSamplingFactor;
      int const n = Random::getRandPoisson(nlambda);
      if (n == 0) return;

      // if (ireso == 17 && temperature > 0.75) { // pion
      //   std::cerr << "start dbg_checkIsotropicPartUCDF" << std::endl;
      //   this->dbg_checkIsotropicPartUCDF();
      //   std::exit(EXIT_FAILURE);
      // }
      // if (ireso == 18) {
      //   std::fprintf(stderr, "jam30x:npos=%g\n", nlambda1);
      // }

      double prob = 1.0;
      if (m_dominating) {
        this->totalIsotropicPartCDF = base::IsotropicPartTotalCDF();
        prob = totalIsotropicPartCDF / totalIntegral;
        if (prob > 1.0) {
          std::cerr << "ParticleSmapleViscous.cxx (SampleResonance): BUG the dominating function actually doesn't dominate (prob=" << prob << ")." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        this->totalIsotropicPartCDF = totalIntegral;
      }

      Particle particle(ireso);
      for (int i = 0; i < n; i++) {
        if (m_dominating && Random::getRand() >= prob) continue;

        if (generateParticleSample(&particle)) {
          particle.e = -1.0; // isospin 等を決めてから mass Pe を決める。
          plist.push_back(new Particle(particle));
        }
      }
    }
  };
}

//=============================================================================
// Implementations

void SampleParticlesC0lrf(
  std::vector<Particle*>& plist,
  HypersurfaceElementC0Lrf const& surface,
  IResonanceList const* rlist,
  double overSamplingFactor,
  bool turnsOffViscousEffect
) {
  SurfaceParticleSampler sampler(&surface);
  sampler.setOverSamplingFactor(overSamplingFactor);
  sampler.setTurnsOffViscousEffect(turnsOffViscousEffect);

  int const iresoN = rlist->numberOfResonances();
  for (int ireso = 0; ireso < iresoN; ireso++)
    sampler.SampleResonance(plist, rlist, ireso);

  // mass の単位は [/fm] で OK
  //   mass の単位は 197.32 [MeV fm] で割ってある
  //   mass [/fm] = M [MeV] / 197.32 [MeV fm] という事。
  // {
  //   for (int ireso = 0; ireso < iresoN; ireso++)
  //     std::cerr << "mass[" << ireso << "]=" << rlist->mass(ireso) * CONST_INVERSE_MEVFM << std::endl;
  //   std::exit(2);
  // }
}

//=============================================================================
//  class OversampledParticleSampleBase
//-----------------------------------------------------------------------------

void OversampledParticleSampleBase::update() {
  // 一括生成済の時
  if (this->pcache.size() > 0) {
    if (++this->indexOfCachedEvents < this->pcache.size())
      return;

    this->pcache.clear();
    this->indexOfCachedEvents = -1;
  }

  // 一括生成要求がある時
  if (this->numberOfExpectedEvents > 0) {
    double const& ncache = this->numberOfExpectedEvents;
    this->updateWithOverSampling(this->m_overSamplingFactor * ncache);
    this->pcache.resize(ncache, std::vector<Particle*>());
    for (std::vector<Particle*>::const_iterator i = this->base::plist.begin(); i != this->base::plist.end(); ++i)
      this->pcache[int(Random::getRand() * ncache)].push_back(*i);

    this->numberOfExpectedEvents = 0;
    this->indexOfCachedEvents = 0;
    return;
  }

  this->updateWithOverSampling(this->m_overSamplingFactor);
}

//=============================================================================
//  sampling for <rfh hypersurface.txt c0lrf>
//-----------------------------------------------------------------------------

ParticleSampleViscous::ParticleSampleViscous(IResonanceList* rlist, std::string const& fname_hypersurface)
  :rlist(rlist),
   fname_hypersurface(fname_hypersurface),
   m_turnsOffViscousEffect(false)
{
  this->m_switchingTemperature = 155.0;
}

void ParticleSampleViscous::updateWithOverSampling(double overSamplingFactor) {
  this->base::clearParticleList();

  std::ifstream ifs(this->fname_hypersurface.c_str());
  if (!ifs) {
    std::cerr << "(ParticleSampleViscous::update): failed to open the hypersurface file(" << fname_hypersurface << ")." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string line;

  double switchingTemperature = this->m_switchingTemperature / CONST_INVERSE_MEVFM; // [/fm]

  // 1行目
  if (std::getline(ifs, line)) {
    std::istringstream s(line);

    // $1-$2 magic
    std::string type, coords;
    s >> type >> coords;
    if (type != "c0lrf" || (coords != "taueta-tilde" && coords != "taueta")) {
      std::cerr << "(ParticleSampleViscous::update):" << fname_hypersurface << ": format not supported." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // $3... Tsw=155/197.32 などを読み取る
    std::string arg;
    while (s >> arg) {
      if (arg.compare(0, 4, "Tsw=", 4) == 0) {
        switchingTemperature = std::atof(arg.c_str() + 4);
      } else {
        std::cerr << "(ParticleSampleViscous::update):" << fname_hypersurface << ":1: unknown argument (arg=" << arg << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  } else {
    std::cerr << "(ParticleSampleViscous::update):" << fname_hypersurface << ": the file is empty." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  double const switchingBeta = 1.0 / switchingTemperature; // [/fm]
  interpolation::TotalCDFInterpolater interp(rlist, switchingBeta);

  bool completed=false;
  while (std::getline(ifs, line)) {
    HypersurfaceElementC0Lrf surface;
    {
      std::istringstream s(line);
      // $1-$13 surface element
      int normalDirection;
      s >> normalDirection;
      if (normalDirection == -1) {
        completed = true;
        break;
      }
      for (int i = 0; i < 4; i++)
        s >> surface.m_ds[i];
      for (int i = 0; i < 4; i++)
        s >> surface.m_pos[i];
      for (int i = 0; i < 4; i++)
        s >> surface.m_dx[i];

      // $14-$18 ideal part
      s >> surface.m_temperature;
      for (int i = 0; i < 4; i++)
        s >> surface.m_velocity[i];

      // $19-$27 stress part
      s >> surface.m_energy;
      s >> surface.m_pressure;
      s >> surface.m_stressMax;
      for (int i = 0; i < 6; i++)
        s >> surface.m_stress[i];
    }

    if (surface.m_temperature < CONST_FreezeoutSkipTemperature) {
      // skip
      continue;
    }

    {
      SurfaceParticleSampler sampler(&surface);
      sampler.setOverSamplingFactor(overSamplingFactor);
      sampler.setTurnsOffViscousEffect(this->m_turnsOffViscousEffect);
      if (std::abs(surface.m_temperature - switchingTemperature) < 1e-5)
        sampler.setCDFInterpolater(&interp);
      else if (surface.m_temperature < switchingTemperature)
        sampler.setDominatingCDFInterpolater(&interp);

      int const iresoN = rlist->numberOfResonances();
      for (int ireso = 0; ireso < iresoN; ireso++)
        sampler.SampleResonance(plist, rlist, ireso);
    }
  }

  if (!completed) {
    std::cerr << "(ParticleSampleViscous::update):" << fname_hypersurface << ": unexpected end of the file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "(ParticleSampleViscous::update):" << fname_hypersurface << ": " << base::plist.size() << " particles were generated," << std::endl;
}

//=============================================================================
//  sampling for <hydrojet freezeout.dat/position.dat>
//-----------------------------------------------------------------------------

ParticleSampleFromHydrojet::ParticleSampleFromHydrojet(
  IResonanceList* rlist,
  std::string const& fname_freezeout,
  std::string const& fname_position
):
  rlist(rlist),
  fname_freezeout(fname_freezeout),
  fname_position(fname_position),
  dx(0.3), dy(0.3), dh(0.3), dtau(0.3)
{
  this->m_switchingTemperature=155.0;
}
ParticleSampleFromHydrojet::ParticleSampleFromHydrojet(
  IResonanceList* rlist,
  std::string const& dname_hydro)
  :rlist(rlist),
   fname_freezeout(dname_hydro + "/freezeout.dat"),
   fname_position(dname_hydro + "/position.dat"),
   dx(0.3), dy(0.3), dh(0.3), dtau(0.3)
{
  this->m_switchingTemperature = 155.0;
}

bool ParticleSampleFromHydrojet::readHypersurfaceElement(HypersurfaceElementC0Lrf& surface, std::ifstream& ifsf, std::ifstream& ifsp) const {
  // isbulk
  int isbulk;
  if (!(ifsf >> isbulk)) {
    std::cerr << "(ParticleSampleFromHydrojet::update):" << fname_freezeout << ": unexpected end of the file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (isbulk < 0) {
    // end
    if (!(ifsp >> isbulk)) {
      std::cerr << "(ParticleSampleFromHydrojet::update):" << fname_position << ": unexpected end of the file." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (isbulk >= 0) {
      std::cerr << "(ParticleSampleFromHydrojet::update):" << fname_position << ": inconsistent number of the entries." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    return false;
  }

  // m_pos
  ifsp >> surface.m_pos[0]; // tau
  ifsp >> surface.m_pos[2]; // x
  ifsp >> surface.m_pos[3]; // y
  ifsp >> surface.m_pos[1]; // eta

  // m_ds (ifsf $2-4)
  {
    double ds;
    ifsf >> ds >> surface.m_ds[2] >> surface.m_ds[3];
    if (isbulk) {
      surface.m_ds[0] = ds;
      surface.m_ds[1] = 0.0;
    } else {
      surface.m_ds[0] = 0.0;
      surface.m_ds[1] = -ds;
    }

    // 基底変換: tilde 座標 → grid 座標
    double const inverseTau = 1.0 / surface.m_pos[0];
    surface.m_ds[0] *= inverseTau;
    surface.m_ds[2] *= inverseTau;
    surface.m_ds[3] *= inverseTau;
  }

  // m_dx
  surface.m_dx[0] = this->dtau;
  surface.m_dx[1] = this->dh;
  surface.m_dx[2] = this->dx;
  surface.m_dx[3] = this->dy;
  if (isbulk)
    surface.m_dx[0] = 0.0;
  else if (surface.m_ds[1])
    surface.m_dx[1] = 0.0;
  else if (surface.m_ds[2])
    surface.m_dx[2] = 0.0;
  else
    surface.m_dx[3] = 0.0;

  // (ifsf $5-7)
  double eta, mu;
  ifsf >> eta >> surface.m_temperature >> mu;
  if (mu != 0.0) {
    std::cout << "ParticleSampleViscous.cxx (ParticleSampleFromHydrojet): a finite baryon chemical potenatial is not supported." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // m_velocity (ifsf $8-10)
  {
    double vxp, vyp, yvp;
    ifsf >> vxp >> vyp >> yvp;
    double const vzp = std::tanh(yvp);
    double const u0 = 1.0 / std::sqrt(1.0 - (vxp * vxp + vyp * vyp + vzp * vzp));
    double const u_lab[4] = {u0, u0 * vzp, u0 * vxp, u0 * vyp};

    double const cosh_ = std::cosh(eta);
    double const sinh_ = std::sinh(eta);
    surface.m_velocity[0] =  u_lab[0] * cosh_ - u_lab[1] * sinh_;
    surface.m_velocity[1] = -u_lab[0] * sinh_ + u_lab[1] * cosh_;
    surface.m_velocity[2] = u_lab[2];
    surface.m_velocity[3] = u_lab[3];
  }

  int iw = 0;
  ifsf >> iw;
  if (iw != 8 && iw != 4) {
    std::cerr << "ParticleSampleViscous.cxx (ParticleSampleFromHydrojet): not supported value: iw=" << iw << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return true;
}

void ParticleSampleFromHydrojet::updateWithOverSampling(double overSamplingFactor) {
  this->base::clearParticleList();

  std::ifstream ifsf(this->fname_freezeout.c_str());
  if (!ifsf) {
    std::cerr << "(ParticleSampleFromHydrojet::update): failed to open the hypersurface element data (" << fname_freezeout << ")." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::ifstream ifsp(this->fname_position.c_str());
  if (!ifsp) {
    std::cerr << "(ParticleSampleFromHydrojet::update): failed to open the hypersurface position data (" << fname_position << ")." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  double switchingTemperature = this->m_switchingTemperature / CONST_INVERSE_MEVFM; // [/fm]

  double const switchingBeta = 1.0 / switchingTemperature; // [/fm]
  std::cout << "ParticleSampleViscous.cxx (interpolation): initializing table..." << std::endl;
  interpolation::TotalCDFInterpolater interp(rlist, switchingBeta);
  std::cout << "ParticleSampleViscous.cxx (interpolation): done" << std::endl;

  HypersurfaceElementC0Lrf surface;

  // dummy values
  surface.m_energy = 1.0;
  surface.m_pressure = 1.0;
  surface.m_stressMax = 0.0;
  for (int i = 0; i < 6; i++)
    surface.m_stress[i] = 0.0;

  int c1 = 0, c2 = 0;
  while (this->readHypersurfaceElement(surface, ifsf, ifsp)) {
    if (surface.m_temperature < CONST_FreezeoutSkipTemperature) continue;

    SurfaceParticleSampler sampler(&surface);
    sampler.setOverSamplingFactor(overSamplingFactor);
    sampler.setTurnsOffViscousEffect(true);

    if (std::abs(surface.m_temperature-switchingTemperature) < 1e-5) {
      sampler.setCDFInterpolater(&interp);
      c1++;
    } else if (surface.m_temperature<switchingTemperature) {
      sampler.setDominatingCDFInterpolater(&interp);
      c1++;
    }
    c2++;

    int const iresoN = rlist->numberOfResonances();
    for (int ireso = 0; ireso < iresoN; ireso++)
      sampler.SampleResonance(plist, rlist, ireso);
  }

  std::cout
    << "ParticleSampleFromHydrojet: done:\n"
    << "  input = " << fname_freezeout << " " << fname_position << "\n"
    << "  interpolation = " << c1 << "/" << c2 << " @ T_sw <= " << this->m_switchingTemperature << "MeV\n"
    << "  number_of_particles = " << base::plist.size() << std::endl;
}

//=============================================================================
// Check codes

namespace {
  class DummyResonanceList: public IResonanceList {
  public:
    virtual int numberOfResonances() const { return 2; }
    virtual double mass(int ir) const { return 135 / CONST_INVERSE_MEVFM; }
    virtual double statisticsSign(int ir) const { return ir == 0 ? +1.0 : -1.0; }

    virtual double chemicalPotential(int ir) const { return 0.0; }
    virtual int numberOfDegrees(int ir) const { return 1; }
  };
}

/// @fn int checkViscousCooperFryeInterpolated(bool debug);
/// @param debug if this argument is false,
///   the program exits with EXIT_FAILURE on the failure of the test.
///   if this argument is true, the function simply returns EXIT_FAILURE
///   on the failure of the test.
///   Otherwise, the function returns EXIT_SUCCESS.
/// @return if the test succeeds, returns EXIT_SUCCESS.
///   Otherwise, returns EXIT_FAILURE or the program exits depending
///   on the argument `debug'.
int checkViscousCooperFryeInterpolated(bool debug) {
  DummyResonanceList reso;

  ViscousCooperFryeIntegrator a;
  double const temperature = 155 / CONST_INVERSE_MEVFM;
  double const beta = 1.0 / temperature;
  a.dsig0           = 0.3;
  a.dsig_abs        = 1.0;
  a.piMax_2enthalpy = 0.5;
  a.vsig            = std::abs(a.dsig0) < a.dsig_abs ? std::abs(a.dsig0) / a.dsig_abs : -1.0;

  a.bmass = beta * reso.mass(0);
  a.sign  = reso.statisticsSign(0);

  interpolation::TotalCDFInterpolater interp(&reso, beta);

  static const int N = 1000;
  for (int it = 0; it < N; it++) {
    double const t = (it + 0.0) / N * SQRT_TANGENT_ASYMPTOTE;
    double const xsig = std::tan(t * t) + a.bmass;
    a.xsig = xsig;
    a.vsig = std::sqrt(1.0 - (a.bmass / a.xsig) * (a.bmass / a.xsig));
    double const oldValue = a.IsotropicPartTotalCDF();
    double const newValue = interp.IsotropicPartTotalCDF(0, a.xsig, a.dsig0, a.dsig_abs, a.piMax_2enthalpy);

    if (debug) {
      std::printf("%lg %lg %lg %lg %lg\n", t, xsig, oldValue, newValue, newValue - oldValue);
    } else if (!(std::abs(newValue - oldValue) < 1e-8 * std::max(oldValue, 1.0))) {
      std::fprintf(
        stderr, "check_interpolation: too large error! (correct=%lg interpolated=%lg delta=%lg)\n",
        oldValue, newValue, newValue - oldValue);
      std::exit(EXIT_FAILURE);
    }
  }

  return EXIT_SUCCESS;
}

}
}
