/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2025, Koichi Murase <myoga.murase at gmail.com>

   SPDX-License-Identifier: GPL-2.0-or-later

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA  */

#include "config.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <utility>
#include <stdexcept>
#include <ksh/integrator.hpp>

#include "args.hpp"
#include "ResonanceList.hpp"

/// @def runjam_cmd_resolist_eos_CheckEnergyFundamentalRelation
/// Perform the integration $e = \int_0^T dT T d^p/dT^2$ and compare the result
/// with the one calculated by $e = sT - p$.  This is output in $6 in the
/// output eos file.  If this macro is not defined, the comparison is skipped
/// and $6 becomes nan.
#undef runjam_cmd_resolist_eos_CheckEnergyFundamentalRelation

namespace idt::runjam {
namespace {

  static const double hbarc_GeVfm = 0.197327053; // hydrojet
  static const double SQRT_TANGENT_ASYMPTOTE = M_SQRT2 / M_2_SQRTPI;
  static const double CONST_ENERGY_MAX_FACTOR = 100.0;

  //! This function is used to determine the temperature corresponding to the
  //! specified energy density $e$.
  //! @param[in] tl,tu [tl, tu)` specifies the search range of the temperature.
  //! @param[in] guess specifies the initial guess.  The initial guess needs to
  //!   be a positive number that is sufficiently close to the correct value.
  //!   If a negative number is specified to `guess`, we automatically
  //!   determine the initial guess based on binary search.
  //! @param[in] func is a function void(double (&buff)[2], double temperature)
  //!   to obtain the enegy density $e$ and its temperature derivative $de/dT$.
  //!   It stores the calculated energy density to `buff[0]`, and the
  //!   derivative of the energy density to `buff[1]`.
  template<typename Derivative>
  double log_log_newton_search(const char* tag, double tl, double tu, double guess, double e, Derivative func) {
    double buff[2];

    int start_newton = 5;
    double tm;
    {
      // check boundary
      func(buff, tl);
      double const el = buff[0];
      if (e < el) {
        tu = tl;
        tl = 0.0;
        start_newton = 999; // Do not use the Newton method
      }

      func(buff, tu);
      double const eu = buff[0];
      if (eu <= e) {
        throw std::range_error(
          idt::util::strprintf(
            "error(%s): energy %g [GeV/fm^3]: out of range [%g, %g)",
            tag,
            e * hbarc_GeVfm,
            el * hbarc_GeVfm,
            eu * hbarc_GeVfm)
        );
      }

      if (guess > 0.0) {
        tm = guess;
        start_newton = 0;
      } else {
        // tm: initial guess

        // Initial guess by the linear interpolation.  This isn't bad but it
        // produces a significantly small estimate ($\times 10^{-18}$ factor
        // difference) as the initial value.
        //tm = tl + (tu - tl) * (e - el) / (eu - el);

        // Initial guess by the linear interpolation in unit of MeV.  This
        // still produces a significantly small initial value with $\times
        // 10^{-5}$.
        //tm = tl + (tu - tl) * (std::pow(e, 0.25) - std::pow(el, 0.25)) / (std::pow(eu, 0.25) - std::pow(el, 0.25));

        // Assume the approximate linear relation in the log-log space,
        // $\ln(e/el) / \ln(eu/el) = \ln(t/tl) / \ln(tu/tl)$.
        double const t1 = 2.0 / hbarc_MeVfm;
        func(buff, t1);
        double const e1 = buff[0];
        tm = std::max(tl, t1 * std::pow(tu / t1, std::log(e / e1) / std::log(eu / e1)));
      }
    }

    int i;
    for (i = 0; i < 50; i++) {
      func(buff, tm);
      double const em = buff[0];
      double const re = (em - e) / std::max(e, em);
      if (std::abs(re) <= 1e-14) break;
      //std::printf("i=%d: %.17g e=%g em=%g re=%g\n", i, tm * hbarc_MeVfm, e, em, re);

      if (e < em) {
        tu = tm;
      } else {
        tl = tm;
      }

      if (i >= start_newton) {
        double const em_T = buff[1];

        // The normal Newton method does not converge well.
        //tm -= (em - e) / em_T;

        // As an improved version of the Newton method, we consider $t^p$ as
        // the horizontal axis, where $p$ is a constant.  The normal Newton
        // method corresponds to $p=1$.  In this case, we solve $(e - em) /
        // (t^p - tm^p) = ded(t^p)$ for the next $t$ to find the update rule.
        // However, this is still slow to converge.  We need to make $p$
        // extremely large such as $p=11$, but it becomes unstable.  constexpr
        // double p = 11.0;
        //tm *= std::pow(1.0 +  p * (e - em) / (tm * em_T), 1.0 / p);

        // As yet another version of the Newton method, we consider it in the
        // log-log space.  We may solve $\ln(e/e_m) / ln(t/t_m) =
        // d\ln(e_m)/d\ln(t_m)$ for the new $t$.  This significantly improves
        // the convergence.  However, it is better to apply this after
        // narrowing down the range to some extent because this is initially
        // unstable.
        tm *= std::pow(e / em, em / (em_T * tm));

        if (tl <= tm && tm < tu) continue;
      }

      tm = 0.5 * (tl + tu);
    }

    //std::printf("i=%d\n", i);
    return tm;
  }

  template<typename GetThermalState>
  void save_eos_table_vs_temperature(idt::runjam::runjam_context& ctx, const char* filename, GetThermalState proc){
    std::FILE* const file = std::fopen(filename, "w");
    if (!file) {
      std::fprintf(stderr, "%s: failed to open the file\n", filename);
      std::exit(1);
    }

#ifdef runjam_cmd_resolist_eos_CheckEnergyFundamentalRelation
    double echeck = 0.0; // fm^{-4}
    double echeck_temp = 0.0;
#endif

    std::fprintf(file, "#T[GeV] e[GeV/fm^3] P[GeV/fm^3] s[fm^{-3}] cs2 (e-3P)/T^4 (e_check-e)/T^4\n");
    static const int itempN = 800;
    double const temp_min =  0.001 / hbarc_GeVfm;
    double const temp_max = 10.000 / hbarc_GeVfm;
    double const dlnT = std::log(temp_max / temp_min) / itempN;
    for (int itemp = 0; itemp <= itempN; itemp++) {
      double const temp = temp_min * std::exp(dlnT * itemp); // fm^{-1}

      double state[4];
      proc(state, temp);
      double const energy_density = state[0];
      double const pressure = state[1];
      double const entropy_density = state[2];
      double const squared_sound_velocity = state[3];

#ifdef runjam_cmd_resolist_eos_CheckEnergyFundamentalRelation
      // e_check += \int_{T_{prev}}^{T} dT T p_TT.
      double integ;
      kashiwa::gauss_legendre_quadrature<32>(1, &integ, echeck_temp, temp, [] (double* integrand, double const T){
        double state[4];
        proc(state, T);

        // de/dT = (de/dp) (dp/dT) = s / cs2
        integrand[0] = state[2] / state[3];
      });
      echeck += integ;
      echeck_temp = temp;
      double const echeck_diff = (echeck - energy_density) / temp4;
#else
      double const echeck_diff = std::numeric_limits<double>::quiet_NaN();
#endif

      double const trace_anomaly = (energy_density - 3.0 * pressure) / std::pow(temp, 4.0);

      std::fprintf(
        file, "%22.16e %22.16e %22.16e %22.16e %22.16e %22.16\n",
        temp * hbarc_GeVfm,
        energy_density * hbarc_GeVfm,
        pressure * hbarc_GeVfm,
        entropy_density,
        squared_sound_velocity,
        trace_anomaly,
        echeck_diff);
    }

    std::fclose(file);
  }

  class EosHRG {
  private:
    struct integrand_for_1d_eos {
      int sign;
      double deg;
      double bmass;
      double bmu;

      /// @param[in] beta inverse temperature in fm^{-1}
      integrand_for_1d_eos(double beta, ResonanceRecord const* reso) {
        sign = -reso->bf;
        bmass = beta * reso->mass;
        bmu = beta * reso->mu;
        deg = reso->deg;
      }

      /*?lwiki
       * @fn void operator()(double* output, double t) const;
       * @param[out] output
       *   The integrand_for_1d_eoss for energy density and pressure
       * @param[in] t
       *   This specifies the integration variable $t$. The energy is calculated
       *   as $x = \beta E = \tau t^2$.
       */
      void operator()(double* output, double t) const {
        // variable transform t -> x
        double const tantt = std::tan(t * t);
        double const jacob = 2 * t * (tantt * tantt + 1);
        double const x = tantt + bmass;

        // d^3p /((2\pi)^3 E) = 4\pi p^2 dp / (8\pi^3 E) = p dE / 2\pi^2
        double const bp2 = x * x - bmass * bmass;
        double const jacob2 = (1.0 / (2.0 * M_PI * M_PI)) * std::sqrt(bp2);

        double f;
        if (x - bmu >= CONST_ENERGY_MAX_FACTOR) {
          f = deg * std::exp(-(x - bmu));
        } else {
          double const exp_ = std::exp(x - bmu);
          f = deg / (exp_ - sign);
        }

        double const w = jacob * jacob2 * f;
        output[0] = w * (x * x);
        output[1] = w * ((1.0 / 3.0) * bp2);
        // output[3] = f * x; // particle number
      }

      void integrand_pressure_TT(double* output, double const t) {
        // variable transform t -> x
        double const tantt = std::tan(t * t);
        double const jacob = 2 * t * (tantt * tantt + 1);
        double const x = tantt + bmass;

        // d^3p /((2\pi)^3 E) = 4\pi p^2 dp / (8\pi^3 E) = p dE / 2\pi^2
        double const bp2 = x * x - bmass * bmass;
        double const jacob2 = (1.0 / (2.0 * M_PI * M_PI)) * std::sqrt(bp2);

        double const xmu = x - bmu;
        double f;
        if (xmu >= CONST_ENERGY_MAX_FACTOR) {
          f = deg * std::exp(-xmu);
        } else {
          f = deg / (std::exp(xmu) - sign);
        }

        double const w = jacob * jacob2 * f;
        double const dlnf = xmu / (1.0 - sign * std::exp(-xmu)); // x * f * exp(x)
        output[0] = w * ((1.0 / 3.0) * bp2);
        output[1] = output[0] * dlnf;
        output[2] = output[1] * (2.0 * dlnf - xmu - 2.0);
      }
    };

  public:
    idt::runjam::ResonanceList rlist;
    EosHRG(idt::runjam::runjam_context& ctx): rlist(ctx) {}

    std::pair<double, double> get_energy_density_and_pressure(
      double temperature //!< [fm^{-1}]
    ) const {
      double energy_density = 0.0; // fm^{-4}
      double pressure = 0.0; // fm^{-4}

      int const iresoN = rlist.size();
      for (int ireso = 0; ireso < iresoN; ireso++) {
        ResonanceRecord const& reso = rlist[ireso];
        integrand_for_1d_eos integ(1.0 / temperature, &reso);
        double result[2];
        kashiwa::gauss_legendre_quadrature<256>(2, &result[0], 0.0, SQRT_TANGENT_ASYMPTOTE, integ);
        energy_density += result[0];
        pressure += result[1];
      }

      double const temp4 = std::pow(temperature, 4.0);
      energy_density *= temp4;
      pressure *= temp4;

      return std::make_pair(energy_density, pressure);
    }

    double pressure(
      double temperature //!< [fm^{-1}]
    ) const {
      double pressure = 0.0; // fm^{-4}

      int const iresoN = rlist.size();
      for (int ireso = 0; ireso < iresoN; ireso++) {
        ResonanceRecord const& reso = rlist[ireso];
        integrand_for_1d_eos integ(1.0 / temperature, &reso);
        double result[2];
        kashiwa::gauss_legendre_quadrature<256>(2, &result[0], 0.0, SQRT_TANGENT_ASYMPTOTE, integ);
        pressure += result[1];
      }
      pressure *= std::pow(temperature, 4.0);

      return pressure;
    }

    void get_pressure_derivatives(
      double* pressure_deriv,
      double temperature //!< [fm^{-1}]
    ) const {
      for (int i = 0; i < 3; i++)
        pressure_deriv[i] = 0.0; // fm^{-4+i}

      int const iresoN = rlist.size();
      for (int ireso = 0; ireso < iresoN; ireso++) {
        ResonanceRecord const& reso = rlist[ireso];
        integrand_for_1d_eos integ(1.0 / temperature, &reso);
        double result[3];
        kashiwa::gauss_legendre_quadrature<256>(3, &result[0], 0.0, SQRT_TANGENT_ASYMPTOTE, [&integ] (double* output, double const t) {
          integ.integrand_pressure_TT(output, t);
        });

        for (int i = 0; i < 3; i++)
          pressure_deriv[i] += result[i];
      }
      for (int i = 0; i < 3; i++)
        pressure_deriv[i] *= std::pow(temperature, 4 - i);
    }

    // state=(e p s cs2)
    void get_thermodynamic_quantities(double* state, double temperature /*!< [fm^{-1}] */) const {
      get_pressure_derivatives(&state[1], temperature);
      state[0] = temperature * state[2] - state[1]; // e = T p_T - p
      state[3] = state[2] / (temperature * state[3]); // cs2 = p_T / (T p_TT)
    }

    double pressure_T(
      double temperature //!< [fm^{-1}]
    ) const {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return pderiv[1];
    }

    double pressure_TT(
      double temperature //!< [fm^{-1}]
    ) const {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return pderiv[2];
    }

    double entropy_density(double temperature) const {
      return pressure_T(temperature);
    }
    double energy_density(double temperature) const {
      auto const [e, _] = get_energy_density_and_pressure(temperature);
      return e;
    }

    // double e_T(double temperature) const {
    //   double pderiv[3];
    //   get_pressure_derivatives(&pderiv[0], temperature);
    //   return temperature * pderiv[2];
    // }

  public:
    double T_vs_e(double e /*!< [fm^{-4}] */, double temperature_guess = -1.0) const {
      if (e <= 0.0) {
        if (e < 0.0) return -T_vs_e(-e); // for regularity
        return 0.0;
      }

      return log_log_newton_search(
        "HRG/T_vs_e", 1e-3 / hbarc_GeVfm, 10.000 / hbarc_GeVfm, temperature_guess, e,
        [this] (double (&buff)[2], double temperature) {
          double pderiv[3];
          get_pressure_derivatives(&pderiv[0], temperature);
          buff[0] = temperature * pderiv[1] - pderiv[0]; // e
          buff[1] = temperature * pderiv[2]; // de/dT
        }
      );
    }

  public:
    static void save_table_vs_temperature(idt::runjam::runjam_context& ctx, const char* filename);
  };

  std::unique_ptr<EosHRG> eosHRG;

  void EosHRG::save_table_vs_temperature(idt::runjam::runjam_context& ctx, const char* filename) {
    save_eos_table_vs_temperature(ctx, filename, [] (double (&state)[4], double temperature) {
      eosHRG->get_thermodynamic_quantities(&state[0], temperature);
    });
  }

  // See Eq. (16) and Table II in HotQCD:2014kol
  class HotQCD2014kol {
    // These coefficients are picked up from HotQCD:2014kol
    static constexpr double tc =  154.00 / hbarc_MeVfm; // [fm^{-1}]
    static constexpr double pi =  95.0 * M_PI * M_PI / 180.0;
    static constexpr double ct =  3.8706;
    static constexpr double an = -8.7704;
    static constexpr double bn =  3.9200;
    static constexpr double cn =  0.0000;
    static constexpr double dn =  0.3419;
    static constexpr double t0 =  0.9761;
    static constexpr double ad = -1.2600;
    static constexpr double bd =  0.8425;
    static constexpr double cd =  0.0000;
    static constexpr double dd = -0.0475;

    // We use the Taylor expansion around $T = 100~\text{MeV}$ for $T <
    // 100~\text{MeV}$.  The above parameters for the Pade approximation in
    // Eq. (16) of HotQCD:2014kol are only valid for the temperature range T =
    // [100, 400] MeV according to the caption in Table II.  In fact, the Pade
    // expansion has a pole around T = 45 MeV, so we cannot use the expression
    // for lower temperature.  We instead use a Taylor expansion of the Pade
    // expansion at T = 100 MeV for a lower temperature T < 100 MeV.
    //
    // Note: I checked different orders of the Taylor expansion, but
    // higher-order expansions are not necessarily stable.  I here decided to
    // use the second-order Taylor expansion.
    //
    // Note: I also checked the behavior of the Taylor expansion at T = 130 MeV
    // since the minimum temperature from the lattice data seems to be around
    // 130 MeV, but I gave up the Taylor expansion around T = 130 MeV because
    // it has a finite gap at T = 100 MeV with the second or third orders, and
    // the higher orders are oscillatory.
    //
    // Note: I also checked the behavior of the Taylor expansion with respect
    // to a different parameter $s = t^4$, whose dimension matches p, but it
    // didn't improve the results (though it's not worse).
    //
    // "tz" is the point where we do the Taylor expansion. "d[0-4]" are the
    // derivatives of the denominator at "tz". "n[0-4]" are of the
    // numerator. "p[1-3]n" are the numerator of the derivatives "p[1-3]".
    static constexpr double tz = 100.0 / 154.0;
    static constexpr double d0 = (((     tz + ad) * tz  + bd) * tz + cd) * tz + dd;
    static constexpr double d1 = ((4.0 * tz + 3.0 * ad) * tz  + 2.0 * bd) * tz + cd;
    static constexpr double d2 = (12.0 * tz + 6.0 * ad) * tz  + 2.0 * bd;
    static constexpr double d3 = 24.0 * tz + 6.0 * ad;
    static constexpr double d4 = 24.0;
    static constexpr double n0 = (((pi * tz + an) * tz  + bn) * tz + cn) * tz + dn;
    static constexpr double n1 = ((4.0 * pi * tz + 3.0 * an) * tz + 2.0 * bn) * tz + cn;
    static constexpr double n2 = (12.0 * pi * tz + 6.0 * an) * tz + 2.0 * bn;
    static constexpr double n3 = 24.0 * pi * tz + 6.0 * an;
    static constexpr double n4 = 24.0 * pi;
    static constexpr double p0 = n0 / d0;
    static constexpr double p1n = n1 * d0 - n0 * d1;
    static constexpr double p2n = (n2 * d0 - n0 * d2) * d0 - 2.0 * p1n * d1;
    static constexpr double p3n =
      (n3 * d0 * d0 - 3 * n1 * d0 * d2 + n0 * (3 * d1 * d2 - d0 * d3)) * d0
      - 3.0 * p2n * d1;
    static constexpr double p1 = p1n / (d0 * d0);
    static constexpr double p2 = p2n / (d0 * d0 * d0);
    static constexpr double p3 = p3n / (d0 * d0 * d0 * d0);

  public:
    static double pressure(
      double temperature //!< [fm^{-1}]
    ) {

      double const t = temperature / tc;
      double expansion;

      if (t < tz) {
        // std::printf("p0 = %21.15e\n", p0);
        // std::printf("p1 = %21.15e\n", p1);
        // std::printf("p2 = %21.15e\n", p2);
        // std::printf("p3 = %21.15e\n", p3);
        // std::exit(1);

        double const dt = t - tz;
        //expansion = p1 * dt + p0; // -> non-monotonic e/T^4
        expansion = (p2 * dt + p1) * dt + p0; // -> fine, but de/dT jumps at T = 100 MeV
        //expansion = ((p3 * dt + p2) * dt + p1) * dt + p0; // -> non-monotonic p/T^4

        // Expansion wrt t^4
        // constexpr double tz4 = tz * tz * tz * tz;
        // constexpr double pb0 = p0;
        // constexpr double pb1 = p1 / (4.0 * tz * tz * tz);
        // constexpr double pb2 = (p2 * tz - 3.0 * p1) / (16.0 * tz * tz * tz * tz4);
        // constexpr double pb3 = ((p3 * tz - 9.0 * p2) * tz + 21.0 * p1) / (64.0 * tz * tz * tz * tz4 * tz4);
        // std::printf("pb0 = %21.15e\n", pb0);
        // std::printf("pb1 = %21.15e\n", pb1);
        // std::printf("pb2 = %21.15e\n", pb2);
        // std::printf("pb3 = %21.15e\n", pb3);
        // std::exit(1);
        // double const dt4 = t * t * t * t - tz4;
        // expansion = pb1 * dt4 + pb0; // -> non-monotonic e/T^4
        // expansion = (pb2 * dt4 + pb1) * dt4 + pb0; // -> non-monotonic e/T^4
        // expansion = ((pb3 * dt4 + pb2) * dt4 + pb1) * dt4 + pb0; // -> non-monotonic p/T^4
      } else {
        double const den = (((     t + ad) * t  + bd) * t + cd) * t + dd;
        double const num = (((pi * t + an) * t  + bn) * t + cn) * t + dn;
        expansion = num / den;
      }

      double const coeff = 0.5 * (1.0 + std::tanh(ct * (t - t0)));
      double const temp4 = std::pow(temperature, 4.0);
      return temp4 * coeff * expansion;
    }

  public:
    /// Calculate the derivatives of the pressure with respect to the
    /// temperature up to the second order, $P, dP/dT, d^2p/dT^2$.
    ///
    /// @param[out] pressure_deriv The obtained pressure derivatives are stored
    /// in this array. pressure_deriv[0] is the pressure $p$, pressure_deriv[1]
    /// is the first-order derivative $dP/dT$, and pressure_deriv[2] is the
    /// second-order derivative $d^2P/dT^2$.
    /// @param[in] temperature The temperature in unit of $\mathrm{fm}^{-1}$.
    static void get_pressure_derivatives(
      double* pressure_deriv,
      double temperature //!< [fm^{-1}]
    ) {
      double const t = temperature / tc;

      //! 0th..2nd order derivatives of the Pade expansion with respect to $t$.
      double ex0;
      double ex1;
      double ex2;

      if (t < tz) {
        double const dt = t - tz;
        ex0 = (p2 * dt + p1) * dt + p0;
        ex1 = 2.0 * p2 * dt + p1;
        ex2 = 2.0 * p2;
      } else {
        //! 0th..2nd order derivatives of the denominator and numerator of the
        //! Pade expansion.
        double const den0 = (((     t + ad) * t  + bd) * t + cd) * t + dd;
        double const den1 = ((4.0 * t + 3.0 * ad) * t  + 2.0 * bd) * t + cd;
        double const den2 = (12.0 * t + 6.0 * ad) * t  + 2.0 * bd;
        double const num0 = (((pi * t + an) * t  + bn) * t + cn) * t + dn;
        double const num1 = ((4.0 * pi * t + 3.0 * an) * t + 2.0 * bn) * t + cn;
        double const num2 = (12.0 * pi * t + 6.0 * an) * t + 2.0 * bn;

        double const p0num = num0;
        double const p1num = num1 * den0 - p0num * den1;
        double const p2num = (num2 * den0 - num0 * den2) * den0 - 2.0 * p1num * den1;
        ex0 = p0num / den0;
        ex1 = p1num / (den0 * den0);
        ex2 = p2num / (den0 * den0 * den0);
      }

      //! 0th..2nd order derivatives of $\tanh$ with respect to $c_t(t-t0)$.
      double const th0 = std::tanh(ct * (t - t0));
      double const th1 = 1.0 - th0 * th0;
      double const th2 = -2.0 * th0 * th1;

      //! 0th..2nd order derivatives of $(1/2)(1+th0)$ with respect to $t$.
      double const coeff0 = 0.5 * (1.0 + th0);
      double const coeff1 = 0.5 * th1 * ct;
      double const coeff2 = 0.5 * th2 * (ct * ct);

      //! 0th..2nd order derivatives of $T^4$ with respect to $t$.
      double const temp40 = std::pow(temperature, 4.0);
      double const temp41 = 4.0 * std::pow(temperature, 3.0) * tc;
      double const temp42 = 12.0 * (temperature * temperature) * (tc * tc);

      pressure_deriv[0] = temp40 * coeff0 * ex0;
      pressure_deriv[1] = (temp40 * (coeff0 * ex1 + coeff1 * ex0) + temp41 * coeff0 * ex0) / tc;
      pressure_deriv[2] = (temp40 * coeff0 * ex2
        + 2.0 * (temp40 * coeff1 + temp41 * coeff0) * ex1
        + (temp40 * coeff2 + 2.0 * temp41 * coeff1 + temp42 * coeff0) * ex0) / (tc * tc);
    }

    static void get_thermodynamic_quantities(double* state, double temperature /*!< [fm^{-1}] */) {
      get_pressure_derivatives(&state[1], temperature);
      state[0] = temperature * state[2] - state[1]; // e = T p_T - p
      state[3] = state[2] / (temperature * state[3]); // cs2 = p_T / (T p_TT)
    }

    static double pressure_T(double temperature /*!< [fm^{-1}] */) {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return pderiv[1];
    }
    static double pressure_TT(double temperature /*!< [fm^{-1}] */) {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return pderiv[2];
    }
    static double entropy_density(double temperature) {
      return pressure_T(temperature);
    }
    static double energy_density(double temperature) {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return temperature * pderiv[1] - pderiv[0];
    }

  public:
    static double T_vs_e(double e /*!< [fm^{-4}] */, double temperature_guess = -1.0) {
      if (e <= 0.0) {
        if (e < 0.0) return -T_vs_e(-e); // for regularity
        return 0.0;
      }

      return log_log_newton_search(
        "HotQCD2014kol/T_vs_e", 1e-10 / hbarc_GeVfm, 10.000 / hbarc_GeVfm, temperature_guess, e,
        [] (double (&buff)[2], double temperature) {
          double pderiv[3];
          get_pressure_derivatives(&pderiv[0], temperature);
          buff[0] = temperature * pderiv[1] - pderiv[0]; // e
          buff[1] = temperature * pderiv[2]; // de/dT
        }
      );
    }

  public:
    static void save_table_vs_temperature(idt::runjam::runjam_context& ctx, const char* filename) {
      save_eos_table_vs_temperature(ctx, filename, [] (double (&state)[4], double temperature) {
        get_thermodynamic_quantities(state, temperature);
      });
    }
  };

  class HRG_QGP {
    static constexpr double delta_Tc = 10.0 / hbarc_MeVfm; // [fm^{-1}]
    static constexpr double Tc = 154.00 / hbarc_MeVfm; // [fm^{-1}]

  public:
    static double pressure(
      double temperature // [fm^{-1}]
    ) {
      double const pHRG = eosHRG->pressure(temperature);
      double const pLattice = HotQCD2014kol::pressure(temperature);
      double const func_T = (temperature - Tc) / delta_Tc;

      double const p_HQ
        = 1.0 / 2.0 * (1.0 - std::tanh(func_T)) * pHRG
        + 1.0 / 2.0 * (1.0 + std::tanh(func_T)) * pLattice;

      return p_HQ;
    }


public:
    static void get_pressure_derivatives(
      double* pressure_deriv,
      double temperature // [fm^{-1}]
    ) {
      double pHRG[3];
      eosHRG->get_pressure_derivatives(&pHRG[0], temperature);

      double pQGP[3];
      HotQCD2014kol::get_pressure_derivatives(&pQGP[0], temperature);

      double const func_T = (temperature - Tc) / delta_Tc;
      double th[3];
      th[0] = std::tanh(func_T);
      th[1] = (1.0 - th[0] * th[0]) / delta_Tc;
      th[2] = -2.0 * th[0] * th[1] / delta_Tc;

      pressure_deriv[0]
        = 0.5 * (1.0 - th[0]) * pHRG[0]
        + 0.5 * (1.0 + th[0]) * pQGP[0];
      pressure_deriv[1]
        = 0.5 * ((1.0 - th[0]) * pHRG[1] - th[1] * pHRG[0])
        + 0.5 * ((1.0 + th[0]) * pQGP[1] + th[1] * pQGP[0]);
      pressure_deriv[2]
        = 0.5 * ((1.0 - th[0]) * pHRG[2] - 2.0 * th[1] * pHRG[1] - th[2] * pHRG[0])
        + 0.5 * ((1.0 + th[0]) * pQGP[2] + 2.0 * th[1] * pQGP[1] + th[2] * pQGP[0]);
    }

    static void get_thermodynamic_quantities(double* state, double temperature /*!< [fm^{-1}] */) {
      get_pressure_derivatives(&state[1], temperature);
      state[0] = temperature * state[2] - state[1]; // e = T p_T - p
      state[3] = state[2] / (temperature * state[3]); // cs2 = p_T / (T p_TT)
    }

    static double pressure_T(
      double temperature // [fm^{-1}]
    ) {
      double p[3];
      get_pressure_derivatives(&p[0], temperature);
      return p[1];
    }
    static double pressure_TT(
      double temperature // [fm^{-1}]
    ) {
      double p[3];
      get_pressure_derivatives(&p[0], temperature);
      return p[2];
    }
    static double entropy_density(double temperature) {
      return pressure_T(temperature);
    }
    static double energy_density(double temperature) {
      double pderiv[3];
      get_pressure_derivatives(&pderiv[0], temperature);
      return temperature * pderiv[1] - pderiv[0];
    }

  public:
    static double T_vs_e(double e /*!< [fm^{-4}] */, double temperature_guess = -1.0) {
      if (e <= 0.0) {
        if (e < 0.0) return -T_vs_e(-e); // for regularity
        return 0.0;
      }

      return log_log_newton_search(
        "HRG_QGP/T_vs_e", 1e-10 / hbarc_GeVfm, 10.000 / hbarc_GeVfm, temperature_guess, e,
        [] (double (&buff)[2], double temperature) {
          double pderiv[3];
          get_pressure_derivatives(&pderiv[0], temperature);
          buff[0] = temperature * pderiv[1] - pderiv[0]; // e
          buff[1] = temperature * pderiv[2]; // de/dT
        }
      );
    }

  public:
    static void save_table_vs_temperature(idt::runjam::runjam_context& ctx, const char* filename) {
      save_eos_table_vs_temperature(ctx, filename, [] (double (&state)[4], double temperature){
        HRG_QGP::get_thermodynamic_quantities(&state[0], temperature);
      });
    }
  };

}


  void debug_T_vs_e_1(double temperature) {
    temperature /= hbarc_MeVfm;
    //double const guess = eosHRG->T_vs_e(eosHRG->energy_density(temperature));
    // double const guess = HotQCD2014kol::T_vs_e(HotQCD2014kol::energy_density(temperature));
    double const guess = HRG_QGP::T_vs_e(HRG_QGP::energy_density(temperature));
    std::printf("(%g) %.17g %g\n", temperature * hbarc_MeVfm, guess * hbarc_MeVfm, guess / temperature - 1.0);
  }

  void debug_T_vs_e() {
    debug_T_vs_e_1(0.0010);
    debug_T_vs_e_1(0.1000);
    debug_T_vs_e_1(0.2000);
    debug_T_vs_e_1(0.3000);
    debug_T_vs_e_1(0.4000);
    debug_T_vs_e_1(0.5000);
    debug_T_vs_e_1(1.0000);
    debug_T_vs_e_1(2.0000);
    debug_T_vs_e_1(3.1415);
    debug_T_vs_e_1(10.000);
    debug_T_vs_e_1(50.000);
    debug_T_vs_e_1(150.000);
    debug_T_vs_e_1(200.000);
    debug_T_vs_e_1(500.000);
    debug_T_vs_e_1(1000.000);
    debug_T_vs_e_1(2000.000);
    debug_T_vs_e_1(5000.000);
  }

  int cmd_resolist_eos(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args) {
    (void) args;

    // initialize
    if (!eosHRG) eosHRG = std::make_unique<EosHRG>(ctx);

    //debug_T_vs_e(); return 0;

    EosHRG::save_table_vs_temperature(ctx, "eos.txt");
    HotQCD2014kol::save_table_vs_temperature(ctx, "eos_lattice.txt");
    HRG_QGP::save_table_vs_temperature(ctx, "eos_QGP_HRG.txt");

    return 0;
  }
}
