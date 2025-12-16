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
#include <ksh/integrator.hpp>

#include "args.hpp"
#include "ResonanceList.hpp"

namespace idt::runjam {
namespace {

  static const double hbarc_GeVfm = 0.197327053; // hydrojet
  static const double SQRT_TANGENT_ASYMPTOTE = M_SQRT2 / M_2_SQRTPI;
  static const double CONST_ENERGY_MAX_FACTOR = 100.0;

  void save_HotQCD2014kol_eos(idt::runjam::runjam_context& ctx);
  void save_HRG_eos(idt::runjam::runjam_context& ctx);
  void save_HRG_QGP_eos(idt::runjam::runjam_context& ctx);

  class eosHRG_t {
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
    };

  public:
    idt::runjam::ResonanceList rlist;
    eosHRG_t(idt::runjam::runjam_context& ctx): rlist(ctx) {}

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
  };

  std::unique_ptr<eosHRG_t> eosHRG;

  void save_HRG_eos(idt::runjam::runjam_context& ctx){
    std::FILE* const file = std::fopen("eos.txt", "w");
    if (!file) {
      std::fprintf(stderr, "eos.txt: failed to open the file\n");
      std::exit(1);
    }

    std::fprintf(file, "#temperature(GeV) energy_density(GeV/fm^3) pressure(GeV/fm^3) (e-3P)/T^4\n");

    static const int itempN = 800;
    double const temp_min =  0.001 / hbarc_GeVfm;
    double const temp_max = 10.000 / hbarc_GeVfm;
    double const dlnT = std::log(temp_max / temp_min) / itempN;
    for (int itemp = 0; itemp <= itempN; itemp++) {
      double const temp = temp_min * std::exp(dlnT * itemp); // fm^{-1}

      auto const [energy_density, pressure] = eosHRG->get_energy_density_and_pressure(temp);

      double const trace_anomaly = (energy_density - 3.0 * pressure) / std::pow(temp, 4.0);

      std::fprintf(file, "%21.15e %21.15e %21.15e %21.15e\n",
        temp * hbarc_GeVfm,
        energy_density * hbarc_GeVfm,
        pressure * hbarc_GeVfm,
        trace_anomaly);
    }

    std::fclose(file);
  }

  // See Eq. (16) and Table II in HotQCD:2014kol
  double pressure_HotQCD2014kol(
    double temperature //!< [fm^{-1}]
  ) {
    constexpr double tc =  154.00 / hbarc_MeVfm; // [fm^{-1}]
    constexpr double pi =  95.0 * M_PI * M_PI / 180.0;
    constexpr double ct =  3.8706;
    constexpr double an = -8.7704;
    constexpr double bn =  3.9200;
    constexpr double cn =  0.0000;
    constexpr double dn =  0.3419;
    constexpr double t0 =  0.9761;
    constexpr double ad = -1.2600;
    constexpr double bd =  0.8425;
    constexpr double cd =  0.0000;
    constexpr double dd = -0.0475;

    double const t = temperature / tc;
    double expansion;

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
    constexpr double tz = 100.0 / 154.0;
    if (t < tz) {
      constexpr double d0 = (((     tz + ad) * tz  + bd) * tz + cd) * tz + dd;
      constexpr double d1 = ((4.0 * tz + 3.0 * ad) * tz  + 2.0 * bd) * tz + cd;
      constexpr double d2 = (12.0 * tz + 6.0 * ad) * tz  + 2.0 * bd;
      constexpr double d3 = 24.0 * tz + 6.0 * ad;
      constexpr double d4 = 24.0;
      constexpr double n0 = (((pi * tz + an) * tz  + bn) * tz + cn) * tz + dn;
      constexpr double n1 = ((4.0 * pi * tz + 3.0 * an) * tz + 2.0 * bn) * tz + cn;
      constexpr double n2 = (12.0 * pi * tz + 6.0 * an) * tz + 2.0 * bn;
      constexpr double n3 = 24.0 * pi * tz + 6.0 * an;
      constexpr double n4 = 24.0 * pi;
      constexpr double p0 = n0 / d0;
      constexpr double p1n = n1 * d0 - n0 * d1;
      constexpr double p2n = (n2 * d0 - n0 * d2) * d0 - 2.0 * p1n * d1;
      constexpr double p3n =
        (n3 * d0 * d0 - 3 * n1 * d0 * d2 + n0 * (3 * d1 * d2 - d0 * d3)) * d0
        - 3.0 * p2n * d1;
      constexpr double p1 = p1n / (d0 * d0);
      constexpr double p2 = p2n / (d0 * d0 * d0);
      constexpr double p3 = p3n / (d0 * d0 * d0 * d0);

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

  void save_HotQCD2014kol_eos(idt::runjam::runjam_context& ctx) {
    static const int itempN = 800;
    double const temp_min =  0.001 / hbarc_GeVfm;
    double const temp_max = 10.000 / hbarc_GeVfm;
    double const dlnT = std::log(temp_max / temp_min) / itempN;

    //lattice QGP
    std::FILE* const file_QGP = std::fopen("eos_lattice.txt", "w");
    if (!file_QGP) {
      std::fprintf(stderr, "eos.txt: failed to open the file\n");
      std::exit(1);
    }

    std::fprintf(file_QGP, "#temperature(GeV) energy_density(GeV/fm^3) pressure(GeV/fm^3) (e-3P)/T^4\n");

    for (int itemp = 0; itemp <= itempN; itemp++) {
      double const temp = temp_min * std::exp(dlnT * itemp); // fm^{-1}
      double energy_density = 0.0; // fm^{-4}

      //dp = s dT
      //de = T ds
      //e = \int de
      //  = \int T ds
      //  = \int T d(dp/dT)  温度0でエネルギー0 なので0からTの積分をする
      //  = \int T dT d^{2}p/dT^{2}

      int intTmax = 1000;
      double const dT = temp / intTmax;
      for (int intT = 0; intT < intTmax; intT++) {
        //積分する
        double T = dT * (intT + 0.5);
        double p_T = pressure_HotQCD2014kol(T);
        double epsilon = dT/10000;
        double p_T_pluss = pressure_HotQCD2014kol(T + epsilon);
        double p_T_minus = pressure_HotQCD2014kol(T - epsilon);
        energy_density += T * dT * ((p_T_pluss - 2*p_T + p_T_minus) / (epsilon *  epsilon));
      }

      double const pressure = pressure_HotQCD2014kol(temp); // fm^{-4}
      double const trace_anomaly = (energy_density - 3.0 * pressure) / std::pow(temp, 4.0);

      std::fprintf(file_QGP, "%21.15e %21.15e %21.15e %21.15e\n",
                   temp * hbarc_GeVfm,
                   energy_density * hbarc_GeVfm,
                   pressure * hbarc_GeVfm,
                   trace_anomaly);
    }
    std::fclose(file_QGP);
  }

  double pressure_HRG_QGP(
    double temperature // [fm^{-1}]
  ) {
    constexpr double delta_Tc = 10.0 / hbarc_MeVfm; // [fm^{-1}]
    constexpr double Tc = 154.00 / hbarc_MeVfm; // [fm^{-1}]

    double const pHRG = eosHRG->pressure(temperature);
    double const pLattice = pressure_HotQCD2014kol(temperature);
    double const func_T = (temperature - Tc)/delta_Tc;

    double p_HQ = 1.0 / 2.0 * (1.0 - std::tanh(func_T)) * pHRG
      + 1.0 / 2.0 * (1.0 + std::tanh(func_T)) * pLattice;

    return p_HQ;
  }

  void save_HRG_QGP_eos(idt::runjam::runjam_context& ctx){
    //output pressure_HRG_QGP

    static const int itempN = 800;
    double const temp_min =  0.001 / hbarc_GeVfm;
    double const temp_max = 10.000 / hbarc_GeVfm;
    double const dlnT = std::log(temp_max / temp_min) / itempN;

    //QGP_HRG
    std::FILE* const file_QGP_HRG = std::fopen("eos_QGP_HRG.txt", "w");
    if (!file_QGP_HRG) {
      std::fprintf(stderr, "eos_QGP_HRG.txt: failed to open the file\n");
      std::exit(1);
    }

    std::fprintf(file_QGP_HRG, "#temperature(GeV) energy_density(GeV/fm^3) pressure(GeV/fm^3) (e-3P)/T^4\n");

    for (int itemp = 0; itemp <= itempN; itemp++) {
      double const temp = temp_min * std::exp(dlnT * itemp); // fm^{-1}
      double energy_density = 0.0; // fm^{-4}

      int intTmax = 10;//QGPの時のように1000にするとすごく時間がかかる
      double const dT = temp / intTmax;

      for (int intT = 0; intT < intTmax; intT++) {
        //積分する
        double T = dT * (intT + 0.5);
        double p_T = pressure_HRG_QGP(T);
        double epsilon = dT/10000;
        double p_T_pluss = pressure_HRG_QGP(T + epsilon);
        double p_T_minus = pressure_HRG_QGP(T - epsilon);
        energy_density += T * dT * ((p_T_pluss - 2*p_T + p_T_minus) / (epsilon *  epsilon));
      }

      double const pressure = pressure_HRG_QGP(temp); // fm^{-4}
      double const trace_anomaly = (energy_density - 3.0 * pressure) / std::pow(temp, 4.0);

      std::fprintf(file_QGP_HRG, "%21.15e %21.15e %21.15e %21.15e\n",
                   temp * hbarc_GeVfm,
                   energy_density * hbarc_GeVfm,
                   pressure * hbarc_GeVfm,
                   trace_anomaly);
    }
    std::fclose(file_QGP_HRG);
  }

}
  int cmd_resolist_eos(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args) {
    (void) args;

    // initialize
    if (!eosHRG) eosHRG = std::make_unique<eosHRG_t>(ctx);

    save_HRG_eos(ctx);
    save_HotQCD2014kol_eos(ctx);
    save_HRG_QGP_eos(ctx);

    return 0;
  }
}
