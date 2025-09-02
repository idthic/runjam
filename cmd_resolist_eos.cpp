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
#include <ksh/integrator.hpp>

#include "args.hpp"
#include "ResonanceList.hpp"

namespace idt::runjam {
namespace {

  static const double hbarc_GeVfm = 0.197327053; // hydrojet
  static const double SQRT_TANGENT_ASYMPTOTE = M_SQRT2 / M_2_SQRTPI;
  static const double CONST_ENERGY_MAX_FACTOR = 100.0;

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
}

  int cmd_resolist_eos(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args) {
    (void) args;

    idt::runjam::ResonanceList rlist(ctx);
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

      double energy_density = 0.0; // fm^{-4}
      double pressure = 0.0; // fm^{-4}

      int const iresoN = rlist.size();
      for (int ireso = 0; ireso < iresoN; ireso++) {
        ResonanceRecord const& reso = rlist[ireso];
        integrand_for_1d_eos integ(1.0 / temp, &reso);
        double result[2];
        kashiwa::gauss_legendre_quadrature<256>(2, &result[0], 0.0, SQRT_TANGENT_ASYMPTOTE, integ);
        energy_density += result[0];
        pressure += result[1];
      }

      double const trace_anomaly = energy_density - 3.0 * pressure;
      double const temp4 = std::pow(temp, 4.0);
      energy_density *= temp4;
      pressure *= temp4;

      std::fprintf(file, "%21.15e %21.15e %21.15e %21.15e\n",
        temp * hbarc_GeVfm,
        energy_density * hbarc_GeVfm,
        pressure * hbarc_GeVfm,
        trace_anomaly);
    }
    std::fclose(file);

    return 0;
  }

}
