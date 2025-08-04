/* This file is a part of runjam <https://github.com/idthic/runjam> and taken
   from libkashiwa <https://github.com/akinomyoga/libkashiwa>.

   Copyright (C) 2011-2020, Koichi Murase <myoga.murase at gmail.com>.

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

#ifndef kashiwa_integrator_hpp
#define kashiwa_integrator_hpp
#include <cstddef>
#include <cmath>
#include <algorithm>
namespace kashiwa {
namespace integrator_detail {
  struct point_weight_pair {
    double t;
    double w;
  };

  template<std::size_t I>
  struct GaussLegendre {
    static const std::size_t order = I;
    static const std::size_t data_size = (I + 1) / 2;
    static point_weight_pair data[data_size];
  };

  template<std::size_t I>
  struct GaussLaguerre {
    static const std::size_t order = I;
    static const std::size_t data_size = I;
    static point_weight_pair data[I];
  };
}

  template<typename F>
  double IntegrateByTrapezoid(const double lower, const double upper, const int iN, const F& f) {
    double const dx = (upper - lower) / iN;

    double s = 0.5 * (f(lower) + f(upper));
    for (int i = 1; i < iN; i++)
      s += f(lower + i * dx);

    return s * dx;
  }

  template<std::size_t I, typename F>
  double IntegrateByGaussLegendre(const double lower, const double upper, const F& f) {
    typedef integrator_detail::GaussLegendre<I> traits_t;
    double const center = 0.5 * (upper + lower);
    double const dxdt   = 0.5 * (upper - lower);

    double s = 0.0;
    for (std::size_t i = 0; i < I / 2; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;
      s += w * (f(center + dxdt * t) + f(center - dxdt * t));
    }

    return s * dxdt;
  }

  template<std::size_t I, typename F>
  double IntegrateByGaussLaguerre(const double lower, const double scale, F f) {
    typedef integrator_detail::GaussLaguerre<I> traits_t;
    double s = 0.0;
    for (std::size_t i = 0; i < traits_t::data_size; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;
      s += w * f(lower + scale * t);
    }
    return s * scale;
  }

  template<std::size_t I>
  class GaussLegendreIntegrator {
    double const xmin;
    double const xmax;
    double const center;
    double const dxdt;
    integrator_detail:: point_weight_pair data[(I + 1) / 2];
  public:
    GaussLegendreIntegrator(double xmin,double xmax):
      xmin(xmin), xmax(xmax),
      center(0.5 * (xmax + xmin)),
      dxdt(0.5 * (xmax - xmin))
    {
      typedef integrator_detail::GaussLegendre<I> traits_t;

      for (int i = 0; i < (I + 1) / 2; i++) {
        double const t = traits_t::data[i].t;
        double const w = traits_t::data[i].w;
        this->data[i].t = dxdt * t;
        this->data[i].w = dxdt * w;
      }
    }
    template<typename F>
    double Integrate(const F& f) {
      double s = 0;
      for (int i = 0; i < I / 2; i++) {
        double const t = this->data[i].t;
        double const w = this->data[i].w;
        s += w * (f(center + t) + f(center - t));
      }

      return s;
    }
  };

  template<int I, typename F>
  void gauss_legendre_quadrature(int N, double* result, double xmin, double xmax, F f) {
    typedef integrator_detail::GaussLegendre<I> traits_t;
    double const center = 0.5 * (xmax + xmin);
    double const dxdt   = 0.5 * (xmax - xmin);

    std::fill(result, result + N, 0.0);

    for (int i = 0; i < I / 2; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;

      double value1[N];
      f(value1, center + dxdt * t);
      double value2[N];
      f(value2, center - dxdt * t);

      for (int k = 0; k < N; k++)
        result[k] += w * (value1[k] + value2[k]);
    }

    for (int k = 0; k < N; k++)
      result[k] *= dxdt;
  }

  template<typename F>
  void gauss_chebyshev_quadrature(int Nd, double* result, int Ni, double xmin, double xmax, F f) {
    double const center = 0.5 * (xmax + xmin);
    double const dxdt   = 0.5 * (xmax - xmin);
    double const dh = M_PI / (2.0 * Ni);

    if (Ni & 1) {
      f(result, center);
    } else {
      std::fill(result, result + Nd, 0.0);
    }

    for (int i = 0, iN = Ni / 2; i < iN; i++) {
      double const t = std::cos((2 * i + 1) * dh);
      double const w = std::sqrt(1.0 - t * t);

      double value1[Nd];
      f(value1, center + dxdt * t);
      double value2[Nd];
      f(value2, center - dxdt * t);

      for (int k = 0; k < Nd; k++)
        result[k] += w * (value1[k] + value2[k]);
    }

    double const factor = dxdt * M_PI / Ni;
    for (int k = 0; k < Nd; k++)
      result[k] *= factor;
  }

}
#endif
