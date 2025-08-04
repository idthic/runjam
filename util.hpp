/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2014-2020, Koichi Murase <myoga.murase at gmail.com>

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

#ifndef runjam_util_hpp
#define runjam_util_hpp
#include <cstddef>
#include <map>
#include <string>

namespace idt {
  // static const double hbarc_MeVfm = 197.32; // hydrojet/src/spectra/ElementReso.h
  // static const double hbarc_MeVfm = 197.327053 ; // 197.327053 hydrojet/src/physicsbase/Const.h
  static const double hbarc_MeVfm = 197.3269718; // 197.3269718(44) idt/rfh/i2/common/def.h
  static const double hbarc_GeVfm = hbarc_MeVfm / 1000.0;
}

namespace idt {
namespace util {
  bool ends_with(const char* s, std::size_t l, const char* suffix, std::size_t len);
  bool ends_with(std::string const& s, std::string const& suffix);
  bool ends_with(std::string const& s, const char* suffix);
  bool ends_with(const char* s, const char* suffix);
}
}

#include <random>

namespace idt {
namespace util {
  std::mt19937& random_engine();
  void set_random_seed(int seed);
  double urand(); // double [0,1)
  double nrand(); // normal distribution N(μ=0, σ^2)
  int irand(int n); // int [0,n)
  std::size_t irand(std::size_t n); // int [0,n)
  int irand_poisson(double lambda);
  int irand_binomial(int trial, double probability);
  std::size_t irand_binomial(std::size_t trial, double probability);
}
}

#endif
