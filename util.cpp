/* Copyright (C) 2014-2020, Koichi Murase @akinomyoga.
   This file is a part of runjam <https://github.com/idthic/runjam>.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA  */

#include "util.hpp"
#include <cstddef>
#include <cstring>
#include <random>

namespace idt::util {

bool ends_with(const char* s, std::size_t l, const char* suffix, std::size_t len) {
  return l >= len && std::memcmp(s + (l - len), suffix, len) == 0;
}
bool ends_with(std::string const& s, std::string const& suffix) {
  return ends_with(s.c_str(), s.size(), suffix.c_str(), suffix.size());
}
bool ends_with(std::string const& s, const char* suffix) {
  return ends_with(s.c_str(), s.size(), suffix, std::strlen(suffix));
}
bool ends_with(const char* s, const char* suffix) {
  return ends_with(s, std::strlen(s), suffix, std::strlen(suffix));
}

//-----------------------------------------------------------------------------

static std::mt19937 eng;

std::mt19937& random_engine() {
  return eng;
}

void set_random_seed(int seed) {
  eng.seed(seed);
}

double urand() {
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(eng);
}

double nrand() {
  static std::normal_distribution<double> dist;
  return dist(eng);
}

int irand(int n) {
  std::uniform_int_distribution<int> dist(0, n - 1);
  return dist(eng);
}
std::size_t irand(std::size_t n) {
  std::uniform_int_distribution<std::size_t> dist(0, n - 1);
  return dist(eng);
}

int irand_poisson(double lambda){
  std::poisson_distribution<int> dist(lambda);
  return dist(eng);
}

int irand_binomial(int trial, double probability) {
  std::binomial_distribution<int> dist(trial, probability);
  return dist(eng);
}
std::size_t irand_binomial(std::size_t trial, double probability) {
  std::binomial_distribution<std::size_t> dist(trial, probability);
  return dist(eng);
}

}
