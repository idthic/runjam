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

#include "util.hpp"
#include <cstddef>
#include <cstring>
#include <random>
#include <string>

#include <cstdarg>

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
// taken from fitmass/common.cpp

std::string vstrprintf(const char* fmt, std::va_list va) {
  std::va_list va1;
  va_copy(va1, va);
  std::size_t const sz = vsnprintf(NULL, 0, fmt, va1);
  va_end(va1);

  char* buff = (char*) std::malloc(sz + 1);
  va_copy(va1, va);
  std::vsprintf(buff, fmt, va1);
  va_end(va1);

  std::string str = buff;
  std::free(buff);
  return str;
}

std::string strprintf(const char* fmt, ...) {
  std::va_list va;
  va_start(va, fmt);
  std::string const str = vstrprintf(fmt, va);
  va_end(va);
  return str;
}

std::string strprintf(std::string const& fmt, ...) {
  std::va_list va;
  va_start(va, fmt);
  std::string const str = vstrprintf(fmt.c_str(), va);
  va_end(va);
  return str;
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
  if (lambda == 0.0) return 0;
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
