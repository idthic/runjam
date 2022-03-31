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
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <random>

namespace idt {
namespace util {

void application_context::set_value(const char* key, int value) {
  char buff[100];
  std::sprintf(buff, "%d", value);
  data[key] = buff;
}

void application_context::set_value(const char* key, double value) {
  char buff[100];
  std::sprintf(buff, "%.15g", value);
  data[key] = buff;
}

const char* application_context::get_value(const char* key) const {
  typedef std::map<std::string, std::string>::const_iterator It;
  It it = data.find(key);
  if (it != data.end())
    return it->second.c_str();
  return std::getenv(key);
}

bool application_context::read_config(std::string& value, const char* key) const {
  const char* str = get_value(key);
  if (str == 0) return false;
  value = str;
  return true;
}

bool application_context::read_config(int& value, const char* key) const {
  const char* str= get_value(key);
  if (str && std::strchr("+-0123456789", *str)) {
    value = std::atoi(str);
    return true;
  } else
    return false;
}

bool application_context::read_config(double& value, const char* key) const {
  const char* str= get_value(key);
  if (str && std::strchr("+-.0123456789", *str)) {
    value = std::atof(str);
    return true;
  } else
    return false;
}

bool application_context::read_config(bool& value, const char* key) const {
  const char* str = get_value(key);
  if (!str || !*str) return false;

  if(isdigit(*str))
    value = std::atoi(str) != 0;
  else
    value = std::strcmp(str, "false") != 0;
  return true;
}


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
}
