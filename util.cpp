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
  if (str && std::isdigit(*str)) {
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


static std::mt19937 random_engine;

void set_random_seed(int seed) {
  random_engine.seed(seed);
}

double urand() {
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(random_engine);
}

int irand(int n) {
  std::uniform_int_distribution<int> dist(0, n - 1);
  return dist(random_engine);
}

int irand_poisson(double lambda){
  std::poisson_distribution<int> dist(lambda);
  return dist(random_engine);
}

}
}
