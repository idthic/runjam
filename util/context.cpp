#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include "context.hpp"

namespace idt {
namespace util {

const char* application_context::get_value(const char* key) const {
  typedef dictionary_type::const_iterator It;
  It it = m_data.find(key);
  if (it != m_data.end())
    return it->second.c_str();
  return std::getenv(key);
}

bool application_context::read_config(const char*& value, const char* key) const {
  const char* str = get_value(key);
  if (str == 0) return false;
  value = str;
  return true;
}
void application_context::set_config(const char* key, const char* value) {
  this->set_value(key, value);
}

bool application_context::read_config(std::string& value, const char* key) const {
  const char* str = get_value(key);
  if (str == 0) return false;
  value = str;
  return true;
}
void application_context::set_config(const char* key, std::string const& value) {
  this->set_value(key, value);
}

bool application_context::read_config(int& value, const char* key) const {
  const char* str= get_value(key);
  if (str && std::strchr("+-0123456789", *str)) {
    value = std::atoi(str);
    return true;
  } else
    return false;
}
void application_context::set_config(const char* key, int value) {
  char buff[100];
  std::sprintf(buff, "%d", value);
  m_data[key] = buff;
}

bool application_context::read_config(double& value, const char* key) const {
  const char* str= get_value(key);
  if (str && std::strchr("+-.0123456789", *str)) {
    value = std::atof(str);
    return true;
  } else
    return false;
}
void application_context::set_config(const char* key, double value) {
  char buff[100];
  std::sprintf(buff, "%.15g", value);
  m_data[key] = buff;
}

bool application_context::read_config(bool& value, const char* key) const {
  const char* str = get_value(key);
  if (!str || !*str) return false;

  if (isdigit(*str))
    value = std::atoi(str) != 0;
  else
    value = std::strcmp(str, "false") != 0;
  return true;
}
void application_context::set_config(const char* key, bool value) {
  set_config(key, value ? 1 : 0);
}

}
}
