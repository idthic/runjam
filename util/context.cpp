#include <cctype>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "context.hpp"

namespace idt {
namespace util {

  static int stricmp(const char* a, const char* b) {
    while (*a && *b) {
      int const diff = (int) std::tolower(*a++) - (int) std::tolower(*b++);
      if (diff) return diff;
    }
    return *a ? 1 : *b ? -1 : 0;
  }

  const char* application_context::get_value(const char* name) const {
    dictionary_type::const_iterator it = m_data.find(name);
    if (it != m_data.end())
      return it->second.c_str();

    return std::getenv(name);
  }

  bool application_context::assert_config(const char* name, const char* value) const {
    return std::strcmp(this->get_config(name, value), value) == 0;
  }

  bool application_context::read_config(const char*& var, const char* name) const {
    const char* str = get_value(name);
    if (str && str[0]) {
      var = str;
      return true;
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, const char* value) {
    this->set_value(name, value);
  }

  bool application_context::read_config(std::string& var, const char* name) const {
    const char* str = get_value(name);
    if (str && str[0]) {
      var = str;
      return true;
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, std::string const& value) {
    this->set_value(name, value);
  }

  bool application_context::read_config(int& var, const char* name) const {
    const char* str = get_value(name);
    if (str && std::strchr("+-0123456789", str[0])) {
      var = std::atoi(str);
      return true;
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, int value) {
    char buff[64];
    std::sprintf(buff, "%d", value);
    this->set_value(name, buff);
  }

  bool application_context::read_config(std::uint32_t& var, const char* name) const {
    const char* str = get_value(name);
    if (str && str[0]) {
      var = (std::uint32_t) strtoll(str, NULL, 10);
      return true;
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, std::uint32_t value) {
    char buff[32];
    std::sprintf(buff, "%" PRIu32 "u", value);
    this->set_value(name, buff);
  }

  bool application_context::read_config(double& var, const char* name) const {
    const char* str = get_value(name);
    if (str && std::strchr("+-.0123456789", str[0])) {
      var = std::atof(str);
      return true;
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, double value) {
    char buff[32];
    std::sprintf(buff, "%.15g", value);
    this->set_value(name, buff);
  }

  bool application_context::read_config(bool& var, const char* name) const {
    const char* str = get_value(name);
    if (str && str[0]) {
      if (stricmp(str, "true") == 0) {
        var = true;
        return true;
      } else if (stricmp(str, "false") == 0) {
        var = false;
        return true;
      } else {
        char* end;
        bool const value = std::strtoll(str, &end, 10);
        if (*end) return false;
        var = value;
        return true;
      }
    } else {
      return false;
    }
  }
  void application_context::set_config(const char* name, bool value) {
    this->set_value(name, value ? "true" : "false");
  }

}
}
