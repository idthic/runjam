// -*- mode: c++; coding: utf-8 -*-
#ifndef idt_util_context_hpp
#define idt_util_context_hpp
#include <cstdint>
#include <string>
#include <utility>
#if __cplusplus >= 201103L
# include <unordered_map>
#else
# include <map>
#endif

namespace idt {
namespace util {

  class application_context {
#if __cplusplus >= 201103L
    typedef std::unordered_map<std::string, std::string> dictionary_type;
#else
    typedef std::map<std::string, std::string> dictionary_type;
#endif
    dictionary_type m_data;

  public:
    const char* get_value(const char* name) const;
    void set_value(const char* name, const char* value) {
      m_data[name] = value;
    }
    void set_value(const char* name, std::string const& value) {
      m_data[name] = value;
    }
#if __cplusplus >= 201103L
    void set_value(const char* name, std::string&& value) {
      m_data[name] = std::move(value);
    }
#endif

  public:
    bool read_config(const char*& var, const char* name) const;
    bool read_config(std::string& var, const char* name) const;
    bool read_config(int& var, const char* name) const;
    bool read_config(std::uint32_t& var, const char* name) const;
    bool read_config(double& var, const char* name) const;
    bool read_config(bool& var, const char* name) const;

  public:
    void set_config(const char* name, const char* value);
    void set_config(const char* name, std::string const& value);
    void set_config(const char* name, int value);
    void set_config(const char* name, std::uint32_t value);
    void set_config(const char* name, double value);
    void set_config(const char* name, bool value);

  public:
    template<typename T>
    bool read_config(T& value, const char* name, T const& default_value) const {
      bool const ret = this->read_config(value, name);
      if (!ret) value = default_value;
      return ret;
    }

    template<typename T>
    T get_config(const char* name, T default_value) const {
      this->read_config(default_value, name);
      return default_value;
    }

    bool assert_config(const char* name, const char* value) const;
  };

}
}
#endif
