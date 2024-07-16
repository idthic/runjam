#ifndef idt_util_context_hpp
#define idt_util_context_hpp
#include <map>
#include <string>

namespace idt {
namespace util {

  class application_context {
    typedef std::map<std::string, std::string> dictionary_type;
    dictionary_type m_data;

  public:
    const char* get_value(const char* key) const;
    void set_value(const char* key, const char* value) {
      m_data[key] = value;
    }
    void set_value(const char* key, std::string const& value) {
      set_value(key, value.c_str());
    }

  public:
    void set_config(const char* key, const char* value);
    void set_config(const char* key, std::string const& value);
    void set_config(const char* key, int value);
    void set_config(const char* key, double value);
    void set_config(const char* key, bool value);

  public:
    bool read_config(const char*& value, const char* key) const;
    bool read_config(std::string& value, const char* key) const;
    bool read_config(int& value, const char* key) const;
    bool read_config(double& value, const char* key) const;
    bool read_config(bool& value, const char* key) const;

  public:
    template<typename T>
    bool read_config(T& value, const char* key, T const& defaultValue) const {
      bool const ret = read_config(value, key);
      if (!ret) value = defaultValue;
      return ret;
    }

    template<typename T>
    T get_config(const char* key, T const& defaultValue) const {
      T value;
      this->read_config<T>(value, key, defaultValue);
      return value;
    }
  };
}
}
#endif
