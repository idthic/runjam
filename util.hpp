#ifndef runjam_util_hpp
#define runjam_util_hpp
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

  class application_context {
    std::map<std::string, std::string> data;

  public:
    void set_value(const char* key, const char* value) {
      data[key] = value;
    }
    void set_value(const char* key, int value);
    void set_value(const char* key, double value);
    void set_value(const char* key, bool value) {
      set_value(key, value ? 1 : 0);
    }
    void set_value(const char* key, std::string const& value) {
      set_value(key, value.c_str());
    }
    const char* get_value(const char* key) const;

  public:
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
      read_config(value, key, defaultValue);
      return value;
    }
  };
}
}

namespace idt {
namespace util {
  void set_random_seed(int seed);
  double urand(); // double [0,1)
  int irand(int n); // int [0,n)
  int irand_poisson(double lambda);
}
}

#endif
