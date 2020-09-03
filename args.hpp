#ifndef hydro2jam_context_hpp
#define hydro2jam_context_hpp
#include <cstring>
#include <cstdio>
#include <string>
#include <map>

namespace idt {
namespace util {

  class application_context {
    std::map<std::string, std::string> data;

  public:
    void set_value(const char* key, const char* value) {
      data[key] = value;
    }
    void set_value(const char* key, int value) {
      char buff[100];
      std::sprintf(buff, "%d", value);
      data[key] = buff;
    }
    void set_value(const char* key, double value) {
      char buff[100];
      std::sprintf(buff, "%.15g", value);
      data[key] = buff;
    }
    void set_value(const char* key, bool value) {
      set_value(key, value ? 1 : 0);
    }
    void set_value(const char* key, std::string const& value) {
      set_value(key, value.c_str());
    }

    const char* get_value(const char* key) const {
      typedef std::map<std::string, std::string>::const_iterator It;
      It it = data.find(key);
      if (it != data.end())
        return it->second.c_str();
      return std::getenv(key);
    }

  public:
    bool read_config(std::string& value, const char* key) const {
      const char* str = get_value(key);
      if (str == 0) return false;
      value = str;
      return true;
    }

    bool read_config(int& value, const char* key) const {
      const char* str= get_value(key);
      if (str && isdigit(*str)) {
        value = std::atoi(str);
        return true;
      } else
        return false;
    }

    bool read_config(double& value, const char* key) const {
      const char* str= get_value(key);
      if (str && std::strchr("+-.0123456789", *str)) {
        value = std::atof(str);
        return true;
      } else
        return false;
    }

    bool read_config(bool& value, const char* key) const {
      const char* str = get_value(key);
      if (!str || !*str) return false;

      if(isdigit(*str))
        value = std::atoi(str) != 0;
      else
        value = std::strcmp(str, "false") != 0;
      return true;
    }

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
namespace hydro2jam {

  class hydro2jam_context: public idt::util::application_context {
  public:
    int seed() const {
      return this->get_config("hydro2jam_seed", 18371);
    }
    int nevent(int defaultValue = 1) const {
      return this->get_config("hydro2jam_nevent", defaultValue);
    }
    std::string outdir() const {
      return this->get_config<std::string>("hydro2jam_output_directory", "out");
    }

    int eospce() const {
      return this->get_config("hydrojet_eospce", 6);
    }
    int kintmp() const {
      return this->get_config("hydrojet_kintmp", 5);
    }
    std::string resodata() const {
      return this->get_config<std::string>("hydro2jam_resodata", "data/ResonanceJam.dat");
    }
    std::string indir() const {
      return this->get_config<std::string>("hydrojet_directory", "test");
    }

  };

}
}

#endif
