#ifndef runjam_context_hpp
#define runjam_context_hpp
#include <cstring>
#include <cstdio>
#include <string>
#include <map>
#include "util.hpp"

namespace idt {
namespace runjam {

  class runjam_context: public idt::util::application_context {
  public:
    int seed() const {
      return this->get_config("runjam_seed", 18371);
    }
    int nevent(int defaultValue = 1) const {
      return this->get_config("runjam_nevent", defaultValue);
    }
    std::string outdir() const {
      std::string ret = this->get_config<std::string>("runjam_output_directory", "out");
      if (ret.size() == 0) ret = ".";
      return ret;
    }

    int eospce() const {
      return this->get_config("hydrojet_eospce", 6);
    }
    int kintmp() const {
      return this->get_config("hydrojet_kintmp", 5);
    }
    std::string resodata() const;
    std::string indir() const {
      return this->get_config<std::string>("hydrojet_directory", "test");
    }
  };

  struct runjam_commandline_arguments {
    std::string subcommand;
    std::string initType;
    std::string initPath;

  public:
    int read(int argc, char** argv, runjam_context& ctx);
  };

}
}

#endif
