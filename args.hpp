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

#ifndef runjam_args_hpp
#define runjam_args_hpp
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

    std::string resodata() const;
    std::string indir() const {
      return this->get_config<std::string>("hydrojet_directory", "test");
    }

    std::string cachedir() const;
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
