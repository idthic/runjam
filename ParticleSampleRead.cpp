/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2014-2020, Koichi Murase <myoga.murase at gmail.com>

   SPDX-License-Identifier: GPL-2.0-or-later

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA  */

#include <cstdlib>
#include <cctype>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ParticleSample.hpp"
#include "ResonanceList.hpp"

namespace idt {
namespace runjam {
namespace {

  static bool is_comment_line(std::string const& line) {
    std::size_t i = 0;
    while (std::isspace(line[i])) i++;
    return line[i] == '#';
  }

  class ParticleSampleRead: public ParticleSampleBase {
    ResonanceList rlist;
    std::string fname_particlesample_dat;
    double m_overSamplingFactor = 1.0;
  public:
    ParticleSampleRead(runjam_context const& ctx, std::string const& fname_particlesample_dat):
      rlist("ResonanceJam.dat"),
      fname_particlesample_dat(fname_particlesample_dat),
      m_overSamplingFactor(ctx.get_config("runjam_oversampling_factor", 1.0))
    {}

  public:
    // この枠組ではファイルに保存されていた粒子集合の oversampling が
    // 分からないので、環境に設定されていた値を信じてそのまま返す。
    virtual double getOverSamplingFactor() const override {
      return m_overSamplingFactor;
    }

  private:
    void readFile() {
      this->clearParticleList();

      int iline = 0;
      {
        std::ifstream ifs(this->fname_particlesample_dat.c_str());
        if (!ifs) goto error_failed_to_open;

        std::string line;
        while (std::getline(ifs, line)) {
          iline++;
          if (is_comment_line(line)) continue;

          std::istringstream is(line);
          double px, py, pz, e, em;
          int ir;
          double tau,rx,ry,eta;
          if (!(is >> px >> py >> pz >> e >> em >> ir >> tau >> rx >> ry >> eta))
            goto error_invalid_format;

          int const pdg = rlist[ir].generatePDGCode();
          this->addParticleMilne(pdg, px, py, pz, em, rx, ry, tau, eta);
        }

        return;
      }

    error_failed_to_open:
      std::cerr
        << "ParticleSampleRead: failed to open the file ("
        << this->fname_particlesample_dat << ")" << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/

    error_invalid_format:
      std::cerr
        << this->fname_particlesample_dat << ":" << iline << ": invalid format (ParticleSampleRead)"
        << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/
    }

    virtual void update() {
      this->readFile();
    }
  };

  class ParticleSampleFactory: ParticleSampleFactoryBase {
    virtual std::unique_ptr<ParticleSampleBase> CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type != "psample") return nullptr;
      return std::make_unique<ParticleSampleRead>(ctx, inputfile);
    }
  } instance;

}
}
}
