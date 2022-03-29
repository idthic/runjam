/* Copyright (C) 2011-2020, Koichi Murase @akinomyoga.
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

#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <util.hpp>
#include <args.hpp>
#include "ResonanceList.hpp"

namespace idt {
namespace runjam {

  int ResonanceRecord::generatePDGCode() const {
    auto const& codes = this->pdg_codes;
    if (codes.size() == 1)
      return codes[0];
    else
      return codes[idt::util::irand(codes.size())];
  }

  void ResonanceList::readFile(std::string const& fn_resodata) {
    std::ifstream ifs(fn_resodata.c_str());
    if (!ifs)  {
      std::cerr << "ResonanceList: failed to open the file '" << fn_resodata << "'." << std::endl;
      std::exit(1);
    }

    this->data.clear();
    this->data.reserve(151);

    std::string line;
    int iline = 0;
    while (std::getline(ifs, line)) {
      iline++;

      // Skip empty lines and comment lines
      std::size_t pos = line.find_first_not_of(" \t");
      if (pos == std::string::npos || line[pos] == '#') continue;

      std::istringstream istr(line);
      ResonanceRecord reso;
      int bftype, pdg;
      std::string name;
      istr >> reso.mass >> reso.deg >> reso.degeff
           >> reso.mu >> bftype >> reso.anti
           >> reso.key >> name;
      while (istr >> pdg)
        reso.pdg_codes.push_back(pdg);
      if (reso.pdg_codes.empty()) {
        std::cerr << fn_resodata << ":" << iline << ": invalid format." << std::endl;
        std::exit(1);
      }

      // // check
      // std::cout << reso.key;
      // for (auto const pdg : reso.pdg_codes) std::cout << " " << pdg;
      // std::cout << std::endl;

      reso.mass /= hbarc_MeVfm; // fm^{-1}
      reso.mu   /= hbarc_MeVfm; // fm^{-1}
      reso.bf   = bftype == 1 ? -1: bftype == 2 ? 1 : bftype; // boson -1, fermion +1
      data.emplace_back(reso);
    }

    if (data.empty()) {
      std::cerr << fn_resodata << ": empty." << std::endl;
      std::exit(1);
    }
  }

}
}
