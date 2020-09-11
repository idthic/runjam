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
      reso.bf   = bftype == 1 ? -1: bftype == 2 ? 1 : bftype;
      data.emplace_back(reso);
    }

    if (data.empty()) {
      std::cerr << fn_resodata << ": empty." << std::endl;
      std::exit(1);
    }
  }

ResonanceRecord ResonanceListPCE::resT[5][21] = {
  //degeff is "effective degree of freedom" which means
  //the d.o.f. for decay into pi^-
  //Note that Lambda bar doesn't decay into pi^-.
  //So just calculate the effective d.o.f. for pbar.

  // Tf = 80 MeV
  {
    {984.8,  3.   ,  1.   ,  329.429,  1,  0},
    {1232.,  16.  ,  5.333,  551.282,  2,  0},
    {1232.,  16.  ,  5.333,  551.282,  2,  1},
    {547. ,  1.   ,  0.23 ,  234.236,  1,  0},
    {957.8,  1.   ,  0.44 ,  276.004,  1,  0},
    {980. ,  1.   ,  1.   ,  190.385,  1,  0},
    {497.7,  2.   ,  0.69 ,  233.670,  1,  0},
    {893.5,  6.   ,  4.   ,  328.863,  1,  0},
    {893.5,  6.   ,  4.   ,  328.863,  1,  0},
    {1115.,  2.   ,  1.28 ,  524.280,  2,  0},
    {1115.,  2.   ,  1.28 ,  524.280,  2,  1},//No decay into pi^-
    {782. ,  3.   ,  2.62 ,  257.020,  1,  0},
    {1020.,  3.   ,  1.47 ,  440.076,  1,  0},
    {768. ,  9.   ,  6.   ,  190.385,  1,  0},
    {800. ,  1.   ,  1.   ,  190.385,  1,  0},
    {1193.,  6.   ,  2.   ,  564.943,  2,  0},
    {1193.,  6.   ,  0.96 ,  564.943,  2,  1},//Sigma bar -> pi^-
    {139. ,  3.   ,  1.   ,  95.1925,  1,  0},
    {493.6,  2.   ,  1.   ,  233.670,  1,  0},
    {939. ,  4.   ,  2.   ,  456.089,  2,  0},
    {939. ,  4.   ,  2.   ,  456.089,  2,  1},
  },

  // Tf = 100 MeV
  {
    {984.8,  3.   ,  1.   ,  256.576,  1,  0},
    {1232.,  16.  ,  5.333,  432.124,  2,  0},
    {1232.,  16.  ,  5.333,  432.124,  2,  1},
    {547. ,  1.   ,  0.23 ,  173.262,  1,  0},
    {957.8,  1.   ,  0.44 ,  220.928,  1,  0},
    {980. ,  1.   ,  1.   ,  166.628,  1,  0},
    {497.7,  2.   ,  0.69 ,  180.805,  1,  0},
    {893.5,  6.   ,  4.   ,  264.120,  1,  0},
    {893.5,  6.   ,  4.   ,  264.120,  1,  0},
    {1115.,  2.   ,  1.28 ,  403.997,  2,  0},
    {1115.,  2.   ,  1.28 ,  403.997,  2,  1},
    {782. ,  3.   ,  2.62 ,  224.948,  1,  0},
    {1020.,  3.   ,  1.47 ,  344.861,  1,  0},
    {768. ,  9.   ,  6.   ,  166.628,  1,  0},
    {800. ,  1.   ,  1.   ,  166.628,  1,  0},
    {1193.,  6.   ,  2.   ,  435.491,  2,  0},
    {1193.,  6.   ,  0.96 ,  435.491,  2,  1},
    {139. ,  3.   ,  1.   ,  83.3141,  1,  0},
    {493.6,  2.   ,  1.   ,  180.805,  1,  0},
    {939. ,  4.   ,  2.   ,  348.810,  2,  0},
    {939. ,  4.   ,  2.   ,  348.810,  2,  1}
  },

  // Tf = 120 MeV
  {
    {984.8,  3.   ,  1.   ,  179.974,  1,  0},
    {1232.,  16.  ,  5.333,  310.546,  2,  0},
    {1232.,  16.  ,  5.333,  310.546,  2,  1},
    {547. ,  1.   ,  0.23 ,  115.292,  1,  0},
    {957.8,  1.   ,  0.44 ,  159.026,  1,  0},
    {980. ,  1.   ,  1.   ,  129.363,  1,  0},
    {497.7,  2.   ,  0.69 ,  128.598,  1,  0},
    {893.5,  6.   ,  4.   ,  193.279,  1,  0},
    {893.5,  6.   ,  4.   ,  193.279,  1,  0},
    {1115.,  2.   ,  1.28 ,  287.104,  2,  0},
    {1115.,  2.   ,  1.28 ,  287.104,  2,  1},
    {782. ,  3.   ,  2.62 ,  174.640,  1,  0},
    {1020.,  3.   ,  1.47 ,  247.722,  1,  0},
    {768. ,  9.   ,  6.   ,  129.363,  1,  0},
    {800. ,  1.   ,  1.   ,  129.363,  1,  0},
    {1193.,  6.   ,  2.   ,  309.536,  2,  0},
    {1193.,  6.   ,  0.96 ,  309.536,  2,  1},
    {139. ,  3.   ,  1.   ,  64.6814,  1,  0},
    {493.6,  2.   ,  1.   ,  128.598,  1,  0},
    {939. ,  4.   ,  2.   ,  245.865,  2,  0},
    {939. ,  4.   ,  2.   ,  245.865,  2,  1}
  },

  // Tf = 140 MeV
  {
    {984.8,  3.   ,  1.   ,  102.415,  1,  0},
    {1232.,  16.  ,  5.333,  186.183,  2,  0},
    {1232.,  16.  ,  5.333,  186.183,  2,  1},
    {547. ,  1.   ,  0.23 ,  61.7501,  1,  0},
    {957.8,  1.   ,  0.44 ,  93.0018,  1,  0},
    {980. ,  1.   ,  1.   ,  81.3296,  1,  0},
    {497.7,  2.   ,  0.69 ,  63.3578,  1,  0},
    {893.5,  6.   ,  4.   ,  117.009,  1,  0},
    {893.5,  6.   ,  4.   ,  117.009,  1,  0},
    {1115.,  2.   ,  1.28 ,  171.004,  2,  0},
    {1115.,  2.   ,  1.28 ,  171.004,  2,  1},
    {782. ,  3.   ,  2.62 ,  109.795,  1,  0},
    {1020.,  3.   ,  1.47 ,  148.084,  1,  0},
    {768. ,  9.   ,  6.   ,  81.3296,  1,  0},
    {800. ,  1.   ,  1.   ,  81.3296,  1,  0},
    {1193.,  6.   ,  2.   ,  184.418,  2,  0},
    {1193.,  6.   ,  0.96 ,  184.418,  2,  1},
    {139. ,  3.   ,  1.   ,  40.6648,  1,  0},
    {493.6,  2.   ,  1.   ,  63.3578,  1,  0},
    {939. ,  4.   ,  2.   ,  145.518,  2,  0},
    {939. ,  4.   ,  2.   ,  145.518,  2,  1}
  },

  // Tf = 160 MeV
  {
    {984.8,  3.   ,  1.   ,  30.6714,  1,  0},
    {1232.,  16.  ,  5.333,  61.4611,  2,  0},
    {1232.,  16.  ,  5.333,  61.4611,  2,  1},
    {547. ,  1.   ,  0.23 ,  16.8848,  1,  0},
    {957.8,  1.   ,  0.44 ,  28.8978,  1,  0},
    {980. ,  1.   ,  1.   ,  27.5734,  1,  0},
    {497.7,  2.   ,  0.69 ,  24.9125,  1,  0},
    {893.5,  6.   ,  4.   ,  38.6992,  1,  0},
    {893.5,  6.   ,  4.   ,  38.6992,  1,  0},
    {1115.,  2.   ,  1.28 ,  56.2548,  2,  0},
    {1115.,  2.   ,  1.28 ,  56.2548,  2,  1},
    {782. ,  3.   ,  2.62 ,  37.2240,  1,  0},
    {1020.,  3.   ,  1.47 ,  48.5552,  1,  0},
    {768. ,  9.   ,  6.   ,  27.5734,  1,  0},
    {800. ,  1.   ,  1.   ,  27.5734,  1,  0},
    {1193.,  6.   ,  2.   ,  60.7121,  2,  0},
    {1193.,  6.   ,  0.96 ,  60.7121,  2,  1},
    {139. ,  3.   ,  1.   ,  13.7867,  1,  0},
    {493.6,  2.   ,  1.   ,  24.9125,  1,  0},
    {939. ,  4.   ,  2.   ,  47.6744,  2,  0},
    {939. ,  4.   ,  2.   ,  47.6744,  2,  1}
  }
};

ResonanceListPCE::ResonanceListPCE(runjam_context const& ctx): base(ctx) {
  int const eospce = ctx.eospce();
  int const kintmp = ctx.kintmp();
  initialize(kintmp, eospce, ctx.resodata());
}
ResonanceListPCE::ResonanceListPCE(int kineticTemp, int eos_pce, std::string const& fn_resodata): base(fn_resodata) {
  initialize(kineticTemp, eos_pce, fn_resodata);
}

void ResonanceListPCE::initialize(int kineticTemp, int eos_pce, std::string const& fn_resodata) {
  int nreso_check;
  switch (eos_pce) {
  case 0: nreso_check = 21; break;
  case 1: nreso_check = 21; break;
  case 2: nreso_check = 151; break;
  case 3: nreso_check = 151; break;
  case 4: nreso_check = 75; break; // EOS-Q
  case 5: nreso_check = 21; break;
  case 6: nreso_check = -1; break; // no-check
  }
  if (nreso_check >= 0 && data.size() != nreso_check) {
    std::cerr << fn_resodata << ": unexpected number of resonances"
              << " [" << data.size() << ", expected: " << nreso_check << ", eospce = " << eos_pce << "]."
              << std::endl;
    std::exit(1);
  }

  if (eos_pce == 1) {
    int const itemp = kineticTemp - 1;
    if (!(0 <= itemp && itemp < 5)) {
      std::cerr << "ResonanceListPCE! unsupported kintmp=" << kineticTemp << "." << std::endl;
      std::exit(1);
    }

    for (int i = 0; i < data.size(); i++) {
      ResonanceRecord& reso1 = data[i];
      ResonanceRecord& reso2 = resT[itemp][i];
      reso1.mu = reso2.mu / hbarc_MeVfm;
      // reso1.mass   = reso2.mass / hbarc_MeVfm;
      // reso1.deg    = reso2.deg;
      // reso1.degeff = reso2.degeff;
      // reso1.bf     = reso2.bf == 1 ? -1 : reso2.bf == 2 ? 1 : reso2.bf;
      // reso1.anti   = reso2.anti;
    }
  }
}

}
}
