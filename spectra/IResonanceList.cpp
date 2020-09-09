#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "IResonanceList.hpp"
#include "../args.hpp"

namespace idt {
namespace runjam {

ResonanceListPCE::resonance ResonanceListPCE::resT[5][21] = {
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

ResonanceListPCE::ResonanceListPCE(runjam_context const& ctx) {
  int const eospce = ctx.eospce();
  int const kintmp = ctx.kintmp();
  std::string const resodata = ctx.resodata();
  initialize(kintmp, eospce, resodata);
}
ResonanceListPCE::ResonanceListPCE(int kineticTemp, int eos_pce, std::string const& fname_rlist) {
  initialize(kineticTemp, eos_pce, fname_rlist);
}

void ResonanceListPCE::initialize(int kineticTemp, int eos_pce, std::string const& fname_rlist) {
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

  if (fname_rlist.size() <= 0) {
    std::cerr << "ResonanceListPCE::.ctor! something is wrong: Resonance table not available " << std::endl;
    std::exit(1);
  }

  std::ifstream fdata(fname_rlist.c_str(), std::ios::in);
  if (!fdata)  {
    std::cerr << "ResonanceListPCE::.ctor! unable to open file " << fname_rlist << std::endl;
    std::exit(1);
  }

  this->data.clear();
  this->data.reserve(151);

  std::string line;
  int iline = 0;
  while (std::getline(fdata, line)) {
    iline++;

    // 空行・コメント行はスキップ
    std::size_t pos = line.find_first_not_of(" \t");
    if (pos == std::string::npos || line[pos] == '#') continue;

    std::istringstream istr(line);
    resonance record;
    int bftype;
    istr >> record.mass
         >> record.deg
         >> record.degeff
         >> record.mu
         >> bftype
         >> record.anti;
    if (!istr) {
      std::cerr << fname_rlist << ":" << iline << ": invalid format." << std::endl;
      std::exit(1);
    }
    record.mass /= hbarc_MeVfm; // fm^{-1}
    record.mu   /= hbarc_MeVfm; // fm^{-1}
    record.bf   = bftype == 1 ? -1: bftype == 2 ? 1 : bftype;
    data.emplace_back(record);
  }
  if (data.empty()) {
    std::cerr << fname_rlist << ": empty." << std::endl;
    std::exit(1);
  } else if (nreso_check >= 0 && data.size() != nreso_check) {
    std::cerr << fname_rlist << ": unexpected number of resonances"
              << " [" << data.size() << ", expected: " << nreso_check << ", eospce = " << eos_pce << "]."
              << std::endl;
    std::exit(1);
  }

	if (eos_pce == 1) {
    int const itemp = kineticTemp - 1;
	  if (itemp >= 0 && itemp < 5) {
	    for (int i = 0; i < data.size(); i++) {
        resonance& recdst = data[i];
        resonance& recsrc = resT[itemp][i];

	      recdst.mu = recsrc.mu;
        // recdst.mass = recsrc.mass;
        // recdst.deg    = recsrc.deg;
        // recdst.degeff = recsrc.degeff;
        // recdst.bf     = recsrc.bf;
        // recdst.anti   = recsrc.anti;

	      recdst.mu /= hbarc_MeVfm;
        // recdst.mass /= hbarc_MeVfm;
        // if (recdst.bf == 1) recdst.bf = -1;
        // if (recdst.bf == 2) recdst.bf = 1;
	    }
	  } else {
      std::cerr << "ResonanceListPCE::.ctor! something is wrong: Resonance table not available " << std::endl;
      std::exit(1);
	  }
	}
}

}
}
