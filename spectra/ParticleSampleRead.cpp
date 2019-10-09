#include <cstdlib>
#include <cctype>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ParticleSampleRead.hpp"

namespace idt {
namespace hydro2jam {

  static bool is_comment_line(std::string const& line) {
    std::size_t i = 0;
    while (std::isspace(line[i])) i++;
    return line[i] == '#';
  }

  void ParticleSampleRead::readFile() {
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
        double px,py,pz,e,em;
        int ir;
        double tau,rx,ry,eta;
        if (!(is >> px >> py >> pz >> e >> em >> ir >> tau >> rx >> ry >> eta))
          goto error_invalid_format;

        this->addParticleTauEta(ir, px, py, pz, em, rx, ry, tau, eta);
      }

      return;
    }

  error_failed_to_open:
    std::cerr
      << "spectra/ParticleSampleRead: failed to open the file ("
      << this->fname_particlesample_dat << ")" << std::endl;
    std::exit(EXIT_FAILURE);
    return; /*NOTREACHED*/

  error_invalid_format:
    std::cerr
      << this->fname_particlesample_dat << ":" << iline << ": invalid format (spectra/ParticleSampleRead)"
      << std::endl;
    std::exit(EXIT_FAILURE);
    return; /*NOTREACHED*/
  }

  void ParticleSampleRead::update() {
    this->readFile();
  }

}
}
