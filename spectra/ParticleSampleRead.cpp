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
    ResonanceListPCE rlist;
    std::string fname_particlesample_dat;
  public:
    ParticleSampleRead(std::string const& fname_particlesample_dat):
      rlist(-1, 6, "ResonanceJam.dat"),
      fname_particlesample_dat(fname_particlesample_dat) {}

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
          double px,py,pz,e,em;
          int ir;
          double tau,rx,ry,eta;
          if (!(is >> px >> py >> pz >> e >> em >> ir >> tau >> rx >> ry >> eta))
            goto error_invalid_format;

          int const pdg = rlist[ir].generatePDGCode();
          this->addParticleTauEta(pdg, px, py, pz, em, rx, ry, tau, eta);
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

    virtual void update() {
      this->readFile();
    }
  };

  class ParticleSampleFactory: ParticleSampleFactoryBase {
    virtual ParticleSampleBase* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type != "psample") return 0;
      return new ParticleSampleRead(inputfile);
    }
  } instance;

}
}
}
