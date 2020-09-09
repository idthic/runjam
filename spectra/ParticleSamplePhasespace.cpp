#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include "IParticleSample.hpp"

namespace idt {
namespace runjam {
namespace {

  class ParticleSampleReadPhasespaceData: public ParticleSampleBase {
    std::string fname_phasespace_dat;
  public:
    ParticleSampleReadPhasespaceData(std::string const& fname_phasespace_dat):
      fname_phasespace_dat(fname_phasespace_dat)
    {
      this->m_currentSampleIndex = -1;
      this->m_numberOfSamples = -1;
    }

  private:
    int m_numberOfSamples;
    int m_currentSampleIndex;
    std::vector<std::vector<Particle*> > pcache;
  public:
    void setNumberOfSamples(int value) {
      this->m_numberOfSamples = value;
    }
    int getNumberOfSamples() const {
      return m_numberOfSamples;
    }

  private:
    virtual ParticleIDType::value_type getParticleIdType() const { return ParticleIDType::PDGCode; }

    virtual std::vector<Particle*> const& getParticleList() const {
      if (m_currentSampleIndex < 0) {
        std::cerr << "ParticleSampleReadPhasespaceData: A phasespace file has not been read." << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (m_currentSampleIndex >= this->pcache.size()) {
        std::cerr << "ParticleSampleReadPhasespaceData: No more samples." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->pcache[m_currentSampleIndex];
    }

    void readPhasespaceDat() {
      this->clearParticleList();

      bool eventCountKnown = m_numberOfSamples >= 0;
      this->pcache.clear();
      if (eventCountKnown)
        this->pcache.reserve(this->m_numberOfSamples);

      int iline;
      {
        std::ifstream ifs(this->fname_phasespace_dat.c_str());
        if (!ifs) goto error_failed_to_open;

        int isample;
        std::string line;
        iline = 0;
        for (isample = 0; ; isample++) {
          iline++;
          if (!std::getline(ifs, line)) goto error_invalid_format;
          int npart, dummy;
          if (2 != std::sscanf(line.c_str(), " %d %d", &npart, &dummy)) break;

          if (isample == this->m_numberOfSamples) goto error_invalid_format;

          if (npart < 0) goto error_invalid_format;

          this->pcache.resize(isample + 1);
          for (int i = 0; i < npart; i++) {
            iline++;
            if (!std::getline(ifs, line)) goto error_invalid_format;
            int kc, kf;
            double px, py, pz, m;
            double x, y, z, t;
            if (10 != std::sscanf(
                line.c_str(), " %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
                &kc, &kf, &px, &py, &pz, &m, &x, &y, &z, &t
              )) goto error_invalid_format;

            this->addParticleMinkowski(kf, px, py, pz, m, x, y, z, t);
            this->pcache[isample].push_back(this->plist.back());
          }
        }

        // 最終行の中身が -999 かどうかチェックする。
        //   条件: line に最終行が入っていること。
        int lastNumber;
        if (1 != std::sscanf(line.c_str(), " %d", &lastNumber) || lastNumber != -999)
          goto error_invalid_format;
        if (eventCountKnown && isample != m_numberOfSamples)
          goto error_invalid_format;
        return;
      }

    error_failed_to_open:
      std::cerr
        << "spectra/ParticleSamplePhasespace.cpp(ParticleSamplePhasespace::update): failed to open the file ("
        << this->fname_phasespace_dat << ")" << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/

    error_invalid_format:
      std::cerr
        << this->fname_phasespace_dat << ":" << iline << ": invalid format (@ParticleSampleReadPhasespaceData::update)"
        << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/
    }

    virtual void update() {
      if (this->m_currentSampleIndex < 0) {
        this->readPhasespaceDat();
        this->m_currentSampleIndex = 0;
      } else {
        this->m_currentSampleIndex++;
      }
    }
  };

  class ParticleSampleReadPhasespaceBinary: public ParticleSampleBase {
    std::string fname_phasespace_bin;
  public:
    ParticleSampleReadPhasespaceBinary(std::string const& fname_phasespace_bin):
      fname_phasespace_bin(fname_phasespace_bin)
    {
      this->m_currentSampleIndex = -1;
      this->m_numberOfSamples = 1000;
    }

  private:
    int m_numberOfSamples;
    int m_currentSampleIndex;
    std::vector<std::vector<Particle*> > pcache;
  public:
    void setNumberOfSamples(int value) {
      this->m_numberOfSamples = value;
    }
    int getNumberOfSamples() const {
      return m_numberOfSamples;
    }

  private:
    virtual ParticleIDType::value_type getParticleIdType() const { return ParticleIDType::PDGCode; }

    virtual std::vector<Particle*> const& getParticleList() const {
      if (m_currentSampleIndex < 0) {
        std::cerr << "ParticleSampleReadPhasespaceBinary: A phasespace file has not been read." << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (m_currentSampleIndex >= this->pcache.size()) {
        std::cerr << "ParticleSampleReadPhasespaceBinary: No more samples." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->pcache[m_currentSampleIndex];
    }

    void readPhasespaceBin() {
      this->clearParticleList();

      bool eventCountKnown = m_numberOfSamples >= 0;
      this->pcache.clear();
      if (eventCountKnown)
        this->pcache.reserve(this->m_numberOfSamples);

      std::ifstream ifs(this->fname_phasespace_bin.c_str(), std::ios::binary);
      int isample = 0;
      if (!ifs) goto error_failed_to_open;
      {
        std::string line;
        std::uint32_t type;
        for (isample = 0; ; isample++) {
          std::ifstream::pos_type const pos_start = ifs.tellg();

          char magic[4];
          if (!ifs.read(magic, sizeof magic)) {
            if (ifs.tellg() == pos_start || ifs.tellg() == (std::ifstream::pos_type) -1) break;
            goto error_invalid_format;
          } else if (std::memcmp(magic, "EvPh", sizeof magic) != 0) {
            goto error_invalid_format;
          }

          std::uint32_t nhadron;
          if (!ifs.read((char*) &nhadron, sizeof nhadron))
            goto error_invalid_format;

          this->pcache.resize(isample + 1);
          for (std::uint32_t i = 0; i < nhadron; i++) {
            std::int32_t ks, kf;
            float px, py, pz, m;
            float x, y, z, t;
            if (!ifs.read((char*) &ks, sizeof ks)) goto error_invalid_format;
            if (!ifs.read((char*) &kf, sizeof kf)) goto error_invalid_format;
            if (!ifs.read((char*) &px, sizeof px)) goto error_invalid_format;
            if (!ifs.read((char*) &py, sizeof py)) goto error_invalid_format;
            if (!ifs.read((char*) &pz, sizeof pz)) goto error_invalid_format;
            if (!ifs.read((char*) &m,  sizeof m )) goto error_invalid_format;
            if (!ifs.read((char*) &x, sizeof x)) goto error_invalid_format;
            if (!ifs.read((char*) &y, sizeof y)) goto error_invalid_format;
            if (!ifs.read((char*) &z, sizeof z)) goto error_invalid_format;
            if (!ifs.read((char*) &t, sizeof t)) goto error_invalid_format;
            this->addParticleMinkowski(kf, px, py, pz, m, x, y, z, t);
            this->pcache[isample].push_back(this->plist.back());
          }
        }

        if (eventCountKnown && isample != m_numberOfSamples)
          goto error_invalid_format;
        return;
      }

    error_failed_to_open:
      std::cerr
        << "spectra/ParticleSamplePhasespace (ParticleSampleReadPhasespaceBinary): failed to open the file ("
        << this->fname_phasespace_bin << ")" << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/

    error_invalid_format:
      std::ostringstream s;
      s << this->fname_phasespace_bin
        << ":#" << std::hex << ifs.tellg()
        << ":isample=" << std::dec << isample;
      std::cerr << s.str() << ": invalid format (ParticleSampleReadPhasespaceBinary)." << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/
    }

    virtual void update() {
      if (this->m_currentSampleIndex < 0) {
        this->readPhasespaceBin();
        this->m_currentSampleIndex = 0;
      } else {
        this->m_currentSampleIndex++;
      }
    }
  };

  class ParticleSampleFactory: ParticleSampleFactoryRegistered {
    virtual IParticleSample* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type == "phase1" || type == "phase") {
        ParticleSampleReadPhasespaceData* psamp = new ParticleSampleReadPhasespaceData(inputfile);
        if (type == "phase1") psamp->setNumberOfSamples(1);
        return psamp;
      } else if (type == "phbin") {
        return new ParticleSampleReadPhasespaceBinary(inputfile);
      }

      return 0;
    }
  } instance;

}
}
}
