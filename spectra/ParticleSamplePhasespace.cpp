#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "IParticleSample.hpp"

namespace idt {
namespace hydro2jam {
namespace {

  class ParticleSampleFromOversampledPhasespace:public ParticleSampleBase{
    std::string fname_phasespace_dat;
  public:
    ParticleSampleFromOversampledPhasespace(std::string const& fname_phasespace_dat)
      :fname_phasespace_dat(fname_phasespace_dat)
    {
      this->m_currentSampleIndex = -1;
      this->m_overSamplingFactor = -1;
    }

  private:
    int m_overSamplingFactor;
    int m_currentSampleIndex;
    std::vector<std::vector<Particle*> > pcache;
  public:
    void setOverSamplingFactor(int value) {
      this->m_overSamplingFactor = value;
    }
    int getOverSamplingFactor() const {
      return m_overSamplingFactor;
    }

  private:
    virtual ParticleIDType::value_type getParticleIdType() const { return ParticleIDType::PDGCode; }

    virtual std::vector<Particle*> const& getParticleList() const{
      if (m_currentSampleIndex < 0) {
        std::cerr<<"ParticleSampleFromOversampledPhasespace: A phasespace file has not been read."<<std::endl;
        std::exit(EXIT_FAILURE);
      } else if (m_currentSampleIndex >= this->pcache.size()) {
        std::cerr<<"ParticleSampleFromOversampledPhasespace: No more samples."<<std::endl;
        std::exit(EXIT_FAILURE);
      }

      return this->pcache[m_currentSampleIndex];
    }

    void readPhasespaceDat() {
      this->clearParticleList();

      bool eventCountKnown = m_overSamplingFactor >= 0;
      this->pcache.clear();
      if (eventCountKnown)
        this->pcache.reserve(this->m_overSamplingFactor);

      int iline;
      {
        std::ifstream ifs(this->fname_phasespace_dat.c_str());
        if (!ifs) goto error_failed_to_open;

        int isample;
        std::string line;
        iline = 0;
        for (isample = 0; ; isample++) {
          iline++;
          if (!std::getline(ifs,line)) goto error_invalid_format;
          int npart,dummy;
          if (2 != std::sscanf(line.c_str()," %d %d",&npart,&dummy)) break;

          if (isample == this->m_overSamplingFactor) goto error_invalid_format;

          if (npart < 0) goto error_invalid_format;

          this->pcache.resize(isample + 1);
          for (int i = 0; i < npart; i++) {
            iline++;
            if (!std::getline(ifs, line)) goto error_invalid_format;
            int kc, kf;
            double px, py, pz, m;
            double x, y, z, t;
            if (10 != std::sscanf(
                line.c_str()," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
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
        if (eventCountKnown && isample != m_overSamplingFactor)
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
        << this->fname_phasespace_dat << ":" << iline << ": invalid format (@ParticleSampleFromOversampledPhasespace::update)"
        << std::endl;
      std::exit(EXIT_FAILURE);
      return; /*NOTREACHED*/
    }

    virtual void update() {
      if (this->m_currentSampleIndex<0){
        this->readPhasespaceDat();
        this->m_currentSampleIndex = 0;
      } else {
        this->m_currentSampleIndex++;
      }
    }
  };

  class ParticleSampleFactory: ParticleSampleFactoryRegistered {
    virtual IParticleSample* CreateInstance(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type == "phase1" || type == "phase") {
        ParticleSampleFromOversampledPhasespace* psamp = new ParticleSampleFromOversampledPhasespace(inputfile);
        if (type == "phase1") psamp->setOverSamplingFactor(1);
        return psamp;
      }

      return 0;
    }
  } instance;

}
}
}
