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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include "ParticleSample.hpp"

namespace idt {
namespace runjam {
namespace {

  class ParticleSampleReadPhasespaceData: public ParticleSampleBase {
    typedef ParticleSampleBase base;
    std::string fname_phasespace_dat;

  private:
    double m_overSamplingFactor; // イベントの一行目の $2 から読み取る
    int    m_numberOfSamples;
    int    m_currentSampleIndex;
  public:
    void setNumberOfSamples(int value) {
      this->m_numberOfSamples = value;
    }
    int getNumberOfSamples() const {
      return m_numberOfSamples;
    }
    virtual double getOverSamplingFactor() const override {
      return m_overSamplingFactor;
    }

  public:
    ParticleSampleReadPhasespaceData(runjam_context const& ctx, std::string const& fname_phasespace_dat):
      fname_phasespace_dat(fname_phasespace_dat),
      m_overSamplingFactor(ctx.get_config("runjam_oversampling_factor", 1.0))
    {
      this->m_currentSampleIndex = -1;
      this->m_numberOfSamples = -1;
    }

  private:
    std::vector<std::size_t> pcache;
    void checkCacheAvailability() {
      if (m_currentSampleIndex < 0) {
        std::cerr << "ParticleSampleReadPhasespaceData: A phasespace file has not been read." << std::endl;
        std::exit(EXIT_FAILURE);
      } else if ((std::size_t) m_currentSampleIndex >= this->pcache.size() - 1) {
        std::cerr << "ParticleSampleReadPhasespaceData: No more samples." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  public:
    virtual Particle* begin() override {
      this->checkCacheAvailability();
      return base::begin() + pcache[m_currentSampleIndex];
    }
    virtual Particle* end() override {
      this->checkCacheAvailability();
      return base::begin() + pcache[m_currentSampleIndex + 1];
    }

  private:
    void readPhasespaceDat() {
      this->clearParticleList();

      bool eventCountKnown = m_numberOfSamples >= 0;
      this->pcache.clear();
      this->pcache.emplace_back(0);
      if (eventCountKnown)
        this->pcache.reserve(this->m_numberOfSamples + 1);

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
          int npart;
          if (2 != std::sscanf(line.c_str(), " %d %lf", &npart, &m_overSamplingFactor)) break;

          if (isample == this->m_numberOfSamples) goto error_invalid_format;

          if (npart < 0) goto error_invalid_format;

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

            this->addParticleCartesian(kf, px, py, pz, m, x, y, z, t);
          }
          this->pcache.emplace_back(base::plist.size());
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
        << "ParticleSamplePhasespace.cpp(ParticleSamplePhasespace::update): failed to open the file ("
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
    typedef ParticleSampleBase base;
    std::string fname_phasespace_bin;

  private:
    double m_overSamplingFactor;
    int   m_numberOfSamples;
    int   m_currentSampleIndex;
  public:
    void setNumberOfSamples(int value) {
      this->m_numberOfSamples = value;
    }
    int getNumberOfSamples() const {
      return m_numberOfSamples;
    }
    // この枠組ではファイルに保存されていた粒子集合の oversampling が
    // 分からないので、環境に設定されていた値を信じてそのまま返す。
    virtual double getOverSamplingFactor() const override {
      return m_overSamplingFactor;
    }

  public:
    ParticleSampleReadPhasespaceBinary(runjam_context const& ctx, std::string const& fname_phasespace_bin):
      fname_phasespace_bin(fname_phasespace_bin),
      m_overSamplingFactor(ctx.get_config("runjam_oversampling_factor", 1.0))
    {
      this->m_currentSampleIndex = -1;
      this->m_numberOfSamples = 1000;
    }

  private:
    std::vector<std::size_t> pcache;
    void checkCacheAvailability() {
      if (m_currentSampleIndex < 0) {
        std::cerr << "ParticleSampleReadPhasespaceBinary: A phasespace file has not been read." << std::endl;
        std::exit(EXIT_FAILURE);
      } else if ((std::size_t) m_currentSampleIndex >= this->pcache.size() - 1) {
        std::cerr << "ParticleSampleReadPhasespaceBinary: No more samples." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  public:
    virtual Particle* begin() override {
      this->checkCacheAvailability();
      return base::begin() + pcache[m_currentSampleIndex];
    }
    virtual Particle* end() override {
      this->checkCacheAvailability();
      return base::begin() + pcache[m_currentSampleIndex + 1];
    }

  private:
    void readPhasespaceBin() {
      this->clearParticleList();

      bool eventCountKnown = m_numberOfSamples >= 0;
      this->pcache.clear();
      this->pcache.emplace_back();
      if (eventCountKnown)
        this->pcache.reserve(this->m_numberOfSamples + 1);

      std::ifstream ifs(this->fname_phasespace_bin.c_str(), std::ios::binary);
      int isample = 0;
      if (!ifs) goto error_failed_to_open;
      {
        std::string line;
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
            this->addParticleCartesian(kf, px, py, pz, m, x, y, z, t);
          }
          this->pcache.emplace_back(base::plist.size());
        }

        if (eventCountKnown && isample != m_numberOfSamples)
          goto error_invalid_format;
        return;
      }

    error_failed_to_open:
      std::cerr
        << "ParticleSamplePhasespace (ParticleSampleReadPhasespaceBinary): failed to open the file ("
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

  class ParticleSampleFactory: ParticleSampleFactoryBase {
    virtual std::unique_ptr<ParticleSampleBase> CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type == "phase1" || type == "phase") {
        ParticleSampleReadPhasespaceData* psamp = new ParticleSampleReadPhasespaceData(ctx, inputfile);
        if (type == "phase1") psamp->setNumberOfSamples(1);
        return std::unique_ptr<ParticleSampleBase>(psamp);
      } else if (type == "phbin") {
        auto psamp = std::make_unique<ParticleSampleReadPhasespaceBinary>(ctx, inputfile);
        psamp->setNumberOfSamples(ctx.get_config("runjam_initial_phbin_size", 1000));
        return psamp;
      }

      return nullptr;
    }
  } instance;

}
}
}
