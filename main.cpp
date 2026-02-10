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
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <chrono>

#include "config.hpp"
#include "args.hpp"
#include "ParticleSample.hpp"
#include "ParticleSampleViscous.hpp"
#include "jamimpl.hpp"

namespace idt::runjam {
  int cmd_resolist_feeddown_factor(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args);
  int cmd_resolist_eos(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args);
}

using namespace idt;
using namespace idt::runjam;

void writePhasespaceData(std::ofstream& ofs, Particle const* begin, Particle const* end, double ntest = 1.0) {
  std::size_t const nv = end - begin;
  ofs << nv << "  " << ntest << std::endl;
  while (begin != end) {
    Particle const& part = *begin++;
    int const ks = runjam::getParticleStableCode(part.pdg);
    int const kf = part.pdg;
    double const px = part.px;
    double const py = part.py;
    double const pz = part.pz;
    double const m  = part.mass;
    double const x  = part.x;
    double const y  = part.y;
    double const z  = part.z;
    double const t  = part.t;
    ofs << std::setw(5)  << ks
        << std::setw(10) << kf
        << std::setw(14) << px
        << std::setw(14) << py
        << std::setw(14) << pz
        << std::setw(14) << m
        << std::setw(14) << x
        << std::setw(14) << y
        << std::setw(14) << z
        << std::setw(14) << t
        << std::endl;
  }
  ofs << std::flush;
}

static bool is_charged(int pdg) {
  return std::abs(pdg) == 211 ||
    std::abs(pdg) == 321 ||
    std::abs(pdg) == 2212;
}

enum phasespace_phbin_type {
  phasespace_phbin_full = 0,
  phasespace_phbin4 = 1,
  phasespace_phbin4_charged = 2,
};

void writePhasespaceBinary(std::ofstream& ofs, Particle const* begin, Particle const* end, phasespace_phbin_type type, double ntest = 1.0) {
  (void) ntest;
  std::uint32_t const nv = end - begin;

  switch (type) {
  case phasespace_phbin4:
    {
      ofs.write("EvPp", 4);
      ofs.write((char*) &nv, 4);
      while (begin != end) {
        Particle const& part = *begin++;
        std::int32_t const pdg = part.pdg;
        float const px = part.px;
        float const py = part.py;
        float const pz = part.pz;
        ofs.write((char*) &pdg, 4);
        ofs.write((char*) &px, 4);
        ofs.write((char*) &py, 4);
        ofs.write((char*) &pz, 4);
      }
      ofs << std::flush;
    }
    break;
  case phasespace_phbin4_charged:
    {
      std::vector<Particle const*> list;
      while (begin != end) {
        Particle const& part = *begin++;
        if (is_charged(part.pdg)) list.push_back(&part);
      }
      std::uint32_t const number_of_particles = list.size();

      ofs.write("EvPp", 4);
      ofs.write((char*) &number_of_particles, 4);
      for (auto const* part: list) {
        std::int32_t const pdg = part->pdg;
        float const px = part->px;
        float const py = part->py;
        float const pz = part->pz;
        ofs.write((char*) &pdg, 4);
        ofs.write((char*) &px, 4);
        ofs.write((char*) &py, 4);
        ofs.write((char*) &pz, 4);
      }
      ofs << std::flush;
    }
    break;
  default:
    {
      ofs.write("EvPh", 4);
      ofs.write((char*) &nv, 4);
      while (begin != end) {
        Particle const& part = *begin++;
        std::size_t const ks = runjam::getParticleStableCode(part.pdg);
        std::size_t const kf = part.pdg;
        float const px = part.px;
        float const py = part.py;
        float const pz = part.pz;
        float const m  = part.mass;
        float const x  = part.x;
        float const y  = part.y;
        float const z  = part.z;
        float const t  = part.t;
        ofs.write((char*) &ks, 4);
        ofs.write((char*) &kf, 4);
        ofs.write((char*) &px, 4);
        ofs.write((char*) &py, 4);
        ofs.write((char*) &pz, 4);
        ofs.write((char*) &m, 4);
        ofs.write((char*) &x, 4);
        ofs.write((char*) &y, 4);
        ofs.write((char*) &z, 4);
        ofs.write((char*) &t, 4);
      }
      ofs << std::flush;
    }
    break;
  }
}

class IObserver {
public:
  virtual void initialize() = 0;
  virtual void process_event(Particle const* begin, Particle const* end, double ntest) = 0;
  virtual void finalize() = 0;
  virtual ~IObserver() {}

public:
  void process_event(std::vector<Particle> const& list, double ntest) {
    return process_event(list.data(), list.data() + list.size(), ntest);
  }
};

class FileWriterBase: public IObserver {
protected:
  std::string filename;
  std::string filename_part;
  std::ios::openmode openmode = std::ios::out;
  std::ofstream ofs;
  FileWriterBase(std::string const& filename, std::ios::openmode openmode = std::ios::out):
    filename(filename), filename_part(filename + ".part"), openmode(openmode) {}
public:
  virtual void initialize() override {
    ofs.open(filename_part.c_str(), openmode);
    if (!ofs) {
      std::cerr << "runjam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(1);
    }
  }
  virtual void finalize() override {
    if (ofs.is_open()) {
      ofs.close();
      if (std::rename(filename_part.c_str(), filename.c_str()) != 0) {
        std::cerr << "runjam: failed to rename '" << filename_part << "' to '" << filename
                  << " (" << std::strerror(errno) << ")" << std::endl;
      }
    }
  }
};

struct PhasespaceDataWriter: public FileWriterBase {
  typedef FileWriterBase base;
  PhasespaceDataWriter(std::string const& filename): base(filename) {}
  virtual void process_event(Particle const* begin, Particle const* end, double ntest) override {
    writePhasespaceData(ofs, begin, end, ntest);
  }
  virtual void finalize() override {
    if (ofs.is_open())
      ofs << -999 << std::endl;
    base::finalize();
  }
};

struct PhasespaceBinaryWriter: public FileWriterBase {
  typedef FileWriterBase base;

private:
  phasespace_phbin_type m_type = phasespace_phbin_full;

public:
  PhasespaceBinaryWriter(std::string const& filename):
    base(filename, std::ios::out | std::ios::binary)
  {
    if (idt::util::ends_with(filename, ".bin4"))
      m_type = phasespace_phbin4;
    else if (idt::util::ends_with(filename, ".bin4ch"))
      m_type = phasespace_phbin4_charged;
  }
  virtual void process_event(Particle const* begin, Particle const* end, double ntest) override {
    writePhasespaceBinary(ofs, begin, end, m_type, ntest);
  }
};

class IndexedPhasespaceDataWriter: public IObserver {
  std::string outdir, suffix;
  std::size_t index;
  std::vector<char> filename;
public:
  IndexedPhasespaceDataWriter(std::string const& outdir, std::string const& suffix, std::size_t startIndex):
    outdir(outdir), suffix(suffix), index(startIndex) {}
  virtual void initialize() override {}
  virtual void process_event(Particle const* begin, Particle const* end, double ntest) {
    filename.resize(outdir.size() + suffix.size() + 50);
    std::sprintf(filename.data(), "%s/dens%06zd%s", outdir.c_str(), index++, suffix.c_str());
    std::ofstream ofs(filename.data());
    writePhasespaceData(ofs, begin, end, ntest);
    ofs << -999 << std::endl;
  }
  virtual void finalize() override {}
};

class Program {
private:
  int cfg_nevent;
  std::string cfg_outdir;
  int cfg_jam_version;


  std::vector<std::unique_ptr<IObserver>> onbefore;
  std::vector<std::unique_ptr<IObserver>> onafter;

public:
  void initialize_parameters(runjam_context const& ctx) {
    cfg_nevent = ctx.nevent(1);
    cfg_outdir = ctx.outdir();
    cfg_jam_version = ctx.get_config("runjam_jam_version", DEFAULT_JAM_VERSION);

    onbefore.clear();
    onafter.clear();
    if (ctx.get_config("runjam_output_phdat0", true)) {
      std::string filename = cfg_outdir + "/phasespace0.dat";
      ctx.read_config(filename, "runjam_fname_phdat0");
      onbefore.emplace_back(std::make_unique<PhasespaceDataWriter>(filename));
    }
    if (ctx.get_config("runjam_output_phdat", true)) {
      std::string filename = cfg_outdir + "/phasespace.dat";
      ctx.read_config(filename, "runjam_fname_phdat");
      onafter.emplace_back(std::make_unique<PhasespaceDataWriter>(filename));
    }
    if (ctx.get_config("runjam_output_phbin0", false)) {
      std::string filename = cfg_outdir + "/phasespace0.bin";
      ctx.read_config(filename, "runjam_fname_phbin0");
      onbefore.emplace_back(std::make_unique<PhasespaceBinaryWriter>(filename));
    }
    if (ctx.get_config("runjam_output_phbin", false)) {
      std::string filename = cfg_outdir + "/phasespace.bin";
      ctx.read_config(filename, "runjam_fname_phbin");
      onafter.emplace_back(std::make_unique<PhasespaceBinaryWriter>(filename));
    }
    if (ctx.get_config("runjam_output_phdat0_indexed", false)) {
      int const startIndex = ctx.get_config("runjam_output_index_start", 0);
      std::string suffix = "_phasespace0.dat";
      onbefore.emplace_back(std::make_unique<IndexedPhasespaceDataWriter>(cfg_outdir, suffix, startIndex));
    }
    if (ctx.get_config("runjam_output_phdat_indexed", false)) {
      int const startIndex = ctx.get_config("runjam_output_index_start", 0);
      std::string suffix = "_phasespace.dat";
      onafter.emplace_back(std::make_unique<IndexedPhasespaceDataWriter>(cfg_outdir, suffix, startIndex));
    }
  }

public:
  void  run(runjam_context const& ctx, ParticleSampleBase& psamp, std::string const& cascadeMode) {
    this->initialize_parameters(ctx);

    int const nprint = 1;

    std::unique_ptr<IJamRunner> runner;
    if (cascadeMode != "sample") {
      runner = create_runner(ctx);
      if (!runner) {
        std::cerr << "runjam: " << cascadeMode << " not supported." << std::endl;
        std::exit(3);
      }
    }

    bool const is_decay = cascadeMode == "decay";
    int const nevent = this->cfg_nevent;
    if (nevent > 0)
      psamp.setAdviceNumberOfExpectedEvents(nevent);

    std::vector<Particle> final_state;
    double averageParticleNumber0 = 0.0;
    double averageParticleNumber = 0.0;

    // initialize
    for (auto const& observer: onbefore)
      observer->initialize();
    if (runner) {
      runner->initialize(ctx, cascadeMode);
      for (auto const& observer: onafter)
        observer->initialize();
    }

    typedef std::chrono::steady_clock clock_t;

    // event loop
    for (int iev = 1; iev <= nevent; iev++) {
      psamp.update();

      if (runner) runner->adjust_mass(psamp);
      //psamp->adjustCenterOfMassByLorentzBoost();

      double const ntest = psamp.getOverSamplingFactor();
      averageParticleNumber0 += (double) psamp.size() / ntest;
      for (auto const& observer: onbefore)
        observer->process_event(psamp.begin(), psamp.end(), ntest);

      if (!runner) continue;

      if (iev % nprint == 0) {
        std::cout
          << "runjam:iev=" << iev << ": " << cascadeMode
          << " start (initial_nparticle=" << psamp.size() << ")" << std::endl;
      }

      clock_t::time_point const time1 = clock_t::now();
      if (is_decay)
        runner->do_decay(psamp, final_state);
      else
        runner->do_cascade(psamp, final_state, iev);
      clock_t::time_point const time2 = clock_t::now();

      if (iev % nprint == 0) {
        double const ata = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() * 0.001;
        std::cout
          << "runjam:iev=" << iev << ": " << cascadeMode
          << " done (final_nparticle=" << final_state.size();
        if (!is_decay)
          std::cout << " average_ncollision=" << runner->get_average_collision_number();
        std::cout << ") " << ata << "sec" << std::endl;
      }

      averageParticleNumber += (double) final_state.size() / ntest;
      for (auto const& observer: onafter)
        observer->process_event(final_state, ntest);
    }

    // finalize
    averageParticleNumber0 /= nevent;
    averageParticleNumber /= nevent;
    std::cout
      << "runjam: the average number of hadrons are "
      << averageParticleNumber0 << " (before JAM) -> "
      << averageParticleNumber << " (after JAM)" << std::endl;
    for (auto const& observer: onbefore)
      observer->finalize();
    if (runner) {
      for (auto const& observer: onafter)
        observer->finalize();
      runner->finalize();
    }
  }
};

int main(int argc, char *argv[]) {
  runjam_context ctx;

  runjam_commandline_arguments args;
  int const ext = args.read(argc, argv, ctx);
  if (ext) return ext;

  std::cout << "runjam [version " << PACKAGE_VERSION << idt::runjam::package_hash;
  for (std::string const& version: runjam::libjam_versions())
    std::cout << ", " << version;
  std::cout << ", seed = " << ctx.seed() << "]" << std::endl;

  idt::util::set_random_seed(ctx.seed());

  if (args.subcommand == "cascade" || args.subcommand == "decay" || args.subcommand == "sample") {
    std::unique_ptr<ParticleSampleBase> psamp = CreateParticleSample(ctx, args.initType, args.initPath);
    if (!psamp) {
      std::cerr << "runjam: failed to initialize ParticleSample of type '" << args.initType << "'." <<  std::endl;
      std::exit(1);
    }

    Program prog;
    prog.run(ctx, *psamp, args.subcommand);
  } else if (args.subcommand == "test-viscous-correction-integration") {
    return checkViscousCooperFryeInterpolated(true);
  } else if (args.subcommand == "resolist-print-ks-and-mass") {
    ResonanceList list(ctx);
    std::cout << std::setprecision(8);
    for (auto const& reso : list) {
      for (int pdg : reso.pdg_codes) {
        if (pdg > 0)
          std::cout << pdg << " " << runjam::getParticleStableCode(pdg) << " " << runjam::getParticleMass(pdg) << std::endl;
      }
    }
    return 0;
  } else if (args.subcommand == "resolist-feeddown-factor") {
    return cmd_resolist_feeddown_factor(ctx, args);
  } else if (args.subcommand == "resolist-eos") {
    return cmd_resolist_eos(ctx, args);
  } else {
    std::cerr << "runjam: unknown subcommand ' " << args.subcommand << "'" << std::endl;
    return 2;
  }

  return 0;
}
