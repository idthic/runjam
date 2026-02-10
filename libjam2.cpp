/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2022-2025, Koichi Murase <myoga.murase at gmail.com>

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

#include "config.hpp"

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>

#include <jam2/JAM.h>
#include "libjam2.hpp"

namespace {

  class initial_condition_adapter: public jam2::UserInitialCondition {
    typedef jam2::UserInitialCondition base;
    static constexpr double DEFAULT_MAX_CM_ENERGY = 20.0;

  public:
    std::vector<idt::runjam::Particle> const* m_particles = nullptr;
    std::unordered_map<int, int> m_new_pdg;

  private:
    // currently unused.
    void determine_max_cm_energy_from_particles(std::vector<idt::runjam::Particle> const& initial_state) {
      // determine maximum invariant mass
      double max_mass2 = 0.0;
      for (std::size_t i = 0; i < initial_state.size(); i++) {
        for (std::size_t j = i + 1; j < initial_state.size(); j++) {
          double const mom0 = initial_state[i].mom[0] + initial_state[j].mom[0];
          double const mom1 = initial_state[i].mom[1] + initial_state[j].mom[1];
          double const mom2 = initial_state[i].mom[2] + initial_state[j].mom[2];
          double const mom3 = initial_state[i].mom[3] + initial_state[j].mom[3];
          double const mass2 = mom0 * mom0 - mom1 * mom1 - mom2 * mom2 - mom3 * mom3;
          if (mass2 > max_mass2) max_mass2 = mass2;
        }
      }

      // [GeV]
      double const max_cm_energy = std::max(DEFAULT_MAX_CM_ENERGY, 5.0 * std::sqrt(max_mass2));
      settings->parm("Beams:eCM", max_cm_energy);
    }

  public:
    initial_condition_adapter(idt::runjam::runjam_context const& ctx, jam2::JAM* jam):
      base(&jam->info, jam->jamParticleData)
    {
      std::string const file = ctx.lookup_data_file("libjam2.pdg_new.txt");
      if (file.empty()) {
        std::fprintf(stderr, "runjam: the data file '%s/libjam2.pdg_new.txt' is not found.\n", ctx.datadir());
        std::exit(EXIT_FAILURE);
      }

      std::ifstream str(file.c_str());
      std::string line;
      std::istringstream istr;
      while (std::getline(str, line)) {
        istr.str(line);
        istr.clear();
        int pdg1, pdg2;
        if (istr >> pdg1 >> pdg2) {
          m_new_pdg.insert(std::make_pair(pdg1, pdg2));
          m_new_pdg.insert(std::make_pair(-pdg1, -pdg2));
        }
      }
    }

    virtual ~initial_condition_adapter() {}

    virtual void init() override {
      // We have been initially determining "Beams:eCM" based on the list of
      // particles using `determine_max_cm_energy_from_particles`, but we now
      // cannot do in that way because we need to initialize Pythia before
      // seeing the initial particles for the multiple events.  We instead just
      // set the fixed value.
      settings->parm("Beams:eCM", DEFAULT_MAX_CM_ENERGY);
    }

    // Convert a PDG Monte-Carlo code used in JAM 1 or the past version of JAM
    // 2 to the latest PDG Monte-Carlo code used in the latest version of JAM
    // 2.7110.
    int get_jam2_pdg(int pdg) const {
      if (auto it = m_new_pdg.find(pdg); it != m_new_pdg.end())
        return it->second;
      return pdg;
    }

    virtual void generate(jam2::Collision* event, int mode = 0) override {
      (void) mode;

#if LIBJAM2_VERSION_ID >= 200080000
      int const optCollisionOrdering = this->fixtime->collisionOrdering();
#endif

      for (idt::runjam::Particle const& particle: *m_particles) {
        int const pdg = get_jam2_pdg(particle.pdg);

        // find this particle in the JAM particle list.
        Pythia8::ParticleDataEntryPtr const pa = jamParticleData->find(pdg);
        Pythia8::Vec4 const r(particle.pos[1], particle.pos[2], particle.pos[3], particle.pos[0]);
        Pythia8::Vec4 const p(particle.mom[1], particle.mom[2], particle.mom[3], particle.mom[0]);

        auto const cp = new jam2::EventParticle(pdg, particle.mass, r, p, pa);
        cp->setPID(jamParticleData->pid(std::abs(pdg)));

#if LIBJAM2_VERSION_ID >= 200080000
        // set evolution time
        if ((10 <= optCollisionOrdering && optCollisionOrdering < 20) || optCollisionOrdering == 21) {
          // If we are in the tau-eta coordinates for jam2, and the particle is
          // outside the light-cone, we ignore the particle.
          double const tau2 = r.e() * r.e() - r.pz() * r.pz();
          if (tau2 < 1e-10) continue;
          double const tau_dummy = std::sqrt(tau2);
          double const lambda_dummy = (p.e() * r.e() - p.pz() * r.pz()) / (p.e() * tau_dummy);
          this->fixtime->initializeEvolutionTime(cp, tau_dummy, lambda_dummy, 100);
        } else if (optCollisionOrdering < 110) {
          // Note: We here assume that pHat = (1, 0, 0, 0).
          this->fixtime->initializeEvolutionTime(cp, r.e(), 1.0 / p.e(), 100);
        } else {
          std::fprintf(stderr, "error(libjam2): unsupported collisionOrdering=%d\n", optCollisionOrdering);
          std::exit(1);
        }

        // set decay time
        cp->sampleLifetime(jamParticleData, this->fixtime);
#else

        // compute decay time if it is resonance.
        double const decay_time = jamParticleData->lifeTime(pa, particle.mass, particle.mom[0]);
        cp->setLifeTime(particle.pos[0] + decay_time);
#endif

        // put this particle into the particle list.
        event->setPList(cp);
      }
    }
  };

}

namespace libjam2 {

  class runner: public irunner {
    jam2::JAM m_jam;
    std::unique_ptr<initial_condition_adapter> m_initial_condition_tmp;
    initial_condition_adapter* m_initial_condition;

  public:
    runner(idt::runjam::runjam_context const& ctx, std::string const& input_filename):
      m_jam(Pythia8_PREFIX "/share/Pythia8/xmldoc", false)
    {
      m_jam.readFile(input_filename);
      m_initial_condition_tmp = std::make_unique<initial_condition_adapter>(ctx, &m_jam);
      m_initial_condition = m_initial_condition_tmp.get();

      Pythia8::Settings* const settings = this->settings();
      settings->mode("Cascade:model", 3);

      // Set dummy value 12 to Cascade:initialCondition. An old JAM2 required
      // to specify a value different from the existing ones to enable the
      // initial condition specified through the argument of "m_jam.init()".
      // The latest version of JAM >= 2.7103 does not have this restriction.
      settings->mode("Cascade:initialCondition", 12);
    }

    Pythia8::Settings* settings() override {
      return m_jam.settings;
    }

    Pythia8::Settings const* settings() const override {
      return m_jam.settings;
    }

    void initialize() override {
      // This moves the ownership into m_jam.
      m_jam.init(m_initial_condition_tmp.release());
    }

    void run(std::vector<idt::runjam::Particle> const& initial_state, std::vector<idt::runjam::Particle>& final_state) override {
      m_initial_condition->m_particles = &initial_state;
      m_jam.next();

      final_state.clear();
      for (jam2::EventParticle* p: m_jam.getEvent()) {
        final_state.emplace_back();
        auto& particle = final_state.back();
        particle.pdg = p->getID();
        particle.mass = p->getMass();
        particle.pos[0] = p->TimeLastColl();
        particle.pos[1] = p->getV(1);
        particle.pos[2] = p->getV(2);
        particle.pos[3] = p->getV(3);
        particle.mom[0] = p->getPe();
        particle.mom[1] = p->getPx();
        particle.mom[2] = p->getPy();
        particle.mom[3] = p->getPz();
      }
    }

    double get_particle_mass(int pdg) const override {
      return m_jam.jamParticleData->find(pdg)->m0();
    }
    int get_particle_stable_code(int pdg) const override {
      // Returns 1 for stable particle and 2 for unstable particle.
      Pythia8::ParticleDataEntryPtr const pd = m_jam.jamParticleData->find(m_initial_condition->get_jam2_pdg(pdg));
      if (!pd) {
        std::cerr << "runjam (libjam2::get_particle_stable_code): particle data for pdg=" << pdg << " not found" << std::endl;
        std::exit(1);
      }
      return pd->mWidth() <= 1e-7 ? 1 : 2;
    }
    double get_event_collision_number() const override {
      return m_jam.getNColl();
    }
  };

  std::unique_ptr<irunner> create_runner(idt::runjam::runjam_context const& ctx, std::string const& input_filename) {
    return std::make_unique<runner>(ctx, input_filename);
  }

  std::string version_string() {
#ifdef LIBJAM2_VERSION
    return "libjam2 " LIBJAM2_VERSION;
#else
    return "libjam2";
#endif
  }
}
