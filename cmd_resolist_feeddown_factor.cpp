#include "config.hpp"
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <fstream>

#include "util.hpp"
#include "args.hpp"
#include "ParticleSample.hpp"
#include "ResonanceList.hpp"
#include "jamimpl.hpp"

namespace idt::runjam {

  int cmd_resolist_feeddown_factor(idt::runjam::runjam_context& ctx, idt::runjam::runjam_commandline_arguments const& args) {
    ctx.set_value("runjam_switch_weak_decay", true);
    ctx.set_value("runjam_phi_decay", true);

    class ParticleSampleResonance: public idt::runjam::ParticleSampleBase {
      std::size_t const m_count;
      idt::runjam::ResonanceRecord const* m_reso = nullptr;
    public:
      ParticleSampleResonance(std::size_t count): m_count(count) {}
      void setResonance(idt::runjam::ResonanceRecord const* reso) { this->m_reso = reso; }
      virtual double getOverSamplingFactor() const override { return 1.0; }
      virtual void update() override {
        if (!m_reso) {
          std::cerr << "FATAL runjam (resolist-feeddown-factor/ParticleSampleResonance): m_reso uinitialized" << std::endl;
          std::exit(9);
        }

        this->clearParticleList();
        for (std::size_t i = 0; i < m_count; i++) {
          int const pdg = m_reso->pdg_codes[idt::util::irand(m_reso->pdg_codes.size())];
          double const rx = 10.0 * idt::util::nrand();
          double const ry = 10.0 * idt::util::nrand();
          double const rz = 10.0 * idt::util::nrand();
          this->addParticleCartesian(pdg, 0.0, 0.0, 0.0, m_reso->mass, rx, ry, rz, 0.0);
        }
      }
    };

    std::unique_ptr<idt::runjam::IJamRunner> runner = idt::runjam::create_runner(ctx);
    if (!runner) {
      std::cerr << "runjam: JAM unsupported." << std::endl;
      std::exit(3);
    }
    runner->initialize(ctx, "decay");

    std::vector<idt::runjam::Particle> final_state;

    static std::size_t const jkbin_size = 1000;
    static std::size_t const jkbin_count = 1000;

    struct jackknife {
      std::size_t count = 0;
      void operator++(int) { count++;}

      std::vector<double> bins;
      void next_bin() {
        bins.push_back((double) count / jkbin_size);
        count = 0;
      }

      double value, error;
      void summarize() {
        double const total = std::accumulate(bins.begin(), bins.end(), 0.0);
        double s1 = 0.0, s2 = 0.0;
        for (auto const& v: bins) {
          double const v1 = (total - v) / (bins.size() - 1);
          s1 += v1;
          s2 += v1 * v1;
        }
        s1 /= bins.size();
        s2 /= bins.size();
        value = total / bins.size();
        error = std::sqrt((bins.size() - 1) * (s2 - s1 * s1));
      }
    };

    ParticleSampleResonance psamp(jkbin_size);
    idt::runjam::ResonanceList rlist(ctx);
    std::ofstream file("feeddown.txt");
    file << "#KEY PDGCODES Nch pi+ pi- pi+-avg K0 Kp Kn p pbar\n";
    for (idt::runjam::ResonanceRecord const& reso: rlist) {
      std::cerr << reso.key << "..." << std::endl;
      psamp.setResonance(&reso);

      jackknife jklist[9];
      int jkindex = 0;
      jackknife& nch   = jklist[jkindex++];
      jackknife& npip  = jklist[jkindex++];
      jackknife& npin  = jklist[jkindex++];
      jackknife& npic  = jklist[jkindex++];
      jackknife& nKp   = jklist[jkindex++];
      jackknife& nKn   = jklist[jkindex++];
      jackknife& nK0   = jklist[jkindex++];
      jackknife& np    = jklist[jkindex++];
      jackknife& npbar = jklist[jkindex++];
      // jackknife& npi0  = jklist[jkindex++];
      // jackknife& nSig  = jklist[jkindex++];
      // jackknife& nLam  = jklist[jkindex++];

      for (std::size_t i = 0; i < jkbin_count; i++) {
        psamp.update();
        runner->do_decay(psamp, final_state);
        for (idt::runjam::Particle& part: final_state) {
          int const pdg = part.pdg;
          switch (std::abs(pdg)) {
          case 211:  nch++; (pdg > 0 ? npip : npin)++; npic++; break;
          case 321:  nch++; (pdg > 0 ? nKp : nKn)++; break;
          case 2212: nch++; (pdg > 0 ? np : npbar)++; break;
          case 311: nK0++; break;
          // case 130: case 310: nK0++; break; // K0S, K0L
          // case 111:  npi0++; break; (decays)
          // case 3212: nLam++; break; // Lambda (decays)
          // case 3112: case 3222: nch++; nSig++; break; // charged Sigma (decays)
          }
        }
        for (auto& jk: jklist) jk.next_bin();
      }
      for (auto& jk: jklist) jk.summarize();

      file << reso.key << " ";
      bool first = true;
      for (int pdg: reso.pdg_codes) {
        if (first) first = false; else file << ",";
        file << pdg;
      }
      file << " ";
      file << reso.deg * nch.value        << "(" << reso.deg * nch.error << ") "        // N_ch
           << reso.deg * npip.value       << "(" << reso.deg * npip.error << ") "       // pi+
           << reso.deg * npin.value       << "(" << reso.deg * npin.error << ") "       // pi-
           << reso.deg * npic.value * 0.5 << "(" << reso.deg * npic.error * 0.5 << ") " // (pi+ + pi-)/2
           << reso.deg * nK0.value        << "(" << reso.deg * nK0.error << ") "        // K0 + K0bar
           << reso.deg * nKp.value        << "(" << reso.deg * nKp.error << ") "        // K+
           << reso.deg * nKn.value        << "(" << reso.deg * nKn.error << ") "        // K-
           << reso.deg * np.value         << "(" << reso.deg * np.error << ") "         // p
           << reso.deg * npbar.value      << "(" << reso.deg * npbar.error << ") "      // pbar
           << std::endl;
    }

    return 0;
  }

}
