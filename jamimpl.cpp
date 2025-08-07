#include "config.hpp"
#include "jamimpl.hpp"
#include <cstdlib>
#include <vector>
#include <string>
#include <memory>

#ifdef USE_LIBJAM2
# include <sstream>
# include <fstream>
# include <Pythia8/Settings.h>
# include <libjam2.hpp>
#endif

#ifdef USE_LIBJAM1
# include "ParticleSample.hpp"
# include "libjam1.hpp"
#endif

namespace idt {
namespace runjam {

#ifdef USE_LIBJAM2
  static libjam2::irunner* runner_instance = nullptr;
#endif

  int getParticleStableCode(int kf) {
#ifdef USE_LIBJAM2
    if (runner_instance)
      return runner_instance->get_particle_stable_code(kf);
#endif

#if defined(USE_LIBJAM1)
    return libjam1::determineStableCode(kf);
#else
    switch (std::abs(kf)) {
    case 111: case 211: case 221: // pions
    case 311: case 321: // kaons
    case 2212: case 2112: // nucleons
    case 3112: case 3122: case 3212: case 3222: // Lambda and Sigma
    case 3312: case 3322: case 3334: // Xi and Omega
      return 1;
    default:
      return 2;
    }
#endif
  }

#if defined(USE_LIBJAM2) || defined(USE_LIBJAM1)
  double getParticleMass(int kf) {
# ifdef USE_LIBJAM2
    if (runner_instance)
      return runner_instance->get_particle_mass(kf);
# endif

# ifdef USE_LIBJAM1
    return libjam1::jamMass(kf);
# else
    std::cerr << "FATAL: libjam2::runner not yet initialized." << std::endl;
    std::exit(3);
# endif
  }
#endif

  class NullRunner: public IJamRunner {
    std::string message;

    void abort() const {
      std::cerr << message << std::endl;
      std::exit(3);
    }

  public:
    NullRunner(const char* message): message(message) {}

    void initialize(runjam_context const& ctx, std::string const& cascadeMode) override {
      (void) ctx;
      (void) cascadeMode;
      this->abort();
    }
    void finalize() override {
      this->abort();
    }
    void adjust_mass(ParticleSampleBase& psample) const override {
      (void) psample;
      this->abort();
    }
    void do_decay(ParticleSampleBase& psamp, std::vector<Particle>& final_state) override {
      (void) psamp;
      (void) final_state;
      this->abort();
    }
    void do_cascade(ParticleSampleBase& psamp, std::vector<Particle>& final_state, int iev = 0) override {
      (void) psamp;
      (void) final_state;
      (void) iev;
      this->abort();
    }
    double get_average_collision_number() const override {
      this->abort();
      return 0.0;
    }
  };

#ifdef USE_LIBJAM2

  class Jam2Runner: public IJamRunner {
    std::unique_ptr<libjam2::irunner> m_runner;
    double m_oversample = 1.0;

  public:
    void initialize(runjam_context const& ctx, std::string const& cascadeMode) override {
      std::string const fconfig = ctx.get_config<std::string>("runjam_jam2_config", "/dev/null");
      m_runner = libjam2::create_runner(ctx, fconfig);
      m_oversample = ctx.get_config("runjam_oversampling_factor", 1.0);
      if (m_oversample != std::floor(m_oversample)) {
        std::cerr << "jam2: 'runjam_oversampling_factor=" << m_oversample << "' needs to be an integer." << std::endl;
        std::exit(2);
      }

      auto const settings = m_runner->settings();

      settings->flag("Random:setSeed", true);
      settings->mode("Random:seed", ctx.seed());
      settings->mode("Cascade:model", cascadeMode == "decay" ? 0 : 3);
      settings->mode("Cascade:overSample", (int) m_oversample);

      std::string const cachedir = ctx.cachedir();
      settings->word("Cascade:dataFile1", cachedir + "/BWintjam2a.dat");
      settings->word("Cascade:dataFile2", cachedir + "/BWintjam2b.dat");
      settings->word("Cascade:bwFileName1", cachedir + "/BWintjam2a.dat");
      settings->word("Cascade:bwFileName2", cachedir + "/BWintjam2b.dat");

      // weak decay
      {
        bool const weakdecay = ctx.get_config("runjam_switch_weak_decay", false);
        std::string const setting = weakdecay ? ":mayDecay = on" : ":mayDecay = off";

        settings->readString("111" + setting); // pi0
        settings->readString("311" + setting); // k0
        settings->readString("-311" + setting); // ak0
        settings->readString("310" + setting); // k0_S
        settings->readString("130" + setting); // k0_L

        settings->readString("411" + setting);  // D+
        settings->readString("421" + setting);  // D0
        settings->readString("221" + setting);  // eta
        settings->readString("331" + setting);  // eta'
        settings->readString("441" + setting);  // xeta_c
        settings->readString("310" + setting);
        settings->readString("431" + setting);
        settings->readString("511" + setting);
        settings->readString("521" + setting);
        settings->readString("531" + setting);
        settings->readString("3122" + setting);
        settings->readString("3112" + setting);
        settings->readString("3212" + setting);
        settings->readString("3222" + setting);
        settings->readString("3312" + setting);
        settings->readString("3322" + setting);
        settings->readString("3334" + setting);
      }

      // phi decay
      settings->flag("333:mayDecay", ctx.get_config("runjam_phi_decays", true));

      m_runner->initialize();
      if (runner_instance == nullptr)
        runner_instance = m_runner.get();
    }

    void finalize() override {
      m_runner.reset();
    }

    void adjust_mass(ParticleSampleBase& psample) const override {
      for (Particle& part: psample) {
        if (part.pdg == 0) continue;
        double const px = part.px;
        double const py = part.py;
        double const pz = part.pz;
        double const e = part.e;
        if (e < 0.0) {
          double const m = m_runner->get_particle_mass(part.pdg);
          part.mass = m;
          part.e = std::sqrt(px * px + py * py + pz * pz + m * m);
        } else {
          part.mass = std::sqrt(e * e - (px * px + py * py + pz * pz));
        }
      }
    }

    void do_cascade(ParticleSampleBase& psamp, std::vector<Particle>& final_state, int iev) override {
      if (psamp.getOverSamplingFactor() != m_oversample) {
        std::cerr << "ParticleSample: mismatching oversampling factor between the initial condition and JAM2 initialization." << std::endl;
        std::exit(1);
      }

      (void) iev;
      std::vector<Particle> initial_state(psamp.begin(), psamp.end());
      m_runner->run(initial_state, final_state);
    }

    void do_decay(ParticleSampleBase& psamp, std::vector<Particle>& final_state) override {
      // settings に指定している。
      this->do_cascade(psamp, final_state, 0);
    }

    double get_average_collision_number() const override {
      return m_runner->get_event_collision_number();
    }
  };

  std::unique_ptr<IJamRunner> create_jam2_runner() {
    return std::make_unique<Jam2Runner>();
  }

#else

  std::unique_ptr<IJamRunner> create_jam2_runner() {
    return std::make_unique<NullRunner>("FATAL: runjam: compiled without the JAM2 support.");
  }

#endif


#ifdef USE_LIBJAM1
  static void storeParticlesInJam1(Particle const* begin, Particle const* end) {
    int nv = 0;
    int nbary = 0;
    int nmeson = 0;

    for (; begin != end; ++begin) {
      Particle const& particle = *begin;
      int const kf = particle.pdg;
      if (kf == 0) continue;

      int const kc = libjam1::jamComp(kf);      // internal particle code.
      int const ks = libjam1::determineStableCode(kf);
      int const ibary = libjam1::getBaryonNumber(kc, kf);  // baryon number
      if (ibary == 0)
        nmeson++;
      else
        nbary++;
      nv++;

      //...Zero the vector.
      libjam1::jamZero(nv);
      libjam1::setK(1,  nv, ks);
      libjam1::setK(2,  nv, kf);
      libjam1::setK(3,  nv, 0);
      libjam1::setK(4,  nv, 0);
      libjam1::setK(5,  nv, -1);
      libjam1::setK(6,  nv, 0);
      libjam1::setK(7,  nv, 1);
      libjam1::setK(8,  nv, 1);
      libjam1::setK(9,  nv, ibary);
      libjam1::setK(10, nv, 0);
      libjam1::setK(11, nv, 0);

      double const x  = particle.x;
      double const y  = particle.y;
      double const z  = particle.z;
      double const t  = particle.t;
      double const px = particle.px  ;
      double const py = particle.py  ;
      double const pz = particle.pz  ;
      double const pe = particle.e   ;
      double const pm = particle.mass;

      libjam1::setP(1, nv, px);
      libjam1::setP(2, nv, py);
      libjam1::setP(3, nv, pz);
      libjam1::setP(4, nv, pe);
      libjam1::setP(5, nv, pm);

      libjam1::setR(1, nv, x);
      libjam1::setR(2, nv, y);
      libjam1::setR(3, nv, z);
      libjam1::setR(4, nv, t);
      libjam1::setR(5, nv, t); // formation time
      libjam1::setV(1, nv, x); // vertex
      libjam1::setV(2, nv, y); // vertex
      libjam1::setV(3, nv, z); // vertex
      libjam1::setV(4, nv, t); // vertex

      // Set resonance decay time.
      double decayTime = 1.0e+35;
      if (ks == 2) decayTime = t + libjam1::jamDecayTime(1, kf, kc, ks, pm, pe);
      libjam1::setV(5, nv, decayTime);
    }

    libjam1::setNV(nv);          // set total number of particles
    libjam1::setNBARY(nbary);    // set total number of baryons
    libjam1::setNMESON(nmeson);  // set total number of mesons
  }

  static int loadParticlesFromJam1(std::vector<Particle>& list) {
    list.clear();
    int const nv = libjam1::getNV();
    int const ntest = libjam1::getMSTC(5);
    for (int i = 1; i <= nv; i++) {
      //if(libjam1::getK(1,i) > 10) continue;
      list.emplace_back();
      Particle& part = list.back();
      part.pdg  = libjam1::getK(2, i);
      part.mass = libjam1::getP(5, i);
      part.px = libjam1::getP(1, i);
      part.py = libjam1::getP(2, i);
      part.pz = libjam1::getP(3, i);
      part.e  = libjam1::getP(4, i);
      part.x  = libjam1::getR(1, i);
      part.y  = libjam1::getR(2, i);
      part.z  = libjam1::getR(3, i);
      part.t  = libjam1::getR(4, i);
    }
    return ntest;
  }

  class Jam1Runner: public IJamRunner {
  public:
    void initialize(runjam_context const& ctx, std::string const& cascadeMode) override {
      (void) cascadeMode;

      int         const cfg_nevent     = ctx.nevent(1);
      std::string const cfg_outdir     = ctx.outdir();
      int         const cfg_jamseed    = ctx.get_config("runjam_jamseed", ctx.seed());
      bool        const cfg_phi_decays = ctx.get_config("runjam_phi_decays", true);
      bool        const cfg_weakdecay  = ctx.get_config("runjam_switch_weak_decay", false);

      // Initialize JAM
      libjam1::setMSTC(  1, cfg_jamseed); // int seed = 1921;
      libjam1::setMSTC(  2, cfg_nevent);  // number of event
      libjam1::setPARC(  6, 5.0); // scale of display
      libjam1::setMSTC(  8,   0); // job mode
      libjam1::setMSTC( 16,   0); // display on/off
      //libjam1::setMSTC( 38,   6); // IO number for jamlist
      //libjam1::setMSTC( 39,   0); //no output fname(4)
      libjam1::setMSTC( 54,   0); // avoid first coll inside the same nucleus off

      // Default settings
      //libjam1::setMSTC(156,   0); // analysis of collision distribution
      //libjam1::setMSTC(161,   0); // no analysis from jam internal subr
      //libjam1::setMSTC(162,   0); // Output collision history
      //libjam1::setMSTC(165,   0); //

      // Additional settings
      libjam1::setMSTC(  4, 100); // user defined frame
      libjam1::setPARC(  7, 1.0); // Output time interval (fm/c)
      //libjam1::setMSTC( 41,   0); // 0:no resonance decay after simulation
      //libjam1::setMSTC( 61,   0); // isotropic resonance decay option
      libjam1::setMSTC( 81,   0); // 1:hard scattering on/off
      libjam1::setMSTC(156,   1); // analysis of collision distribution
      libjam1::setMSTC(161,   0); // no analysis from jam internal subr
      libjam1::setMSTC(162,   1); // Output collision history
      libjam1::setMSTC(165,   1); //

      // JAM*.DAT files
      if (cfg_outdir.size() > 0) {
        std::string const dir = cfg_outdir + "/";
        libjam1::setFNAME(2, dir + libjam1::getFNAME(2));
        libjam1::setFNAME(3, dir + libjam1::getFNAME(3));
        libjam1::setFNAME(4, dir + libjam1::getFNAME(4));
        libjam1::setFNAME(8, dir);
      }

      //libjam1::setMDCY(libjam1::jamComp(111) , 1, 0); // no pi0 decay
      //libjam1::setMDCY(libjam1::jamComp(3122), 1, 1); // Lambda decay
      //libjam1::setMDCY(libjam1::jamComp(3222), 1, 1); // Sigma- decay
      //libjam1::setMDCY(libjam1::jamComp(3212), 1, 1); // Sigma0 decay
      //libjam1::setMDCY(libjam1::jamComp(3112), 1, 1); // Sigma+ decay
      if (!cfg_phi_decays)
        libjam1::setMDCY(libjam1::jamComp(333), 1, 0); // no phi decay

      if (cfg_weakdecay) libjam1::setMSTC(42, 0); // allow weak decay

      // JAMINIT
      std::string frame = "user"; // comp. frame in this case, user defined
      double dt = 100.0;          // collision time(fm/c)
      int nstep = 1;              // time step (i.e. no time step)
      double bmin = 0.0;         // minimum impact parameter (dummy)
      double bmax = 0.0;         // maximum impact parameter (dummy)
      libjam1::setPARD(16, 10.0); // user defined frame.
      libjam1::jamInit(cfg_nevent, bmin, bmax, dt, nstep, frame.c_str(), "p ", "p ", "2gev");
    }

    void finalize() override {
      libjam1::jamFin();
    }

    void adjust_mass(ParticleSampleBase& psample) const override {
      for (Particle& part: psample) {
        if (part.pdg == 0) continue;
        double const px = part.px;
        double const py = part.py;
        double const pz = part.pz;
        double const e = part.e;
        if (e < 0.0) {
          double const m = runjam::getParticleMass(part.pdg);
          part.mass = m;
          part.e = std::sqrt(px * px + py * py + pz * pz + m * m);
        } else {
          part.mass = std::sqrt(e * e - (px * px + py * py + pz * pz));
        }
      }
    }

    void do_decay(ParticleSampleBase& psamp, std::vector<Particle>& final_state) override {
      double const ntest = psamp.getOverSamplingFactor();
      libjam1::setMSTC(5, ntest);
      storeParticlesInJam1(psamp.begin(), psamp.end());
      libjam1::finalResonanceDecay();
      loadParticlesFromJam1(final_state);
    }

    void do_cascade(ParticleSampleBase& psamp, std::vector<Particle>& final_state, int iev = 0) override {
      double const ntest = psamp.getOverSamplingFactor();
      libjam1::setMSTC(5, ntest);
      storeParticlesInJam1(psamp.begin(), psamp.end());
      libjam1::jamEvt(iev);
      loadParticlesFromJam1(final_state);
    }

    double get_average_collision_number() const override {
      int const ntest = libjam1::getMSTC(5);
      return (libjam1::getMSTD(41) + libjam1::getMSTD(42)) / ntest;
    }
  };

  std::unique_ptr<IJamRunner> create_jam1_runner() {
    return std::make_unique<Jam1Runner>();
  }

#else

  std::unique_ptr<IJamRunner> create_jam1_runner() {
    return std::make_unique<NullRunner>("FATAL: runjam: compiled without the JAM1 support.");
  }

#endif

  std::unique_ptr<IJamRunner> create_runner(runjam_context const& ctx) {
#if defined(USE_LIBJAM1) || defined(USE_LIBJAM2)
    int const libjam_version = ctx.get_config("runjam_jam_version", DEFAULT_JAM_VERSION);
    switch (libjam_version) {
    case 1:
#ifdef USE_LIBJAM1
      return runjam::create_jam1_runner();
#else
      std::cerr << "runjam was compiled without the support for JAM1" << std::endl;
      break;
#endif
    case 2:
#ifdef USE_LIBJAM2
      return runjam::create_jam2_runner();
#else
      std::cerr << "runjam was compiled without the support for JAM2" << std::endl;
      break;
#endif
    default:
      std::cerr << "runjam: invalid value runjam_jam_version='" << libjam_version << "'" << std::endl;
      break;
    }
#else
    std::cerr << "  runjam was compiled without the support for JAM" << std::endl;
#endif
    return std::unique_ptr<IJamRunner>();
  }

}
}
