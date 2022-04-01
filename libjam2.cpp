#include "config.hpp"
#include "libjam2.hpp"
#include <jam2/JAM.h>

namespace {

  class initial_condition_adapter: public jam2::UserInitialCondition {
    typedef jam2::UserInitialCondition base;
    static constexpr double DEFAULT_MAX_CM_ENERGY = 20.0;

  public:
    std::vector<idt::runjam::Particle> const* m_particles = nullptr;

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
    initial_condition_adapter(jam2::JAM* jam):
      base(jam->settings, jam->jamParticleData, jam->rndm) {}

    virtual ~initial_condition_adapter() {}

    virtual void init() override {
      // We have been initially determining "Beams:eCM" based on the list of
      // particles using `determine_max_cm_energy_from_particles`, but we now
      // cannot do in that way because we need to initialize Pythia before
      // seeing the initial particles for the multiple events.  We instead just
      // set the fixed value.
      settings->parm("Beams:eCM", DEFAULT_MAX_CM_ENERGY);
    }

    virtual void generate(jam2::Collision* event, int mode = 0) override {
      (void) mode;

      for (idt::runjam::Particle const& particle: *m_particles) {
        // find this particle in the JAM particle list.
        Pythia8::ParticleDataEntry* const pa = jamParticleData->find(particle.pdg);
        Pythia8::Vec4 const r(particle.pos[1], particle.pos[2], particle.pos[3], particle.pos[0]);
        Pythia8::Vec4 const p(particle.mom[1], particle.mom[2], particle.mom[3], particle.mom[0]);

        auto const cp = new jam2::EventParticle(particle.pdg, particle.mass, r, p, pa);
        cp->setPID(jamParticleData->pid(std::abs(particle.pdg)));

        // compute decay time if it is resonance.
        double const decay_time = jamParticleData->lifeTime(pa, particle.mass, particle.mom[0]);
        cp->setLifeTime(particle.pos[0] + decay_time);

        // put this particle into the particle list.
        event->setPList(cp);
      }
    }
  };

}

namespace libjam2 {

  class runner: public irunner {
    jam2::JAM m_jam;
    initial_condition_adapter* m_initial_condition;

  public:
    runner(std::string const& input_filename):
      m_jam(input_filename, Pythia8_PREFIX "/share/Pythia8/xmldoc", false)
    {
      m_initial_condition = new initial_condition_adapter(&m_jam);

      Pythia8::Settings* const settings = this->settings();
      settings->mode("Cascade:model", 3);
      settings->mode("Cascade:initialCondition", 3);
    }

    Pythia8::Settings* settings() override {
      return m_jam.settings;
    }

    Pythia8::Settings const* settings() const override {
      return m_jam.settings;
    }

    void initialize() override {
      m_jam.init(m_initial_condition);
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
      // stable particle = 1
      // unstable particle = 2
      return m_jam.jamParticleData->find(pdg)->mWidth() <= 1e-7 ? 1 : 2;
    }
    double get_event_collision_number() const override {
      return m_jam.getNColl();
    }
  };

  std::unique_ptr<irunner> create_runner(std::string const& input_filename) {
    return std::make_unique<runner>(input_filename);
  }

}
