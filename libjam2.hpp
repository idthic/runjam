#include <memory>
#include <vector>
#include <string>
#include <Pythia8/Settings.h>

namespace jam { struct JAM; }

namespace libjam2 {

  struct phasespace_particle {
    int pdg_code;
    double mass;
    double pos[4];
    double mom[4];
  };

  class irunner {
  public:
    virtual ~irunner() {}

    virtual Pythia8::Settings* settings() = 0;
    virtual Pythia8::Settings const* settings() const = 0;
    virtual void initialize() = 0;
    virtual void run(std::vector<phasespace_particle> const& initial_state, std::vector<phasespace_particle>& final_state) = 0;

    virtual double get_particle_mass(int pdg) const = 0;
    virtual int get_particle_stable_code(int pdg) const = 0;

    //! The expected number of collisions in the last event, i.e., ncoll / oversample
    virtual double get_event_collision_number() const = 0;
  };

  std::unique_ptr<irunner> create_runner(std::string const& input_filename = "/dev/null");
}
