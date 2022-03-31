#include "config.hpp"
#include <string>
#include <memory>
#include "ParticleSample.hpp"
#include "args.hpp"

#if defined(USE_LIBJAM2)
namespace libjam2 { class runner; }
#endif

namespace idt {
namespace runjam {

#if defined(USE_LIBJAM2) || defined(USE_LIBJAM1)
  int getParticleStableCode(int kf);
  double getParticleMass(int kf);
#endif

  class IJamRunner {
  public:
    virtual void initialize(runjam_context const& ctx, std::string const& cascadeMode) = 0;
    virtual void finalize() = 0;
    virtual void adjust_mass(ParticleSampleBase& psample) const = 0;
    virtual void do_decay(ParticleSampleBase& psamp, std::vector<Particle>& final_state) = 0;
    virtual void do_cascade(ParticleSampleBase& psamp, std::vector<Particle>& final_state, int iev = 0) = 0;
    virtual double get_average_collision_number() const = 0;
    virtual ~IJamRunner() {}
  };

  std::unique_ptr<IJamRunner> create_jam1_runner();
  std::unique_ptr<IJamRunner> create_jam2_runner();
  std::unique_ptr<IJamRunner> create_runner(runjam_context const& ctx);
}
}
