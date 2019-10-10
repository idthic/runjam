#ifndef hydro2jam_util_PyRand_hpp
#define hydro2jam_util_PyRand_hpp

#include "Random.hpp"

namespace idt {
namespace hydro2jam {
//The random number generator from Pythia6.

class PyRand: public Random {
private:
  int    mrpy1,mrpy2,mrpy3,mrpy4,mrpy5,mrpy6;
  double rrpy[101];
  double* rrpy98;
  double* rrpy99;
  double* rrpy00;
public:
  PyRand():Random() { init();}
  PyRand(int i) : Random(i) { init(); }
  virtual ~PyRand() { }
  virtual double rand();

protected:
  void init();
};

}
}

#endif /* hydro2jam_util_PyRand_hpp */
