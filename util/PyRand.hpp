#ifndef runjam_util_PyRand_hpp
#define runjam_util_PyRand_hpp

#include "Random.hpp"

namespace idt {
namespace runjam {
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

#endif /* runjam_util_PyRand_hpp */
