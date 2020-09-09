// -*- mode: c++ -*-
#ifndef runjam_util_Random_hpp
#define runjam_util_Random_hpp

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

namespace idt {
namespace runjam {

class Random {
protected:
  int seed;
public:
  Random()      { setSeed(); }
  Random(int s) { setSeed(s); }
  virtual ~Random() { }
  virtual void   setSeed()      { setSeed(::time(0)); }
  virtual void   setSeed(int s) { srand48(s); seed = s; }
  virtual double rand()         { return drand48(); }
  virtual unsigned long lrand() { return lrand48(); }
  double Gauss(double mean, double sigma);

private:
  static Random* srand;
  //static bool     first;

public:
  static void setRandom(Random* r) {
    if(srand) {
	    std::cerr << "warning: setRandom should be called once!!!"
           << std::endl;
	    return;
    }
    srand = r;
  }

  static Random* getRandom() { return srand; }

  /// @fn getRand
  /// generates a random number following an uniform distribution with the interval [0,1).
  static double getRand() { return srand->rand(); }

  // ToDo: this not efficient. improve the implementation (KM)
  /// @fn getRandomPoisson
  /// generates a random number following an Poisson distribution with the specified mean value.
  /// @param value the mean value of a Poisson distribution
  static int getRandPoisson(double value);

  static long getUniformLong() { return srand->lrand(); }

};

//extern Random *gRandomX;

struct RandomURGB {
  typedef unsigned long result_type;
  static result_type min() { return 0; }
  static result_type max() { return ((result_type) 1 << 31) - 1; }
  result_type operator()() const { return Random::getUniformLong(); }
};

}
}

#endif
