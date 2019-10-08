// -*- mode:c++ -*-
#ifndef Random_h
#define Random_h

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

class Random{
protected:
  int seed;
public:
  Random()      { seed = ::time(0); /* srand=0; */ }
  Random(int s) { seed = s;         /* srand=0; */ }
  virtual ~Random() { }
  virtual void   setSeed()      {seed = ::time(0); srand48(seed);}
  virtual void   setSeed(int s) {srand48(s); seed=s;}
  virtual double rand()         { return drand48();}
  double Gauss(double mean, double sigma);

private:
  static Random*  srand;
  //static bool     first;

public:
  static void setRandom(Random* r) {
    if(srand) {
	    cerr << "warning: setRandom should be called once!!!"
           << endl;
	    return;
    }
    srand = r;
  }

  static Random* getRandom(){return srand;}

  /// @fn getRand
  /// generates a random number following an uniform distribution with the interval [0,1).
  static double getRand(){return srand->rand();}

  // ToDo: this not efficient. improve the implementation (KM)
  /// @fn getRandomPoisson
  /// generates a random number following an Poisson distribution with the specified mean value.
  /// @param value the mean value of a Poisson distribution
  static int getRandPoisson(double value);

};

//extern Random *gRandomX;

#endif
