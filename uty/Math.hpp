// -*- mode:c++;indent-tabs-mode:nil -*-
#ifndef Math_h
#define Math_h
#include <cmath>
#include <cstdlib>
#include <algorithm>

namespace idt {
namespace hydro2jam {

class Math
{
public:
  // ordinary Gamma function Gamma(x) for positive, real arguments.
  static double tgamma(double x);

  // Euler's beta function, requires ordinary Gamma function
  static double eulBet(double x,double y)
  { return tgamma(x)*tgamma(y)/tgamma(x+y);}

  // log(n!).
  static double logFactorial(int n);
  static double sign(double a, double b);
  static int    sign(int a, int b);
  static double sq(const double x) {return x*x;}

  static double BesselJ0(double x);
  static double BesselJ1(double x);
  static double BesselJ(int n, double x);

  /* Bessel I_0(x) function in double precision */
  static double BesselI0(double x);
};

inline double Math::sign(double a, double b) {
  return b >= 0.0 ? std::abs(a) : -std::abs(a);
}

inline int Math::sign(int a, int b) {
  return b >= 0 ? std::abs(a) : -std::abs(a);
}

template<class T> inline const T& min(const T& a, const T& b, const T& c)
{ return std::min(std::min(a, b), c); }

template<class T> inline const T& max(const T& a, const T& b, const T& c)
{ return std::max(std::max(a, b), c); }

}
}

#endif // Math_h
