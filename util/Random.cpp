#define _USE_MATH_DEFINES
#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif
#include <cmath>

#include "Random.hpp"
#include "Math.hpp"

namespace idt {
namespace runjam {

Random* Random::srand = 0;
//bool EventBase::first = true;

double Random::Gauss(double mean, double sigma)
{
  double  x, y, z, result;
  do {
    y = this->rand();
  } while (!y);
  z = this->rand();
  x = z * 6.283185;
  result = mean + sigma * std::sin(x) * std::sqrt(-2 * std::log(y));
  return result;
}

int Random::getRandPoisson(double lambda) {
  if (lambda < 30) {
    double const exp_ = std::exp(-lambda);
    int n;
    double xp = 1.0;
    for (n = 0; (xp *= Random::getRand()) > exp_;n++);
    return n;
  } else {
    // From http://www.johndcook.com/blog/cpp_random_number_generation/
    //      http://www.johndcook.com/SimpleRNG.cpp

    // "Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
    // Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
    // The article is on pages 29-35. The algorithm given here is on page 32.

    double const c        = 0.767 - 3.36 / lambda;
    double const ibeta    = std::sqrt(3.0 * lambda) * (1.0 / M_PI);
    double const k        = std::log(c * ibeta) - lambda;
    double const lnLambda = std::log(lambda);

    for (;;) {
      double const u = Random::getRand();
      double const y = std::log(1.0 / u - 1.0);
      int    const n = (int) std::floor(lambda - y*ibeta + 0.5);
      if (n < 0) continue;
      double const lhs = y + std::log(Random::getRand()*u*u);
      double const rhs = k + n * lnLambda - Math::logFactorial(n);
      if (lhs <= rhs) return n;
    }
  }
}

}
}
