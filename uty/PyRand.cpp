#include "PyRand.hpp"

namespace idt {
namespace hydro2jam {

//...Generates random numbers uniformly distributed between
//...0 and 1, excluding the endpoints.

void PyRand::init()
{

  mrpy1 = seed;
  //...initial values for the random number generator.
  //mrpy1 = 19780503;

  //...Initialize generation from given seed.
  //if(mrpy2 != 0) return;

  int ij = (mrpy1/30082) % 31329;
  int kl = mrpy1 % 30082;
  int i = (ij/177) % 177 + 2;
  int j = (ij % 177) + 2;
  int k = (kl/169) % 178+1;
  int l = kl % 169;
  for (int ii=1; ii<= 97; ii++) {
    double s=0.0;
    double t=0.5;
    for(int jj=1; jj<= 48; jj++) {
      int m = ((i*j) % 179 *k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1)%169;
      if((l*m)%64 >= 32) s += t;
      t *=0.5;
    }
    rrpy[ii]=s;
  }

  double twom24 = 1.0;
  for (int i24 = 1; i24 <= 24; i24++) {
    twom24 *= 0.5;
  }
  rrpy[98] = 362436.0*twom24;
  rrpy[99] = 7654321.0*twom24;
  rrpy[100] = 16777213.0*twom24;

  mrpy2 = 1;
  mrpy3 = 0;
  mrpy4 = 97;
  mrpy5 = 33;
  mrpy6 = 0;
  rrpy98 = &rrpy[98];
  rrpy99 = &rrpy[99];
  rrpy00 = &rrpy[100];

}

//...Generate next random number.
double PyRand::rand()
{
  double runi;

  while(1) {
    runi = rrpy[mrpy4]-rrpy[mrpy5];
    if(runi < 0.0) runi += 1.0;
    rrpy[mrpy4] = runi;
    mrpy4--;
    if(mrpy4 == 0) mrpy4 = 97;
    mrpy5--;
    if(mrpy5 == 0) mrpy5 = 97;
    *rrpy98 -= *rrpy99;
    if(*rrpy98 < 0.0) *rrpy98 += *rrpy00;
    runi -= *rrpy98;
    if(runi < 0.0) runi += 1.0;
    if(runi > 0.0 && runi < 1.0) break;
  }

  //...Update counters. Random number to output.
  mrpy3++;
  if(mrpy3 == 1000000000) {
    mrpy2++;
    mrpy3=0;
  }

  return runi;
}

}
}
