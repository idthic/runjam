#ifndef PyRand_h
#define PyRand_h

#include "Random.h"

//The random number generator from Pythia6.

class PyRand: public Random
{
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


#endif /* base_rand_PyRand_h */
