#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef USE_JAM
// complie example for gcc:   g++ Jam1.cxx  -L/home/ynara/lib -ljam -lg2c
// complie example for Alpha: cxx Jam1.cxx  -L/home/ynara/lib -ljam -lm -lfor
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "Jam1.h"

using namespace std;

#ifdef HYDROJET_FORTRAN_NOEXTNAME
#  define jaminit_ jaminit
#  define jamevt_  jamevt
#  define jamfin_  jamfin
#  define jamzero_ jamzero
#  define jamcomp_ jamcomp
#  define jamdtim_ jamdtim
#  define jamlist_ jamlist
#  define jamchge_ jamchge
#  define jamfdec_ jamfdec
#  define pjmass_  pjmass
#endif

extern "C" void   jaminit_(int *mev, double *bmin, double *bmax, double *dt,
                           int *nstp, const char *frame, const char *proj,
                           const char *targ, const char *cwin,
                           int lf, int lp, int lt, int lc);
extern "C" void   jamevt_(int *iev);
extern "C" void   jamfin_();
extern "C" void   jamzero_(int *ip);
extern "C" int    jamcomp_(int *kf);
extern "C" double jamdtim_(int *im,int *kf,int *kc,int *ks,double *em,double *ee);
extern "C" void    jamlist_(int *imode);
extern "C" int     jamchge_(int *kf);
extern "C" void    jamfdec_();
extern "C" double  pjmass_(int *kf);

Jam1::Jam1()
{
}

Jam1::~Jam1()
{
}

void Jam1::jamInit(int mev, double bmin, double bmax, double dt, int nstp,
      const char *frame, const char *proj, const char *targ, const char *cwin)
{
    int  s1    = strlen(frame);
    int  s2    = strlen(proj);
    int  s3    = strlen(targ);
    int  s4    = strlen(cwin);
    jaminit_(&mev, &bmin,&bmax, &dt, &nstp,frame,proj,targ,cwin,s1,s2,s3,s4);

}

void Jam1::print(int i)
{
    cout << i
	<< " KF= " << getK(2,i)
	<< endl
	<< " K= "   << getK(1,i)
	<< setw(5) << getK(2,i)
	<< setw(5) << getK(3,i)
	<< setw(5) << getK(4,i)
	<< setw(5) << getK(5,i)
	<< setw(5) << getK(6,i)
	<< setw(5) << getK(7,i)
	<< setw(5) << getK(8,i)
	<< setw(5) << getK(9,i)
	<< setw(5) << getK(10,i)
	<< setw(5) << getK(11,i)
	<< endl
	<< " R= "   << getR(1,i)
	<< setw(10) << getR(2,i)
	<< setw(10) << getR(3,i)
	<< setw(10) << getR(4,i)
	<< setw(10) << getR(5,i)
	<< endl
	<< " P= "   << getP(1,i)
	<< setw(15) << getP(2,i)
	<< setw(15) << getP(3,i)
	<< setw(15) << getP(4,i)
	<< setw(15) << getP(5,i)
	<< endl
	<< " V= "   << getV(1,i)
	<< setw(10) << getV(2,i)
	<< setw(10) << getV(3,i)
	<< setw(10) << getV(4,i)
	<< setw(10) << getV(5,i)
	<< endl;
}

void Jam1::jamEvt(int iev)
{
    jamevt_(&iev);
}

void Jam1::jamFin()
{
    jamfin_();
}

int Jam1::jamComp(int kf)
{
    return jamcomp_(&kf);
}

void Jam1::jamZero(int i)
{
    jamzero_(&i);
}

double Jam1::jamDecayTime(int im,int kf, int kc, int ks, double em,double ee)
{
    return jamdtim_(&im,&kf,&kc,&ks,&em,&ee);
}

int Jam1::jamCharge(int kf) {
    return jamchge_(&kf);
}

double Jam1::jamMass(int kf) {
    return pjmass_(&kf);
}

void Jam1::jamList(int i) {
    jamlist_(&i);
}

char* Jam1::getFNAME(int i) const
{
    static char buf[81]="";
    strncpy(buf,jamdat3_.FNAME[i-1],80);
    buf[80]=0;
    return buf;
    //return jamdat3_.FNAME[i-1][strlen(jamdat3_.FNAME[i-1])];
}

void  Jam1::setFNAME(int i, const char* v)
{
    strncpy(jamdat3_.FNAME[i-1],v,strlen(v));
}

void  Jam1::finalResonanceDecay()
{
    jamfdec_();
}

//#define MAIN
#ifdef MAIN
#include <iostream>
using namespace std;
int main() {

    Jam1 *jam = new Jam1();

    //for(int i=0;i<200;i++) cout << jam->getMSTC(i) << endl;

    int mev=1;
    double bmin=0.0;
    double bmax=0.0;
    double dt=100;
    int nstep=1;
    jam->jamInit(mev,bmin,bmax,dt,nstep,"nn","p ","p ","200.0gev");
    cout << "after jamini" << endl;
    jam->jamEvt(1);
    cout << "after jamevt" << endl;
    jam->jamFin();
    delete jam;

    return 0;
}
#endif

#endif
