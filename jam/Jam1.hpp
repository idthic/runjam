// -*- mode:c++ -*-
#ifndef Jam1_h
#define Jam1_h
#include <cstdlib>
#include "Jam1Common.hpp"

namespace idt {
namespace hydro2jam {

class Jam1 {
public:
  Jam1();
  ~Jam1();

  int    getNV()                       { return jamevnt2_.NV; }
  int    getNBARY()                    { return jamevnt2_.NBARY; }
  int    getNMESON()                   { return jamevnt2_.NMESON; }
  int    getK(int i, int ip)           { return jamevnt1_.K[ip-1][i-1]; }
  double getR(int i, int ip)           { return jamevnt1_.R[ip-1][i-1]; }
  double getP(int i, int ip)           { return jamevnt1_.P[ip-1][i-1]; }
  double getV(int i, int ip)           { return jamevnt1_.V[ip-1][i-1]; }

  void   setNV(int n)                  { jamevnt2_.NV = n;    }
  void   setNBARY(int n)               { jamevnt2_.NBARY = n; }
  void   setNMESON(int n)              { jamevnt2_.NMESON = n; }
  void   setK(int i, int ip, int k)    { jamevnt1_.K[ip-1][i-1] = k; }
  void   setR(int i, int ip, double r) { jamevnt1_.R[ip-1][i-1] = r; }
  void   setP(int i, int ip, double p) { jamevnt1_.P[ip-1][i-1] = p; }
  void   setV(int i, int ip, double v) { jamevnt1_.V[ip-1][i-1] = v; }

  static int         getMSTC(int i) { return jamdat1_.MSTC[i-1]; }
  static double      getPARC(int i) { return jamdat1_.PARC[i-1]; }
  static int         getMSTD(int i) { return jamdat1_.MSTD[i-1]; }
  static double      getPARD(int i) { return jamdat1_.PARD[i-1]; }

  void        setMSTC(int i, int m   ) { jamdat1_.MSTC[i-1] = m; }
  void        setPARC(int i, double p) { jamdat1_.PARC[i-1] = p; }
  void        setMSTD(int i, int m   ) { jamdat1_.MSTD[i-1] = m; }
  void        setPARD(int i, double p) { jamdat1_.PARD[i-1] = p; }
  void        addPARD(int i, double p) { jamdat1_.PARD[i-1] += p; }

  static int         getKCHG(int ip, int i) { return jydat2_.KCHG[i-1][ip-1]; }
  static double      getPMAS(int ip, int i) { return jydat2_.PMAS[i-1][ip-1]; }
  static double      getPARF        (int i) { return jydat2_.PARF[i-1]; }
  static double      getVCKM(int i,  int j) { return jydat2_.VCKM[j-1][i-1]; }

  void        setKCHG(int ip, int i, int k   ) { jydat2_.KCHG[i-1][ip-1] = k; }
  void        setPMAS(int ip, int i, double m) { jydat2_.PMAS[i-1][ip-1] = m; }
  void        setPARF        (int i, double p) { jydat2_.PARF[i-1]       = p; }
  void        setVCKM (int i, int j, double v) { jydat2_.VCKM[j-1][i-1]  = v; }


  static int         getMDCY(int i, int j) { return jydat3_.MDCY[j-1][i-1]; }
  static int         getMDME(int i, int j) { return jydat3_.MDME[j-1][i-1]; }
  static double      getBRAT       (int i) { return jydat3_.BRAT[i-1]; }
  static int         getKFDP(int i, int j) { return jydat3_.KFDP[j-1][i-1]; }

  void        setMDCY(int i, int j, int m) { jydat3_.MDCY[j-1][i-1] = m; }
  void        setMDME(int i, int j, int m) { jydat3_.MDME[j-1][i-1] = m; }
  void        setBRAT(int i, double b)     { jydat3_.BRAT[i-1]      = b; }
  void        setKFDP(int i, int j, int k) { jydat3_.KFDP[j-1][i-1] = k; }

  char*  getFNAME(int i) const;
  void   setFNAME(int i, const char* v);

  static int isign(const int a, const int b) { return  b >= 0 ? std::abs(a): -std::abs(a); }

  int  getBaryonNumber(int kc, int kf) {return isign(jydat2_.KCHG[5][kc-1],kf);}

  void  jamInit(int mev, double bmin, double bmax, double dt, int nstp,
  const char* frame, const char* proj, const char* targ, const char* cwin);
  void    jamEvt(int iev);
  void    jamFin();
  static int     jamComp(int kf);
  void    jamZero(int i);
  double  jamDecayTime(int im,int kf, int kc, int ks, double em,double ee);
  void           jamList(int i);
  static double  jamMass(int kf);
  static int     jamCharge(int kf);
  void           finalResonanceDecay();
  void         print(int i);

};

}
}

#endif
