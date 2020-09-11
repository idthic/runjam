#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "libjam.hpp"

#include "config.hpp"
// This number should be the same as "mxv" in jam1.inc
#ifndef CONFIG_JAM_MXV
# define CONFIG_JAM_MXV  200000
//#define CONFIG_JAM_MXV  100000
//#define CONFIG_JAM_MXV  50000
//#define CONFIG_JAM_MXV  30000
#endif

#ifdef CONFIG_FORTRAN_NOEXTNAME
#  define jamevnt1_ jamevnt1
#  define jamevnt2_ jamevnt2
#  define jamdat1_  jamdat1
#  define jamdat2_  jamdat2
#  define jamdat3_  jamdat3
#  define jydat2_   jydat2
#  define jydat3_   jydat3
#  define jydat4_   jydat4
#
#  define jaminit_  jaminit
#  define jamevt_   jamevt
#  define jamfin_   jamfin
#  define jamzero_  jamzero
#  define jamcomp_  jamcomp
#  define jamdtim_  jamdtim
#  define jamlist_  jamlist
#  define jamchge_  jamchge
#  define jamfdec_  jamfdec
#  define pjmass_   pjmass
#endif


namespace idt {
namespace libjam {

  static std::size_t strlen_fortran(const char* buff, int bufferSize = -1) {
    std::size_t ret = 0, run = 0;
    const char* end = buff + bufferSize;
    for (; buff != end && *buff; buff++) {
      if (*buff == ' ') {
        run++;
      } else {
        ret += run + 1;
        run = 0;
      }
    }
    return ret;
  }

  static char* strcpy_f2c(char* dst, std::size_t szdst, const char* src, std::size_t szsrc) {
    std::size_t len = std::min(szdst, strlen_fortran(src, szsrc));
    std::memcpy(dst, src, len);
    if (len < szdst) dst[len++] = '\0';
    return dst + len;
  }

  static bool strcpy_c2f(char* dst, std::size_t szdst, const char* src) {
    char* end = dst + szdst;
    while (dst != end && *src)
      *dst++ = *src++;
    while (dst != end)
      *dst++ = ' ';
    return !*src;
  }

  static int isign(const int a, const int b) {
    return  b >= 0 ? std::abs(a): -std::abs(a);
  }

  //---------------------------------------------------------------------------
  // JAM data section

  struct Jamevnt1 {
    double R[CONFIG_JAM_MXV][5];
    double P[CONFIG_JAM_MXV][5];
    double V[CONFIG_JAM_MXV][5];
    int    K[CONFIG_JAM_MXV][11];
  };
  extern "C" Jamevnt1 jamevnt1_;

  int    getK(int i, int ip)           { return jamevnt1_.K[ip-1][i-1]; }
  double getR(int i, int ip)           { return jamevnt1_.R[ip-1][i-1]; }
  double getP(int i, int ip)           { return jamevnt1_.P[ip-1][i-1]; }
  double getV(int i, int ip)           { return jamevnt1_.V[ip-1][i-1]; }
  void   setK(int i, int ip, int k)    { jamevnt1_.K[ip-1][i-1] = k; }
  void   setR(int i, int ip, double r) { jamevnt1_.R[ip-1][i-1] = r; }
  void   setP(int i, int ip, double p) { jamevnt1_.P[ip-1][i-1] = p; }
  void   setV(int i, int ip, double v) { jamevnt1_.V[ip-1][i-1] = v; }

  struct Jamevnt2 {
    int NV;
    int NBARY;
    int NMESON;
  };
  extern "C" Jamevnt2 jamevnt2_;

  int    getNV()                       { return jamevnt2_.NV; }
  int    getNBARY()                    { return jamevnt2_.NBARY; }
  int    getNMESON()                   { return jamevnt2_.NMESON; }
  void   setNV(int n)                  { jamevnt2_.NV = n;    }
  void   setNBARY(int n)               { jamevnt2_.NBARY = n; }
  void   setNMESON(int n)              { jamevnt2_.NMESON = n; }

  struct Jamdat1 {
    int    MSTC[200];
    double PARC[200];
    int    MSTD[200];
    double PARD[200];
  };
  extern "C" Jamdat1 jamdat1_;

  int    getMSTC(int i) { return jamdat1_.MSTC[i-1]; }
  double getPARC(int i) { return jamdat1_.PARC[i-1]; }
  int    getMSTD(int i) { return jamdat1_.MSTD[i-1]; }
  double getPARD(int i) { return jamdat1_.PARD[i-1]; }
  void   setMSTC(int i, int m   ) { jamdat1_.MSTC[i-1] = m; }
  void   setPARC(int i, double p) { jamdat1_.PARC[i-1] = p; }
  void   setMSTD(int i, int m   ) { jamdat1_.MSTD[i-1] = m; }
  void   setPARD(int i, double p) { jamdat1_.PARD[i-1] = p; }
  void   addPARD(int i, double p) { jamdat1_.PARD[i-1] += p; }

  struct Jamdat2 {
    int    MSTE[200];
    double PARE[200];
  };
  extern "C" Jamdat2 jamdat2_;

  struct Jamdat3 {
    char  FNAME[8][80];
  };
  extern "C" Jamdat3 jamdat3_;

  char*  getFNAME(int i) {
    static char buf[81];
    strcpy_f2c(buf, 80, jamdat3_.FNAME[i - 1], 80);
    //strncpy(buf, jamdat3_.FNAME[i - 1], 80);
    buf[80] = '\0';
    return buf;
  }
  void   setFNAME(int i, const char* v) {
    if (!strcpy_c2f(jamdat3_.FNAME[i - 1], 80, v)) {
      std::cerr << "runjam: filename \"" << v << "\" is too long" << std::endl;
      std::exit(1);
    }
  }
  void   setFNAME(int i, std::string const& v) { setFNAME(i, v.c_str()); }

  struct Jydat2 {
    int    KCHG[7][500];
    double PMAS[4][500];
    double PARF[2000];
    double VCKM[4][4];
  };
  extern "C" Jydat2 jydat2_;

  int    getKCHG(int ip, int i) { return jydat2_.KCHG[i-1][ip-1]; }
  double getPMAS(int ip, int i) { return jydat2_.PMAS[i-1][ip-1]; }
  double getPARF        (int i) { return jydat2_.PARF[i-1]; }
  double getVCKM(int i,  int j) { return jydat2_.VCKM[j-1][i-1]; }
  void   setKCHG(int ip, int i, int k   ) { jydat2_.KCHG[i-1][ip-1] = k; }
  void   setPMAS(int ip, int i, double m) { jydat2_.PMAS[i-1][ip-1] = m; }
  void   setPARF        (int i, double p) { jydat2_.PARF[i-1]       = p; }
  void   setVCKM (int i, int j, double v) { jydat2_.VCKM[j-1][i-1]  = v; }
  int    getBaryonNumber(int kc, int kf) { return isign(jydat2_.KCHG[5][kc-1],kf); }


  static const int KNDCAY1 = 4000; //should be 8000 for pythia62

  struct Jydat3 {
    int    MDCY[3][500];
    int    MDME[3][KNDCAY1];
    double BRAT[KNDCAY1];
    int    KFDP[5][KNDCAY1];
  };
  extern "C" Jydat3 jydat3_;

  int    getMDCY(int i, int j) { return jydat3_.MDCY[j-1][i-1]; }
  int    getMDME(int i, int j) { return jydat3_.MDME[j-1][i-1]; }
  double getBRAT       (int i) { return jydat3_.BRAT[i-1]; }
  int    getKFDP(int i, int j) { return jydat3_.KFDP[j-1][i-1]; }
  void   setMDCY(int i, int j, int m) { jydat3_.MDCY[j-1][i-1] = m; }
  void   setMDME(int i, int j, int m) { jydat3_.MDME[j-1][i-1] = m; }
  void   setBRAT(int i, double b)     { jydat3_.BRAT[i-1]      = b; }
  void   setKFDP(int i, int j, int k) { jydat3_.KFDP[j-1][i-1] = k; }

  struct Jydat4 {
    char  CHAF[2][500][16];	// here I needed manual intervention
  };
  extern "C" Jydat4 jydat4_;

  //---------------------------------------------------------------------------
  // JAM routines

  extern "C" void   jaminit_(
    int *mev, double *bmin, double *bmax, double *dt,
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

  void jamInit(int mev, double bmin, double bmax, double dt, int nstp,
    const char *frame, const char *proj, const char *targ, const char *cwin) {
    int s1 = std::strlen(frame);
    int s2 = std::strlen(proj);
    int s3 = std::strlen(targ);
    int s4 = std::strlen(cwin);
    jaminit_(&mev, &bmin,&bmax, &dt, &nstp, frame, proj, targ, cwin, s1, s2, s3, s4);
  }

  void jamEvt(int iev) {
    jamevt_(&iev);
  }

  void jamFin() {
    jamfin_();
  }

  int jamComp(int kf) {
    return jamcomp_(&kf);
  }

  void jamZero(int i) {
    jamzero_(&i);
  }

  double jamDecayTime(int im, int kf, int kc, int ks, double em, double ee) {
    return jamdtim_(&im, &kf, &kc, &ks, &em, &ee);
  }

  int jamCharge(int kf) {
    return jamchge_(&kf);
  }

  double jamMass(int kf) {
    return pjmass_(&kf);
  }

  void jamList(int i) {
    jamlist_(&i);
  }

  void  finalResonanceDecay() {
    jamfdec_();
  }

  //---------------------------------------------------------------------------
  // Utilities

  void printParticleInformation(int i) {
    std::cout
      << i
      << " KF= " << getK(2, i)
      << std::endl
      << " K= "   << getK(1, i)
      << std::setw(5) << getK(2, i)
      << std::setw(5) << getK(3, i)
      << std::setw(5) << getK(4, i)
      << std::setw(5) << getK(5, i)
      << std::setw(5) << getK(6, i)
      << std::setw(5) << getK(7, i)
      << std::setw(5) << getK(8, i)
      << std::setw(5) << getK(9, i)
      << std::setw(5) << getK(10, i)
      << std::setw(5) << getK(11, i)
      << std::endl
      << " R= "   << getR(1, i)
      << std::setw(10) << getR(2, i)
      << std::setw(10) << getR(3, i)
      << std::setw(10) << getR(4, i)
      << std::setw(10) << getR(5, i)
      << std::endl
      << " P= "   << getP(1, i)
      << std::setw(15) << getP(2, i)
      << std::setw(15) << getP(3, i)
      << std::setw(15) << getP(4, i)
      << std::setw(15) << getP(5, i)
      << std::endl
      << " V= "   << getV(1, i)
      << std::setw(10) << getV(2, i)
      << std::setw(10) << getV(3, i)
      << std::setw(10) << getV(4, i)
      << std::setw(10) << getV(5, i)
      << std::endl;
  }

}
}
