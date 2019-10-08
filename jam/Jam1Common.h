#ifndef Jam1Common_h
#define Jam1Common_h

// This number should be the same as "mxv" in jam1.inc
#define JAM_MXV  200000
//#define JAM_MXV  100000
//#define JAM_MXV  50000
//#define JAM_MXV  30000

#ifdef HYDROJET_FORTRAN_NOEXTNAME
#  define jamevnt1_ jamevnt1
#  define jamevnt2_ jamevnt2
#  define jamdat1_  jamdat1
#  define jamdat2_  jamdat2
#  define jamdat3_  jamdat3
#  define jydat2_   jydat2
#  define jydat3_   jydat3
#  define jydat4_   jydat4
#endif

struct Jamevnt1{
  double R[JAM_MXV][5];
  double P[JAM_MXV][5];
  double V[JAM_MXV][5];
  int    K[JAM_MXV][11];
};
extern Jamevnt1 jamevnt1_;

struct Jamevnt2{
    int NV;
    int NBARY;
    int NMESON;
};
extern Jamevnt2 jamevnt2_;

struct Jamdat1 {
  int    MSTC[200];
  double PARC[200];
  int    MSTD[200];
  double PARD[200];
};
extern Jamdat1 jamdat1_;

struct Jamdat2 {
  int    MSTE[200];
  double PARE[200];
};
extern Jamdat2 jamdat2_;

struct Jamdat3 {
    char  FNAME[8][80];
};
extern Jamdat3 jamdat3_;

struct Jydat2 {
  int    KCHG[7][500];
  double PMAS[4][500];
  double PARF[2000];
  double VCKM[4][4];
};
extern Jydat2 jydat2_;

int   const KNDCAY1  =  4000; //should be 8000 for pythia62

struct Jydat3 {
  int    MDCY[3][500];
  int    MDME[3][KNDCAY1];
  double BRAT[KNDCAY1];
  int    KFDP[5][KNDCAY1];
};
extern Jydat3 jydat3_;

struct Jydat4 {
  char  CHAF[2][500][16];	// here I needed manual intervention
};
extern Jydat4 jydat4_;

#endif
