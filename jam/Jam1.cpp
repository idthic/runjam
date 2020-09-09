#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "Jam1.hpp"

namespace idt {
namespace runjam {

Jam1::Jam1() {}

Jam1::~Jam1() {}

void Jam1::jamInit(int mev, double bmin, double bmax, double dt, int nstp,
  const char *frame, const char *proj, const char *targ, const char *cwin) {
  int s1 = std::strlen(frame);
  int s2 = std::strlen(proj);
  int s3 = std::strlen(targ);
  int s4 = std::strlen(cwin);
  jaminit_(&mev, &bmin,&bmax, &dt, &nstp, frame, proj, targ, cwin, s1, s2, s3, s4);
}

void Jam1::print(int i) {
  std::cout << i
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

void Jam1::jamEvt(int iev) {
  jamevt_(&iev);
}

void Jam1::jamFin() {
  jamfin_();
}

int Jam1::jamComp(int kf) {
  return jamcomp_(&kf);
}

void Jam1::jamZero(int i) {
  jamzero_(&i);
}

double Jam1::jamDecayTime(int im, int kf, int kc, int ks, double em, double ee) {
  return jamdtim_(&im, &kf, &kc, &ks, &em, &ee);
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


static std::size_t fort_strlen(const char* buff, int bufferSize = -1) {
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
  std::size_t len = std::min(szdst, fort_strlen(src, szsrc));
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

char* Jam1::getFNAME(int i) const {
  static char buf[81];
  strcpy_f2c(buf, 80, jamdat3_.FNAME[i - 1], 80);
  //strncpy(buf, jamdat3_.FNAME[i - 1], 80);
  buf[80] = '\0';
  return buf;
}

void  Jam1::setFNAME(int i, const char* v) {
  if (!strcpy_c2f(jamdat3_.FNAME[i - 1], 80, v)) {
    std::cerr << "runjam: filename \"" << v << "\" is too long" << std::endl;
    std::exit(1);
  }
  //strncpy(jamdat3_.FNAME[i - 1], v, std::strlen(v));
}

void  Jam1::finalResonanceDecay() {
  jamfdec_();
}

}
}
