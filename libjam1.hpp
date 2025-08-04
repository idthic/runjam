/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2006, Yasushi Nara, 2011-2020, Koichi Murase.

   SPDX-License-Identifier: GPL-2.0-or-later

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA  */

#ifndef runjam_libjam1_hpp
#define runjam_libjam1_hpp
#include <string>

namespace idt {
namespace libjam1 {

  int    getK(int i, int ip);
  double getR(int i, int ip);
  double getP(int i, int ip);
  double getV(int i, int ip);
  void   setK(int i, int ip, int k);
  void   setR(int i, int ip, double r);
  void   setP(int i, int ip, double p);
  void   setV(int i, int ip, double v);

  int    getNV();
  int    getNBARY();
  int    getNMESON();
  void   setNV(int n);
  void   setNBARY(int n);
  void   setNMESON(int n);

  int    getMSTC(int i);
  double getPARC(int i);
  int    getMSTD(int i);
  double getPARD(int i);
  void   setMSTC(int i, int m   );
  void   setPARC(int i, double p);
  void   setMSTD(int i, int m   );
  void   setPARD(int i, double p);
  void   addPARD(int i, double p);

  char*  getFNAME(int i);
  void   setFNAME(int i, const char* v);
  void   setFNAME(int i, std::string const& v);

  int    getKCHG(int ip, int i);
  double getPMAS(int ip, int i);
  double getPARF        (int i);
  double getVCKM(int i,  int j);
  void   setKCHG(int ip, int i, int k   );
  void   setPMAS(int ip, int i, double m);
  void   setPARF        (int i, double p);
  void   setVCKM (int i, int j, double v);
  int    getBaryonNumber(int kc, int kf);

  int    getMDCY(int i, int j);
  int    getMDME(int i, int j);
  double getBRAT       (int i);
  int    getKFDP(int i, int j);
  void   setMDCY(int i, int j, int m);
  void   setMDME(int i, int j, int m);
  void   setBRAT(int i, double b);
  void   setKFDP(int i, int j, int k);

  void   jamInit(int mev, double bmin, double bmax, double dt, int nstp,
    const char* frame, const char* proj, const char* targ, const char* cwin);
  void   jamEvt(int iev);
  void   jamFin();
  int    jamComp(int kf);
  void   jamZero(int i);
  double jamDecayTime(int im,int kf, int kc, int ks, double em,double ee);
  void   jamList(int i);
  double jamMass(int kf);
  int    jamCharge(int kf);
  void   finalResonanceDecay();

  void printParticleInformation(int i);
  int determineStableCode(int kf);

}
}

#endif
