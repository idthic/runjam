/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2011-2020, Koichi Murase <myoga.murase at gmail.com>

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

#ifndef runjam_ResonanceList_hpp
#define runjam_ResonanceList_hpp
#include <vector>
#include <string>
#include <args.hpp>

namespace idt {
namespace runjam {

  struct ResonanceRecord {
    double mass;   //!< 共鳴の質量。単位は [fm^{-1}]
    double deg;    //!< 共鳴の自由度 (isospin, spin, etc.) の数 g
    double degeff; //!< 崩壊させた後の pion の数 * g。
    double mu;     //!< 共鳴の化学ポテンシャル
    int bf;        //!< 共鳴の統計符号。Boson に対して -1、fermion に対して +1
    int anti;      //!< 反粒子の時に 1。それ以外の時 0
    std::string key;
    std::vector<int> pdg_codes;

  public:
    //! PDG Monte-Carlo code を生成します。
    int generatePDGCode() const;
  };

  class ResonanceList {
  protected:
    std::vector<ResonanceRecord> data;
  public:
    ResonanceRecord const& operator[](int ireso) const {
      return this->data[ireso];
    }
    ResonanceRecord& operator[](int ireso) {
      return this->data[ireso];
    }
    std::vector<ResonanceRecord>::iterator begin() { return data.begin(); }
    std::vector<ResonanceRecord>::iterator end() { return data.end(); }
    std::vector<ResonanceRecord>::const_iterator begin() const { return data.begin(); }
    std::vector<ResonanceRecord>::const_iterator end() const { return data.end(); }
    std::size_t size() const{ return this->data.size(); }

  public:
    void readFile(std::string const& fn_resodata);

  public:
    ResonanceList() {}
    ResonanceList(std::string const& fn_resodata) { readFile(fn_resodata); }
    ResonanceList(runjam_context const& ctx) { readFile(ctx.resodata()); }
  };

}
}

#endif
