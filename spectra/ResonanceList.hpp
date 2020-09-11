// -*- mode: c++ -*-
#ifndef runjam_spectra_IResonanceData_hpp
#define runjam_spectra_IResonanceData_hpp
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
    virtual ~ResonanceList() {}
  };

  class ResonanceListPCE: public ResonanceList {
    typedef ResonanceList base;
  private:
    static ResonanceRecord resT[5][21];

    void initialize(int kineticTemp, int eos_pce,std::string const& fn_resodata);
  public:
    ResonanceListPCE(runjam_context const& ctx);
    ResonanceListPCE(int kineticTemp, int eos_pce,std::string const& fn_resodata);
  };

}
}

#endif
