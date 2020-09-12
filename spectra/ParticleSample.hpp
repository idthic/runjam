// -*- mode: c++ -*-
#ifndef runjam_spectra_ParticleSample_hpp
#define runjam_spectra_ParticleSample_hpp
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <args.hpp>

namespace idt {
namespace runjam {

  struct vector4 {
    double data[4];

  public:
    double& operator[](std::size_t index) { return data[index]; }
    double const& operator[](std::size_t index) const { return data[index]; }

    vector4 operator+() const { return *this; }
    vector4 operator-() const { return {-data[0], -data[1], -data[2], -data[3]}; }
    vector4 operator~() const { return {data[0], -data[1], -data[2], -data[3]}; }

    double operator*(const vector4& rhs) const {
      return data[0] * rhs.data[0] - data[1] * rhs.data[1] - data[2] * rhs.data[2] - data[3] * rhs.data[3];
    }

    vector4& operator*=(double scalar) {
      data[0] *= scalar;
      data[1] *= scalar;
      data[2] *= scalar;
      data[3] *= scalar;
      return *this;
    }
    vector4& operator/=(double scalar) {
      return *this *= 1.0 / scalar;
    }

    vector4& operator+=(const vector4& rhs) {
      data[0] += rhs.data[0];
      data[1] += rhs.data[1];
      data[2] += rhs.data[2];
      data[3] += rhs.data[3];
      return *this;
    }
    vector4& operator-=(const vector4& rhs) {
      data[0] -= rhs.data[0];
      data[1] -= rhs.data[1];
      data[2] -= rhs.data[2];
      data[3] -= rhs.data[3];
      return *this;
    }
    vector4 operator+(vector4 const& rhs) const {
      vector4 ret(*this);
      ret += rhs;
      return ret;
    }
    vector4 operator-(vector4 const& rhs) const {
      vector4 ret(*this);
      ret -= rhs;
      return ret;
    }

    void boost(vector4 const& u) {
      vector4& p = *this;
      double const udotp = u[1] * p[1] + u[2] * p[2] + u[3] * p[3];
      double const coeff = udotp / (u[0] + 1) + p[0];
      p[0] = u[0] * p[0] + udotp;
      p[1] += u[1] * coeff;
      p[2] += u[2] * coeff;
      p[3] += u[3] * coeff;
    }
  };

  struct Particle {
    int pdg;
    double mass;

    union {
      vector4 pos;
      struct {
        double t;
        double x;
        double y;
        double z;
      };
    };

    union {
      vector4 mom;
      struct {
        double e ;
        double px;
        double py;
        double pz;
      };
    };
  };

  class ParticleSampleBase {
  protected:
    std::vector<Particle> plist;

  protected:
    void clearParticleList();
    void shuffleParticleList();

    /// @param[in] px [GeV/c] momentum in the x-direction
    /// @param[in] py [GeV/c] momentum in the y-direction
    /// @param[in] pz [GeV/c] momentum in the z-direction
    /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
    void addParticleCartesian(int iReso, double px, double py, double pz, double m, double x, double y, double z, double t);

    /// @param[in] px [GeV/c] momentum in the x-direction
    /// @param[in] py [GeV/c] momentum in the y-direction
    /// @param[in] pz [GeV/c] momentum in the z-direction
    /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
    void addParticleMilne(int iReso, double px, double py, double pz, double m, double x, double y, double tau, double eta);

  public:
    ParticleSampleBase() {}

    virtual ~ParticleSampleBase() { this->clearParticleList(); }

    /// @fn void setAdviceNumberOfExpectedEvents(int nEvents);
    /// このインスタンスを用いて生成すると予想されるイベントの数を指定します。
    virtual void setAdviceNumberOfExpectedEvents(int nEvents) {}

    virtual double getOverSamplingFactor() const = 0;

    /// @fn Particle* begin();
    /// @fn Particle* end();
    /// @fn Particle const* begin() const;
    /// @fn Particle const* end() const;
    /// \~en retrieves the generated particle iterators.
    /// \~ja 生成した粒子集合のイテレータを取得します。
    virtual Particle* begin() { return &plist[0]; }
    virtual Particle* end() { return &plist[0] + plist.size(); }
    Particle const* begin() const { return const_cast<ParticleSampleBase*>(this)->begin(); }
    Particle const* end() const { return const_cast<ParticleSampleBase*>(this)->end(); }
    std::size_t size() const { return end() - begin(); }

    /// @fn void update();
    /// \~en generates a resonance distribution
    /// \~ja 粒子分布の生成を実行します。
    virtual void update() = 0;

    void adjustCenterOfMassByGalileiBoost();
    void adjustCenterOfMassByLorentzBoost();
  };

  class OversampledParticleSampleBase: public ParticleSampleBase {
    typedef ParticleSampleBase base;

    // 実装: 複数事象一括生成について。
    //
    // 1 numberOfExpectedEvents が有限の値に設定されている時、一括生成
    //   が要求されている事を意味する。numberOfExpectedEvents は
    //   setAdviceNumberOfExpectedEvents を通して設定できる。一括生成
    //   が要求されている時に update が呼ばれると一括生成が実行され、
    //   numberOfExpectedEvents は 0 にクリアされる。
    //

    // 2 一括生成された粒子は base::plist に保持される。事象 #ievent
    //   の粒子一覧は base::plist 内の一区間として表現される。その範囲
    //   は [prange[ievent], prange[ievent]) になる。"pcached.size() >
    //   0" の時、一括生成された事象が未だ残っている事を表す。
    //   ("pcache.size() > 0" の間 base::plist には一括生成された全粒
    //   子が格納されている事になる。)

    //
    // 3 pcache の事象を使い切ると pcache はクリアされる。この場合は通
    //   常の 1 事象の生成が行われる。その過程で、今迄一括生成の全粒子
    //   plist も解放・クリアされる。

  private:
    double m_overSamplingFactor;
  public:
    void setOverSamplingFactor(double value) {
      this->m_overSamplingFactor = value;
    }
    virtual double getOverSamplingFactor() const override {
      return this->m_overSamplingFactor;
    }

  private:
    int numberOfExpectedEvents;
  public:
    virtual void setAdviceNumberOfExpectedEvents(int nEvents) override {
      this->numberOfExpectedEvents = nEvents;
      this->indexOfCachedEvents = -1;
    }

  private:
    int indexOfCachedEvents;
    std::vector<std::size_t> pcache;
  public:
    virtual Particle* begin() override {
      if (this->indexOfCachedEvents >= 0)
        return base::begin() + pcache[indexOfCachedEvents];
      else
        return base::begin();
    }
    virtual Particle* end() override {
      if (this->indexOfCachedEvents >= 0)
        return base::begin() + pcache[indexOfCachedEvents + 1];
      else
        return base::end();
    }

  public:
    OversampledParticleSampleBase(runjam_context const& ctx):
      m_overSamplingFactor(ctx.get_config("runjam_oversampling_factor", 1.0))
    {
      this->setAdviceNumberOfExpectedEvents(0);
    }

  public:
    virtual void updateWithOverSampling(double overSamplingFactor) = 0;
    virtual void update() override;
  };

  class ParticleSampleFactoryBase {
  public:
    virtual std::unique_ptr<ParticleSampleBase> CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) = 0;
    virtual ~ParticleSampleFactoryBase() {}
  protected:
    ParticleSampleFactoryBase() { Register(this); }
    void Register(ParticleSampleFactoryBase* factory);
  };

  std::unique_ptr<ParticleSampleBase> CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile);

}
}

#endif
