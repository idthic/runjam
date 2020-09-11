// -*- mode: c++ -*-
#ifndef runjam_spectra_ParticleSample_hpp
#define runjam_spectra_ParticleSample_hpp
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <args.hpp>

namespace idt {
namespace runjam {

  struct Particle {
    int pdg;
    double mass;

    double x;
    double y;
    double z;
    double t;

    double e ;
    double px;
    double py;
    double pz;
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
    void addParticleMinkowski(int iReso, double px, double py, double pz, double m, double x, double y, double z, double t);

    /// @param[in] px [GeV/c] momentum in the x-direction
    /// @param[in] py [GeV/c] momentum in the y-direction
    /// @param[in] pz [GeV/c] momentum in the z-direction
    /// @param[in] m  [GeV/c] mass 粒子番号から計算する場合は -1.0 を指定する。
    void addParticleTauEta(int iReso, double px, double py, double pz, double m, double x, double y, double tau, double eta);

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

    /// @fn void update();
    /// \~en generates a resonance distribution
    /// \~ja 粒子分布の生成を実行します。
    virtual void update() = 0;

  private:
    // コピー禁止
    ParticleSampleBase(ParticleSampleBase const&) {
      std::cerr << "ParticleSampleBase(copy ctor): copy not supported!" << std::endl;
      std::exit(1);
    }
    ParticleSampleBase& operator=(ParticleSampleBase const&) {
      std::cerr << "ParticleSampleBase(copy assign): copy not supported!" << std::endl;
      std::exit(1);
    }
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
    OversampledParticleSampleBase(runjam_context const& ctx) {
      this->m_overSamplingFactor = ctx.get_config("runjam_oversampling_factor", 1.0);
      this->setAdviceNumberOfExpectedEvents(0);
    }

  public:
    virtual void updateWithOverSampling(double overSamplingFactor) = 0;

    virtual void update() {
      // 一括生成済の時
      if (this->pcache.size() > 1) {
        if (++this->indexOfCachedEvents < this->pcache.size() - 1)
          return;

        this->pcache.clear();
        this->indexOfCachedEvents = -1;
      }

      // 一括生成要求がある時
      if (this->numberOfExpectedEvents > 0) {
        int ncache = this->numberOfExpectedEvents;
        this->updateWithOverSampling(this->m_overSamplingFactor * ncache);

        // 二項分布で各イベントの粒子数を決定する
        this->pcache.clear();
        this->pcache.reserve(ncache + 1);
        this->pcache.emplace_back(0);
        std::size_t pos = 0, size = base::plist.size();
        for (; ncache > 1; ncache--) {
          pos += idt::util::irand_binomial(size - pos, 1.0 / ncache);
          this->pcache.emplace_back(pos);
        }
        this->pcache.emplace_back(size);

        this->shuffleParticleList();
        this->numberOfExpectedEvents = 0;
        this->indexOfCachedEvents = 0;
        return;
      }

      this->updateWithOverSampling(this->m_overSamplingFactor);
    }
  };

  class ParticleSampleFactoryBase {
  public:
    virtual ParticleSampleBase* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) = 0;
    virtual ~ParticleSampleFactoryBase() {}
  protected:
    ParticleSampleFactoryBase() { Register(this); }
    void Register(ParticleSampleFactoryBase* factory);
  };

  ParticleSampleBase* CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile);

}
}

#endif
