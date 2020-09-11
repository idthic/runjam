// -*- mode: c++ -*-
#ifndef runjam_spectra_IParticleSample_hpp
#define runjam_spectra_IParticleSample_hpp
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

  class IParticleSample {
  public:
    virtual ~IParticleSample() {}

    virtual void setAdviceNumberOfExpectedEvents(int nEvents) {}

    /// @fn void update();
    /// \~en generates a resonance distribution
    /// \~ja 粒子分布の生成を実行します。
    virtual void update() = 0;

    /// @fn std::vector<Particle*> const& getParticleList() const;
    /// \~en retrieves the generated resonance distribution.
    /// \~ja 生成した粒子分布を取得します。
    virtual std::vector<Particle*> const& getParticleList() const = 0;
  };

  class ParticleSampleBase: public IParticleSample {
  protected:
    std::vector<Particle*> plist;

  protected:
    void clearParticleList();

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
    virtual std::vector<Particle*> const& getParticleList() const { return this->plist; }

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

  private:
    double m_overSamplingFactor;
  public:
    void setOverSamplingFactor(double value) {
      this->m_overSamplingFactor = value;
    }
    double getOverSamplingFactor() const {
      return this->m_overSamplingFactor;
    }

    // 実装: 複数事象一括生成について。
    //
    // 1 numberOfExpectedEvents が有限の値に設定されている時、一括生成
    //   が要求されている事を意味する。numberOfExpectedEvents は
    //   setAdviceNumberOfExpectedEvents を通して設定できる。一括生成
    //   が要求されている時に update が呼ばれると一括生成が実行され、
    //   numberOfExpectedEvents は 0 にクリアされる。
    //
    // 2 一括生成された粒子は base::plist に保持され寿命が管理される。
    //   同時に、事象 #ievent の粒子一覧は pcache[ievent] に記録される。
    //   pcache.size()>0 の時、一括生成された事象が未だ残っている事を
    //   表す。(pcache.size()>0 の間 base::plist には一括生成された全
    //   粒子が格納されている事になる。)
    //
    // 3 pcache の事象を使い切ると pcache はクリアされる。この場合は通
    //   常の 1 事象の生成が行われる。その過程で、今迄一括生成の全粒子
    //   plist も解放・クリアされる。

  private:
    int numberOfExpectedEvents;
    int indexOfCachedEvents;
    std::vector<std::vector<Particle*> > pcache;

  public:
    virtual std::vector<Particle*> const& getParticleList() const {
      if (this->indexOfCachedEvents >= 0)
        return this->pcache[indexOfCachedEvents];
      else
        return this->plist;
    }

  public:
    virtual void setAdviceNumberOfExpectedEvents(int nEvents) override {
      this->numberOfExpectedEvents = nEvents;
      this->indexOfCachedEvents = -1;
    }

  public:
    OversampledParticleSampleBase(runjam_context const& ctx) {
      this->m_overSamplingFactor = ctx.get_config("runjam_oversampling_factor", 1.0);
      this->setAdviceNumberOfExpectedEvents(0);
    }

    virtual void updateWithOverSampling(double overSamplingFactor) = 0;

    virtual void update() {
      // 一括生成済の時
      if (this->pcache.size() > 0) {
        if (++this->indexOfCachedEvents < this->pcache.size())
          return;

        this->pcache.clear();
        this->indexOfCachedEvents = -1;
      }

      // 一括生成要求がある時
      if (this->numberOfExpectedEvents > 0) {
        int ncache = this->numberOfExpectedEvents;
        this->updateWithOverSampling(this->m_overSamplingFactor * ncache);
        this->pcache.resize(ncache, std::vector<Particle*>());
        for (std::vector<Particle*>::const_iterator i = this->base::plist.begin(); i != this->base::plist.end(); ++i)
          this->pcache[idt::util::irand(ncache)].push_back(*i);

        this->numberOfExpectedEvents = 0;
        this->indexOfCachedEvents = 0;
        return;
      }

      this->updateWithOverSampling(this->m_overSamplingFactor);
    }
  };

  class IParticleSampleFactory {
  public:
    virtual IParticleSample* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) = 0;
    virtual ~IParticleSampleFactory() {}
  protected:
    void Register(IParticleSampleFactory* factory);
  };

  class ParticleSampleFactoryRegistered: public IParticleSampleFactory {
  protected:
    ParticleSampleFactoryRegistered() { Register(this); }
  };

  IParticleSample* CreateParticleSample(runjam_context const& ctx, std::string const& type, std::string const& inputfile);

}
}

#endif
