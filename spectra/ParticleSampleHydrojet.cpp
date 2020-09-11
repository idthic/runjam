#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef _OPENMP
# include <omp.h>
#endif

#include <util.hpp>
#include <ksh/integrator.hpp>
#include <ksh/phys/Minkowski.hpp>
#include <spectra/IResonanceList.hpp>
#include <spectra/IParticleSample.hpp>
#include <spectra/IntegratedCooperFrye.hpp>

#include "HydroSpectrum.hpp"

using namespace idt;
using namespace idt::runjam;

namespace {

  static const std::size_t HYDRO2JAM_ITERATION_MAX = 20;
  static const double HYDRO2JAM_TEMPERATURE_MIN = 0.01; // 2014-07-30 (unit: [/fm])

  // static const double HYDRO2JAM_FACRANMAX = 1.6;
  // static const double HYDRO2JAM_FACRANMAX = 1.7;// 2010-06-28, lower switching T, larger radial flow
  // static const double HYDRO2JAM_FACRANMAX = 1.8;// 2019-08-27, LHC, larger radial flow
  static const double HYDRO2JAM_FACRANMAX = 2.0; // 2014-07-30 for RFH

  class ParticleSampleHydrojet: public OversampledParticleSampleBase {
    typedef OversampledParticleSampleBase base;
  private:
    ResonanceListPCE rlist;
    HydroSpectrum m_hf;

  private:
    bool flag_negative_contribution = false;
    bool cfg_reverse_particles;
    bool cfg_shuffle_particles;

  private:
    std::string fn_freezeout_dat;
    std::string fn_position_dat;
  public:
    void setHypersurfaceFilenames(std::string const& fn_freezeout_dat, std::string const& fn_position_dat) {
      this->fn_freezeout_dat = fn_freezeout_dat;
      this->fn_position_dat = fn_position_dat;
    }

  private:
    int    cfg_baryon_disable;
    double cfg_baryon_tmpf;  // [unit: fm^{-1}]
    double cfg_baryon_mubf;  // [unit: fm^{-1}]
    double cfg_baryon_meanf;
  public:
    void setBaryonFree(int i) { cfg_baryon_disable = i; }
    void setTMPF(double t) { cfg_baryon_tmpf = t / hbarc_GeVfm; }
    void setMUBF(double m) { cfg_baryon_mubf = m / hbarc_GeVfm; }

  private:
    double dx, dy, dh, dtau;
  public:
    void setDx(double d) { dx = d; }
    void setDy(double d) { dy = d; }
    void setDh(double d) { dh = d; }
    void setDtau(double d) { dtau = d; }
    double getDtau() { return dtau; }
    double getDx() { return dx; }
    double getDy() { return dy; }
    double getDh() { return dh; }

  private:
    bool cache_available;
    std::vector<std::string>   cache_fname;
    std::vector<std::ifstream> cache_ifsPos;
    std::vector<std::ifstream> cache_ifsNeg;
    bool openCooperFryeCacheForRead();
    void createCooperFryeCache(int ireso);
    void createCooperFryeCache();

  public:
    ParticleSampleHydrojet(runjam_context const& ctx, std::string const& cachedir, std::string suffix,
      int kin, int eos_pce, std::string const& fn_resodata);

  private:
    void initialize();
    void analyze(double oversamplingFactor);
    void finalize();
    void generateParticle(
      double vx, double vy, double vz,
      double ds0, double dsx, double dsy,
      double dsz, int ir, int ipos,
      double tau, double xx, double yy, double eta, int reflection = 0);
    void updateWithOverSampling(double oversamplingFactor) {
      this->base::clearParticleList();
      this->initialize();
      this->analyze(oversamplingFactor);
      this->finalize();
    }
    virtual void update() override {
      this->updateWithOverSampling(1.0);
    }
  };

  ParticleSampleHydrojet::ParticleSampleHydrojet(
    runjam_context const& ctx, std::string const& cachedir, std::string suffix,
    int kintmp, int eos_pce, std::string const& fn_resodata
  ): base(ctx), rlist(kintmp, eos_pce, fn_resodata), m_hf(kintmp, eos_pce) {
    int const nreso_loop = this->rlist.numberOfResonances();
    this->cache_fname.resize(nreso_loop);
    for (int i = 0; i < nreso_loop; i++) {
      std::ostringstream sstr;
      if (cachedir.size() > 0) sstr << cachedir << "/";
      sstr << "ELEMENT." << rlist[i].key << suffix;
      cache_fname[i] = sstr.str();
    }

    cache_available = false;

    // constants
  	cfg_baryon_tmpf = 0.16 / hbarc_GeVfm;
  	cfg_baryon_mubf = 1.6 / hbarc_GeVfm;
  	cfg_baryon_meanf = 0.45 / hbarc_GeVfm;

    // 2013/04/23, KM, reverse z axis
    this->cfg_reverse_particles = ctx.get_config("hydrojet_reverse_particles", false);
    if (this->cfg_reverse_particles)
      std::cout << "ParticleSampleHydrojet: ReverseParticleList mode enabled!" << std::endl;

    // 2013/04/30, KM, shuffle the particle list
    this->cfg_shuffle_particles = ctx.get_config("hydrojet_shuffle_particles", false);
    if (this->cfg_shuffle_particles)
      std::cout << "ParticleSampleHydrojet: ShuffleParticleList mode enabled!" << std::endl;
  }

  // 一度に 151x2 のファイルを開いて書き込むとディスクに悪いので共鳴毎に処理する。
  void ParticleSampleHydrojet::createCooperFryeCache(int ireso) {
    //---------------------------------------------------------------------------
    // initialize files

    // open output files
    std::ofstream ostr_elm;
    if (ireso < 21) {
      ostr_elm.open((cache_fname[ireso]).c_str());
      if (!ostr_elm) {
        std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << cache_fname[ireso] << std::endl;
        std::exit(1);
      }
    }
    std::ofstream ostr_pos((cache_fname[ireso] + ".POS").c_str());
    if (!ostr_pos) {
      std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << cache_fname[ireso] << ".POS" << std::endl;
      std::exit(1);
    }
    std::ofstream ostr_neg((cache_fname[ireso] + ".NEG").c_str());
    if (!ostr_neg) {
      std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << cache_fname[ireso] << ".NEG" << std::endl;
      std::exit(1);
    }

    ResonanceListPCE::resonance& recreso = this->rlist[ireso];

    m_hf.openFDataFile(fn_freezeout_dat);

    //---------------------------------------------------------------------------
    while (!m_hf.readFData()) {
      if (!cfg_baryon_disable) {
        recreso.mu = 0.0;
        if (recreso.bf == 1) {
          double const mub = cfg_baryon_mubf * sqrt(1.0 - m_hf.tf * m_hf.tf / cfg_baryon_tmpf / cfg_baryon_tmpf);
          double const mean_field = cfg_baryon_meanf * m_hf.nbf;
          recreso.mu   = mub - mean_field;
          if (recreso.anti) recreso.mu = -mub + mean_field;
          // if (recreso.anti) recreso.mu = -mub - mean_field;
        }
      }

      if (m_hf.tf == 0.0) {
        std::cerr << "(ParticleSampleHydrojet::createCooperFryeCache) TF=ZERO" << std::endl;
        //write(24, 9000)0.0
        continue;
      }

      //-------------------------------------------------------------------------
      // integration

      double npos = 0.0;
      double nneg = 0.0;

      double const gamma = 1.0 / sqrt(1.0 - m_hf.vx * m_hf.vx - m_hf.vy * m_hf.vy - m_hf.vz * m_hf.vz);
      kashiwa::phys::vector4 u(gamma, m_hf.vx * gamma, m_hf.vy * gamma, m_hf.vz * gamma);
      kashiwa::phys::vector4 ds(m_hf.ds0, -m_hf.dsx, -m_hf.dsy, -m_hf.dsz);
      double const beta = 1.0 / m_hf.tf;

      if (recreso.bf == -1) {
        idt::runjam::IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
      } else {
        idt::runjam::IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
      }

      double const n = (nneg + npos) * recreso.degeff;
      npos *= recreso.deg;
      nneg *= recreso.deg;
      //-------------------------------------------------------------------------

      if (ireso < 21) ostr_elm << n << "\n";
      ostr_pos << npos << "\n";
      ostr_neg << nneg << "\n";
    }

    //  close files.
    if (ireso < 21) ostr_elm.close();
    ostr_pos.close();
    ostr_neg.close();

    m_hf.closeFDataFile();
  }

  void ParticleSampleHydrojet::createCooperFryeCache() {
    int const nreso_loop = this->rlist.numberOfResonances();
    for(int ireso = 0; ireso < nreso_loop; ireso++)
      this->createCooperFryeCache(ireso);
  }

  bool ParticleSampleHydrojet::openCooperFryeCacheForRead() {
    int const nreso_loop = rlist.numberOfResonances();

    cache_ifsPos.clear();
    for (int i = 0; i < nreso_loop; i++) {
      std::string fnpos = cache_fname[i] + ".POS";
      cache_ifsPos.emplace_back(fnpos.c_str());
      if (!cache_ifsPos.back()) goto failed;
    }

    cache_ifsNeg.clear();
    for (int i = 0; i < nreso_loop; i++) {
      std::string fnneg = cache_fname[i] + ".NEG";
      cache_ifsNeg.emplace_back(fnneg.c_str());
      if (!cache_ifsNeg.back()) goto failed;
    }
    return true;

  failed:
    cache_ifsPos.clear();
    cache_ifsNeg.clear();
    return false;
  }

  void ParticleSampleHydrojet::initialize() {
    // Open files for input.
    m_hf.openFDataFile(this->fn_freezeout_dat);
    m_hf.openPDataFile(this->fn_position_dat);

    std::cout << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): checking Cooper-Frye cache files (.POS/.NEG)... " << std::flush;
    if (this->openCooperFryeCacheForRead()) {
      cache_available = true;
      std::cout << "yes" << std::endl;
      return;
    } else {
      std::cout << "no (incomplete).\n";
      std::cout << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): entering delayed Cooper-Frye evaluation mode." << std::endl;
    }
  }

  void ParticleSampleHydrojet::analyze(double oversamplingFactor) {
    int const nreso_loop = rlist.numberOfResonances();

    double ran;

    while (!m_hf.readFData()) {
      m_hf.readPData();

      if (!cfg_baryon_disable) {
        for (int i = 0; i < nreso_loop; i++) {
          rlist[i].mu = 0.0;
          if (rlist[i].bf == 1) {
            double const mub = cfg_baryon_mubf * sqrt(1.0 - m_hf.tf * m_hf.tf / cfg_baryon_tmpf / cfg_baryon_tmpf);
            double const mean_field = cfg_baryon_meanf * m_hf.nbf;
            rlist[i].mu   = mub - mean_field;
            if (rlist[i].anti) rlist[i].mu = -mub + mean_field;
            // if (rlist[i].anti) rlist[i].mu = -mub - mean_field;
          }
        }
      }

      if (m_hf.tf < HYDRO2JAM_TEMPERATURE_MIN) {
        if (m_hf.tf == 0.0)
          std::cerr << "ParticleSampleHydrojet! (warning) zero temperature surface." << std::endl;
        continue;
      }

      // Loop over all particles.
      for (int ir = 0; ir < nreso_loop; ir++) {
        double numResPos, numResNeg;
        if (cache_available) {
          cache_ifsPos[ir] >> numResPos;
          if (numResPos > 1.0)
            std::cout << "Suspicious fluid element! ireso =" << ir
                      << " " << numResPos << std::endl;

          cache_ifsNeg[ir] >> numResNeg;
          if (numResNeg > 1.0)
            std::cout << "Suspicious fluid element! ireso =" << ir << std::endl;
        } else {
          double npos = 0.0;
          double nneg = 0.0;

          double gamma = 1.0 / sqrt(1.0 - m_hf.vx * m_hf.vx - m_hf.vy * m_hf.vy - m_hf.vz * m_hf.vz);
          kashiwa::phys::vector4 u(gamma, m_hf.vx * gamma, m_hf.vy * gamma, m_hf.vz * gamma);
          kashiwa::phys::vector4 ds(m_hf.ds0, -m_hf.dsx, -m_hf.dsy, -m_hf.dsz);
          double beta = 1.0 / m_hf.tf;

          if (rlist[ir].bf == -1) {
            idt::runjam::IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, rlist[ir].mass, rlist[ir].mu);
          } else {
            idt::runjam::IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, rlist[ir].mass, rlist[ir].mu);
          }

          double n = (nneg + npos) * rlist[ir].degeff;
          numResPos = npos * rlist[ir].deg;
          numResNeg = nneg * rlist[ir].deg;
        }
        numResPos *= oversamplingFactor;
        numResNeg *= oversamplingFactor;

        int reflection_count, reflection_step;
        switch (m_hf.iw % 4) {
        case 1: // iw = 1, 5
          // reflection = 0, 1, 2, 3
          reflection_count = 4;
          reflection_step = 1;
          break;
        case 2: // iw = 2, 6
          // reflection = 0, 2
          reflection_count = 2;
          reflection_step = 2;
          break;
        case 3: // iw = 3, 7
          // reflection = 0, 1
          reflection_count = 2;
          reflection_step = 1;
          break;
        default: // iw = 4, 8
          // reflection = 0
          reflection_count = 1;
          reflection_step = 1;
          break;
        }

        if (m_hf.bulk == 1) {
          m_hf.ds0 = m_hf.dss * std::cosh(m_hf.hh);
          m_hf.dsz = -m_hf.dss * std::sinh(m_hf.hh);
        } else {
          m_hf.ds0 = m_hf.dss * std::sinh(m_hf.hh);
          m_hf.dsz = -m_hf.dss * std::cosh(m_hf.hh);
        }

        // positive contribution
        {
          int const ipos = 1;
          int const n = idt::util::irand_poisson(numResPos * reflection_count);
          for (int i = 0; i < n; i++) {
            int const reflection = idt::util::irand(reflection_count) * reflection_step;
            generateParticle(m_hf.vx, m_hf.vy, m_hf.yv, m_hf.ds0, m_hf.dsx, m_hf.dsy, m_hf.dsz, ir, ipos, m_hf.tau, m_hf.xx, m_hf.yy, m_hf.eta, reflection);
          }
        }

        // negative contribution
        if (flag_negative_contribution) {
          int const ipos = 0;
          int const n = idt::util::irand_poisson(numResNeg * reflection_count);
          for (int i = 0; i < n; i++) {
            int const reflection = idt::util::irand(reflection_count) * reflection_step;
            generateParticle(m_hf.vx, m_hf.vy, m_hf.yv, m_hf.ds0, m_hf.dsx, m_hf.dsy, m_hf.dsz, ir, ipos, m_hf.tau, m_hf.xx, m_hf.yy, m_hf.eta, reflection);
          }
        }
      }
    }

    // 2013/04/30, KM, shuffle the particle list
    if (this->cfg_shuffle_particles) {
#if __cplusplus >= 201703L
      std::shuffle(this->plist.begin(), this->plist.end(), RandomURGB());
#else
      std::random_shuffle(this->plist.begin(), this->plist.end());
#endif
    }
  }

  void ParticleSampleHydrojet::finalize() {
    m_hf.closeFDataFile();
    m_hf.closePDataFile();
    if (cache_available) {
      cache_ifsPos.clear();
      cache_ifsNeg.clear();
    }
  }

  void ParticleSampleHydrojet::generateParticle(double vx, double vy, double yv,
  	  double ds0, double dsx, double dsy, double dsz, int ir, int ipos,
  	  double tau, double x0, double y0, double eta0, int reflection)
  {
    // Note: surface element
    //   dsx = tau*dy*deta*dtau
    //   dsy = tau*dx*deta*dtau
    //   dss = dtau*dx*dy
    if (reflection & 1) {
      vy = -vy;
      dsy = -dsy;
      y0 = -y0;
    }
    if (reflection & 2) {
      vx = -vx;
      yv = -yv;
      dsx = -dsx;
      dsz = -dsz;
      x0 = -x0;
      eta0 = -eta0;
    }

    double p[58], pw[58];

    double const ptmid = 1e3 / hbarc_MeVfm;
    double const dx = getDx();
    double const dy = getDy();
    double const dh = getDh();
    double const dtau = getDtau();
    double const vz = tanh(yv);
    double const gamma =  cosh(yv) / sqrt(1.0 - (vx * vx + vy * vy) * cosh(yv) * cosh(yv));
    double const beta = 1.0 / m_hf.tf;
    double const mres = rlist[ir].mass;
    double const mres2 = mres*mres;
    double prds, prx, pry, prz, er, pu;

    double const mu = rlist[ir].mu;
    double const sgn = rlist[ir].bf;
    auto integrand = [mres2, beta, mu, sgn] (double const p) {
      double const energy = std::sqrt(p * p + mres2);
      double const x = (energy - mu) * beta;
      if (x >= 30.0) return 0.0;
      return p * p / (std::exp(x) + sgn);
    };
    double const fm
      = kashiwa::IntegrateByGaussLegendre<38>(0.0, ptmid, integrand)
      + kashiwa::IntegrateByGaussLaguerre<20>(ptmid, 1.0, integrand);

    double ranmax = dx * dy * dh * tau * HYDRO2JAM_FACRANMAX;
    if (m_hf.bulk == 0) {
      if (dsx != 0.0 || dsy != 0.0) {
        ranmax = dx * dh * tau * dtau * HYDRO2JAM_FACRANMAX;
      } else {
        ranmax = dtau * dx * dy * HYDRO2JAM_FACRANMAX;
      }
    }

    double ranemis;
    do {
      do {
        // Generate momentum [0:6GeV/c] according to Bose/Fermi
        // distribution in local rest frame using bisection method
        double r1 = idt::util::urand() * fm;
        double pmax = 6000.0 / hbarc_MeVfm;
        double pmin = 0.0;
        double ppp = (pmax + pmin) * 0.5;

        double const mu = rlist[ir].mu;
        double const sgn = rlist[ir].bf;
        for (int id = 0; id < HYDRO2JAM_ITERATION_MAX; id++) {
          ppp = (pmax + pmin) * 0.5;
          double const fp = kashiwa::IntegrateByGaussLegendre<12>(0.0, ppp,
            [beta, mu, mres2, sgn] (double p) {
              double energy = std::sqrt(p * p + mres2);
              double aaa = (energy - mu) * beta;
              return aaa < 100.0 ? p * p / (std::exp(aaa) + sgn) : 0.0;
            });
          if (fp > r1)
            pmax = ppp;
          else
            pmin = ppp;
        }

        // random variable on unit sphere
        double r2 = -2 * idt::util::urand() + 1;
        double theta = std::acos(r2);
        double phi = 2 * M_PI * idt::util::urand();

        // uniform random number on surface
        double prxd = ppp * std::sin(theta) * std::cos(phi);
        double pryd = ppp * std::sin(theta) * std::sin(phi);
        double przd = ppp * std::cos(theta);
        double erd = std::sqrt(ppp * ppp + mres2);

        //Lorentz transformation by flow velocity
        double ddd = gamma * (erd + (prxd * vx + pryd * vy + przd * vz) * gamma/(1. + gamma));

        //Momentum in lab. frame
        prx = prxd + vx * ddd;
        pry = pryd + vy * ddd;
        prz = przd + vz * ddd;

        double prt = std::sqrt(prx * prx + pry * pry);
        er = std::sqrt(prt * prt + prz * prz + mres2);

        //	double mrt = sqrt(prt * prt + mres2);
        //	double yr = log((er + prz) / (er - prz)) * 0.5;

        pu = gamma * (er - vx * prx - vy * pry - vz * prz); //pu = erd
        prds = er * ds0 + prx * dsx + pry * dsy + prz * dsz;
        if (!ipos) prds = -(er * ds0 + prx * dsx + pry * dsy + prz * dsz);

      } while (prds < 0.0);

      if (prds / pu / gamma > ranmax) {
        std::cout
          << "Warning: prds/pu/gamma is greater than maximum random number. "
          << "Please increase 'HYDRO2JAM_FACRANMAX'"
          << " at least "
          << prds / pu / gamma / ranmax * HYDRO2JAM_FACRANMAX
          << " [at ParticleSampleHydrojet::generateParticle]"
          << std::endl;
      }

      ranemis = ranmax * idt::util::urand();

    } while (ranemis > prds / pu / gamma);

    //↓ は bulk emission の時だけしか正しく無い気がする by KM
    //Uniformly distributed in a fluid element in coordinate space
    double ran1 = idt::util::urand();
    double ran2 = idt::util::urand();
    double ran3 = idt::util::urand();
    double xx = x0 + dx * (ran1 - 0.5);
    double yy = y0 + dy * (ran2 - 0.5);
    double eta= eta0 + dh * (ran3 - 0.5);

    // 2013/04/23, KM, reverse z axis
    if (this->cfg_reverse_particles) {
      prz = -prz;
      eta = -eta;
    }

    // Note: ipos=1 の時が positive contribution による粒子。ipos=0 の時
    //   は negative contribution による粒子。元々の hydrojet では
    //   ipos=0 の粒子もファイルに出力する機能があった。
    if (ipos) {
      int const pdg = rlist.generatePDGCode(ir);
      double const px_GeV = prx * hbarc_GeVfm;
      double const py_GeV = pry * hbarc_GeVfm;
      double const pz_GeV = prz * hbarc_GeVfm;

      // use jamMass (JAM初期化時に jam->jamMass() で mass と e を自動決定)
      double const mass_GeV = -1.0; // mres * hbarc_GeVfm;

      base::addParticleTauEta(pdg, px_GeV, py_GeV, pz_GeV, mass_GeV, xx, yy, eta, tau);
    }
  }

  class ParticleSampleFactory: ParticleSampleFactoryRegistered {
    virtual IParticleSample* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type != "hydrojet.original") return 0;

      std::string const cachedir = ctx.indir();
      int const kintmp = ctx.kintmp();
      int const eospce = ctx.eospce();
      std::string const resodata = ctx.resodata();
      std::string const fn_freezeout_dat = inputfile + "/freezeout.dat";
      std::string const fn_position_dat = inputfile + "/position.dat";

      ParticleSampleHydrojet* psamp
        = new ParticleSampleHydrojet(ctx, cachedir, ".PC170", kintmp, eospce, resodata);
      psamp->setDtau(ctx.get_config("hydrojet_deltat", 0.3));
      psamp->setDx(ctx.get_config("hydrojet_deltax", 0.3));
      psamp->setDy(ctx.get_config("hydrojet_deltay", 0.3));
      psamp->setDh(ctx.get_config("hydrojet_deltah", 0.3));

      psamp->setBaryonFree(ctx.get_config("hydrojet_baryonfree", 1));
      psamp->setHypersurfaceFilenames(fn_freezeout_dat, fn_position_dat);
      return psamp;
    }
  } instance;

}
