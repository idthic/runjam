#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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

  class ParticleSampleHydrojet: public HydroSpectrum, public IParticleSample {
  private:
    ResonanceListPCE rlist;
    std::vector<Particle*> plist;

    bool flag_negative_contribution = false;
    int    baryonfree;
    double tmpf;
    double mubf;
    double meanf;
    bool cfg_reverse_particles;
    bool cfg_shuffle_particles;

    bool mode_delayed_cooperfrye;
    std::vector<std::ifstream> resDataPos;
    std::vector<std::ifstream> resDataNeg;

    std::vector<std::string>   elemFile;

  public:
    ParticleSampleHydrojet(runjam_context const& ctx, std::string const& cachedir, std::string* outf, int kin, int eos_pce, std::string const& fname);
    ~ParticleSampleHydrojet();

    void setBaryonFree(int i) { baryonfree = i; }
    void setTMPF(double t) { tmpf = t / hbarc_GeVfm; }
    void setMUBF(double m) { mubf = m / hbarc_GeVfm; }

  private:
    bool openCooperFryeCacheForRead();
    void initialize();
    void analyze(double oversamplingFactor);
    void finalize();

  private:
    void createCooperFryeCache(int ireso);
    void createCooperFryeCache();

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

    virtual std::vector<Particle*> const& getParticleList() const override { return plist; }

  private:
    std::string fn_freezeout_dat;
    std::string fn_position_dat;
  public:
    void setHypersurfaceFilenames(std::string const& fn_freezeout_dat, std::string const& fn_position_dat) {
      this->fn_freezeout_dat = fn_freezeout_dat;
      this->fn_position_dat = fn_position_dat;
    }
    void updateWithOverSampling(double oversamplingFactor) {
      this->initialize();
      this->analyze(oversamplingFactor);
      this->finalize();
    }
    virtual void update() override {
      this->updateWithOverSampling(1.0);
    }

  private:
    void generateParticle(
      double vx, double vy, double vz,
      double ds0, double dsx, double dsy,
      double dsz, int ir, int ipos,
      double tau, double xx, double yy, double eta, int reflection = 0);
    void putParticle(
      double px, double py, double pz,
      double e, double m, int ir, double tau, double x,
      double y, double eta, int ipos);
    void outputData(
      double prx, double pry, double prz,
      double er, double mres, int ir, double tau, double xx,
      double yy, double eta, int ipos);
  };

ParticleSampleHydrojet::ParticleSampleHydrojet(runjam_context const& ctx, std::string const& cachedir, std::string* outf, int kintmp, int eos_pce, std::string const& fname):
  HydroSpectrum(kintmp, eos_pce), rlist(kintmp, eos_pce, fname)
{
	plist.clear();

  int const nreso_loop = this->rlist.numberOfResonances();
  this->elemFile.resize(nreso_loop);
  for (int i = 0; i < nreso_loop; i++) {
    if (cachedir.size() >0)
      elemFile[i] = cachedir + "/" + outf[i];
    else
      elemFile[i] = outf[i];
  }

  mode_delayed_cooperfrye = false;

  // constants
	tmpf = 0.16 / hbarc_GeVfm;
	mubf = 1.6 / hbarc_GeVfm;
	meanf = 0.45 / hbarc_GeVfm;

  // 2013/04/23, KM, reverse z axis
  this->cfg_reverse_particles = ctx.get_config("hydrojet_reverse_particles", false);
  if (this->cfg_reverse_particles)
    std::cout << "ParticleSampleHydrojet: ReverseParticleList mode enabled!" << std::endl;

  // 2013/04/30, KM, shuffle the particle list
  this->cfg_shuffle_particles = ctx.get_config("hydrojet_shuffle_particles", false);
  if (this->cfg_shuffle_particles)
    std::cout << "ParticleSampleHydrojet: ShuffleParticleList mode enabled!" << std::endl;
}

ParticleSampleHydrojet::~ParticleSampleHydrojet() {
  if (plist.size() > 0) {
    std::vector<Particle*>::iterator cp;
    for (cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }
}

  // 一度に 151x2 のファイルを開いて書き込むとディスクに悪いので共鳴毎に処理する。
  void ParticleSampleHydrojet::createCooperFryeCache(int ireso) {
    //---------------------------------------------------------------------------
    // initialize files

    // open output files
    std::ofstream ostr_elm;
    if (ireso < 21) {
      ostr_elm.open((elemFile[ireso]).c_str());
      if (!ostr_elm) {
        std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << elemFile[ireso] << std::endl;
        std::exit(1);
      }
    }
    std::ofstream ostr_pos((elemFile[ireso] + ".POS").c_str());
    if (!ostr_pos) {
      std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << elemFile[ireso] << ".POS" << std::endl;
      std::exit(1);
    }
    std::ofstream ostr_neg((elemFile[ireso] + ".NEG").c_str());
    if (!ostr_neg) {
      std::cerr << "ParticleSampleHydrojet::createCooperFryeCache! failed to create file " << elemFile[ireso] << ".NEG" << std::endl;
      std::exit(1);
    }

    ResonanceListPCE::resonance& recreso = this->rlist[ireso];

    openFDataFile(fn_freezeout_dat);

    //---------------------------------------------------------------------------
    while (!readFData()) {
      if (!baryonfree) {
        recreso.mu = 0.0;
        if (recreso.bf == 1) recreso.mu = mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
        if (recreso.anti) recreso.mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) + meanf * nbf;
        //if (recreso.anti) recreso.mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
        //if (recreso.bf==1) recreso.mu = mub;
        //if (recreso.anti) recreso.mu = -mub;
      }

      if (tf == 0.0) {
        std::cerr << "(ParticleSampleHydrojet::createCooperFryeCache) TF=ZERO" << std::endl;
        //write(24, 9000)0.0
        continue;
      }

      //-------------------------------------------------------------------------
      // integration

      double npos = 0.0;
      double nneg = 0.0;

      double const gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
      kashiwa::phys::vector4 u(gamma, vx * gamma, vy * gamma, vz * gamma);
      kashiwa::phys::vector4 ds(ds0, -dsx, -dsy, -dsz);
      double const beta = 1.0 / tf;

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

    closeFDataFile();
  }

  void ParticleSampleHydrojet::createCooperFryeCache() {
    int const nreso_loop = this->rlist.numberOfResonances();
    for(int ireso = 0; ireso < nreso_loop; ireso++)
      this->createCooperFryeCache(ireso);
  }

bool ParticleSampleHydrojet::openCooperFryeCacheForRead() {
  int const nreso_loop = rlist.numberOfResonances();

  resDataPos.clear();
  for (int i = 0; i < nreso_loop; i++) {
    std::string fnpos = elemFile[i] + ".POS";
    resDataPos.emplace_back(fnpos.c_str());
    if (!resDataPos.back()) goto failed;
  }

  resDataNeg.clear();
  for (int i = 0; i < nreso_loop; i++) {
    std::string fnneg = elemFile[i] + ".NEG";
    resDataNeg.emplace_back(fnneg.c_str());
    if (!resDataNeg.back()) goto failed;
  }
  return true;

failed:
  resDataPos.clear();
  resDataNeg.clear();
  return false;
}

void ParticleSampleHydrojet::initialize() {
  // Open files for input.
  openFDataFile(this->fn_freezeout_dat);
  openPDataFile(this->fn_position_dat);

  std::cout << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): checking Cooper-Frye cache files (.POS/.NEG)... " << std::flush;
  if (this->openCooperFryeCacheForRead()) {
    std::cout << "yes" << std::endl;
    return;
  } else {
    std::cout << "no (incomplete).\n";
    std::cout << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): entering delayed Cooper-Frye evaluation mode." << std::endl;
    mode_delayed_cooperfrye = true;
  }
}

void ParticleSampleHydrojet::analyze(double oversamplingFactor) {
  int const nreso_loop = rlist.numberOfResonances();

  if (plist.size() > 0) {
    std::vector<Particle*>::iterator cp;
    for(cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }

  double ran;

  while (!readFData()) {
    if (!baryonfree) {
      for (int i = 0; i < nreso_loop; i++) {
        rlist[i].mu = 0.0;
        if (rlist[i].bf == 1) rlist[i].mu = mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
        if (rlist[i].anti) rlist[i].mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) + meanf * nbf;
        //if (rlist[i].anti) rlist[i].mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
        //if (rlist[i].bf == 1) rlist[i].mu = mub;
        //if (rlist[i].anti) rlist[i].mu = -mub;
      }
    }

    readPData();

    if (tf < HYDRO2JAM_TEMPERATURE_MIN) {
      if (tf == 0.0)
        std::cerr << "ParticleSampleHydrojet! (warning) zero temperature surface." << std::endl;
      continue;
    }

    // Loop over all particles.
    for (int ir = 0; ir < nreso_loop; ir++) {
      double numResPos, numResNeg;
      if (!mode_delayed_cooperfrye) {
        resDataPos[ir] >> numResPos;
        if (numResPos > 1.0)
          std::cout << "Suspicious fluid element! ireso =" << ir
                    << " " << numResPos << std::endl;

        resDataNeg[ir] >> numResNeg;
        if (numResNeg > 1.0)
          std::cout << "Suspicious fluid element! ireso =" << ir << std::endl;
      } else {
        double npos = 0.0;
        double nneg = 0.0;

        double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
        kashiwa::phys::vector4 u(gamma, vx * gamma, vy * gamma, vz * gamma);
        kashiwa::phys::vector4 ds(ds0, -dsx, -dsy, -dsz);
        double beta = 1.0 / tf;

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
      switch (iw % 4) {
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

      if (bulk == 1) {
        ds0 = dss * std::cosh(hh);
        dsz = -dss * std::sinh(hh);
      } else {
        ds0 = dss * std::sinh(hh);
        dsz = -dss * std::cosh(hh);
      }

      // positive contribution
      {
        int const ipos = 1;
        int const n = idt::util::irand_poisson(numResPos * reflection_count);
        for (int i = 0; i < n; i++) {
          int const reflection = idt::util::irand(reflection_count) * reflection_step;
          generateParticle(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta, reflection);
        }
      }

      // negative contribution
      if (flag_negative_contribution) {
        int const ipos = 0;
        int const n = idt::util::irand_poisson(numResNeg * reflection_count);
        for (int i = 0; i < n; i++) {
          int const reflection = idt::util::irand(reflection_count) * reflection_step;
          generateParticle(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta, reflection);
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
  closeFDataFile();
  closePDataFile();
  if (!mode_delayed_cooperfrye) {
    resDataPos.clear();
    resDataNeg.clear();
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
  double const beta= 1./tf;
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
  if (bulk == 0) {
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

  putParticle(prx, pry, prz, er, mres, ir, tau, xx, yy, eta, ipos);
}

void ParticleSampleHydrojet::putParticle(double px, double py, double pz,
	double e, double m, int ir, double tau, double x,
	double y, double eta, int ipos)
{
  // Note: ipos=1 の時が positive contribution による粒子。ipos=0 の時
  //   は negative contribution による粒子。元々の hydrojet では
  //   ipos=0 の粒子もファイルに出力する機能があった。
  if (!ipos) return;

  Particle* part = new Particle(rlist.generatePDGCode(ir));
  part->px = px * hbarc_GeVfm;
  part->py = py * hbarc_GeVfm;
  part->pz = pz * hbarc_GeVfm;
  //part->setPe(std::sqrt(m*m+px*px+py*py+pz*pz)*hbarc_GeVfm);
  part->e = -1.0; // onshell (JAM初期化時に jam->jamMass() で自動決定させる)
  part->x = x;
  part->y = y;
  part->t = tau * std::cosh(eta);
  part->z = tau * std::sinh(eta);
  this->plist.push_back(part);
}

static std::string elementOutputFilenames[151] = {
  "ELEMENT.A0.PC170",
  "ELEMENT.DELTA.PC170",
  "ELEMENT.DELTABAR.PC170",
  "ELEMENT.ETA.PC170",
  "ELEMENT.ETAP.PC170",
  "ELEMENT.F0.PC170",
  "ELEMENT.KBAR.PC170",//previously "ELEMENT.K0S.PC170"
  "ELEMENT.KSTAR.PC170",
  "ELEMENT.KSTARBAR.PC170",
  "ELEMENT.LAMBDA.PC170",
  "ELEMENT.LAMBDABAR.PC170",
  "ELEMENT.OMEGA.PC170",
  "ELEMENT.PHI.PC170",
  "ELEMENT.RHO.PC170",
  "ELEMENT.SIGMA.PC170",
  "ELEMENT.SIGMAB.PC170",
  "ELEMENT.SIGMABBAR.PC170",
  "ELEMENT.PI.PC170",
  "ELEMENT.K.PC170",
  "ELEMENT.PRO.PC170",
  "ELEMENT.PBAR.PC170",
  "ELEMENT.A1_1260.PC170",
  "ELEMENT.A2_1320.PC170",
  "ELEMENT.B1_1235.PC170",
  "ELEMENT.PI_1300.PC170",
  "ELEMENT.PI2_1670.PC170",
  "ELEMENT.RHO_1465.PC170",
  "ELEMENT.F2_1270.PC170",
  "ELEMENT.F2P_1525.PC170",
  "ELEMENT.H1_1170.PC170",
  "ELEMENT.F0P.PC170",
  "ELEMENT.H1P.PC170",
  "ELEMENT.ETA_1295.PC170",
  "ELEMENT.F1_1285.PC170",
  "ELEMENT.F1_1420.PC170",
  "ELEMENT.F0_1300.PC170",
  "ELEMENT.OMEGA_1420.PC170",
  "ELEMENT.F1_1510.PC170",
  "ELEMENT.OMEGA_1600.PC170",
  "ELEMENT.K1_1270.PC170",
  "ELEMENT.K1BAR_1270.PC170",
  "ELEMENT.K2S_1430.PC170",
  "ELEMENT.K2SBAR_1430.PC170",
  "ELEMENT.K0S_1430.PC170",
  "ELEMENT.K0SBAR_1430.PC170",
  "ELEMENT.K1_1400.PC170",
  "ELEMENT.K1BAR_1400.PC170",
  "ELEMENT.KSTAR_1410.PC170",
  "ELEMENT.KSTARBAR_1410.PC170",
  "ELEMENT.N_1440.PC170",
  "ELEMENT.NBAR_1440.PC170",
  "ELEMENT.N_1520.PC170",
  "ELEMENT.NBAR_1520.PC170",
  "ELEMENT.N_1535.PC170",
  "ELEMENT.NBAR_1535.PC170",
  "ELEMENT.N_1650.PC170",
  "ELEMENT.NBAR_1650.PC170",
  "ELEMENT.DELTA_1600.PC170",
  "ELEMENT.DELTABAR_1600.PC170",
  "ELEMENT.DELTA_1620.PC170",
  "ELEMENT.DELTABAR_1620.PC170",
  "ELEMENT.LAMBDA_1405.PC170",
  "ELEMENT.LAMBDABAR_1405.PC170",
  "ELEMENT.LAMBDA_1520.PC170",
  "ELEMENT.LAMBDABAR_1520.PC170",
  "ELEMENT.LAMBDA_1600.PC170",
  "ELEMENT.LAMBDABAR_1600.PC170",
  "ELEMENT.LAMBDA_1670.PC170",
  "ELEMENT.LAMBDABAR_1670.PC170",
  "ELEMENT.SIGMAB_1385.PC170",
  "ELEMENT.SIGMABBAR_1385.PC170",
  "ELEMENT.SIGMAB_1660.PC170",
  "ELEMENT.SIGMABBAR_1660.PC170",
  "ELEMENT.SIGMAB_1670.PC170",
  "ELEMENT.SIGMABBAR_1670.PC170",
  "ELEMENT.XI.PC170",
  "ELEMENT.XIBAR.PC170",
  "ELEMENT.XI_1530.PC170",
  "ELEMENT.XIBAR_1530.PC170",
  "ELEMENT.OMEGAB.PC170",
  "ELEMENT.OMEGABBAR.PC170",
  "ELEMENT.PHI_1680.PC170",
  "ELEMENT.RHO_1700.PC170",
  "ELEMENT.KSTAR_1680.PC170",
  "ELEMENT.KSTARBAR_1680.PC170",
  "ELEMENT.K_3_1780.PC170",
  "ELEMENT.K_3BAR_1780.PC170",
  "ELEMENT.K_2_1770.PC170",
  "ELEMENT.K_2BAR_1770.PC170",
  "ELEMENT.K_2_1820.PC170",
  "ELEMENT.K_2BAR_1820.PC170",
  "ELEMENT.N_1675.PC170",
  "ELEMENT.NBAR_1675.PC170",
  "ELEMENT.N_1680.PC170",
  "ELEMENT.NBAR_1680.PC170",
  "ELEMENT.N_1700.PC170",
  "ELEMENT.NBAR_1700.PC170",
  "ELEMENT.N_1710.PC170",
  "ELEMENT.NBAR_1710.PC170",
  "ELEMENT.N_1720.PC170",
  "ELEMENT.NBAR_1720.PC170",
  "ELEMENT.N_1990.PC170",
  "ELEMENT.NBAR_1990.PC170",
  "ELEMENT.DELTA_1700.PC170",
  "ELEMENT.DELTABAR_1700.PC170",
  "ELEMENT.DELTA_1900.PC170",
  "ELEMENT.DELTABAR_1900.PC170",
  "ELEMENT.DELTA_1905.PC170",
  "ELEMENT.DELTABAR_1905.PC170",
  "ELEMENT.DELTA_1910.PC170",
  "ELEMENT.DELTABAR_1910.PC170",
  "ELEMENT.DELTA_1920.PC170",
  "ELEMENT.DELTABAR_1920.PC170",
  "ELEMENT.DELTA_1930.PC170",
  "ELEMENT.DELTABAR_1930.PC170",
  "ELEMENT.DELTA_1950.PC170",
  "ELEMENT.DELTABAR_1950.PC170",
  "ELEMENT.LAMBDA_1690.PC170",
  "ELEMENT.LAMBDABAR_1690.PC170",
  "ELEMENT.LAMBDA_1800.PC170",
  "ELEMENT.LAMBDABAR_1800.PC170",
  "ELEMENT.LAMBDA_1810.PC170",
  "ELEMENT.LAMBDABAR_1810.PC170",
  "ELEMENT.LAMBDA_1820.PC170",
  "ELEMENT.LAMBDABAR_1820.PC170",
  "ELEMENT.LAMBDA_1830.PC170",
  "ELEMENT.LAMBDABAR_1830.PC170",
  "ELEMENT.LAMBDA_1890.PC170",
  "ELEMENT.LAMBDABAR_1890.PC170",
  "ELEMENT.LAMBDA_2100.PC170",
  "ELEMENT.LAMBDABAR_2100.PC170",
  "ELEMENT.LAMBDA_2110.PC170",
  "ELEMENT.LAMBDABAR_2110.PC170",
  "ELEMENT.SIGMAB_1750.PC170",
  "ELEMENT.SIGMABBAR_1750.PC170",
  "ELEMENT.SIGMAB_1775.PC170",
  "ELEMENT.SIGMABBAR_1775.PC170",
  "ELEMENT.SIGMAB_1915.PC170",
  "ELEMENT.SIGMABBAR_1915.PC170",
  "ELEMENT.SIGMAB_1940.PC170",
  "ELEMENT.SIGMABBAR_1940.PC170",
  "ELEMENT.SIGMAB_2030.PC170",
  "ELEMENT.SIGMABBAR_2030.PC170",
  "ELEMENT.XI_1690.PC170",
  "ELEMENT.XIBAR_1690.PC170",
  "ELEMENT.XI_1820.PC170",
  "ELEMENT.XIBAR_1820.PC170",
  "ELEMENT.XI_1950.PC170",
  "ELEMENT.XIBAR_1950.PC170",
  "ELEMENT.XI_2030.PC170",
  "ELEMENT.XIBAR_2030.PC170",
};

class ParticleSampleFactory: ParticleSampleFactoryRegistered {
  virtual IParticleSample* CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
    if (type != "hydrojet.original") return 0;

    std::string const indir = ctx.indir();
    int const kintmp = ctx.kintmp();
    int const eospce = ctx.eospce();
    std::string const resodata = ctx.resodata();
    std::string const fn_freezeout_dat = inputfile + "/freezeout.dat";
    std::string const fn_position_dat = inputfile + "/position.dat";

    ParticleSampleHydrojet* psamp
      = new ParticleSampleHydrojet(ctx, indir, elementOutputFilenames, kintmp, eospce, resodata);
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
