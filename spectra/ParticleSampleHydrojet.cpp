#include <cmath>
#include <cstdlib>
#include <new>
#include <algorithm>
#include <ksh/phys/Minkowski.hpp>
#include <util/Constants.hpp>
#include <util/Random.hpp>

#include "IParticleSample.hpp"
#include "IntegratedCooperFrye.hpp"
#include "ElementReso.hpp"

namespace idt {
namespace hydro2jam {
namespace {

  static const double FreezeoutSkipTemperature = 0.01; // unit: [/fm]

  class ParticleSampleHydrojet: public ElementReso, public IParticleSample {
  private:
    std::vector<std::ifstream> resDataPos;
    std::vector<std::ifstream> resDataNeg;
    int di; // iteration number in bisection method

    std::vector<Particle*> plist;
    int    baryonfree;
    double tmpf;
    double mubf;
    double meanf;

    bool mode_delayed_cooperfrye;
    bool fReverseParticleList;
    bool fShuffleParticleList;

  public:
    ParticleSampleHydrojet(std::string const& dir, std::string* outf, int kin, int eos_pce, std::string const& fname);
    ~ParticleSampleHydrojet();
    void setBaryonFree(int i) { baryonfree = i; }
    void setTMPF(double t) { tmpf = t / hbarc_MeVfm * 1000.0; }
    void setMUBF(double m) { mubf = m / hbarc_MeVfm * 1000.0; }

  private:
    bool tryOpenCooperFryeCache();
  public:
    void initialize(std::string const& fn, std::string const& fn_p);
    void analyze(std::string fn, std::string fn_p);
    void finish();

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

    std::vector<Particle*> const& getParticleList() const { return plist; }

  private:
    std::string fn_freezeout_dat;
    std::string fn_position_dat;
  public:
    void setHypersurfaceFilenames(std::string const& fn_freezeout_dat, std::string const& fn_position_dat) {
      this->fn_freezeout_dat = fn_freezeout_dat;
      this->fn_position_dat = fn_position_dat;
    }
    void update() {
      this->analyze(this->fn_freezeout_dat, this->fn_position_dat);
    }

  private:
    void getSample(
      double vx, double vy, double vz,
      double ds0, double dsx, double dsy,
      double dsz, int ir, int ipos,
      double tau, double xx, double yy, double eta);
    void putParticle(
      double px, double py, double pz,
      double e, double m, int ir, double tau, double x,
      double y, double eta, int ipos);
    void outputData(
      double prx, double pry, double prz,
      double er, double mres, int ir, double tau, double xx,
      double yy, double eta, int ipos);
  };

ParticleSampleHydrojet::ParticleSampleHydrojet(std::string const& dir, std::string* outf, int kin, int eos_pce, std::string const& fname):
  ElementReso(dir, outf, kin, eos_pce, fname)
{
	//	di = 10;
	di = 20; // iteration number in bisection method
	plist.clear();

  mode_delayed_cooperfrye = false;

  // constants
	tmpf = 0.16 / hbarc_MeVfm * 1000.0;
	mubf = 1.6 / hbarc_MeVfm * 1000.0;
	meanf = 0.45 / hbarc_MeVfm * 1000.0;

  // 2013/04/23, KM, reverse z axis
  {
    const char* env = std::getenv("ParticleSample_ReverseParticleList");
    this->fReverseParticleList = env && std::atoi(env);
    if (this->fReverseParticleList)
      std::cout << "ParticleSampleHydrojet: ReverseParticleList mode enabled!" << std::endl;
  }

  // 2013/04/30, KM, shuffle the particle list
  {
    const char* env = std::getenv("ParticleSample_ShuffleParticleList");
    this->fShuffleParticleList = env && std::atoi(env);
    if (this->fShuffleParticleList)
      std::cout << "ParticleSampleHydrojet: ShuffleParticleList mode enabled!" << std::endl;
  }
}

ParticleSampleHydrojet::~ParticleSampleHydrojet() {
  if (plist.size() > 0) {
    std::vector<Particle*>::iterator cp;
    for (cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }
}

bool ParticleSampleHydrojet::tryOpenCooperFryeCache() {
  int const nreso_loop = rlist.numberOfResonances();
  // if(baryonfree)nreso_loop = 20;

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
    resDataNeg[i].open(fnneg.c_str(), std::ios::in);
    if (!resDataNeg.back()) goto failed;
  }
  return true;

failed:
  resDataPos.clear();
  return false;
}

void ParticleSampleHydrojet::initialize(std::string const& fn_freezeout_dat, std::string const& fn_position_dat) {
  // Open files for input.
  openFDataFile(fn_freezeout_dat);
  openPDataFile(fn_position_dat);

  std::cout << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): checking Cooper-Frye cache files (.POS/.NEG)... " << std::flush;
  if (this->tryOpenCooperFryeCache()) {
    std::cout << "yes" << std::endl;
    return;
  }

  std::cout
    << "no(incomplete).\n"
    << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::initialize): entering delayed Cooper-Frye evaluation mode." << std::endl;
  mode_delayed_cooperfrye = true;
}

//void ParticleSampleHydrojet::analyze(std::string fn_freezeout_dat, std::string fn_position_dat, std::string fn_ecc)
void ParticleSampleHydrojet::analyze(std::string fn_freezeout_dat, std::string fn_position_dat) {
  int const nreso_loop = rlist.numberOfResonances();

  // Open files.
  initialize(fn_freezeout_dat, fn_position_dat);

  if (plist.size() > 0) {
    std::vector<Particle*>::iterator cp;
    for(cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }

  double ran;
  int nsamp = 1;
  int ipos = 1;
  double numResPos;
  double numResNeg;

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

    // 2014-07-30
    if (tf < FreezeoutSkipTemperature) {
      if (tf == 0.0)
        std::cerr << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::analyze)! TF=ZERO" << std::endl;
      //std::cerr << "ParticleSampleHydrojet.cpp(ParticleSampleHydrojet::analyze): skipped lowT surface." << std::endl;
      continue;
    }

    // Loop over all particles.
    for (int ir = 0; ir < nreso_loop; ir++) {
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
          idt::hydro2jam::IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, rlist[ir].mass, rlist[ir].mu);
        } else {
          idt::hydro2jam::IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, rlist[ir].mass, rlist[ir].mu);
        }

        double n = (nneg + npos) * rlist[ir].degeff;
        numResPos = npos * rlist[ir].deg;
        numResNeg = nneg * rlist[ir].deg;
      }

      if (bulk == 1) {
        ds0 = dss * cosh(hh);
        dsz = -dss * sinh(hh);
      } else {
        ds0 = dss * sinh(hh);
        dsz = -dss * cosh(hh);
      }

      for (int isamp = 0; isamp < nsamp; isamp++) {
        if ((iw == 1)  ||  (iw == 5)) {
          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, -vy, yv, ds0, dsx, -dsy, dsz, ir, ipos, tau, xx, -yy, eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx, vy, -yv, ds0, -dsx, dsy, -dsz, ir, ipos, tau, -xx, yy, -eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx, -vy, -yv, ds0, -dsx, -dsy, -dsz, ir, ipos, tau, -xx, -yy, -eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, -vy, yv, ds0, dsx, -dsy, dsz, ir, ipos, tau, xx, -yy, eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx, vy, -yv, ds0, -dsx, dsy, -dsz, ir, ipos, tau, -xx, yy, -eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx, -vy, -yv, ds0, -dsx, -dsy, -dsz, ir, ipos, tau, -xx, -yy, -eta);

        } else if((iw == 3)  ||  (iw == 7)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx, vy, -yv, ds0, -dsx, dsy, -dsz, ir, ipos, tau, -xx, yy, -eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx, vy, -yv, ds0, -dsx, dsy, -dsz, ir, ipos, tau, -xx, yy, -eta);

        } else if((iw == 2)  ||  (iw == 6)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, -vy, yv, ds0, dsx, -dsy, dsz, ir, ipos, tau, xx, -yy, eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, -vy, yv, ds0, dsx, -dsy, dsz, ir, ipos, tau, xx, -yy, eta);

        } else if((iw == 4)  ||  (iw == 8)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx, vy, yv, ds0, dsx, dsy, dsz, ir, ipos, tau, xx, yy, eta);

        }
      }
    }
  }

  // 2013/04/30, KM, shuffle the particle list
  if (this->fShuffleParticleList) {
#if __cplusplus >= 201703L
    std::shuffle(this->plist.begin(), this->plist.end(), RandomURGB());
#else
    std::random_shuffle(this->plist.begin(), this->plist.end());
#endif
  }

  finish();
}

void ParticleSampleHydrojet::finish()
{
  closeFDataFile();
  closePDataFile();
  if (!mode_delayed_cooperfrye) {
    resDataPos.clear();
    resDataNeg.clear();
  }
}

void ParticleSampleHydrojet::
getSample(double vx, double vy, double yv,
	  double ds0, double dsx, double dsy, double dsz, int ir, int ipos,
	  double tau, double x0, double y0, double eta0)
{

  double p[58], pw[58];
  double p1[12], pw1[12];
  // double  p1[38], pw1[38];

  double ptmid=1e3 / hbarc_MeVfm;
  double dx = getDx();
  double dy = getDy();
  double dh = getDh();
  double dtau = getDtau();
  double vz = tanh(yv);
  double gamma =  cosh(yv) / sqrt(1.0 - (vx * vx + vy * vy) * cosh(yv) * cosh(yv));
  double beta= 1./tf;
  double mres = rlist[ir].mass;
  double mres2 = mres*mres;
  double prds, prx, pry, prz, er, pu;

  GauLag(0.0, ptmid, p, pw);
  double fm = 0.0;
  for (int ip = 0; ip < 58; ip++) {
    double eee = sqrt(p[ip] * p[ip] + mres2);
    double aaa = (eee - rlist[ir].mu) * beta;
    if (aaa < 30.0) {
      aaa = exp(aaa);
      fm += p[ip] * p [ip] * pw[ip] / (aaa + rlist[ir].bf);
    }
  }

  double ranemis;
  // double facranmax = 1.6;
  // double facranmax = 1.7;//06/28/2010, lower switching T, larger radial flow
  // double facranmax = 1.8;//08/27/2019, LHC, larger radial flow
  double facranmax = 2.0; // 2014-07-30 for RFH
  double ranmax = dx * dy * dh * tau * facranmax;

  // surface:
  //dsx = tau*dy*deta*dtau
  //dsy = tau*dx*deta*dtau
  //dss = dtau*dx*dy
  if (bulk == 0) {
    if (dsx != 0.0 || dsy != 0.0) {
      ranmax = dx * dh * tau * dtau * facranmax;
      // Assuming dx = dy
      //	ranmax = dx*dx*dx*2.0*tau*facranmax;
      //                ^^^Coming from MABIKI/3
      //                      in FreezeOutHyperSurface.cxx

    } else {
      ranmax = dtau * dx * dy * facranmax;
      //	ranmax = dx*dx*dx*2.0*facranmax;
      //                ^^^Coming from MABIKI/3
      //                      in FreezeOutHyperSurface.cxx
    }
  }
  do {
    do {
      //Generate momentum [0:6GeV/c]
      //according to Bose/Fermi distribution
      //in local rest frame
      //using bisection method
      double r1 = Random::getRand() * fm;
      double pmax = 6000.0 / hbarc_MeVfm;
      double pmin = 0.0;
      double ppp = (pmax + pmin) * 0.5;

      for (int id = 0; id < di; id++) {
        ppp = (pmax + pmin) * 0.5;
        Gauss12(0.0, ppp, p1, pw1);
        //	    Gauss38(0.0, ppp, p1, pw1);
        double fp = 0.0;
        for (int ip = 0; ip < 12;ip++) {
          // for(int ip=0;ip<38;ip++) {
          double eee = sqrt(p1[ip] * p1[ip] + mres2);
          double aaa = (eee - rlist[ir].mu) * beta;
          if (aaa < 100.0) {
            aaa = exp(aaa);
            fp += p1[ip] * p1[ip] * pw1[ip]/(aaa + rlist[ir].bf);
          }
        }
        double f = fp - r1;
        if (f > 0.0) pmax = ppp;
        else pmin = ppp;
      }
      // random variable on unit sphere
      double r2 = -2 * Random::getRand() + 1;
      double theta = acos(r2);
      double phd = 2 * M_PI * Random::getRand();

      // uniform random number on surface
      double prxd = ppp * sin(theta) * cos(phd);
      double pryd = ppp * sin(theta) * sin(phd);
      double przd = ppp * cos(theta);
      double erd = sqrt(ppp * ppp + mres2);

      //Lorentz transformation by flow velocity
      double ddd = gamma * (erd + (prxd * vx + pryd * vy + przd * vz) * gamma/(1. + gamma));

      //Momentum in lab. frame
      prx = prxd + vx * ddd;
      pry = pryd + vy * ddd;
      prz = przd + vz * ddd;

      double prt = sqrt(prx * prx + pry * pry);
      er = sqrt(prt * prt + prz * prz + mres2);

      //	double mrt = sqrt(prt * prt + mres2);
      //	double yr = log((er + prz) / (er - prz)) * 0.5;

      pu = gamma * (er - vx * prx - vy * pry - vz * prz); //pu = erd
      prds = er * ds0 + prx * dsx + pry * dsy + prz * dsz;
      if (!ipos) prds = -(er * ds0 + prx * dsx + pry * dsy + prz * dsz);

    } while (prds < 0.0);

    if (prds / pu / gamma > ranmax) {
      std::cout
        << "Warning: prds/pu/gamma is greater than maximum random number. "
        << "Please increase 'facranmax'"
        << " at least "
        << prds / pu / gamma / ranmax * facranmax
        << " [at ParticleSampleHydrojet::getSample]"
        << std::endl;
    }

    ranemis = ranmax * Random::getRand();

  } while (ranemis > prds / pu / gamma);

  //↓ は bulk emission の時だけしか正しく無い気がする by KM
  //Uniformly distributed in a fluid element in coordinate space
  double ran1 = Random::getRand();
  double ran2 = Random::getRand();
  double ran3 = Random::getRand();
  double xx = x0 + dx * (ran1 - 0.5);
  double yy = y0 + dy * (ran2 - 0.5);
  double eta= eta0 + dh * (ran3 - 0.5);

  // 2013/04/23, KM, reverse z axis
  if (this->fReverseParticleList) {
    prz = -prz;
    eta = -eta;
  }

  putParticle(prx, pry, prz, er, mres, ir, tau, xx, yy, eta, ipos);
}

void ParticleSampleHydrojet::putParticle(double px, double py, double pz,
	double e, double m, int ir, double tau, double x,
	double y, double eta, int ipos)
{
  // Note: ipos=1 の時が positive contribution による粒子。
  //   ipos=0 の時は negative contribution による粒子。
  //   元々の hydrojet では ipos=0 の粒子もファイルに出力する機能があった。
  if (!ipos) return;

  Particle* jp = new(std::nothrow) Particle(ir);
  if (!jp) {
    std::cerr << "(ParticleSampleHydrojet::putPartile:) No more memory" << std::endl;
    exit(1);
  }

  jp->px = px * hbarc_GeVfm;
  jp->py = py * hbarc_GeVfm;
  jp->pz = pz * hbarc_GeVfm;
  //jp->setPe(std::sqrt(m*m+px*px+py*py+pz*pz)*hbarc_GeVfm);
  jp->e = -1.0; // 自動で計算する様にする
  jp->x = x;
  jp->y = y;
  jp->t = tau * std::cosh(eta);
  jp->z = tau * std::sinh(eta);
  this->plist.push_back(jp);
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
  virtual IParticleSample* CreateInstance(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile) {
    if (type != "hydrojet.original") return 0;

    std::string const indir = ctx.indir();
    int const kintmp = ctx.kintmp();
    int const eospce = ctx.eospce();
    std::string const resodata = ctx.resodata();
    std::string const fn_freezeout_dat = inputfile + "/freezeout.dat";
    std::string const fn_position_dat = inputfile + "/position.dat";

    ParticleSampleHydrojet* psamp =
      new ParticleSampleHydrojet(indir, elementOutputFilenames, kintmp, eospce, resodata);
    psamp->setDtau(ctx.get_config("hydro2jam_deltat", 0.3));
    psamp->setDx(ctx.get_config("hydro2jam_deltax", 0.3));
    psamp->setDy(ctx.get_config("hydro2jam_deltay", 0.3));
    psamp->setDh(ctx.get_config("hydro2jam_deltah", 0.3));

    psamp->setBaryonFree(ctx.get_config("hydrojet_baryonfree", 1));
    psamp->setHypersurfaceFilenames(fn_freezeout_dat, fn_position_dat);
    return psamp;
  }
} instance;

}
}
}
