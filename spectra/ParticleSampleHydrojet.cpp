#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#ifdef HYDROJET_ALPHA
# include <strstream>
#else
# include <sstream>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#include <util.hpp>
#include <args.hpp>
#include <ksh/integrator.hpp>
#include <spectra/ResonanceList.hpp>
#include <spectra/ParticleSample.hpp>

using namespace idt;
using namespace idt::runjam;

namespace {

  static const std::size_t HYDROJET_ITERATION_MAX = 20;
  static const double HYDROJET_TEMPERATURE_MIN = 0.01; // 2014-07-30 (unit: [/fm])

  // static const double HYDROJET_FACRANMAX = 1.6;
  // static const double HYDROJET_FACRANMAX = 1.7;// 2010-06-28, lower switching T, larger radial flow
  // static const double HYDROJET_FACRANMAX = 1.8;// 2019-08-27, LHC, larger radial flow
  static const double HYDROJET_FACRANMAX = 2.0; // 2014-07-30 for RFH

  //===========================================================================
  // @fn IntegrateBosonCooperFrye
  // @fn IntegrateFermionCooperFrye

  namespace CooperFryeIntegration {
    static const double sqrtTangentLowerBound = 0;
    static const double sqrtTangentUpperBound = std::sqrt(M_PI/2);

    template<int BF>
    double Integral2(double xsig, double bmu) {
      return kashiwa::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound, sqrtTangentUpperBound, [=](double t) -> double {
        double tantt = std::tan(t * t);
        double jacob = 2*t*(tantt * tantt + 1);
        double x = tantt + xsig;
        double fE2 = x * x / (std::exp(x - bmu) - BF);
        return jacob * fE2;
      });
    }
    template<int BF>
    double IntegralP(double xsig, double bmu, double bmass) {
      return kashiwa::IntegrateByGaussLegendre<100>(sqrtTangentLowerBound, sqrtTangentUpperBound, [=](double t) -> double {
        double tantt = std::tan(t * t);
        double jacob = 2 * t * (tantt * tantt + 1);
        double x = tantt + xsig;
        double fEP = x * std::sqrt(x * x - bmass * bmass) / (std::exp(x - bmu) - BF);
        return jacob * fEP;
      });
    }
    template<int BF>
    double Integral0(double xsig, double bmu) {
      return -BF * std::log(1 - BF * std::exp(bmu - xsig));
    }

    template<int BF>
    void IntegrateCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
      dNPos = 0;
      dNNeg = 0;

      double const pi_beta3 = M_PI / (beta * beta * beta) / (8.0 * M_PI * M_PI * M_PI);
      double const bmu = beta * mu;
      double const bmass = beta * mass;
      double const dS0_ = u * ds;
      double const dS0 = std::abs(dS0_);
      (dS0_ >= 0 ? dNPos : dNNeg) = 4 * pi_beta3 * dS0 * IntegralP<BF>(bmass, bmu, bmass);

      double const dsds = ds * ds;
      if (dsds >= 0) return; // if(timelike)return;
      double const dSz = std::sqrt(dS0 * dS0 - dsds);
      double const vsig = dS0/dSz;
      double const xsig = bmass / std::sqrt(1 - vsig * vsig);
      double const dNSec = pi_beta3 * (
        dSz * ((vsig * vsig + 1) * Integral2<BF>(xsig, bmu) - bmass * bmass * Integral0<BF>(xsig, bmu))
        - 2 * dS0 * IntegralP<BF>(xsig, bmu, bmass));

      dNPos += dNSec;
      dNNeg += dNSec;
    }

    void IntegrateBosonCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
      return IntegrateCooperFrye<+1>(dNPos, dNNeg, u, ds, beta, mass, mu);
    }
    void IntegrateFermionCooperFrye(double& dNPos, double& dNNeg, const vector4& u, const vector4& ds, double beta, double mass, double mu) {
      return IntegrateCooperFrye<-1>(dNPos, dNNeg, u, ds, beta, mass, mu);
    }
  }

  using CooperFryeIntegration::IntegrateBosonCooperFrye;
  using CooperFryeIntegration::IntegrateFermionCooperFrye;

  //===========================================================================
  // HypersurfaceReader

  class HypersurfaceReader {
#ifdef HYDROJET_ALPHA
    static double clamp(double value, double lower, double upper) {
      return value < lower ? lower : value > uppper ? upper : value;
    }
    // 1 / (exp(x) + sgn)
    static double thermal_distribution(double x, double sgn) {
      double const value = 1.0 / (std::exp(clamp(x, -50.0, 50.0)) + sgn);
      return value < 0.0 ? 0.0 : value;
    }
#else
    // 1 / (exp(x) + sgn)
    static double thermal_distribution(double x, double sgn) {
      double const value = x < 30.0 ? 1.0 / (std::exp(x) + sgn) : 0.0;
      return value < 0.0 ? 0.0 : value;
    }
#endif

  public:
    int degpi, degk, degp; // degree of freedom.
    double  mpi, mk, mpro; // masses of pions, kaons and protons.
    double  mupi, muk, mup; // chemical potentials.
    int co; //=1: v1  =2:v2  =3: v3

    int fin;

    // FData
    int bulk, iw;
    double dss, ds0, dsx, dsy, dsz, tf, nbf, vx, vy, vz, yv, hh, beta;

    // PData
    double tau, xx, yy, eta;

    // EData
    double hh_ecc[100], coe_ecc[100], ecc_ecc[100];
    double area_ecc[100], aveene_ecc[100], eccp_ecc[100];

    std::string fdata_fname;
    std::string pdata_fname;
    std::string eccdata_fname;
    std::ifstream fdata;
    std::ifstream pdata;
    std::ifstream eccdata;

  private:
    bool cfg_rotate_freezeout;

  private:
    static constexpr double mPion = 139.0;
    //static constexpr double mPion  = 1020.0;//For phi meson
    //static constexpr double mPion  = 3096.9;//For J/psi meson
    static constexpr double mKaon    = 493.6;
    static constexpr double mProton  = 939.0;

  public:
    HypersurfaceReader(runjam_context const& ctx);

    void openFDataFile(std::string const& fn_freezeout_dat);
    void openPDataFile(std::string const& fn_position_dat);
    void openEDataFile(std::string const& fn_ecc);
    void closeFDataFile() { fdata.clear(); fdata.close(); }
    void closePDataFile() { pdata.clear(); pdata.close(); }
    void closeEDataFile() { eccdata.clear();eccdata.close(); }
    int  readFData();
    int  readPData();
    int  readEData();

    void setCoef(int c) {co = c;}

  protected:
    double pti, cosp, sinp, gam;

    /// @fn double thermaldist(double ee, double pz, double mu, int iw, int sgn);
    /// @param[in] ee
    /// @param[in] pz
    /// @param[in] mu
    /// @param[in] iw
    /// @param[in] sgn 1: fermion  -1:boson
    /// @remarks pti, cosp, sinp, and gam need to be set before call of
    ///   this function.
    double thermaldist(double ee, double pz, double mu, int iw, int sgn);
    //double thermaldist(double mt, double yy, double mu, int iw, int sgn);
  };

  HypersurfaceReader::HypersurfaceReader(runjam_context const& ctx) {
    fin  = -1;

    int const eospce = ctx.eospce();

    //...degree of freedom  [pi+ or pi- or pi0]
    degpi = 1;
    //      degpi  = 3;//For  phi meson
    degk  = 1;
    degp  = 2;

    //...pion mass [fm^-1]
    mpi   = 0.0;
    mk    = 0.0;
    mpro  = 0.0;
    if (eospce != 10) {
      mpi   = mPion / hbarc_MeVfm;
      mk    = mKaon / hbarc_MeVfm;
      mpro  = mProton / hbarc_MeVfm;
    }

    mupi = 0.0; //!< chemical potential at f.o. (pion)   from pi to delta(1232)
    muk = 0.0;  //!< chemical potential at f.o. (kaon)   from pi to delta(1232)
    mup = 0.0;  //!< chemical potential at f.o. (proton) from pi to delta(1232)

    if (eospce == 1) {
      int kintmp = ctx.kintmp();
      switch (kintmp) {
        //   (tf=80mev)
      case 1:
        mupi = 0.951925e+02 / hbarc_MeVfm;
        muk = 0.233670e+03 / hbarc_MeVfm;
        mup = 0.456089e+03 / hbarc_MeVfm;
        std::cout << "HypersurfaceReader Tf=80MeV" << std::endl;
        break;

        //     (tf=100mev)
      case 2:
        mupi = 0.833141e+02 / hbarc_MeVfm;
        //    mupi = 0.356941e+03 / hbarc_MeVfm;//For phi-meson
        //    mupi = 0.0 / hbarc_MeVfm;//For J/psi
        muk = 0.180805e+03 / hbarc_MeVfm;
        mup = 0.348810e+03 / hbarc_MeVfm;
        std::cout << "HypersurfaceReader Tf=100MeV" << std::endl;
        break;

        //     (tf=120mev)
      case 3:
        mupi = 0.646814e+02 / hbarc_MeVfm;
        muk = 0.128598e+03 / hbarc_MeVfm;
        mup = 0.245865e+03 / hbarc_MeVfm;
        std::cout << "HypersurfaceReader Tf=120MeV" << std::endl;
        break;

        //     (tf=140mev)
      case 4:
        mupi = 0.406648e+02 / hbarc_MeVfm;
        muk = 0.633578e+02 / hbarc_MeVfm;
        mup = 0.145518e+03 / hbarc_MeVfm;
        std::cout << "HypersurfaceReader Tf=140MeV" << std::endl;
        break;

        //     (tf=160mev)
      case 5:
        mupi = 0.137867e+02 / hbarc_MeVfm;
        muk = 0.249125e+02 / hbarc_MeVfm;
        mup = 0.476744e+02 / hbarc_MeVfm;
        std::cout << "HypersurfaceReader Tf=160MeV" << std::endl;
        break;

      default:
        std::cerr << "HypersurfaceReader::  Not avaiable sorry" << std::endl;
        std::cerr << " kinetic temperature = " << kintmp << std::endl;
        std::exit(1);
      }
    }

    //      co = 1;
    //     ^^^^^^^---> V1
    co = 2;
    //     ^^^^^^^---> V2
    //      co = 3;
    //     ^^^^^^^---> V3
    //      co = 4;
    //     ^^^^^^^---> V3

    // KM, 2013/04/20, initializes the flag to rotate the data in freezeout.dat
    this->cfg_rotate_freezeout = ctx.get_config("hydrojet_rotate_freezeout", false);
    if (cfg_rotate_freezeout)
      std::cout << "HypersurfaceReader: RotateFreezeoutData mode enabled!" << std::endl;
  }

  double HypersurfaceReader::thermaldist(double ee, double pz, double mu, int iw, int sgn) {
    if (iw == 1 || iw == 5) {
      double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
      double pds4 = ee * ds0 + pti * cosp * dsx - pti * sinp * dsy + pz * dsz;
      double pds6 = ee * ds0 - pti * cosp * dsx + pti * sinp * dsy - pz * dsz;
      double pds7 = ee * ds0 - pti * cosp * dsx - pti * sinp * dsy - pz * dsz;
      double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
      double pu4 = gam * (ee - pti * cosp * vx + pti * sinp * vy - pz * vz);
      double pu6 = gam * (ee + pti * cosp * vx - pti * sinp * vy + pz * vz);
      double pu7 = gam * (ee + pti * cosp * vx + pti * sinp * vy + pz * vz);
      double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
      double bose4 = pds4 * thermal_distribution(beta * (pu4 - mu), sgn);
      double bose6 = pds6 * thermal_distribution(beta * (pu6 - mu), sgn);
      double bose7 = pds7 * thermal_distribution(beta * (pu7 - mu), sgn);
      return bose1 + bose4 + bose6 + bose7;
    } else if (iw == 3 || iw == 7) {
      double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
      double pds6 = ee * ds0 - pti * cosp * dsx + pti * sinp * dsy - pz * dsz;
      double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
      double pu6 = gam * (ee + pti * cosp * vx - pti * sinp * vy + pz * vz);
      double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
      double bose6 = pds6 * thermal_distribution(beta * (pu6 - mu), sgn);
      return bose1 + bose6;
    } else if (iw == 2 || iw == 6) {
      double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
      double pds4 = ee * ds0 + pti * cosp * dsx - pti * sinp * dsy + pz * dsz;
      double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
      double pu4 = gam * (ee - pti * cosp * vx + pti * sinp * vy - pz * vz);
      double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
      double bose4 = pds4 * thermal_distribution(beta * (pu4 - mu), sgn);
      return bose1 + bose4;
    } else if (iw == 4 || iw == 8) {
      double pds1 = ee * ds0 + pti * cosp * dsx + pti * sinp * dsy + pz * dsz;
      double pu1 = gam * (ee - pti * cosp * vx - pti * sinp * vy - pz * vz);
      double bose1 = pds1 * thermal_distribution(beta * (pu1 - mu), sgn);
      return bose1;
    } else {
      std::cerr << "HydroSpec::thermaldist funny iw " << iw << std::endl;
      std::exit(1);
    }
  }

  void HypersurfaceReader::openFDataFile(std::string const& fn_freezeout_dat) {
    if (fdata.is_open()) {
      std::cerr << "HypersurfaceReader::openFDataFile! file already opened '" << fn_freezeout_dat << "'" << std::endl;
      std::exit(1);
    }
    fdata_fname = fn_freezeout_dat;
    fdata.open(fn_freezeout_dat.c_str());
    if (!fdata) {
      std::cerr << "HypersurfaceReader::openFDataFile! unable to open file '" << fn_freezeout_dat << "'" << std::endl;
      std::exit(1);
    }
  }
  void HypersurfaceReader::openPDataFile(std::string const& fn_position_dat) {
    if (pdata.is_open()) {
      std::cerr << "HypersurfaceReader::openPDataFile! file already opened '" << fn_position_dat << "'" << std::endl;
      std::exit(1);
    }
    pdata_fname = fn_position_dat;
    pdata.open(fn_position_dat.c_str());
    if (!pdata)  {
      std::cerr << "HypersurfaceReader::openPDataFile! unable to open file '" << fn_position_dat << "'" << std::endl;
      std::exit(1);
    }
  }
  void HypersurfaceReader::openEDataFile(std::string const& fn_ecc) {
    if (eccdata.is_open()) {
      std::cerr << "HypersurfaceReader::openEDataFile! file already opened '" << fn_ecc << "'" << std::endl;
      std::exit(1);
    }
    eccdata_fname = fn_ecc;
    eccdata.open(fn_ecc.c_str());
    if (!eccdata)  {
      std::cerr << "HypersurfaceReader::openEDataFile! unable to open file '" << fn_ecc << "'" << std::endl;
      std::exit(1);
    }
  }

  // input data of hydrodynamic flow
  int HypersurfaceReader::readFData() {
    if (!fdata) {
      std::cerr << "HypersurfaceReader::readFData: unexpected end of freezeout.dat file" << std::endl;
      std::exit(1);
    }

    fdata >> bulk;
    if (bulk == fin) return 1;

    fdata >> dss; // $2 time component of the surface
    fdata >> dsx; // $3 x component of the surface
    fdata >> dsy; // $4 y component of the surface
    // fdata >> dsz;
    fdata >> hh;  // $5 eta of the cell?
    fdata >> tf;  // $6 temperature of the fluid
    fdata >> nbf; // $7 baryon number density of the fluid
    fdata >> vx;  // $8 velocity along the x axis
    fdata >> vy;  // $9 velocity along the y axis

    fdata >> yv;  // $10 rapidity
    // vz = tanh(yv);

    fdata >> iw;

    // 2013/04/20, KM, reverse data
    if (this->cfg_rotate_freezeout) {
      // hh = -hh;
      // if (bulk != 1)
      //   dss = -dss;

      if (bulk == 1)
        hh = -hh;
      yv = -yv;
      if (iw != 8) {
        std::cerr << "HypersurfaceReader_RotateFreezeoutData: not supported for the case that iw==" << iw << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // if (tf * 197.32 > 100.0) cout << tf * 197.32 << std::endl;
    beta = 1.0 / tf;

    vz = tanh(yv);

    if (bulk == 1) {
      ds0 = dss * cosh(hh);
      dsz = -dss * sinh(hh);
    } else {
      ds0 = dss * sinh(hh);
      dsz = -dss * cosh(hh);
    }

    return 0;
  }
  int HypersurfaceReader::readPData() {
    pdata >> tau;
    pdata >> xx;
    pdata >> yy;
    pdata >> eta;

    // 2013/04/20, KM, reverse data
    if (this->cfg_rotate_freezeout) {
      eta=-eta;
    }

    return 0;
  }
  int HypersurfaceReader::readEData() {
    char buffer[1280];
    std::string line;
    std::size_t index = 0, iline = 0;
    while (std::getline(eccdata, line)) {
      iline++;
      if (index >= 100) {
        std::cerr << eccdata_fname << ":" << iline << ":"
                  << " too many data" << std::endl;
        std::exit(1);
      }

      int r = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf",
        &hh_ecc[index],
        &coe_ecc[index],
        &ecc_ecc[index],
        &area_ecc[index],
        &aveene_ecc[index],
        &eccp_ecc[index]);
      if (r != 6) {
        std::cerr << eccdata_fname << ":" << iline << ":"
                  << " invalid format" << std::endl;
        std::exit(1);
      }
      index++;
    }

    return index;
  }

  //===========================================================================
  // ParticleSampleHydrojet

  class ParticleSampleHydrojet: public OversampledParticleSampleBase {
    typedef OversampledParticleSampleBase base;
  private:
    ResonanceList rlist;
    HypersurfaceReader m_hf;

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
    ParticleSampleHydrojet(runjam_context const& ctx, std::string const& cachedir, std::string suffix);

  private:
    void initialize();
    void analyze(double oversamplingFactor);
    void finalize();
    void generateParticle(
      double vx, double vy, double vz,
      double ds0, double dsx, double dsy,
      double dsz, ResonanceRecord const& reso, int ipos,
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
    runjam_context const& ctx, std::string const& cachedir, std::string suffix
  ): base(ctx), rlist(ctx), m_hf(ctx) {
    int const nreso_loop = this->rlist.size();
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

    ResonanceRecord& recreso = this->rlist[ireso];
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
      vector4 u {gamma, m_hf.vx * gamma, m_hf.vy * gamma, m_hf.vz * gamma};
      vector4 ds {m_hf.ds0, -m_hf.dsx, -m_hf.dsy, -m_hf.dsz};
      double const beta = 1.0 / m_hf.tf;

      if (recreso.bf == -1) {
        IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
      } else {
        IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
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
    int const nreso_loop = this->rlist.size();
    for(int ireso = 0; ireso < nreso_loop; ireso++)
      this->createCooperFryeCache(ireso);
  }

  bool ParticleSampleHydrojet::openCooperFryeCacheForRead() {
    int const nreso_loop = this->rlist.size();

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

    std::cout << "ParticleSampleHydrojet: checking Cooper-Frye cache files (.POS/.NEG)... " << std::flush;
    if (this->openCooperFryeCacheForRead()) {
      cache_available = true;
      std::cout << "yes" << std::endl;
      return;
    } else {
      std::cout << "no (incomplete).\n";
      std::cout << "ParticleSampleHydrojet: entering delayed Cooper-Frye evaluation mode." << std::endl;
    }
  }

  void ParticleSampleHydrojet::analyze(double oversamplingFactor) {
    int const nreso_loop = this->rlist.size();

    double ran;

    while (!m_hf.readFData()) {
      m_hf.readPData();

      if (!cfg_baryon_disable) {
        for (ResonanceRecord& reso: rlist) {
          reso.mu = 0.0;
          if (reso.bf == 1) {
            double const mub = cfg_baryon_mubf * sqrt(1.0 - m_hf.tf * m_hf.tf / cfg_baryon_tmpf / cfg_baryon_tmpf);
            double const mean_field = cfg_baryon_meanf * m_hf.nbf;
            reso.mu   = mub - mean_field;
            if (reso.anti) reso.mu = -mub + mean_field;
            // if (reso.anti) reso.mu = -mub - mean_field;
          }
        }
      }

      if (m_hf.tf < HYDROJET_TEMPERATURE_MIN) {
        if (m_hf.tf == 0.0)
          std::cerr << "ParticleSampleHydrojet! (warning) zero temperature surface." << std::endl;
        continue;
      }

      // Loop over all particles.
      for (int ireso = 0; ireso < nreso_loop; ireso++) {
        ResonanceRecord const& reso = rlist[ireso];

        double numResPos, numResNeg;
        if (cache_available) {
          cache_ifsPos[ireso] >> numResPos;
          if (numResPos > 1.0)
            std::cout << "Suspicious fluid element! ireso =" << ireso
                      << " " << numResPos << std::endl;

          cache_ifsNeg[ireso] >> numResNeg;
          if (numResNeg > 1.0)
            std::cout << "Suspicious fluid element! ireso =" << ireso << std::endl;
        } else {
          double npos = 0.0;
          double nneg = 0.0;

          double gamma = 1.0 / sqrt(1.0 - m_hf.vx * m_hf.vx - m_hf.vy * m_hf.vy - m_hf.vz * m_hf.vz);
          vector4 u {gamma, m_hf.vx * gamma, m_hf.vy * gamma, m_hf.vz * gamma};
          vector4 ds {m_hf.ds0, -m_hf.dsx, -m_hf.dsy, -m_hf.dsz};
          double beta = 1.0 / m_hf.tf;

          if (reso.bf == -1) {
            IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, reso.mass, reso.mu);
          } else {
            IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, reso.mass, reso.mu);
          }

          double n = (nneg + npos) * reso.degeff;
          numResPos = npos * reso.deg;
          numResNeg = nneg * reso.deg;
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
            generateParticle(m_hf.vx, m_hf.vy, m_hf.yv, m_hf.ds0, m_hf.dsx, m_hf.dsy, m_hf.dsz, reso, ipos, m_hf.tau, m_hf.xx, m_hf.yy, m_hf.eta, reflection);
          }
        }

        // negative contribution
        if (flag_negative_contribution) {
          int const ipos = 0;
          int const n = idt::util::irand_poisson(numResNeg * reflection_count);
          for (int i = 0; i < n; i++) {
            int const reflection = idt::util::irand(reflection_count) * reflection_step;
            generateParticle(m_hf.vx, m_hf.vy, m_hf.yv, m_hf.ds0, m_hf.dsx, m_hf.dsy, m_hf.dsz, reso, ipos, m_hf.tau, m_hf.xx, m_hf.yy, m_hf.eta, reflection);
          }
        }
      }
    }

    // 2013/04/30, KM, shuffle the particle list
    if (this->cfg_shuffle_particles) this->shuffleParticleList();
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
      double ds0, double dsx, double dsy, double dsz, ResonanceRecord const& reso, int ipos,
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
    double const vz = std::tanh(yv);
    double const gamma =  std::cosh(yv) / std::sqrt(1.0 - (vx * vx + vy * vy) * std::cosh(yv) * std::cosh(yv));
    double const beta = 1.0 / m_hf.tf;
    double const mres = reso.mass;
    double const mres2 = mres*mres;
    double prds, prx, pry, prz, er, pu;

    double const mu = reso.mu;
    double const sgn = reso.bf;
    auto integrand = [mres2, beta, mu, sgn] (double const p) {
      double const energy = std::sqrt(p * p + mres2);
      double const x = (energy - mu) * beta;
      if (x >= 30.0) return 0.0;
      return p * p / (std::exp(x) + sgn);
    };
    double const fm
      = kashiwa::IntegrateByGaussLegendre<38>(0.0, ptmid, integrand)
      + kashiwa::IntegrateByGaussLaguerre<20>(ptmid, 1.0, integrand);

    double ranmax = dx * dy * dh * tau * HYDROJET_FACRANMAX;
    if (m_hf.bulk == 0) {
      if (dsx != 0.0 || dsy != 0.0) {
        ranmax = dx * dh * tau * dtau * HYDROJET_FACRANMAX;
      } else {
        ranmax = dtau * dx * dy * HYDROJET_FACRANMAX;
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

        double const mu = reso.mu;
        double const sgn = reso.bf;
        for (int id = 0; id < HYDROJET_ITERATION_MAX; id++) {
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

        //  double mrt = sqrt(prt * prt + mres2);
        //  double yr = log((er + prz) / (er - prz)) * 0.5;

        pu = gamma * (er - vx * prx - vy * pry - vz * prz); //pu = erd
        prds = er * ds0 + prx * dsx + pry * dsy + prz * dsz;
        if (!ipos) prds = -(er * ds0 + prx * dsx + pry * dsy + prz * dsz);

      } while (prds < 0.0);

      if (prds / pu / gamma > ranmax) {
        std::cout
          << "Warning: prds/pu/gamma is greater than maximum random number. "
          << "Please increase 'HYDROJET_FACRANMAX'"
          << " at least "
          << prds / pu / gamma / ranmax * HYDROJET_FACRANMAX
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
      int const pdg = reso.generatePDGCode();
      double const px_GeV = prx * hbarc_GeVfm;
      double const py_GeV = pry * hbarc_GeVfm;
      double const pz_GeV = prz * hbarc_GeVfm;
      double const mass_GeV = mres * hbarc_GeVfm;
      base::addParticleMilne(pdg, px_GeV, py_GeV, pz_GeV, mass_GeV, xx, yy, eta, tau);
    }
  }

  class ParticleSampleFactory: ParticleSampleFactoryBase {
    virtual std::unique_ptr<ParticleSampleBase> CreateInstance(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
      if (type != "hydrojet.original") return nullptr;

      std::string const cachedir = ctx.indir();
      std::string const fn_freezeout_dat = inputfile + "/freezeout.dat";
      std::string const fn_position_dat = inputfile + "/position.dat";

      ParticleSampleHydrojet* psamp = new ParticleSampleHydrojet(ctx, cachedir, ".PC170");
      psamp->setDtau(ctx.get_config("hydrojet_deltat", 0.3));
      psamp->setDx(ctx.get_config("hydrojet_deltax", 0.3));
      psamp->setDy(ctx.get_config("hydrojet_deltay", 0.3));
      psamp->setDh(ctx.get_config("hydrojet_deltah", 0.3));

      psamp->setBaryonFree(ctx.get_config("hydrojet_baryonfree", 1));
      psamp->setHypersurfaceFilenames(fn_freezeout_dat, fn_position_dat);
      return std::unique_ptr<ParticleSampleBase>(psamp);
    }
  } instance;

}
