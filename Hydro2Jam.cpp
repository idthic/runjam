//...A main program to use the initial condition of hadronic cascade
//...from hydro
#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Hydro2Jam.hpp"
#include "args.hpp"
#include "util/Random.hpp"
#include "spectra/IParticleSample.hpp"

namespace idt {
namespace hydro2jam {

//=========Set input values and switches ========================
Hydro2Jam::Hydro2Jam(hydro2jam_context const& ctx) {
  this->initialize(ctx);
}

void Hydro2Jam::initialize(hydro2jam_context const& ctx) {
  this->nevent = ctx.nevent(1);
  ctx.read_config(this->dumpPhaseSpaceData, "hydro2jam_phasespace_enabled", 1);

  std::string const outdir = ctx.outdir();

  int const seed = ctx.seed();

  jam = new Jam1();

  //....Initialize JAM
  jam->setMSTC(1, ctx.get_config("hydro2jam_jamseed", seed)); // int seed = 1921;
  jam->setMSTC(2, this->nevent); // number of event.
  //jam->setMSTC(38,6);          // io number for jamlist.
  jam->setMSTC(8,0);             // job mode.
  jam->setMSTC(16,0);            // display on/off.
  jam->setPARC(6,5.0);           // scale of display
  jam->setMSTC(54,0);            //avoid first coll inside the same nucleus off
  //jam->setMSTC(39,0);          //no output fname(4)

  //....Switch on some analysis.
  jam->setMSTC(156,0);           // analysis of collision distribution
  jam->setMSTC(161,0);           // no analysis from jam internal subr.
  jam->setMSTC(162,0);           // Output collision histroy
  jam->setMSTC(165,0);           //
  jam->setPARC(7, 1.0);          // Output time interval (fm/c)
  jam->setMSTC(81,0);            // 1:hard scattering on/off
  jam->setMSTC(4,100);           // user defined frame.

  //jam->setMSTC(61,0);          // isotropic resonance decay option

  //....Initial setting for JAM.
  std::string frame = "user"; // comp. frame in this case, user defined
  double dt = 100.0;          // collision time(fm/c)
  int nstep = 1;              // time step (i.e. no time step)

  //...dummy in this case.
  double bmin = 0.0;         // minimum impact parameter (dummy)
  double bmax = 0.0;         // maximum impact parameter (dummy)

  jam->setPARD(16,10.0);         // user defined frame.

  std::string fn2 = jam->getFNAME(2);
  std::string fn3 = jam->getFNAME(3);
  std::string fn4 = jam->getFNAME(4);
  if (outdir.size() > 0) {
    std::string const dir = outdir + "/";
    fn2 = dir + fn2;
    fn3 = dir + fn3;
    fn4 = dir + fn4;

    jam->setFNAME(2, fn2.c_str());
    jam->setFNAME(3, fn3.c_str());
    jam->setFNAME(4, fn4.c_str());
    jam->setFNAME(8, dir.c_str());
  }

  //...Initialize jam.
  jam->jamInit(this->nevent, bmin, bmax, dt, nstep, "user", "p ", "p ", "2gev");

  numberTestParticle = 1;

  if (dumpPhaseSpaceData) {
    std::string fname_phasespace;
    std::string fname_phasespace0;
    ctx.read_config<std::string>(fname_phasespace, "hydro2jam_phasespace_fname", "phasespace.dat");
    ctx.read_config<std::string>(fname_phasespace0, "hydro2jam_phasespace_fname0", "phasespace0.dat");
    if (outdir.size() > 0) {
      if (fname_phasespace[0] != '/')
        fname_phasespace = outdir + "/" + fname_phasespace;
      if (fname_phasespace0[0] != '/')
        fname_phasespace0 = outdir + "/" + fname_phasespace0;
    }

    ofs.open(fname_phasespace.c_str());
    if (!ofs) {
      std::cerr << "hydro2jam: failed to open '" << fname_phasespace << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs0.open(fname_phasespace0.c_str());
    if (!ofs0) {
      std::cerr << "hydro2jam: failed to open '" << fname_phasespace0 << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (ctx.get_config("hydro2jam_output_phbin", false)) {
    std::string filename = outdir + "/phasespace.bin";
    ofs_bin.open(filename, std::ios::binary);
    if (!ofs_bin) {
      std::cerr << "hydro2jam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (ctx.get_config("hydro2jam_output_phbin0", false)) {
    std::string filename = outdir + "/phasespace0.bin";
    ofs_bin0.open(filename, std::ios::binary);
    if (!ofs_bin0) {
      std::cerr << "hydro2jam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  //jam->setMDCY(jam->jamComp(111) ,1,0);   // no pi0 decay
  //jam->setMDCY(jam->jamComp(3122),1,1);   // Lambda decay
  //jam->setMDCY(jam->jamComp(3222),1,1);   // Sigma- decay
  //jam->setMDCY(jam->jamComp(3212),1,1);   // Sigma0 decay
  //jam->setMDCY(jam->jamComp(3112),1,1);   // Sigma+ decay

  if (!ctx.get_config("hydro2jam_phi_decays", true))
    jam->setMDCY(jam->jamComp(333), 1, 0); // no phi decay
}

Hydro2Jam::~Hydro2Jam() {
  if (dumpPhaseSpaceData) {
    ofs << -999 << std::endl;
    ofs.close();
    ofs0 << -999 << std::endl;
    ofs0.close();
  }

  delete jam;
}

void outputPhaseSpaceBinary(Jam1* jam, std::ofstream& ofs) {
  std::uint32_t const nv = jam->getNV();
  ofs.write("EvPh", 4);
  ofs.write((char*) &nv, 4);
  for (std::uint32_t i = 1; i <= nv; i++) {
    std::int32_t const ks = jam->getK(1,i);
    std::int32_t const kf = jam->getK(2,i);
    float const px = jam->getP(1,i);
    float const py = jam->getP(2,i);
    float const pz = jam->getP(3,i);
    float const m = jam->getP(5,i);
    float const x = jam->getR(1,i);
    float const y = jam->getR(2,i);
    float const z = jam->getR(3,i);
    float const t = jam->getR(4,i);
    ofs.write((char*) &ks, 4);
    ofs.write((char*) &kf, 4);
    ofs.write((char*) &px, 4);
    ofs.write((char*) &py, 4);
    ofs.write((char*) &pz, 4);
    ofs.write((char*) &m, 4);
    ofs.write((char*) &x, 4);
    ofs.write((char*) &y, 4);
    ofs.write((char*) &z, 4);
    ofs.write((char*) &t, 4);
  }
}

void Hydro2Jam::generateEvent(IParticleSample* psamp, std::string const& cascadeMode) {
  int const nevent = this->nevent;
  aveNumberPart1 = 0.0;
  aveNumberPart2 = 0.0;
  jam->setMSTC(5,numberTestParticle);

  if (nevent * numberTestParticle > 0)
    psamp->setAdviceNumberOfExpectedEvents(nevent * numberTestParticle);

  bool const flagSampleOnly = cascadeMode == "sample";
  bool const flagDecayOnly = cascadeMode == "decay";

  int nprint = 1;

  //...Simulation start.
  for (int iev = 1; iev <= nevent; iev++) {
    //...Set initial particle momentum and coordinate in JAM.
    nv=0;
    nbary=0;
    nmeson=0;

    // Sample particles from hydro output.
    for (int i = 0; i < numberTestParticle; i++) {
      psamp->update();
      initJam(psamp);
    }

    if (iev % nprint == 0) {
      std::cout << "hydro2jam:iev=" << iev << ": "
                << "sampling done. The number of initial test particles is nv=" << nv << "." << std::endl;
    }

    aveNumberPart1 += (double)nv/(nevent*numberTestParticle);

    //...C.M.correction.
    //cmCorrection(nv);

    if (ofs0.is_open())
      printPhaseSpaceData(ofs0); // output the distribution to phasespace0.dat
    if (ofs_bin0.is_open())
      outputPhaseSpaceBinary(jam, ofs_bin0);

    if (flagSampleOnly) continue;

    if (flagDecayOnly) {
      // perform resonance decay (no rescatterings will be performed.)
      jam->finalResonanceDecay();
      if (iev % nprint == 0) {
        std::cout
          << "hydro2jam:iev=" << iev << ": "
          << "decay done." << std::endl;
      }
    } else {
      jam->jamEvt(iev);
      if (iev % nprint == 0) {
        std::cout
          << "hydro2jam:iev=" << iev << ": "
          << "cascade done. The average number of collisions is "
          << (jam->getMSTD(41) + jam->getMSTD(42)) / numberTestParticle << "." << std::endl;
      }
    }

    if (ofs.is_open())
      printPhaseSpaceData(ofs); // output the distribution to phasespace0.dat
    if (ofs_bin.is_open())
      outputPhaseSpaceBinary(jam, ofs_bin);
  }

  jam->jamFin();
}

void Hydro2Jam::printPhaseSpaceData(std::ofstream& output) {
  int nv = jam->getNV();
  output << nv << "  " << numberTestParticle << std::endl;
  for (int i = 1; i <= nv; i++) {
    //if(jam->getK(1,i) > 10) continue;
    int ks=jam->getK(1,i);
    int kf=jam->getK(2,i);
    double px=jam->getP(1,i);
    double py=jam->getP(2,i);
    double pz=jam->getP(3,i);
    //double pe=jam->getP(4,i);
    double m = jam->getP(5,i);
    double x=jam->getR(1,i);
    double y=jam->getR(2,i);
    double z=jam->getR(3,i);
    double t=jam->getR(4,i);
    output << std::setw(5) << ks
           << std::setw(10) << kf
           << std::setw(14) << px
           << std::setw(14) << py
           << std::setw(14) << pz
           << std::setw(14) << m
           << std::setw(14) << x
           << std::setw(14) << y
           << std::setw(14) << z
           << std::setw(14) << t
           << std::endl;
  }
}

void Hydro2Jam::initJam(IParticleSample* psamp) {
  std::vector<Particle*> const& plist = psamp->getParticleList();

  std::vector<Particle*>::const_iterator mp;
  for (mp = plist.begin(); mp != plist.end(); mp++) {
    Particle const* const particle = *mp;

    int kf; // PDG particle code.
    switch(psamp->getParticleIdType()){
    case ParticleIDType::HydroParticleID:
      kf = sampleJamID(particle->id + 1);
      break;
    case ParticleIDType::PDGCode:
      kf = particle->id;
      break;
    default:
      std::cerr<<"Hydro2Jam::initJam: invalid value of psamp->getParticleIdType()."<<std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (kf == 0) continue;

    int kc = jam->jamComp(kf);      // internal particle code.
    int ibary = jam->getBaryonNumber(kc, kf);  // baryon number
    if (ibary == 0)
      nmeson++;
    else
      nbary++;

    nv++;

    //...Zero the vector.
    jam->jamZero(nv);

    int ks = 2;
    if (jam->getPMAS(kc, 2) <= 1e-7 || jam->getMDCY(kc, 1) == 0
       || jam->getMDCY(kc, 2) == 0 || jam->getMDCY(kc, 3) == 0) ks = 1;
    jam->setK(1,  nv, ks);
    jam->setK(2,  nv, kf);
    jam->setK(3,  nv, 0);
    jam->setK(4,  nv, 0);
    jam->setK(5,  nv, -1);
    jam->setK(6,  nv, 0);
    jam->setK(7,  nv, 1);
    jam->setK(8,  nv, 1);
    jam->setK(9,  nv, ibary);
    jam->setK(10, nv, 0);
    jam->setK(11, nv, 0);

    double const x  = particle->x;
    double const y  = particle->y;
    double const z  = particle->z;
    double const t  = particle->t;
    double const px = particle->px;
    double const py = particle->py;
    double const pz = particle->pz;
    double pe = particle->e;
    double pm;
    if (pe < 0.0) {
      pm = jam->jamMass(kf);
      pe = std::sqrt(px * px + py * py + pz * pz + pm * pm);
    } else {
      pm = std::sqrt(pe * pe - (px * px + py * py + pz * pz));
    }

    jam->setP(1, nv, px);
    jam->setP(2, nv, py);
    jam->setP(3, nv, pz);
    jam->setP(4, nv, pe);
    jam->setP(5, nv, pm);

    jam->setR(1, nv, x);
    jam->setR(2, nv, y);
    jam->setR(3, nv, z);
    jam->setR(4, nv, t);
    jam->setR(5, nv, t); // formation time

    //...Vertex
    jam->setV(1, nv, x);
    jam->setV(2, nv, y);
    jam->setV(3, nv, z);
    jam->setV(4, nv, t);

    //.....Set resonance decay time.
    jam->setV(5, nv, 1.0e+35);
    if (jam->getK(1, nv) == 2)
      jam->setV(5, nv, t + jam->jamDecayTime(1, kf, kc, ks, pm, pe));
  }

  jam->setNV(nv);          // set total number of particles.
  jam->setNBARY(nbary);    // set total number of baryon.
  jam->setNMESON(nmeson);  // set total number of mesons.
}

void Hydro2Jam::cmCorrection() {
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double s = 0.0;
  for (int i = 1; i <= nv; i++) {
    px += jam->getP(1, i);
    py += jam->getP(2, i);
    pz += jam->getP(3, i);
    cx += jam->getR(1, i) * jam->getP(5, i);
    cy += jam->getR(2, i) * jam->getP(5, i);
    cz += jam->getR(3, i) * jam->getP(5, i);
    s  += jam->getP(5, i);
  }
  cx = -cx / s;
  cy = -cy / s;
  cz = -cz / s;
  px = -px / nv;
  py = -py / nv;
  pz = -pz / nv;

  //cout << "cx= " << cx << " cy= " << cy << " cz= " << cz << endl;
  //cout << "px= " << px << " py= " << py << " pz= " << pz << endl;
  //cin.get();

  for (int i = 1;i <= nv; i++) {
    jam->setR(1, i, jam->getR(1, i) + cx);
    jam->setR(2, i, jam->getR(2, i) + cy);
    jam->setR(3, i, jam->getR(3, i) + cz);
    jam->setV(1, i, jam->getR(1, i));
    jam->setV(2, i, jam->getR(2, i));
    jam->setV(3, i, jam->getR(3, i));
    double const p1 = jam->getP(1,i) + px;
    double const p2 = jam->getP(2,i) + py;
    double const p3 = jam->getP(3,i) + pz;
    double const m  = jam->getP(5,i);
    double const e = sqrt(m*m + p1*p1 + p2*p2 + p3*p3);
    jam->setP(1, i, p1);
    jam->setP(2, i, p2);
    jam->setP(3, i, p3);
    jam->setP(4, i, e);
  }
}

class HydroParticleCodeTable{
  static const int nresonance=151;
  std::vector<int> table[nresonance];

  void registerMapping(int hydroParticleCode, int jamCode1, int jamCode2 = 0, int jamCode3 = 0, int jamCode4 = 0){
    if (1<=hydroParticleCode&&hydroParticleCode<=nresonance){
      std::vector<int>& list=table[hydroParticleCode-1];
      list.clear();
      list.push_back(jamCode1);
      if (jamCode2)list.push_back(jamCode2);
      if (jamCode3)list.push_back(jamCode3);
      if (jamCode4)list.push_back(jamCode4);
    }else{
      std::cerr<<"Hydro2Jam.cxx(HydroParticleCodeTable::registerMapping): invalid initialization"<<std::endl;
      std::exit(1);
    }
  }

public:
  HydroParticleCodeTable(){
    this->registerMapping(  1, 10111, 10211,-10211);           // a0
    this->registerMapping(  2,  1114,  2114,  2214,  2224);    // Delta
    this->registerMapping(  3, -1114, -2114, -2214, -2224);    // Delta-bar
    this->registerMapping(  4,   221);                         // eta
    this->registerMapping(  5,   331);                         // eta'
    this->registerMapping(  6, 10221);                         // f_0(980)
    this->registerMapping(  7,  -311,  -321);                  // Kbar
    this->registerMapping(  8,   313,   323);                  // K*(892)
    this->registerMapping(  9,  -313,  -323);                  // K*(892)_bar
    this->registerMapping( 10,  3122);                         // Lambda
    this->registerMapping( 11, -3122);                         // Lambda-bar
    this->registerMapping( 12,   223);                         // omega
    this->registerMapping( 13,   333);                         // phi(1020)
    this->registerMapping( 14,   113,   213,  -213);           // rho(770)
    this->registerMapping( 15, 10220);                         // f_0(600) i.e. sigma meson
    this->registerMapping( 16,  3112,  3212,  3222);           // Sigma
    this->registerMapping( 17, -3112, -3212, -3222);           // Sigma-bar
    this->registerMapping( 18,   111,   211,  -211);           // pion
    this->registerMapping( 19,   321,   311);                  // K=(K_0,K+)
    this->registerMapping( 20,  2112,  2212);                  // nucleon
    this->registerMapping( 21, -2112, -2212);                  // anti-nucleon
    this->registerMapping( 22, 20113, 20213,-20213);           // a_1(1260)
    this->registerMapping( 23,   115,   215,  -215);           // a_2(1320)
    this->registerMapping( 24, 10113, 10213,-10213);           // b_1(1235)
    this->registerMapping( 25, 20111, 20211,-20211);           // pi(1300)
    this->registerMapping( 26, 10115, 10215,-10215);           // pi_2(1670)
    this->registerMapping( 27, 30113, 30213,-30213);           // rho(1465)
    this->registerMapping( 28,   225);                         // f_2(1270)
    this->registerMapping( 29,   335);                         // f_2'(1525)
    this->registerMapping( 30, 10223);                         // h_1(1170)
    this->registerMapping( 31, 10331);                         // f'_0
    this->registerMapping( 32, 10333);                         // h'_1
    this->registerMapping( 33, 20221);                         // eta(1295)
    this->registerMapping( 34, 20223);                         // f_1(1285)
    this->registerMapping( 35, 20333);                         // f_1(1420)
    this->registerMapping( 36, 30221);                         // f_0(1300)
    this->registerMapping( 37, 30223);                         // omega(1420)
    this->registerMapping( 38, 50223);                         // f_1(1510)
    this->registerMapping( 39, 60223);                         // omega(1600)
    this->registerMapping( 40, 10313, 10323);                  // K_1(1270)
    this->registerMapping( 41,-10313,-10323);                  // K_1(1270)_bar
    this->registerMapping( 42,   315,   325);                  // K_2*(1430)
    this->registerMapping( 43,  -315,  -325);                  // K_2*(1430)_bar
    this->registerMapping( 44, 10311, 10321);                  // K_0*(1430)
    this->registerMapping( 45,-10311,-10321);                  // K_0*(1430)_bar
    this->registerMapping( 46, 20313, 20323);                  // K_1(1400)
    this->registerMapping( 47,-20313,-20323);                  // K_1(1400)_bar
    this->registerMapping( 48, 30313, 30323);                  // K*(1410)
    this->registerMapping( 49,-30313,-30323);                  // K*(1410)_bar
    this->registerMapping( 50, 12112, 12212);                  // N(1440)
    this->registerMapping( 51,-12112,-12212);                  // N(1440) bar
    this->registerMapping( 52,  1214,  2124);                  // N(1520)
    this->registerMapping( 53, -1214, -2124);                  // N(1520) bar
    this->registerMapping( 54, 22112, 22212);                  // N(1535)
    this->registerMapping( 55,-22112,-22212);                  // N(1535) bar
    this->registerMapping( 56, 32112, 32212);                  // N(1650)
    this->registerMapping( 57,-32112,-32212);                  // N(1650) bar
    this->registerMapping( 58, 31114, 32114, 32214, 32224);    // Delta(1600)
    this->registerMapping( 59,-31114,-32114,-32214,-32224);    // Delta(1600) bar
    this->registerMapping( 60,  1112,  1212,  2122,  2222);    // Delta(1620)
    this->registerMapping( 61, -1112, -1212, -2122, -2222);    // Delta(1620) bar
    this->registerMapping( 62, 13122);                         // Lambda(1405)
    this->registerMapping( 63,-13122);                         // Lambda(1405) bar
    this->registerMapping( 64,  3124);                         // Lambda(1520)
    this->registerMapping( 65, -3124);                         // Lambda(1520) bar
    this->registerMapping( 66, 23122);                         // Lambda(1600)
    this->registerMapping( 67,-23122);                         // Lambda(1600) bar
    this->registerMapping( 68, 33122);                         // Lambda(1670)
    this->registerMapping( 69,-33122);                         // Lambda(1670) bar
    this->registerMapping( 70,  3114,  3214,  3224);           // Sigma(1385)
    this->registerMapping( 71, -3114, -3214, -3224);           // Sigma(1385) bar
    this->registerMapping( 72, 13112, 13212, 13222);           // Sigma(1660)
    this->registerMapping( 73,-13112,-13212,-13222);           // Sigma(1660) bar
    this->registerMapping( 74, 13114, 13214, 13224);           // Sigma(1670)
    this->registerMapping( 75,-13114,-13214,-13224);           // Sigma(1670) bar
    this->registerMapping( 76,  3322,  3312);                  // Xi
    this->registerMapping( 77, -3322, -3312);                  // Xi bar
    this->registerMapping( 78,  3314,  3324);                  // Xi(1530)
    this->registerMapping( 79, -3314, -3324);                  // Xi(1530) bar
    this->registerMapping( 80,  3334);                         // Omega
    this->registerMapping( 81, -3334);                         // Omega bar
    this->registerMapping( 82, 30333);                         // phi(1680)
    this->registerMapping( 83, 40113, 40213,-40213);           // rho(1700)
    this->registerMapping( 84, 40313, 40323);                  // K*(1680)
    this->registerMapping( 85,-40313,-40323);                  // K*(1680)_bar
    this->registerMapping( 86,   317,   327);                  // K_3(1780)
    this->registerMapping( 87,  -317,  -327);                  // K_3(1780)_bar
    this->registerMapping( 88, 10315, 10325);                  // K_2(1770)
    this->registerMapping( 89,-10315,-10325);                  // K_2(1770)_bar
    this->registerMapping( 90, 20315, 20325);                  // K_2(1820)
    this->registerMapping( 91,-20315,-20325);                  // K_2(1820)_bar
    this->registerMapping( 92,  2116,  2216);                  // N(1675)
    this->registerMapping( 93, -2116, -2216);                  // N(1675)_bar
    this->registerMapping( 94, 12116, 12216);                  // N(1680)
    this->registerMapping( 95,-12116,-12216);                  // N(1680)_bar
    this->registerMapping( 96, 21214, 22124);                  // N(1700)
    this->registerMapping( 97,-21214,-22124);                  // N(1700)_bar
    this->registerMapping( 98, 42112, 42212);                  // N(1710)
    this->registerMapping( 99,-42112,-42212);                  // N(1710)_bar
    this->registerMapping(100, 31214, 32124);                  // N(1720)
    this->registerMapping(101,-31214,-32124);                  // N(1720)_bar
    this->registerMapping(102, 11218, 12128);                  // N(1990)
    this->registerMapping(103,-11218,-12128);                  // N(1990)_bar
    this->registerMapping(104, 11114, 12114, 12214, 12224);    // D(1700)
    this->registerMapping(105,-11114,-12114,-12214,-12224);    // D(1700)_bar
    this->registerMapping(106, 11112, 11212, 12122, 12222);    // D(1900)
    this->registerMapping(107,-11112,-11212,-12122,-12222);    // D(1900)_bar
    this->registerMapping(108,  1116,  1216,  2126,  2226);    // D(1905)
    this->registerMapping(109, -1116, -1216, -2126, -2226);    // D(1905)_bar
    this->registerMapping(110, 21112, 21212, 22122, 22222);    // D(1910)
    this->registerMapping(111,-21112,-21212,-22122,-22222);    // D(1910)_bar
    this->registerMapping(112, 21114, 22114, 22214, 22224);    // D(1920)
    this->registerMapping(113,-21114,-22114,-22214,-22224);    // D(1920)_bar
    this->registerMapping(114, 11116, 11216, 12126, 12226);    // D(1930)
    this->registerMapping(115,-11116,-11216,-12126,-12226);    // D(1930)_bar
    this->registerMapping(116,  1118,  2118,  2218,  2228);    // D(1950)
    this->registerMapping(117, -1118, -2118, -2218, -2228);    // D(1950)_bar
    this->registerMapping(118, 13124);                         // L(1690)0
    this->registerMapping(119,-13124);                         // L(1690)0_bar
    this->registerMapping(120, 43122);                         // L(1800)0
    this->registerMapping(121,-43122);                         // L(1800)0_bar
    this->registerMapping(122, 53122);                         // L(1810)0
    this->registerMapping(123,-53122);                         // L(1810)0_bar
    this->registerMapping(124,  3126);                         // L(1820)0
    this->registerMapping(125, -3126);                         // L(1820)0_bar
    this->registerMapping(126, 13126);                         // L(1830)0
    this->registerMapping(127,-13126);                         // L(1830)0_bar
    this->registerMapping(128, 23124);                         // L(1890)0
    this->registerMapping(129,-23124);                         // L(1890)0_bar
    this->registerMapping(130,  3128);                         // L(2100)0
    this->registerMapping(131, -3128);                         // L(2100)0_bar
    this->registerMapping(132, 23126);                         // L(2110)0
    this->registerMapping(133,-23126);                         // L(2110)0_bar
    this->registerMapping(134, 23112, 23212, 23222);           // S(1750)
    this->registerMapping(135,-23112,-23212,-23222);           // S(1750)_bar
    this->registerMapping(136,  3116,  3216,  3226);           // S(1775)
    this->registerMapping(137, -3116, -3216, -3226);           // S(1775)_bar
    this->registerMapping(138, 13116, 13216, 13226);           // S(1915)
    this->registerMapping(139,-13116,-13216,-13226);           // S(1915)_bar
    this->registerMapping(140, 23114, 23214, 23224);           // S(1940)
    this->registerMapping(141,-23114,-23214,-23224);           // S(1940)_bar
    this->registerMapping(142,  3118,  3218,  3228);           // S(2030)
    this->registerMapping(143, -3118, -3218, -3228);           // S(2030)_bar
    this->registerMapping(144, 13312, 13322);                  // X(1690)
    this->registerMapping(145,-13312,-13322);                  // X(1690)_bar
    this->registerMapping(146, 13314, 13324);                  // X(1820)
    this->registerMapping(147,-13314,-13324);                  // X(1820)_bar
    this->registerMapping(148, 23312, 23322);                  // X(1950)
    this->registerMapping(149,-23312,-23322);                  // X(1950)_bar
    this->registerMapping(150,  3316,  3326);                  // X(2030)
    this->registerMapping(151, -3316, -3326);                  // X(2030)_bar
  }

  int generateJamID(int hydroParticleCode){
    if (1 <= hydroParticleCode && hydroParticleCode <= nresonance) {
      std::vector<int>& list = table[hydroParticleCode - 1];
      if (list.size() == 1)
        return list[0];
      else if (list.size() > 0) {
        int index = int(Random::getRand()*list.size());
        if (index == list.size()) index--;
        return list[index];
      }
    }

    std::cerr << "Hydro2Jam.cxx: unexpected hydroParticleCode=" << hydroParticleCode << std::endl;
    return 0;
  }
};

int Hydro2Jam::sampleJamID(int irshift){
  static HydroParticleCodeTable table;
  return table.generateJamID(irshift);
}

}
}
