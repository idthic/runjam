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
#include "RunJam.hpp"
#include "args.hpp"
#include "util/Random.hpp"
#include "spectra/IParticleSample.hpp"

namespace idt {
namespace runjam {

//=========Set input values and switches ========================
RunJam::RunJam(runjam_context const& ctx) {
  this->initialize(ctx);
}

void RunJam::initialize(runjam_context const& ctx) {
  this->nevent = ctx.nevent(1);

  std::string const outdir = ctx.outdir();

  int const seed = ctx.seed();

  jam = new Jam1();

  //....Initialize JAM
  jam->setMSTC(1, ctx.get_config("runjam_jamseed", seed)); // int seed = 1921;
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

  if (ctx.get_config("runjam_phasespace_enabled", true)) {
    std::string fname_phasespace;
    std::string fname_phasespace0;
    ctx.read_config<std::string>(fname_phasespace, "runjam_phasespace_fname", "phasespace.dat");
    ctx.read_config<std::string>(fname_phasespace0, "runjam_phasespace_fname0", "phasespace0.dat");
    if (outdir.size() > 0) {
      if (fname_phasespace[0] != '/')
        fname_phasespace = outdir + "/" + fname_phasespace;
      if (fname_phasespace0[0] != '/')
        fname_phasespace0 = outdir + "/" + fname_phasespace0;
    }

    ofs.open(fname_phasespace.c_str());
    if (!ofs) {
      std::cerr << "runjam: failed to open '" << fname_phasespace << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs0.open(fname_phasespace0.c_str());
    if (!ofs0) {
      std::cerr << "runjam: failed to open '" << fname_phasespace0 << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (ctx.get_config("runjam_output_phbin", false)) {
    std::string filename = outdir + "/phasespace.bin";
    ofs_bin.open(filename, std::ios::binary);
    if (!ofs_bin) {
      std::cerr << "runjam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (ctx.get_config("runjam_output_phbin0", false)) {
    std::string filename = outdir + "/phasespace0.bin";
    ofs_bin0.open(filename, std::ios::binary);
    if (!ofs_bin0) {
      std::cerr << "runjam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  //jam->setMDCY(jam->jamComp(111) ,1,0);   // no pi0 decay
  //jam->setMDCY(jam->jamComp(3122),1,1);   // Lambda decay
  //jam->setMDCY(jam->jamComp(3222),1,1);   // Sigma- decay
  //jam->setMDCY(jam->jamComp(3212),1,1);   // Sigma0 decay
  //jam->setMDCY(jam->jamComp(3112),1,1);   // Sigma+ decay

  if (!ctx.get_config("runjam_phi_decays", true))
    jam->setMDCY(jam->jamComp(333), 1, 0); // no phi decay
}

RunJam::~RunJam() {
  if (ofs.is_open()) {
    ofs << -999 << std::endl;
    ofs.close();
  }
  if (ofs0.is_open()) {
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

void RunJam::printPhaseSpaceData(std::ofstream& output) {
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

void RunJam::initJam(IParticleSample* psamp) {
  std::vector<Particle*> const& plist = psamp->getParticleList();

  std::vector<Particle*>::const_iterator mp;
  for (mp = plist.begin(); mp != plist.end(); mp++) {
    Particle const* const particle = *mp;

    int kf = particle->pdg;
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
      pm = Jam1::jamMass(kf);
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

void RunJam::cmCorrection() {
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

void RunJam::generateEvent(IParticleSample* psamp, std::string const& cascadeMode) {
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
      std::cout << "runjam:iev=" << iev << ": "
                << "sampling done. The number of initial test particles is nv=" << nv << "." << std::endl;
    }

    aveNumberPart1 += (double) nv / (nevent * numberTestParticle);

    //...C.M.correction.
    //cmCorrection(nv);

    if (ofs0.is_open())
      printPhaseSpaceData(ofs0); // output the distribution to phasespace0.dat
    if (ofs_bin0.is_open())
      outputPhaseSpaceBinary(jam, ofs_bin0);

    if (flagSampleOnly) continue;

    if (flagDecayOnly) {
      jam->finalResonanceDecay();
      if (iev % nprint == 0) {
        std::cout
          << "runjam:iev=" << iev << ": "
          << "decay done." << std::endl;
      }
    } else {
      jam->jamEvt(iev);
      if (iev % nprint == 0) {
        std::cout
          << "runjam:iev=" << iev << ": "
          << "cascade done. The average number of collisions is "
          << (jam->getMSTD(41) + jam->getMSTD(42)) / numberTestParticle << "." << std::endl;
      }
    }

    aveNumberPart1 += (double) jam->getNV() / (nevent * numberTestParticle);

    if (ofs.is_open())
      printPhaseSpaceData(ofs); // output the distribution to phasespace0.dat
    if (ofs_bin.is_open())
      outputPhaseSpaceBinary(jam, ofs_bin);
  }

  jam->jamFin();
}

}
}
