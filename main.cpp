#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "config.hpp"
#include "args.hpp"
#include "spectra/ParticleSample.hpp"
#include "spectra/ParticleSampleViscous.hpp"
#include "libjam.hpp"

using namespace idt;
using namespace idt::runjam;

void writePhasespaceData(std::FILE* f, Particle* begin, Particle* end, int ntest = 1) {
  std::size_t const nhadron = end - begin;
  std::fprintf(f, "%-4zd %d\n", nhadron, ntest);
  for (; begin != end; begin++) {
    Particle& particle = *begin;

    // (1) kf ... PDG particle code
    int const kf = particle.pdg;
    if (kf == 0) {
      std::cerr << "runjam: invalid value of kf (PDG particle code)." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int const ks = libjam::determineStableCode(kf);

    // (3) px,py,pz,m
    double const px = particle.px;
    double const py = particle.py;
    double const pz = particle.pz;
    double const pe = particle.e;
    double const m = particle.mass;

    // (3) x,y,z,t
    double const x = particle.x;
    double const y = particle.y;
    double const z = particle.z;
    double const t = particle.t;

    std::fprintf(
      f, "%5d %9d %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g\n",
      ks, kf, px, py, pz, m, x, y, z, t);
  }
  std::fprintf(f, "-999\n");
}

void writePhasespaceData(std::ofstream& ofs) {
  int const nv = libjam::getNV();
  int const numberTestParticle = libjam::getMSTC(5);
  ofs << nv << "  " << numberTestParticle << std::endl;
  for (int i = 1; i <= nv; i++) {
    //if(libjam::getK(1,i) > 10) continue;
    int ks = libjam::getK(1, i);
    int kf = libjam::getK(2, i);
    double px = libjam::getP(1, i);
    double py = libjam::getP(2, i);
    double pz = libjam::getP(3, i);
    //double pe=libjam::getP(4, i);
    double m = libjam::getP(5, i);
    double x = libjam::getR(1, i);
    double y = libjam::getR(2, i);
    double z = libjam::getR(3, i);
    double t = libjam::getR(4, i);
    ofs << std::setw(5) << ks
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

class IObserver {
public:
  virtual void initialize() = 0;
  virtual void process_event() = 0;
  virtual void finalize() = 0;
  virtual ~IObserver() {}
};

class FileWriterBase: public IObserver {
protected:
  std::string filename;
  std::ios::openmode openmode = std::ios::out;
  std::ofstream ofs;
  FileWriterBase(std::string const& filename, std::ios::openmode openmode = std::ios::out):
    filename(filename), openmode(openmode) {}
public:
  virtual void initialize() override {
    ofs.open(filename.c_str(), openmode);
    if (!ofs) {
      std::cerr << "runjam: failed to open '" << filename << "' for write." << std::endl;
      std::exit(1);
    }
  }
  virtual void finalize() override {
    if (ofs.is_open())
      ofs.close();
  }
};

struct PhasespaceDataWriter: public FileWriterBase {
  typedef FileWriterBase base;
  PhasespaceDataWriter(std::string const& filename): base(filename) {}
  virtual void process_event() override {
    writePhasespaceData(ofs);
  }
  virtual void finalize() override {
    if (ofs.is_open())
      ofs << -999 << std::endl;
    base::finalize();
  }
};

struct PhasespaceBinaryWriter: public FileWriterBase {
  typedef FileWriterBase base;
  PhasespaceBinaryWriter(std::string const& filename):
    base(filename, std::ios::out | std::ios::binary) {}
  virtual void process_event() override {
    std::uint32_t const nv = libjam::getNV();
    ofs.write("EvPh", 4);
    ofs.write((char*) &nv, 4);
    for (std::uint32_t i = 1; i <= nv; i++) {
      std::int32_t const ks = libjam::getK(1,i);
      std::int32_t const kf = libjam::getK(2,i);
      float const px = libjam::getP(1,i);
      float const py = libjam::getP(2,i);
      float const pz = libjam::getP(3,i);
      float const m = libjam::getP(5,i);
      float const x = libjam::getR(1,i);
      float const y = libjam::getR(2,i);
      float const z = libjam::getR(3,i);
      float const t = libjam::getR(4,i);
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
};

class IndexedPhasespaceDataWriter: public IObserver {
  std::string outdir, suffix;
  std::size_t index;
  std::vector<char> filename;
public:
  IndexedPhasespaceDataWriter(std::string const& outdir, std::string const& suffix, std::size_t startIndex):
    outdir(outdir), suffix(suffix), index(startIndex) {}
  virtual void initialize() override {}
  virtual void process_event() {
    filename.resize(outdir.size() + suffix.size() + 50);
    std::sprintf(&filename[0], "%s/dens%06zd%s", outdir.c_str(), index++, suffix.c_str());
    std::ofstream ofs(&filename[0]);
    writePhasespaceData(ofs);
    ofs << -999 << std::endl;
  }
  virtual void finalize() override {}
};

void forceJamMass(ParticleSampleBase& psample) {
  for (Particle& part: psample) {
    if (part.pdg == 0) continue;
    double const px = part.px;
    double const py = part.py;
    double const pz = part.pz;
    double const e = part.e;
    if (e < 0.0) {
      double const m = libjam::jamMass(part.pdg);
      part.mass = m;
      part.e = std::sqrt(px * px + py * py + pz * pz + m * m);
    } else {
      part.mass = std::sqrt(e * e - (px * px + py * py + pz * pz));
    }
  }
}

class Program {
private:
  int cfg_nevent;
  std::string cfg_outdir;
  int cfg_seed;
  int cfg_jamseed;
  bool cfg_phi_decays;
  bool cfg_weakdecay;

  std::vector<std::unique_ptr<IObserver>> onbefore;
  std::vector<std::unique_ptr<IObserver>> onafter;

public:
  Program(runjam_context const& ctx) {
    cfg_nevent     = ctx.nevent(1);
    cfg_outdir     = ctx.outdir();
    cfg_seed       = ctx.seed();
    cfg_jamseed    = ctx.get_config("runjam_jamseed", cfg_seed);
    cfg_phi_decays = ctx.get_config("runjam_phi_decays", true);
    cfg_weakdecay  = ctx.get_config("runjam_switch_weak_decay", false);

    if (ctx.get_config("runjam_output_phdat0", true)) {
      std::string filename = cfg_outdir + "/phasespace0.dat";
      ctx.read_config(filename, "runjam_fname_phdat0");
      onbefore.emplace_back(new PhasespaceDataWriter(filename));
    }
    if (ctx.get_config("runjam_output_phdat", true)) {
      std::string filename = cfg_outdir + "/phasespace.dat";
      ctx.read_config(filename, "runjam_fname_phdat");
      onafter.emplace_back(new PhasespaceDataWriter(filename));
    }
    if (ctx.get_config("runjam_output_phbin0", false)) {
      std::string filename = cfg_outdir + "/phasespace0.bin";
      ctx.read_config(filename, "runjam_fname_phbin0");
      onbefore.emplace_back(new PhasespaceBinaryWriter(filename));
    }
    if (ctx.get_config("runjam_output_phbin", false)) {
      std::string filename = cfg_outdir + "/phasespace.bin";
      ctx.read_config(filename, "runjam_fname_phbin");
      onafter.emplace_back(new PhasespaceBinaryWriter(filename));
    }
    if (ctx.get_config("runjam_output_phdat0_indexed", false)) {
      int const startIndex = ctx.get_config("runjam_output_index_start", 0);
      std::string suffix = "_phasespace0.dat";
      onbefore.emplace_back(new IndexedPhasespaceDataWriter(cfg_outdir, suffix, startIndex));
    }
    if (ctx.get_config("runjam_output_phdat_indexed", false)) {
      int const startIndex = ctx.get_config("runjam_output_index_start", 0);
      std::string suffix = "_phasespace.dat";
      onafter.emplace_back(new IndexedPhasespaceDataWriter(cfg_outdir, suffix, startIndex));
    }
  }

  void initializeJam() {
    // Initialize JAM
    libjam::setMSTC(1, cfg_jamseed); // int seed = 1921;
    libjam::setMSTC(2, cfg_nevent);  // number of event.
    //libjam::setMSTC(38,6); // io number for jamlist.
    libjam::setMSTC(8,0);    // job mode.
    libjam::setMSTC(16,0);   // display on/off.
    libjam::setPARC(6,5.0);  // scale of display
    libjam::setMSTC(54,0);   // avoid first coll inside the same nucleus off
    //libjam::setMSTC(39,0); //no output fname(4)

    // Switch on some analysis.
    libjam::setMSTC(156,0);  // analysis of collision distribution
    libjam::setMSTC(161,0);  // no analysis from jam internal subr.
    libjam::setMSTC(162,0);  // Output collision history
    libjam::setMSTC(165,0);  //
    libjam::setPARC(7, 1.0); // Output time interval (fm/c)
    libjam::setMSTC(81,0);   // 1:hard scattering on/off
    libjam::setMSTC(4,100);  // user defined frame.
    //libjam::setMSTC(61,0); // isotropic resonance decay option

    // Additional settings
    libjam::setMSTC(156, 1); // analysis of collision distribution
    libjam::setMSTC(161, 0); // no analysis from jam internal subr.
    libjam::setMSTC(162, 1); // Output collision history
    libjam::setMSTC(165, 1); //
    //libjam::setMSTC(41,0); // 0:no resonance decay after simulation.

    // JAM*.DAT files
    if (cfg_outdir.size() > 0) {
      std::string const dir = cfg_outdir + "/";
      libjam::setFNAME(2, dir + libjam::getFNAME(2));
      libjam::setFNAME(3, dir + libjam::getFNAME(3));
      libjam::setFNAME(4, dir + libjam::getFNAME(4));
      libjam::setFNAME(8, dir);
    }

    //libjam::setMDCY(libjam::jamComp(111) ,1,0);   // no pi0 decay
    //libjam::setMDCY(libjam::jamComp(3122),1,1);   // Lambda decay
    //libjam::setMDCY(libjam::jamComp(3222),1,1);   // Sigma- decay
    //libjam::setMDCY(libjam::jamComp(3212),1,1);   // Sigma0 decay
    //libjam::setMDCY(libjam::jamComp(3112),1,1);   // Sigma+ decay
    if (!cfg_phi_decays)
      libjam::setMDCY(libjam::jamComp(333), 1, 0); // no phi decay

    if (cfg_weakdecay) libjam::setMSTC(42, 0); // allow weak decay

    // JAMINIT
    std::string frame = "user"; // comp. frame in this case, user defined
    double dt = 100.0;          // collision time(fm/c)
    int nstep = 1;              // time step (i.e. no time step)
    double bmin = 0.0;         // minimum impact parameter (dummy)
    double bmax = 0.0;         // maximum impact parameter (dummy)
    libjam::setPARD(16, 10.0); // user defined frame.
    libjam::jamInit(this->cfg_nevent, bmin, bmax, dt, nstep, frame.c_str(), "p ", "p ", "2gev");
  }
  void finalizeJam() {
    libjam::jamFin();
  }

public:
  void  generateEvent(ParticleSampleBase* psamp, std::string const& cascadeMode) {
    int const nevent = this->cfg_nevent;

    if (nevent > 0)
      psamp->setAdviceNumberOfExpectedEvents(nevent);

    bool const flagSampleOnly = cascadeMode == "sample";
    bool const flagDecayOnly = cascadeMode == "decay";

    int nprint = 1;

    for (auto const& observer: onbefore)
      observer->initialize();
    for (auto const& observer: onafter)
      observer->initialize();

    double averageParticleNumber0 = 0.0;
    double averageParticleNumber = 0.0;
    for (int iev = 1; iev <= nevent; iev++) {
      psamp->update();
      forceJamMass(*psamp);
      //psamp->adjustCenterOfMassByLorentzBoost();
      double const ntest = psamp->getOverSamplingFactor();
      libjam::setMSTC(5, ntest);
      storeParticlesInJam(psamp->begin(), psamp->end());

      if (iev % nprint == 0) {
        std::cout
          << "runjam:iev=" << iev << ": "
          << "sampling done. The initial test-particle number is nv=" << libjam::getNV() << "." << std::endl;
      }

      averageParticleNumber0 += (double) libjam::getNV() / ntest;

      for (auto const& observer: onbefore)
        observer->process_event();

      if (flagSampleOnly) continue;

      if (flagDecayOnly) {
        libjam::finalResonanceDecay();
        if (iev % nprint == 0) {
          std::cout
            << "runjam:iev=" << iev << ": "
            << "decay done." << std::endl;
        }
      } else {
        libjam::jamEvt(iev);
        if (iev % nprint == 0) {
          std::cout
            << "runjam:iev=" << iev << ": "
            << "cascade done. The average collision number is "
            << (libjam::getMSTD(41) + libjam::getMSTD(42)) / ntest << "." << std::endl;
        }
      }

      averageParticleNumber += (double) libjam::getNV() / ntest;

      for (auto const& observer: onafter)
        observer->process_event();
    }

    for (auto const& observer: onbefore)
      observer->finalize();
    for (auto const& observer: onafter)
      observer->finalize();

    averageParticleNumber0 /= nevent;
    averageParticleNumber /= nevent;
    std::cout
      << "runjam: the average number of hadrons are "
      << averageParticleNumber0 << " (before JAM) -> "
      << averageParticleNumber << " (after JAM)" << std::endl;
  }

  static void storeParticlesInJam(Particle const* begin, Particle const* end) {
    int nv = 0;
    int nbary = 0;
    int nmeson = 0;

    for (; begin != end; ++begin) {
      Particle const& particle = *begin;
      int const kf = particle.pdg;
      if (kf == 0) continue;

      int const kc = libjam::jamComp(kf);      // internal particle code.
      int const ks = libjam::determineStableCode(kf);
      int const ibary = libjam::getBaryonNumber(kc, kf);  // baryon number
      if (ibary == 0)
        nmeson++;
      else
        nbary++;
      nv++;

      //...Zero the vector.
      libjam::jamZero(nv);
      libjam::setK(1,  nv, ks);
      libjam::setK(2,  nv, kf);
      libjam::setK(3,  nv, 0);
      libjam::setK(4,  nv, 0);
      libjam::setK(5,  nv, -1);
      libjam::setK(6,  nv, 0);
      libjam::setK(7,  nv, 1);
      libjam::setK(8,  nv, 1);
      libjam::setK(9,  nv, ibary);
      libjam::setK(10, nv, 0);
      libjam::setK(11, nv, 0);

      double const x  = particle.x;
      double const y  = particle.y;
      double const z  = particle.z;
      double const t  = particle.t;
      double const px = particle.px  ;
      double const py = particle.py  ;
      double const pz = particle.pz  ;
      double const pe = particle.e   ;
      double const pm = particle.mass;

      libjam::setP(1, nv, px);
      libjam::setP(2, nv, py);
      libjam::setP(3, nv, pz);
      libjam::setP(4, nv, pe);
      libjam::setP(5, nv, pm);

      libjam::setR(1, nv, x);
      libjam::setR(2, nv, y);
      libjam::setR(3, nv, z);
      libjam::setR(4, nv, t);
      libjam::setR(5, nv, t); // formation time
      libjam::setV(1, nv, x); // vertex
      libjam::setV(2, nv, y); // vertex
      libjam::setV(3, nv, z); // vertex
      libjam::setV(4, nv, t); // vertex

      // Set resonance decay time.
      double decayTime = 1.0e+35;
      if (ks == 2) decayTime = t + libjam::jamDecayTime(1, kf, kc, ks, pm, pe);
      libjam::setV(5, nv, decayTime);
    }

    libjam::setNV(nv);          // set total number of particles
    libjam::setNBARY(nbary);    // set total number of baryons
    libjam::setNMESON(nmeson);  // set total number of mesons
  }
};

void doCascade(runjam_context const& ctx, std::string const& type, std::string const& inputfile, std::string const& cascadeMode) {
  ParticleSampleBase* psamp = CreateParticleSample(ctx, type, inputfile);
  if (!psamp) {
    std::cerr << "runjam: failed to initialize ParticleSample of type '" << type << "'." <<  std::endl;
    std::exit(1);
  }

  Program prog(ctx);
  std::cout << "runjam: JAM hadronic cascade start" << std::endl;
  prog.initializeJam();
  prog.generateEvent(psamp, cascadeMode);
  prog.finalizeJam();
  std::cout << "runjam: JAM hadronic cascade end" << std::endl;
  delete psamp;
}

void doGeneratePhasespace0(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
  ParticleSampleBase* psamp = CreateParticleSample(ctx, type, inputfile);
  if (!psamp) {
    std::cerr << "runjam: failed to initialize ParticleSample of type '" << type << "'." <<  std::endl;
    std::exit(1);
  }

  int const nevent = ctx.nevent(1000);
  std::string const outdir = ctx.outdir();
  int const ibase = ctx.get_config("runjam_output_index_start", 0);
  if (nevent > 0)
    psamp->setAdviceNumberOfExpectedEvents(nevent);

  for (int i = 0; i < nevent; i++) {
    psamp->update();
    forceJamMass(*psamp);

    std::vector<char> buff(outdir.size() + 50);
    char* pbuff = &buff[0];
    std::sprintf(pbuff, "%s/dens%06d_phasespace0.dat", outdir.c_str(), ibase + i);
    std::FILE* f = std::fopen(pbuff, "w");
    if (!f) {
      std::cerr << "runjam: failed to open the file '" << pbuff << "'" << std::endl;
      std::exit(1);
    }
    writePhasespaceData(f, psamp->begin(), psamp->end());
    std::fclose(f);
  }
  delete psamp;
}

//-----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  runjam_context ctx;

  runjam_commandline_arguments args;
  int const ext = args.read(argc, argv, ctx);
  if (ext) return ext;

  std::cout << "runjam [version " << PACKAGE_VERSION << PACKAGE_HASH << ", seed = " << ctx.seed() << "]" << std::endl;

  idt::util::set_random_seed(ctx.seed());

  if (args.subcommand == "generate-phasespace0") {
    // test 用
    doGeneratePhasespace0(ctx, args.initType, args.initPath);
  } else if (args.subcommand == "test-viscous-correction-integration") {
    return checkViscousCooperFryeInterpolated(true);
  } else if (args.subcommand == "cascade") {
    std::string mode = ctx.get_config<std::string>("runjam_cascade_mode", "cascade");
    if (mode == "cascade" && ctx.get_config("runjam_decay_only", false)) mode = "decay";
    doCascade(ctx, args.initType, args.initPath, mode);
  } else if (args.subcommand == "decay" || args.subcommand == "sample") {
    doCascade(ctx, args.initType, args.initPath, args.subcommand);
  } else {
    std::cerr << "runjam: unknown subcommand ' " << args.subcommand << "'" << std::endl;
    return 2;
  }

  return 0;
}
