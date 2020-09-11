#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "config.hpp"
#include "RunJam.hpp"
#include "spectra/ParticleSampleViscous.hpp"
#include "jam/Jam1.hpp"

using namespace idt::runjam;

void doCascade(runjam_context const& ctx, std::string const& type, std::string const& inputfile, std::string const& cascadeMode) {
  std::cout << "JAM hadronic cascade start" << std::endl;

  RunJam runjam(ctx);
  runjam.setMSTC(156, 1); // analysis of collision distribution
  runjam.setMSTC(161, 0); // no analysis from jam internal subr.
  runjam.setMSTC(162, 1); // Output collision histroy
  runjam.setMSTC(165, 1); //
  //runjam.setMSTC(41,0); // 0:no resonance decay after simulation.

  int const ntest = ctx.get_config("runjam_oversampling_factor", 1);
  bool const sw_weakdecay = ctx.get_config("runjam_switch_weak_decay", false);
  runjam.setNumberOfTestParticle(ntest);
  if (sw_weakdecay) runjam.setWeakDecay();

  std::cout << "jam event generation start" << std::endl;

  ParticleSampleBase* psamp = CreateParticleSample(ctx, type, inputfile);
  if (!psamp) {
    std::cerr << "runjam: failed to initialize ParticleSample of type '" << type << "'." <<  std::endl;
    std::exit(1);
  }

  runjam.generateEvent(psamp, cascadeMode);
  delete psamp;

  std::cout
    << "runjam: the average number of hadrons are "
    << runjam.getIniAverageParticleNumber1() << " (before JAM) -> "
    << " nhadron_after_cascade=" << runjam.getIniAverageParticleNumber2() << " (after JAM)"
    << std::endl;
}

void savePhasespaceData(std::string fname, std::vector<Particle*> plist) {
  std::FILE* f = std::fopen(fname.c_str(), "w");
  if (!f) {
    std::cerr << "runjam: failed to open the file '" << fname << "'" << std::endl;
    return;
  }

  int const nhadron = plist.size();
  std::fprintf(f, "%-4d 1\n", nhadron);
  for (std::vector<Particle*>::iterator i = plist.begin(); i != plist.end(); ++i) {
    Particle* particle = *i;

    // (1) kf ... PDG particle code
    int kf = particle->pdg;
    if (kf == 0) {
      std::cerr << "runjam: invalid value of kf (PDG particle code)." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // (2) ks ... 安定粒子なら 1, 不安定粒子なら 2.
    int ks = 2;
    int const kc = Jam1::jamComp(kf); // jam internal particle code.
    if (Jam1::getPMAS(kc,2) <= 1e-7 || Jam1::getMDCY(kc,1) == 0
      || Jam1::getMDCY(kc,2) == 0 || Jam1::getMDCY(kc,3) == 0) ks = 1;

    // (3) px,py,pz,m
    double const px = particle->px;
    double const py = particle->py;
    double const pz = particle->pz;
    double pe = particle->e;
    double m;
    if (pe < 0.0) {
      // onshell jamMass を用いて決定
      m = Jam1::jamMass(kf);
      pe = std::sqrt(px*px+py*py+pz*pz+m*m);
    } else {
      m = std::sqrt(pe*pe-(px*px+py*py+pz*pz));
    }

    // (3) x,y,z,t
    double const x = particle->x;
    double const y = particle->y;
    double const z = particle->z;
    double const t = particle->t;

    std::fprintf(
      f, "%5d %9d %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g\n",
      ks, kf, px, py, pz, m, x, y, z, t);
  }
  std::fprintf(f, "-999\n");
  std::fclose(f);
}

void doGeneratePhasespace0(runjam_context const& ctx, std::string const& type, std::string const& inputfile) {
  int const nevent = ctx.nevent(1000);
  int const ibase = ctx.get_config("runjam_ievent_begin", 0);
  double const ntest = ctx.get_config("runjam_oversampling_factor", 1.0);
  std::string const outdir = ctx.outdir();

  runjam_context ctx1(ctx);
  ctx1.set_value("runjam_nevent", 1);
  ctx1.set_value("runjam_oversampling_factor", nevent * ntest);

  ParticleSampleBase* psamp = CreateParticleSample(ctx1, type, inputfile);
  psamp->update();

  // 振り分け
  std::vector<Particle*> const& plist = psamp->getParticleList();
  std::vector<std::vector<Particle*> > phases((std::size_t) nevent);
  for (std::vector<Particle*>::const_iterator i = plist.begin(); i != plist.end(); ++i)
    phases[idt::util::irand(nevent)].push_back(*i);

  // 保存
  for (int i = 0; i < nevent; i++) {
    std::vector<char> fn(outdir.size() + 50);
    std::sprintf(&fn[0], "%s/dens%06d_phasespace0.dat", outdir.c_str(), ibase + i);
    savePhasespaceData(&fn[0], phases[i]);
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
