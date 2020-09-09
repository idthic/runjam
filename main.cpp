#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "config.hpp"
#include "util/Random.hpp"
#include "util/PyRand.hpp"
#include "RunJam.hpp"
#include "spectra/ParticleSampleViscous.hpp"
#include "jam/Jam1.hpp"

using namespace idt::runjam;

//-----------------------------------------------------------------------------
// parameters

struct RunjamCommandlineArguments {
  std::string subcommand;
  std::string initType;
  std::string initPath;

public:
  RunjamCommandlineArguments() {
    this->subcommand = "cascade";
    this->initType = "c0lrf";
    this->initPath = "hypersurface_v1.txt";

    // this->initType = "hydrojet.original";
    // this->initPath = "test";
  }

private:
  static void cmd_version() {
    std::cout << "runjam (idt) " << PACKAGE_VERSION << PACKAGE_HASH << std::endl;
  }
  static void cmd_help() {
    std::printf(
      "usage: runjam [SUBCOMMAND] [OPTIONS|VAR=VALUE]\n"
      "\n"
      "SUBCOMMAND\n"
      "  cascade (default)\n"
      "  decay\n"
      "  sample\n"
      "\n"
      "OPTIONS and VARIABLES\n"
      "      runjam_cascade_mode=SUBCOMMAND\n"
      "  -n, runjam_nevent=INT [1]              number of events to process\n"
      "  -s, runjam_seed=INT [18371]            seed for random numbers\n"
      "  -t, runjam_oversampling_factor=NUM [1] ntest\n"
      "  -w, runjam_switch_weak_decay=INT [0]   sw_weakdecay {0 | 1}\n"
      "      runjam_phi_decays=BOOL [true]\n"
      "      runjam_decay_only=BOOL [false]\n"
      "  --resodata, runjam_resodata=FILE [ResonanceJam.dat] resodata\n"
      "\n"
      " Output options\n"
      "  -o,        runjam_output_directory=DIR [out]   directory of output files\n"
      "  -d,        runjam_phasespace_enabled=INT [1]   dump phasespace\n"
      "  --fphase,  runjam_phasespace_fname=FILE [phasespace.dat]   output filename\n"
      "  --fphase0, runjam_phasespace_fname0=FILE [phasespace0.dat] output filename\n"
      "             runjam_output_phbin=BOOL [false]    output binary phasespace\n"
      "             runjam_output_phbin0=BOOL [false]   output binary phasespace0\n"
      "\n"
      " Initialiation options\n"
      "  -i ICSPEC       specify initial condition\n"
      "    c0lrf:FILE    sample particles from the hypersurface data from\n"
      "                  rfh c0lrf format \"hypersurface_v1.txt\"\n"
      "    hydrojet:DIR  sample particles using the hypersurface data from\n"
      "                  hydrojet (DIR/freezeout.dat, DIR/position.dat)\n"
      "    phase:FILE    load particle lists from the text format \"phasespace.dat\".\n"
      "    phase1:FILE   load a particle list from FILE. This performs an additional\n"
      "                  check to require that FILE contains only a single event.\n"
      "    phbin:FILE    load particle lists from the binary format \"ph000k.bin\".\n"
      "    psample:FILE  read a particle list from the file in the format of\n"
      "                  runjam \"particlesample_pos.dat\"\n"
      "\n"
      " Options for hydrojet hypersurface\n"
      "  --hydrojet-ftemp, hydrojet_kintmp=INT [5]       freezeout temperature type\n"
      "  --hydrojet-pce,   hydrojet_eospce=INT [6]       eos_pce\n"
      "  --hydrojet-bfree, hydrojet_baryonfree=INT [1]   baryonfree\n"
      "  --hydrojet-dir,   hydrojet_directory=DIR [test] directory of freezeout.dat\n"
      "  --hydrojet-dt,    hydrojet_deltat=NUM [0.3]     delta tau\n"
      "  --hydrojet-dx,    hydrojet_deltax=NUM [0.3]     delta x\n"
      "  --hydrojet-dy,    hydrojet_deltay=NUM [0.3]     delta y\n"
      "  --hydrojet-dh,    hydrojet_deltah=NUM [0.3]     delta eta\n"
      "\n"
      " Options for c0lrf sampler\n"
      "  --switching-temperature, runjam_switching_temperature=TEMP [155]\n"
      "          an advice to switching temperature in MeV\n"
      "  runjam_turnsOffViscousEffect=INT\n"
      "\n"
      " Other options\n"
      "  --help          show this help\n"
      "  --version       show version information\n"
      "\n"
      "SAMPLE\n"
      "\n"
      "$ ./runjam cascade -s 12345 -o jam --hydrojet-dir hydro\n"
      "$ ./runjam cascade -s 12345 -o jam -i phase:phasespace0.in\n"
      "$ ./runjam decay -s 12345 -o jam -i phase:phasespace0.in\n"
      "\n"
    );
  }

private:
  int argc;
  char** argv;
  runjam_context* ctx;

  int i;
  const char* arg;
  int arg_index;
  const char* arg_optarg;
  bool flag_error;
  bool flag_version;
  bool flag_help;

  const char* get_optarg() {
    if (arg_optarg) {
      const char* ret = arg_optarg;
      arg_optarg = nullptr;
      return ret;
    } else if (i < argc) {
      return argv[i++];
    } else {
      std::cerr << "runjam:$" << arg_index << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
      return NULL;
    }
  }

  void assign_optarg(const char* key) {
    if (const char* optarg = get_optarg())
      ctx->set_value(key, optarg);
  }

  void assign_optarg_double(const char* key) {
    if (const char* optarg = get_optarg()) {
      if (*optarg) {
        ctx->set_value(key, std::atof(optarg));
      } else {
        std::cerr << "runjam:option(" << arg << "): the argument of the option is empty." << std::endl;
        flag_error = true;
      }
    }
  }

  void assign_optarg_int(const char* key) {
    if (const char* optarg = get_optarg())
      ctx->set_value(key, std::atoi(optarg));
  }

  void assign_optarg_input() {
    const char* optarg = get_optarg();
    if (!optarg) return;

    std::string spec = optarg;

    std::size_t const pos = spec.find(':');
    if (pos == std::string::npos) {
      std::cerr << "runjam: unrecognized input '-i " << optarg << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->initType = spec.substr(0, pos);
    this->initPath = spec.substr(pos + 1);
  }

private:
  void read_longname_option() {
    std::string longname = arg + 2;
    std::size_t ia;
    if ((ia = longname.find('=', 0)) != std::string::npos) {
      arg_optarg = &arg[2 + ia + 1]; // = の右辺
      longname = longname.substr(0, ia);
    }

    if (longname == "help") {
      flag_help = true;
    } else if (longname == "version") {
      flag_version = true;
    } else if (longname == "resodata") {
      assign_optarg("runjam_resodata");
    } else if (longname == "fphase") {
      assign_optarg("runjam_phasespace_fname");
    } else if (longname == "fphase0") {
      assign_optarg("runjam_phasespace_fname0");

    } else if (longname == "hydrojet-dir") {
      assign_optarg("hydrojet_directory");
    } else if (longname == "hydrojet-ftemp")  {
      assign_optarg_int("hydrojet_kintmp");
    } else if (longname == "hydrojet-pce") {
      assign_optarg_int("hydrojet_eospce");
    } else if (longname == "hydrojet-bfree") {
      assign_optarg_int("hydrojet_baryonfree");
    } else if (longname == "hydrojet-dt") {
      assign_optarg_double("hydrojet_deltat");
    } else if (longname == "hydrojet-dx") {
      assign_optarg_double("hydrojet_deltax");
    } else if (longname == "hydrojet-dy") {
      assign_optarg_double("hydrojet_deltay");
    } else if (longname == "hydrojet-dh") {
      assign_optarg_double("hydrojet_deltah");

    } else if (longname == "switching-temperature") {
      assign_optarg_double("runjam_switching_temperature");
    } else {
      std::cerr << "unknown option '" << arg << "'" << std::endl;
      flag_error = true;
    }
  }

  bool process_option(int optchar) {
    switch (optchar) {
    case 's': assign_optarg_int("runjam_seed"); break;
    case 'n': assign_optarg_int("runjam_nevent"); break;
    case 't': assign_optarg_int("runjam_oversampling_factor"); break;
    case 'o': assign_optarg("runjam_output_directory"); break;
    case 'd': assign_optarg_int("runjam_phasespace_enabled"); break;
    case 'w': assign_optarg_int("runjam_swtich_weak_decay"); break;
    case 'i': assign_optarg_input(); break;
    default: return false;
    }
    return true;
  }
  void read_option() {
    char c;
    for (const char* p = arg + 1; p && (c = *p); p++) {
      if (p[1]) arg_optarg = p + 1;
      if (!process_option(c)) {
        std::cerr << "runjam:$" << arg_index << ": unknown option '-" << c << "'" << std::endl;
        flag_error = true;
      }
    }
  }

  bool read_assign() {
    const char* p = arg;
    while (std::isalnum(*p) && *p == '_') p++;
    if (p == arg || *p != '=') return false;

    std::string const name(arg, p);
    const char* const value = p + 1;
    ctx->set_value(name.c_str(), value);
    return true;
  }
public:
  int read(int argc, char** argv, runjam_context& ctx) {
    this->argc = argc;
    this->argv = argv;
    this->ctx = &ctx;
    this->flag_error = false;
    this->flag_help = false;
    this->flag_version = false;

    for (i = 1; i < argc; ) {
      arg = argv[arg_index = i++];
      arg_optarg = nullptr;
      if (arg[0] == '-') {
        if (arg[1] == '-') {
          this->read_longname_option();
        } else {
          this->read_option();
        }
      } else {
        // option 以外の文字列
        if (arg_index == 1) {
          this->subcommand = arg;
        } else if (!read_assign()) {
          std::cerr << "runjam:$" << arg_index << ": unrecognized argument '" << arg << "'" << std::endl;
          flag_error = true;
        }
      }
    }

    bool flag_exit = false;
    int exit_status = EXIT_SUCCESS;
    if (flag_error) exit_status = 2;

    if (flag_version) {
      cmd_version();
      flag_exit = 1;
    }
    if (flag_help) {
      cmd_help();
      flag_exit = 1;
    }
    if (flag_exit)
      std::exit(exit_status);

    return exit_status;
  }
};

void doCascade(runjam_context const& ctx, std::string const& type, std::string const& inputfile, std::string const& cascadeMode) {
  std::cout << "JAM hadronic cascade start" << std::endl;

  RunJam jam(ctx);
  jam.setMSTC(156, 1); // analysis of collision distribution
  jam.setMSTC(161, 0); // no analysis from jam internal subr.
  jam.setMSTC(162, 1); // Output collision histroy
  jam.setMSTC(165, 1); //
  //jam.setMSTC(41,0); // 0:no resonance decay after simulation.

  int const ntest = ctx.get_config("runjam_oversampling_factor", 1);
  bool const sw_weakdecay = ctx.get_config("runjam_swtich_weak_decay", 0);
  jam.setNumberOfTestParticle(ntest);
  if (sw_weakdecay) jam.setWeakDecay();

  std::cout << "jam event generation start" << std::endl;

  IParticleSample* psamp = CreateParticleSample(ctx, type, inputfile);
  if (!psamp) {
    std::cerr << "runjam: failed to initialize ParticleSample of type '" << type << "'." <<  std::endl;
    std::exit(1);
  }

  jam.generateEvent(psamp, cascadeMode);
  delete psamp;

  std::cout
    << "Average initial particle number from hydrojet:"
    << " before decay= " << jam.getIniAverageParticleNumber1()
    << " after decay= " << jam.getIniAverageParticleNumber2()
    << std::endl;
}

void savePhasespaceData(std::string fname, std::vector<Particle*> plist, ParticleIDType::value_type idtype) {
  std::FILE* f = std::fopen(fname.c_str(), "w");
  if (!f) {
    std::cerr << "runjam: failed to open the file '" << fname << "'" << std::endl;
    return;
  }

  int const nhadron = plist.size();
  std::fprintf(f, "%-4d 1\n", nhadron);
  for (std::vector<Particle*>::iterator i = plist.begin(); i != plist.end(); ++i) {
    Particle* particle = *i;

    //---------------------------------
    // (1) kf ... PDG particle code
    int kf;
    switch (idtype) {
    case ParticleIDType::HydroParticleID:
      {
        int ir = particle->id;
        kf = RunJam::sampleJamID(ir + 1);
      }
      break;
    case ParticleIDType::PDGCode:
      kf = particle->id;
      break;
    default:
      std::cerr << "runjam: invalid value of psamp->getParticleIdType()." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (kf == 0) {
      std::cerr << "runjam: invalid value of kf (PDG particle code)." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    //---------------------------------
    // (2) ks ... 安定粒子なら 1, 不安定粒子なら 2.
    int ks = 2;
    int const kc = Jam1::jamComp(kf); // jam internal particle code.
    if (Jam1::getPMAS(kc,2) <= 1e-7 || Jam1::getMDCY(kc,1) == 0
      || Jam1::getMDCY(kc,2) == 0 || Jam1::getMDCY(kc,3) == 0) ks = 1;

    //---------------------------------
    // (3) px,py,pz,m
    double const px = particle->px;
    double const py = particle->py;
    double const pz = particle->pz;
    double pe = particle->e;
    double m;
    if (pe < 0.0) {
      m = Jam1::jamMass(kf);
      pe = std::sqrt(px*px+py*py+pz*pz+m*m);
    } else {
      m = std::sqrt(pe*pe-(px*px+py*py+pz*pz));
    }

    //---------------------------------
    // (3) x,y,z,t
    double const x = particle->x;
    double const y = particle->y;
    double const z = particle->z;
    double const t = particle->t;

    //---------------------------------
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

  IParticleSample* psamp = CreateParticleSample(ctx1, type, inputfile);
  psamp->update();

  // 振り分け
  std::vector<Particle*> const& plist = psamp->getParticleList();
  std::vector<std::vector<Particle*> > phases((std::size_t) nevent);
  for (std::vector<Particle*>::const_iterator i = plist.begin(); i != plist.end(); ++i)
    phases[std::min(int(Random::getRand() * nevent), nevent - 1)].push_back(*i);

  // 保存
  for (int i = 0; i < nevent; i++) {
    std::vector<char> fn(outdir.size() + 50);
    std::sprintf(&fn[0], "%s/dens%06d_phasespace0.dat", outdir.c_str(), ibase + i);
    savePhasespaceData(&fn[0], phases[i], psamp->getParticleIdType());
  }

  delete psamp;
}

//-----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  runjam_context ctx;

  RunjamCommandlineArguments args;
  int const ext = args.read(argc, argv, ctx);
  if (ext) return ext;

  std::cout << "runjam [version " << PACKAGE_VERSION << PACKAGE_HASH << ", seed = " << ctx.seed() << "]" << std::endl;
  Random rand(ctx.seed());
  // Random rand(ctx.seed());
  Random::setRandom(&rand);

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
