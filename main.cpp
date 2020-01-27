#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "util/Random.hpp"
#include "util/PyRand.hpp"
#include "Hydro2Jam.hpp"
#include "spectra/ParticleSampleViscous.hpp"
#include "jam/Jam1.hpp"

#define PACKAGE_VERSION "0.1a"

using namespace idt::hydro2jam;

//-----------------------------------------------------------------------------
// parameters

struct Hydro2jamCommandlineArguments {
  std::string subcommand;
  std::string initType;
  std::string initPath;

public:
  Hydro2jamCommandlineArguments(){
    this->subcommand = "cascade";
    this->initType = "c0lrf";
    this->initPath = "hypersurface_v1.txt";

    // this->initType = "hydrojet.original";
    // this->initPath = "test";
  }

private:
  static void cmd_help(){
    std::printf(
      "usage: hydro2jam [subcommand] [options]\n"
      "\n"
      "SUBCOMMAND\n"
      "  cascade (default)\n"
      "  decay\n"
      "  sample\n"
      "\n"
      "OPTIONS and VARIABLES\n"
      "      hydro2jam_cascade_mode=SUBCOMMAND\n"
      "  -n, hydro2jam_nevent=INT [1]              number of events to process\n"
      "  -s, hydro2jam_seed=INT [18371]            seed for random numbers\n"
      "  -t, hydro2jam_oversampling_factor=NUM [1] ntest\n"
      "  -w, hydro2jam_switch_weak_decay=INT [0]   sw_weakdecay {0 | 1}\n"
      "      hydro2jam_phi_decays=BOOL [true]\n"
      "      hydro2jam_decay_only=BOOL [false]\n"
      "  -resodata, hydro2jam_resodata=FILE [data/ResonanceJam.dat] resodata\n"
      "\n"
      " Output options\n"
      "  -dirJAM, hydro2jam_output_directory=DIR [jam] directory of output files\n"
      "  -d,  hydro2jam_phasespace_enabled=INT   [1] dump phasespace\n"
      "  -f,  hydro2jam_phasespace_fname=FILE    [phasespace.dat] output filename\n"
      "  -f0, hydro2jam_phasespace_fname0=FILE   [phasespace0.dat] output filename\n"
      "       hydro2jam_output_phbin=BOOL        [false] output binary phasespace\n"
      "       hydro2jam_output_phbin0=BOOL       [false] output binary phasespace0\n"
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
      "                  hydro2jam \"particlesample_pos.dat\"\n"
      "\n"
      " Options for hydrojet hypersurface\n"
      "  -ftemp, hydrojet_kintmp=INT     [5] freezeout temperature type\n"
      "  -pce,   hydrojet_eospce=INT     [6] eos_pce\n"
      "  -bfree, hydrojet_baryonfree=INT [1] baryonfree\n"
      "  -dir,   hydrojet_directory=DIR  [test] directory of freezeout.dat\n"
      "  -dt,    hydrojet_deltat=NUM     [0.3] delta tau\n"
      "  -dx,    hydrojet_deltax=NUM     [0.3] delta x\n"
      "  -dy,    hydrojet_deltay=NUM     [0.3] delta y\n"
      "  -dh,    hydrojet_deltah=NUM     [0.3] delta eta\n"
      "\n"
      " Options for c0lrf sampler\n"
      "  --switching-temperature, hydro2jam_switching_temperature=TEMP [155]\n"
      "          an advice to switching temperature in MeV\n"
      "  hydro2jam_turnsOffViscousEffect=INT\n"
      "\n"
      " Other options\n"
      "  --help          show this help\n"
      "\n"
      "SAMPLE\n"
      "\n"
      "$ ./hydro2jam -s 12345 -dir hydro -dirJAM jam\n"
      "$ ./hydro2jam -s 12345 -i phase:phasespace0.in -dirJAM jam\n"
      "$ ./hydro2jam decay -s 12345 -i phase:phasespace0.in -dirJAM jam\n"
      "\n"
    );
  }

private:
  int argc;
  char** argv;
  hydro2jam_context* ctx;

  int i;
  char* arg;
  bool flag_error;

  const char* get_optarg() {
    if (i < argc) {
      return argv[i++];
    } else {
      std::cerr << "hydro2jam:$" << i-1 << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
      return NULL;
    }
  }

  void assign_optarg(const char* key) {
    if (const char* optarg = get_optarg())
      ctx->set_value(key, optarg);
  }

  void assign_optarg_double(const char* key, const char* value) {
    if (*value) {
      ctx->set_value(key, std::atof(value));
    } else {
      std::cerr << "hydro2jam:option(" << arg << "): the argument of the option is empty." << std::endl;
      flag_error = true;
    }
  }

  void assign_optarg_double(const char* key) {
    if (const char* optarg = get_optarg())
      assign_optarg_double(key, optarg);
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
      std::cerr << "unrecognized input '-i " << arg << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->initType = spec.substr(0, pos);
    this->initPath = spec.substr(pos + 1);
  }

private:
  void read_longname_option() {
    std::string longname = arg + 2;
    std::string a;
    std::size_t ia;
    if ((ia = longname.find('=', 0)) != std::string::npos) {
      a = longname.substr(ia + 1);
      longname = longname.substr(0, ia);
    }

    if (longname == "help") {
      cmd_help();
      std::exit(EXIT_SUCCESS);
    } else if (longname == "switching-temperature") {
      if (ia != std::string::npos)
        assign_optarg_double("hydro2jam_switching_temperature", a.c_str());
      else
        assign_optarg_double("hydro2jam_switching_temperature");
    } else {
      std::cerr << "unknown option '" << arg << "'" << std::endl;
      flag_error = true;
    }
  }

  void read_option() {
    if (std::strcmp(arg, "-s") == 0) {
      assign_optarg_int("hydro2jam_seed");
    } else if (std::strcmp(arg, "-n") == 0) {
      assign_optarg_int("hydro2jam_nevent");
    } else if (std::strcmp(arg, "-t") == 0) {
      assign_optarg_int("hydro2jam_oversampling_factor");
    } else if (std::strcmp(arg, "-dirJAM") == 0) {
      assign_optarg("hydro2jam_output_directory");
    } else if (std::strcmp(arg, "-w") == 0) {
      assign_optarg_int("hydro2jam_swtich_weak_decay");
    } else if (std::strcmp(arg, "-resodata") == 0) {
      assign_optarg("hydro2jam_resodata");
    } else if (std::strcmp(arg, "-d") == 0) {
      assign_optarg_int("hydro2jam_phasespace_enabled");
    } else if (std::strcmp(arg, "-f") == 0) {
      assign_optarg("hydro2jam_phasespace_fname");
    } else if (std::strcmp(arg, "-f0") == 0) {
      assign_optarg("hydro2jam_phasespace_fname0");
    } else if (std::strcmp(arg, "-dir") == 0) {
      assign_optarg("hydrojet_directory");
    } else if (std::strcmp(arg, "-ftemp") == 0)  {
      assign_optarg_int("hydrojet_kintmp");
    } else if (std::strcmp(arg, "-pce") == 0) {
      assign_optarg_int("hydrojet_eospce");
    } else if (std::strcmp(arg, "-bfree") == 0) {
      assign_optarg_int("hydrojet_baryonfree");
    } else if (std::strcmp(arg, "-dt") == 0) {
      assign_optarg_double("hydrojet_deltat");
    } else if (std::strcmp(arg, "-dx") == 0) {
      assign_optarg_double("hydrojet_deltax");
    } else if (std::strcmp(arg, "-dy") == 0) {
      assign_optarg_double("hydrojet_deltay");
    } else if (std::strcmp(arg, "-dh") == 0) {
      assign_optarg_double("hydrojet_deltah");
    } else if (std::strcmp(arg, "-i") == 0) {
      assign_optarg_input();
    } else {
      std::cerr << "hydro2jam:$" << i-1 << ": unknown option '" << arg << "'" << std::endl;
      flag_error = true;
    }
  }
public:
  int read(int argc, char** argv, hydro2jam_context& ctx) {
    this->argc = argc;
    this->argv = argv;
    this->ctx = &ctx;
    this->flag_error = false;

    for (i = 1; i < argc; ) {
      arg = argv[i++];
      if (arg[0] == '-') {
        if (arg[1] == '-') {
          this->read_longname_option();
        } else {
          this->read_option();
        }
      } else {
        // option 以外の文字列
        if (i-1 == 1) {
          this->subcommand = arg;
        } else {
          std::cerr << "hydro2jam:$" << i-1 << ": unrecognized argument '" << arg << "'" << std::endl;
          flag_error = true;
        }
      }
    }

    if (flag_error) return 2;
    return 0;
  }
};

void doCascade(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile, std::string const& cascadeMode) {
  std::cout << "JAM hadronic cascade start" << std::endl;

  Hydro2Jam jam(ctx);
  jam.setMSTC(156, 1); // analysis of collision distribution
  jam.setMSTC(161, 0); // no analysis from jam internal subr.
  jam.setMSTC(162, 1); // Output collision histroy
  jam.setMSTC(165, 1); //
  //jam.setMSTC(41,0); // 0:no resonance decay after simulation.

  int const ntest = ctx.get_config("hydro2jam_oversampling_factor", 1);
  bool const sw_weakdecay = ctx.get_config("hydro2jam_swtich_weak_decay", 0);
  jam.setNumberOfTestParticle(ntest);
  if (sw_weakdecay) jam.setWeakDecay();

  std::cout << "jam event generation start" << std::endl;

  IParticleSample* psamp = CreateParticleSample(ctx, type, inputfile);
  if (!psamp) {
    std::cerr << "hydro2jam: failed to initialize ParticleSample of type '" << type << "'." <<  std::endl;
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
    std::cerr << "hydro2jam: failed to open the file '" << fname << "'" << std::endl;
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
        kf = Hydro2Jam::sampleJamID(ir + 1);
      }
      break;
    case ParticleIDType::PDGCode:
      kf = particle->id;
      break;
    default:
      std::cerr << "hydro2jam: invalid value of psamp->getParticleIdType()." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (kf == 0) {
      std::cerr << "hydro2jam: invalid value of kf (PDG particle code)." << std::endl;
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

void doGeneratePhasespace0(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile) {
  int const nevent = ctx.nevent(1000);
  int const ibase = ctx.get_config("hydro2jam_ievent_begin", 0);
  double const ntest = ctx.get_config("hydro2jam_oversampling_factor", 1.0);
  std::string const outdir = ctx.outdir();

  hydro2jam_context ctx1(ctx);
  ctx1.set_value("hydro2jam_nevent", 1);
  ctx1.set_value("hydro2jam_oversampling_factor", nevent * ntest);

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
  hydro2jam_context ctx;

  Hydro2jamCommandlineArguments args;
  int const ext = args.read(argc, argv, ctx);
  if (ext) return ext;

  std::cout << "hydro2jam [version " << PACKAGE_VERSION << ", seed = " << ctx.seed() << "]" << std::endl;
  Random rand(ctx.seed());
  // Random rand(ctx.seed());
  Random::setRandom(&rand);

  if (args.subcommand == "generate-phasespace0") {
    // test 用
    doGeneratePhasespace0(ctx, args.initType, args.initPath);
  } else if (args.subcommand == "test-viscous-correction-integration") {
    return checkViscousCooperFryeInterpolated(true);
  } else if (args.subcommand == "cascade") {
    std::string mode = ctx.get_config<std::string>("hydro2jam_cascade_mode", "cascade");
    if (mode == "cascade" && ctx.get_config("hydro2jam_decay_only", false)) mode = "decay";
    doCascade(ctx, args.initType, args.initPath, mode);
  } else if (args.subcommand == "decay" || args.subcommand == "sample") {
    doCascade(ctx, args.initType, args.initPath, args.subcommand);
  } else {
    std::cerr << "hydro2jam: unknown subcommand ' " << args.subcommand << "'" << std::endl;
    return 2;
  }

  return 0;
}
