#include <cstdlib>
#include <iostream>
#include <cstring>

#include "util/Random.hpp"
#include "util/PyRand.hpp"
#include "Hydro2Jam.hpp"
#include "spectra/ParticleSampleViscous.hpp"

#define PACKAGE_VERSION "0.1a"

using namespace std;
using namespace idt::hydro2jam;

void generatePhasespace0(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile);

//-----------------------------------------------------------------------------
// parameters

int ntest = 1;
int baryonfree = 1;
int sw_weakdecay = 0; // this does not work now, sorry. do not put =1.

struct Hydro2jamCommandlineArguments {
  std::string initType;
  std::string initPath;

public:
  Hydro2jamCommandlineArguments(){
    this->initType = "c0lrf";
    this->initPath = "hypersurface_v1.txt";

    // this->initType = "hydrojet.original";
    // this->initPath = "test";
  }

private:
  static void cmd_help(){
    std::printf(
      "usage: hydro2jam [options]\n"
      "\n"
      "OPTIONS\n"
      "  -n INT          number of events to process [default: 1]\n"
      "  -s INT          seed\n"
      "  -t INT          ntest\n"
      "  -f FILE         output phasespace.dat\n"
      "  -ftemp FLOAT    freezeout temperature (kintmp)\n"
      "  -f0 FILE        output phasespace0.dat\n"
      "  -dir DIR        directory of freezeout.dat\n"
      "  -dirJAM DIR     directory of output files\n"
      "  -d INT          dumpPhaseSpace = 0 | 1\n"
      "  -w INT          sw_weakdecay = 0 | 1\n"
      "  -pce INT        eos_pce\n"
      "  -bfree INT      baryonfree\n"
      "  -resodata STR   resodata\n"
      "  -dx FLOAT       deltax\n"
      "  -dy FLOAT       deltay\n"
      "  -dh FLOAT       deltah\n"
      "\n"
      "  -i ICSPEC       specify initial condition\n"
      "    c0lrf:FILENAME    sample particles using the rfh output\n"
      "                      \"hypersurface_v1.txt\"\n"
      "    hydrojet:DIR      load IC from hydrojet output\n"
      "                      (DIR/freezeout.dat, DIR/position.dat)\n"
      "    phase:FILENAME    load ICs from \"phasespace.dat\" format file.\n"
      "    phase1:FILENAME   load a single IC from PHASESPACE. This performs an\n"
      "                      additional check that PHASESPACE has only one event.\n"
      "    phbin:FILENAME    load ICs from \"ph000k.bin\" format file.\n"
      "    psample:FILENAME  read a particle list from the file in the format of\n"
      "                      hydro2jam \"particlesample_pos.dat\"\n"
      "\n"
      "  --help          show this help\n"
      "  --switching-temperature=TEMP [MeV]\n"
      "\n"
      "ENVIRONMENT VARIABLES\n"
      "  hydro2jam_phi_decays  [true]\n"
      "  hydro2jam_decay_only  [false]\n"
      "\n"
      "SAMPLE\n"
      "$ ./hydro2jam -s 12345 -dir hydro -dirJAM jam\n"
      "$ ./hydro2jam -s 12345 -i phase:phasespace0.in -dirJAM jam\n"
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

  void assign_optarg(const char* key) {
    if (i < argc) {
      ctx->set_value(key, argv[i++]);
    } else {
      std::cerr << "hydro2jam:arg#" << i << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
    }
  }

  void assign_optarg_double(const char* key) {
    if (i < argc) {
      ctx->set_value(key, std::atof(argv[i++]));
    } else if (!*argv[i]) {
      std::cerr << "hydro2jam:option(" << arg << "): the argument of the option is empty." << std::endl;
      flag_error = true;
    } else {
      std::cerr << "hydro2jam:arg#" << i << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
    }
  }

  void assign_optarg_int(const char* key) {
    if (i < argc) {
      ctx->set_value(key, std::atoi(argv[i++]));
    } else {
      std::cerr << "hydro2jam:arg#" << i << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
    }
  }

  void assign_optarg_input() {
    if (i >= argc) {
      std::cerr << "hydro2jam:arg#" << i << " (" << arg << "): missing optional argument." << std::endl;
      flag_error = true;
      return;
    }

    std::string spec = argv[i++];

    std::size_t const pos = spec.find(':');
    if (pos == std::string::npos) {
      std::cerr << "unrecognized input '-i " << arg << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->initType = spec.substr(0, pos);
    this->initPath = spec.substr(pos + 1);
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
          std::string longname = arg + 2;
          std::string a;
          std::size_t ia;
          if ((ia = longname.find('=', 0)) != std::string::npos) {
            a = longname.substr(ia + 1);
            longname = longname.substr(0, ia);
          }

          if (longname == "debug20150102") {
            std::exit(checkViscousCooperFryeInterpolated(true));
          } else if (longname == "check") {
            checkViscousCooperFryeInterpolated(false);
            std::exit(EXIT_SUCCESS);
          } else if (longname == "help") {
            cmd_help();
            std::exit(EXIT_SUCCESS);
          } else if (longname == "switching-temperature") {
            assign_optarg_double("hydro2jam_switching_temperature");
          } else {
            std::cerr << "unknowon option '" << arg << "'" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        } else {
          if (!std::strcmp(arg, "-s")) {
            assign_optarg_int("hydro2jam_seed");
            std::cout << "hydro2jam: randomSeed is set to '" << ctx.get_value("hydro2jam_seed") << "'"
                      << " [hydro2jam version " << PACKAGE_VERSION << "]" << std::endl;
          } else if (!std::strcmp(arg, "-t")) {
            ntest = std::atoi(argv[++i]);
          } else if (!std::strcmp(arg, "-n")) {
            assign_optarg_int("hydro2jam_nevent");
          } else if (!std::strcmp(arg, "-f")) {
            assign_optarg("hydro2jam_phasespace_fname");
          } else if (!std::strcmp(arg, "-f0")) {
            assign_optarg("hydro2jam_phasespace_fname0");
          } else if (!std::strcmp(arg, "-ftemp"))  {
            assign_optarg_int("hydrojet_kintmp");
          } else if (!std::strcmp(arg, "-dir")) {
            assign_optarg("hydrojet_directory");
          } else if (!std::strcmp(arg, "-dirJAM")) {
            assign_optarg("hydro2jam_output_directory");
          } else if (!std::strcmp(arg, "-d")) {
            assign_optarg_int("hydro2jam_phasespace_enabled");
          } else if (!std::strcmp(arg, "-w")) {
            sw_weakdecay = std::atoi(argv[++i]);
          } else if (!std::strcmp(arg, "-pce")) {
            assign_optarg_int("hydrojet_eospce");
          } else if (!std::strcmp(arg, "-bfree")) {
            baryonfree = std::atoi(argv[++i]);
          } else if (!std::strcmp(arg, "-resodata")) {
            assign_optarg("hydro2jam_resodata");
          } else if (!std::strcmp(arg, "-dt")) {
            assign_optarg_double("hydrojet_deltat");
          } else if (!std::strcmp(arg, "-dx")) {
            assign_optarg_double("hydrojet_deltax");
          } else if (!std::strcmp(arg, "-dy")) {
            assign_optarg_double("hydrojet_deltay");
          } else if (!std::strcmp(arg, "-dh")) {
            assign_optarg_double("hydrojet_deltah");
          } else if (!std::strcmp(arg, "-i")) {
            assign_optarg_input();
          } else {
            std::cerr << "hydro2jam:#" << i << ": unknown option '" << arg << "'" << std::endl;
            flag_error = true;
          }
        }
      } else {
        // option 以外の文字列
        if (arg == "generate-phasespace0") {
          // test 用
          Random* rand = new PyRand(ctx.seed());
          Random::setRandom(rand);
          generatePhasespace0(ctx, this->initType, this->initPath);
          delete rand;
          std::exit(EXIT_SUCCESS);
        } else {
          std::cerr << "hydro2jam:#" << i << ": unknown argument '" << arg << "'" << std::endl;
          flag_error = true;
        }
      }
    }

    if (flag_error) return 2;
    return 0;
  }
} args;

//-----------------------------------------------------------------------------
// routines

#include <cmath>
#include "jam/Jam1.hpp"

void savePhasespaceData(std::string fname, std::vector<Particle*> plist, ParticleIDType::value_type idtype) {
  std::FILE* f = std::fopen(fname.c_str(), "w");
  if (!f) {
    std::cerr << "hydro2jam(savePhasespaceData): failed to open the file '" << fname << "'" << std::endl;
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
      std::cerr << "hydro2jam.cxx(savePhasespaceData): invalid value of psamp->getParticleIdType()." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (kf == 0) {
      std::cerr << "hydro2jam.cxx(savePhasespaceData): invalid value of kf (PDG particle code)." << std::endl;
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

void hadronicCascade(hydro2jam_context const& ctx) {
  std::cout << "JAM hadronic cascade start" << std::endl;

  Hydro2Jam* jam = new Hydro2Jam(ctx);
  jam->setMSTC(156, 1);        // analysis of collision distribution
  jam->setMSTC(161, 0);        // no analysis from jam internal subr.
  jam->setMSTC(162, 1);        // Output collision histroy
  jam->setMSTC(165, 1);        //
  //jam->setMSTC(41,0);         // 0:no resonance decay after simulation.

  jam->setNumberOfTestParticle(ntest);
  if (sw_weakdecay) jam->setWeakDecay(); //allow weak decays

  cout << "jam event generation start" << endl;

  if (IParticleSample* psamp = CreateParticleSample(ctx, args.initType, args.initPath)) {
    jam->generateEvent(psamp);
    delete psamp;
  } else {
    std::cerr << "hydro2jam: failed to initialize ParticleSample of type '" << args.initType << "'." <<  std::endl;
    std::exit(1);
  }

  std::cout << "Average initial particle number from hydrojet:"
            << " before decay= " << jam->getIniAverageParticleNumber1()
            << " after decay= " << jam->getIniAverageParticleNumber2()
            << std::endl;

  delete jam;
}

void generatePhasespace0(hydro2jam_context const& ctx, std::string const& type, std::string const& inputfile) {
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

  int const ext = args.read(argc, argv, ctx);
  if (ext) return ext;

  // Set Random number generator.
  Random* rand = new PyRand(ctx.seed());
  // Random* rand = new Random(randomSeed);
  Random::setRandom(rand);

  hadronicCascade(ctx);

  delete rand;

  return 0;
}
