#include <cstdlib>
#include <iostream>
#include <cstring>

#include "ksh/util.hpp"

#include "util/Random.hpp"
#include "util/PyRand.hpp"
#include "Hydro2Jam.hpp"

#include "spectra/ParticleSamplePhasespace.hpp"
#include "spectra/ParticleSampleViscous.hpp"
#include "spectra/ParticleSampleRead.hpp"
#include "spectra/ParticleSampleHydrojet.hpp"

#define PACKAGE_VERSION "0.1a"

using namespace std;
using namespace idt::hydro2jam;

void generatePhasespace0(hydro2jam_context const& ctx, std::string const& inputfile);

//-----------------------------------------------------------------------------
// parameters

//int seed = 1921;

int ntest = 1;
int baryonfree = 1;
int sw_weakdecay = 0; // this does not work now, sorry. do not put =1.
double deltat = 0.3;
double deltax = 0.3;
double deltay = 0.3;
double deltah = 0.3;

enum InitialType {
  InitialType_None = 0,
  InitialType_PHASE1,
  InitialType_C0LRF,
  InitialType_HYDRO,
  InitialType_PHASE,
  InitialType_PSAMPLE,
};

struct Hydro2jamCommandlineArguments {
  InitialType jamInitType;
  std::string fnameInitialPhasespaceData;

  double switchingTemperature;
public:
  Hydro2jamCommandlineArguments(){
    this->jamInitType = InitialType_None;
    this->switchingTemperature = -1.0;
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
      "    phase:FILENAME    load ICs from PHASESPACE (phasespace.dat)\n"
      "    phase1:FILENAME   load a single IC from PHASESPACE. This performs an\n"
      "                      additional check that PHASESPACE has only one event.\n"
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

public:
  int read(int argc, char** argv, hydro2jam_context& ctx) {
    for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        if (argv[i][1] == '-') {
          std::string longname = argv[i] + 2;
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
            if (ia == std::string::npos) {
              std::cerr << "hydro2jam:option(" << argv[i] << "): the argument of the option is not specified." << std::endl;
              std::exit(EXIT_SUCCESS);
            }
            if (!a.size()) {
              std::cerr << "hydro2jam:option(" << argv[i] << "): the argument of the option is empty." << std::endl;
              std::exit(EXIT_SUCCESS);
            }
            switchingTemperature = std::atof(a.c_str());
          } else {
            std::cerr << "unknowon option '" << argv[i] << "'" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        } else {
          if (!strcmp(argv[i], "-s")) {
            int const randomSeed = atoi(argv[++i]);
            ctx.set_value("hydro2jam_seed", randomSeed);
            std::cout << "hydro2jam: randomSeed is set to '" << randomSeed << "'"
                      << " [hydro2jam version " << PACKAGE_VERSION << "]" << std::endl;
          } else if (!strcmp(argv[i], "-t")) {
            ntest = atoi(argv[++i]);
          } else if (!strcmp(argv[i], "-n")) {
            ctx.set_value("hydro2jam_nevent", atoi(argv[++i]));
          } else if (!strcmp(argv[i], "-f")) {
            ctx.set_value("hydro2jam_phasespace_fname", argv[++i]);
          } else if (!strcmp(argv[i], "-f0")) {
            ctx.set_value("hydro2jam_phasespace_fname0", argv[++i]);
          } else if (!strcmp(argv[i], "-ftemp"))  {
            ctx.set_value("hydrojet_kintmp", atoi(argv[++i]));
          } else if (!strcmp(argv[i], "-dir")) {
            ctx.set_value("hydrojet_directory", argv[++i]);
          } else if (!strcmp(argv[i], "-dirJAM")) {
            ctx.set_value("hydrojet_output_directory", argv[++i]);
          } else if (!strcmp(argv[i], "-d")) {
            ctx.set_value("hydro2jam_phasespace_enabled", atoi(argv[++i]));
          } else if (!strcmp(argv[i], "-w")) {
            sw_weakdecay = atoi(argv[++i]);
          } else if (!strcmp(argv[i], "-pce")) {
            ctx.set_value("hydrojet_eospce", atoi(argv[++i]));
          } else if (!strcmp(argv[i], "-bfree")) {
            baryonfree = atoi(argv[++i]);
          } else if (!strcmp(argv[i], "-resodata")) {
            ctx.set_value("hydrojet_resodata", argv[++i]);
          } else if (!strcmp(argv[i], "-dx")) {
            deltax = atof(argv[++i]);
          } else if (!strcmp(argv[i], "-dy")) {
            deltay = atof(argv[++i]);
          } else if (!strcmp(argv[i], "-dh")) {
            deltah = atof(argv[++i]);
          } else if (!strcmp(argv[i], "-i")) {
            if (argv[i + 1] == 0) {
              std::cerr << "option `-i': missing an argument" << std::endl;
              std::exit(EXIT_FAILURE);
            }

            std::string spec = argv[++i];

            std::size_t const pos = spec.find(':');
            if (pos == std::string::npos) {
              std::cerr << "unrecognized input '-i " << argv[i] << "'" << std::endl;
              std::exit(EXIT_FAILURE);
            }

            std::string type(spec, 0, pos);
            this->fnameInitialPhasespaceData = spec.substr(pos + 1);

            if (type == "phase" || type == "kawaguchi")
              this->jamInitType = InitialType_PHASE;
            else if (type == "phase1")
              this->jamInitType = InitialType_PHASE1;
            else if (type == "c0lrf")
              this->jamInitType = InitialType_C0LRF;
            else if (type == "hydrojet")
              this->jamInitType = InitialType_HYDRO;
            else if (type == "psample")
              this->jamInitType = InitialType_PSAMPLE;
            else {
              std::cerr << "unrecognized input '-i " << argv[i] << "'" << std::endl;
              std::exit(EXIT_FAILURE);
            }

          } else {
            if (argv[i][0] == '-')
              std::cerr << "unknown option '" << argv[i] << "'" << std::endl;
            else
              std::cerr << "unknown argument '" << argv[i] << "'" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
      } else {
        // option 以外の文字列
        std::string const arg = argv[i];
        if (arg == "generate-phasespace0") {
          // test 用
          Random* rand = new PyRand(ctx.seed());
          Random::setRandom(rand);
          generatePhasespace0(ctx, this->fnameInitialPhasespaceData);
          delete rand;
          std::exit(EXIT_SUCCESS);
        }
      }
    }

    return 0;
  }
} args;

//-----------------------------------------------------------------------------
// routines

#include <cmath>
#include "jam/Jam1.hpp"

void savePhasespaceData(std::string fname, std::vector<Particle*> plist, ParticleIDType::value_type idtype) {
  std::FILE* f=std::fopen(fname.c_str(), "w");
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
    if (
      Jam1::getPMAS(kc,2) <= 1e-7 || Jam1::getMDCY(kc,1) == 0
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

#ifdef USE_JAM
void hadronicCascadeC0Lrf(hydro2jam_context const& ctx, Hydro2Jam* jam, std::string const& fname) {
  ResonanceListPCE reso(ctx);

  ParticleSampleViscous* psamp = new ParticleSampleViscous(&reso, fname);
  if (args.switchingTemperature > 0.0)
    psamp->setSwitchingTemperature(args.switchingTemperature);
  jam->generateEvent(psamp);
  delete psamp;
}

void hadronicCascadeHydrojet(hydro2jam_context const& ctx, Hydro2Jam* jam, std::string const& dname) {
  ResonanceListPCE reso(ctx);

  ParticleSampleFromHydrojet* psamp = new ParticleSampleFromHydrojet(&reso, dname);
  psamp->setDtau(deltat);
  psamp->setDh(deltah);
  psamp->setDx(deltax);
  psamp->setDy(deltay);
  if (args.switchingTemperature > 0.0)
    psamp->setSwitchingTemperature(args.switchingTemperature);

  jam->generateEvent(psamp);
  delete psamp;
}
#endif

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
  switch (args.jamInitType) {
  case InitialType_C0LRF:
    hadronicCascadeC0Lrf(ctx, jam, args.fnameInitialPhasespaceData);
    break;
  case InitialType_HYDRO:
    hadronicCascadeHydrojet(ctx, jam, args.fnameInitialPhasespaceData);
    break;
  case InitialType_PHASE1:
  case InitialType_PHASE:
    {
      ParticleSampleFromOversampledPhasespace* psamp = new ParticleSampleFromOversampledPhasespace(args.fnameInitialPhasespaceData);
      if (args.jamInitType == InitialType_PHASE1)
        psamp->setOverSamplingFactor(1);
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  case InitialType_PSAMPLE:
    {
      ParticleSampleRead* psamp = new ParticleSampleRead(args.fnameInitialPhasespaceData);
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  default: // default dir = "test"
    {
      IParticleSample* psamp = CreateParticleSampleHydrojet(ctx);
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  }

  std::cout << "Average initial particle number from hydrojet:"
            << " before decay= " << jam->getIniAverageParticleNumber1()
            << " after decay= " << jam->getIniAverageParticleNumber2()
            << std::endl;

  delete jam;
}

void generatePhasespace0(hydro2jam_context const& ctx, std::string const& inputfile) {
  int const nevent = ctx.nevent(1000);
  int const ibase = ctx.get_config("hydro2jam_ievent_begin", 0);
  std::string const outdir = ctx.outdir();

  ResonanceListPCE reso(ctx);

  ParticleSampleViscous* psamp = new ParticleSampleViscous(&reso, inputfile);
  psamp->setOverSamplingFactor(nevent);
  if (ctx.get_config("hydro2jam_turnsOffViscousEffect", false))
    psamp->setTurnsOffViscousEffect(true);
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
