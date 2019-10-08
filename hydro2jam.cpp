//
// example main program for the jam afterburner
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "user/HydroSpec.h"
#include "uty/Random.h"
#include "uty/PyRand.h"
#include "user/Hydro2Jam.h"

#include "spectra/ElementReso.h"
#include "spectra/ParticleSampleDebugJamAsymmetry.h"
#include "spectra/ParticleSampleFromPhasespaceDat.h"
#include "spectra/ParticleSampleFromOversampledPhasespace.h"
#include "spectra/ParticleSampleViscous.h"

#include "ksh/util.h"

using namespace std;

void generatePhasespaceData20141020(int ibase, std::string dirJAM);

//-----------------------------------------------------------------------------
// parameters

#define DATADIR "dict"
int seed = 1921;
int nhistout = 10;   // output histgram every "nhistout" event.
std::string dir = "test";
std::string dirJAM = "jam";
int randomSeed = 18371;
int eos_pce = 6;
int freezeoutTemp = 5;// 1:80MeV 2:100MeV 3:120MeV 4:140MeV 5:160MeV, This option works only for eos_pce = 1.
std::string resodata = DATADIR "/ResonanceJam.dat";
int baryonfree = 1;
int ntest = 1;
int sw_weakdecay = 0; // this does not work now, sorry. do not put =1.
int dumpPhaseSpace = 1; // =1: output phase space data in JAM.
std::string fnamePS = "phasespace.dat";
std::string fnamePS0 = "phasespace0.dat";
double deltat = 0.3;
double deltax = 0.3;
double deltay = 0.3;
double deltah = 0.3;

bool outputHistogramFlag = false;
enum InitialType {
  InitialType_None = 0,
  InitialType_PHASE,
  InitialType_C0LRF,
  InitialType_HYDRO,
  InitialType_Debug201304,
  InitialType_KAWAGUCHI,
};
struct Hydro2jamCommandlineArguments {
  InitialType jamInitType;
  std::string fnameInitialPhasespaceData;

  int nev;

  int makeElem;

  double switchingTemperature;
public:
  Hydro2jamCommandlineArguments(){
    // default values

    this->jamInitType=InitialType_None;

    // 2014-06-18 KM: changed from 10 to 1 since 1 is more reasonable in the e-by-e picture
    this->nev=1;

    this->makeElem=0;
    this->switchingTemperature=-1.0;
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
      "  -elm INT        makeElem\n"
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
      "  -ho INT         nhistout\n"
      "  -dx FLOAT       deltax\n"
      "  -dy FLOAT       deltay\n"
      "  -dh FLOAT       deltah\n"
      "\n"
      "  -i ICSPEC       specify initial condition\n"
      "      debug201304          IC for debugging (see ParticleSampleDebugJamAsymmetry.cxx)\n"
      "      phase:PHASESPACE     load IC from PHASESPACE (phasespace.dat)\n"
      "      c0lrf:HYPERSURFACE   load IC from rfh output (hypersurface.txt)\n"
      "      hydrojet:DIR         load IC from hydrojet output\n"
      "                           (DIR/freezeout.dat, DIR/position.dat)\n"
      "      kawaguchi:PHASESPACE 暫定的機能\n"
      "\n"
      "  --disable-hist  do not output *.hist files [default]\n"
      "  --enable-hist   output *.hist files\n"
      "  --help          show this help\n"
      "  --switching-temperature=TEMP [MeV]\n"
      "\n"
      "ENVIRONMENT VARIABLES\n"
      "  hydro2jam_phi_decays  [true]\n"
      "  hydro2jam_decay_only  [false]\n"
      "\n"
      "SAMPLE\n"
      "$ ./hydro2jam -s 12345 -dir hydro -dirJAM jam\n"
      "$ ./hydro2jam -s 12345 -i phase:phasespace0.in -dirJAM jam --disable-hist\n"
      "\n"
    );
  }

public:
  int read(int argc, char** argv) {
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

          if (longname == "disable-hist") {
            outputHistogramFlag = false;
          } else if (longname == "enable-hist") {
            outputHistogramFlag = true;
          } else if (longname == "debug201304") {
            this->jamInitType = InitialType_Debug201304;
          } else if (longname == "debug20150102") {
            int checkViscousCooperFryeInterpolated(bool debug);
            std::exit(checkViscousCooperFryeInterpolated(true));
          } else if (longname == "check") {
            int checkViscousCooperFryeInterpolated(bool debug);
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
            randomSeed = atoi(argv[++i]);
            std::cout << "hydro2jam: randomSeed is set to '" << randomSeed << "' (hydrojet version = " << PACKAGE_VERSION << ")" << std::endl;
          } else if (!strcmp(argv[i], "-t"))      ntest = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-n"))        this->nev = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-elm"))      this->makeElem = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-f")) {
            fnamePS = argv[++i];
          }else if (!strcmp(argv[i], "-ftemp")) {
            // ? バグ? (本当は別のオプションを割り当てるべきではないか?)
            freezeoutTemp = atoi(argv[++i]);
          }else if (!strcmp(argv[i], "-f0"))      fnamePS0 = argv[++i];
          else if (!strcmp(argv[i], "-dir"))      dir = argv[++i];
          else if (!strcmp(argv[i], "-dirJAM"))   dirJAM = argv[++i];
          else if (!strcmp(argv[i], "-d"))        dumpPhaseSpace = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-w"))        sw_weakdecay = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-pce"))      eos_pce = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-bfree"))    baryonfree = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-resodata")) resodata = argv[++i];
          else if (!strcmp(argv[i], "-ho"))       nhistout = atoi(argv[++i]);
          else if (!strcmp(argv[i], "-dx"))       deltax = atof(argv[++i]);
          else if (!strcmp(argv[i], "-dy"))       deltay = atof(argv[++i]);
          else if (!strcmp(argv[i], "-dh"))       deltah = atof(argv[++i]);
          else if (!strcmp(argv[i], "-i")) {
            if (argv[i + 1] == 0) {
              std::cerr << "option `-i': missing an argument" << std::endl;
              std::exit(EXIT_FAILURE);
            }

            std::string spec = argv[++i];
            if(spec == "debug201304") {
              this->jamInitType = InitialType_Debug201304;
            } else if (0 == spec.compare(0,6,"phase:",6)) {
              this->jamInitType=InitialType_PHASE;
              this->fnameInitialPhasespaceData = spec.substr(6);
            } else if (0 == spec.compare(0,6,"c0lrf:",6)) {
              this->jamInitType=InitialType_C0LRF;
              this->fnameInitialPhasespaceData = spec.substr(6);
            } else if (0 == spec.compare(0,9,"hydrojet:",9)) {
              this->jamInitType=InitialType_HYDRO;
              this->fnameInitialPhasespaceData = spec.substr(9);
            } else if (0 == spec.compare(0,10,"kawaguchi:",10)) {
              this->jamInitType=InitialType_KAWAGUCHI;
              this->fnameInitialPhasespaceData = spec.substr(10);
            } else {
              std::cerr << "unrecognized option '-i " << argv[i] << "'" << std::endl;
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
        if (arg == "generate_phasespace0") {
          // test 用
          Random* rand = new PyRand(randomSeed);
          Random::setRandom(rand);
          {
            int const iEv_begin = kashiwa::getenvAsInt("iEv_begin", 50000);
            generatePhasespaceData20141020(iEv_begin, dirJAM);
          }
          delete rand;
          std::exit(EXIT_SUCCESS);
        }
      }
    }

    return 0;
  }
} args;


#undef DATADIR

//-----------------------------------------------------------------------------
// routines

#include <cmath>
#ifdef USE_JAM
# include "jam/Jam1.h"
#endif

void savePhasespaceData(std::string fname, std::vector<HydroParticleCF*> plist, ParticleIDType::value_type idtype) {
  std::FILE* f=std::fopen(fname.c_str(), "w");
  if (!f) {
    std::cerr << "hydro2jam(savePhasespaceData): failed to open the file '" << fname << "'" << std::endl;
    return;
  }

  int const nhadron = plist.size();
  std::fprintf(f, "%-4d 1\n", nhadron);
  for (std::vector<HydroParticleCF*>::iterator i = plist.begin(); i != plist.end(); ++i) {
    HydroParticleCF* particle = *i;

    //---------------------------------
    // (1) kf ... PDG particle code
    int kf;
    switch (idtype) {
    case ParticleIDType::HydroParticleID:
      {
        int ir = particle->getID();
        kf = Hydro2Jam::getJamID(ir + 1);
      }
      break;
    case ParticleIDType::PDGCode:
      kf = particle->getID();
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
#ifdef USE_JAM
    int const kc = Jam1::jamComp(kf); // jam internal particle code.
    if (
      Jam1::getPMAS(kc,2) <= 1e-7 || Jam1::getMDCY(kc,1) == 0
      || Jam1::getMDCY(kc,2) == 0 || Jam1::getMDCY(kc,3) == 0) ks = 1;
#endif
    //---------------------------------
    // (3) px,py,pz,m
    double const px = particle->getPx();
    double const py = particle->getPy();
    double const pz = particle->getPz();
    double pe = particle->getPe();
    double m;
    if (pe < 0.0) {
#ifdef USE_JAM
      m = Jam1::jamMass(kf);
      pe = std::sqrt(px*px+py*py+pz*pz+m*m);
#else
      m = -1.0;
#endif
    } else {
      m = std::sqrt(pe*pe-(px*px+py*py+pz*pz));
    }
    //---------------------------------
    // (3) x,y,z,t
    double const x = particle->getX();
    double const y = particle->getY();
    double const z = particle->getZ();
    double const t = particle->getT();
    //---------------------------------

    std::fprintf(
      f, "%5d %9d %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g\n",
      ks, kf, px, py, pz, m, x, y, z, t);
  }
  std::fprintf(f, "-999\n");
  std::fclose(f);
}

void getHydroSpectra() {
  HydroSpec *spec = new HydroSpec();
  spec->setDirectory("data");  // name of directory for output spectra.
  spec->setDirectoryResonance(dir);
  spec->setFreezeoutData(dir + "/freezeout.dat");
  spec->setPositionData(dir + "/position.dat");
  spec->setFileEccentricity(dir + "/eccentricity.dat");
  spec->setKineticFreezeOutTemp(freezeoutTemp);
  spec->setEoSOpt(eos_pce);
  spec->setResoData(resodata);
  spec->setDirectSpectra(0);
  spec->setDirectSpectraPt(0);
  spec->setResonanceSpectra(0);
  spec->setResonanceElement(1);
  spec->setSeed(120419);
  spec->setRapidity(0);
  spec->setRapidityPt(2);

  spec->makeResonanceElement();

  delete spec;
}

#ifdef USE_JAM
void hadronicCascadeC0Lrf(Hydro2JamInitParams& iparam, Hydro2Jam* jam, std::string const& fname) {
  ResonanceListPCE reso(iparam.kintmp, iparam.eos_pce, ::resodata);

  ParticleSampleViscous* psamp = new ParticleSampleViscous(&reso, fname);
  if (args.switchingTemperature > 0.0)
    psamp->setSwitchingTemperature(args.switchingTemperature);
  jam->generateEvent(psamp);
  delete psamp;
}

void hadronicCascadeHydrojet(Hydro2JamInitParams& iparam, Hydro2Jam* jam, std::string const& dname) {
  ResonanceListPCE reso(iparam.kintmp, iparam.eos_pce, ::resodata);

  ParticleSampleFromHydrojet* psamp = new ParticleSampleFromHydrojet(&reso,dname);
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

void hadronicCascade(int mevent) {
#ifdef USE_JAM
  cout << "JAM hadronic cascade start" << endl;
  //int seed = (int)Random::getRand()
  Hydro2Jam* jam;

  Hydro2JamInitParams iparam;
  iparam.mevent              = mevent;
  iparam.seed                = randomSeed;
  iparam.dir_reso            = dir;
  iparam.kintmp              = freezeoutTemp;
  iparam.eos_pce             = eos_pce;
  iparam.dir                 = dirJAM;
  iparam.dpd                 = dumpPhaseSpace;
  iparam.fnamePS             = fnamePS;
  iparam.fnamePS0            = fnamePS0;
  iparam.outputHistogramFlag = outputHistogramFlag;
  jam = new Hydro2Jam(iparam);

  jam->setNumberOfHistgramOutput(nhistout);
  jam->setResoData(resodata);
  jam->setIsFile(0);
  jam->setMSTC(156,1);        // analysis of collision distribution
  jam->setMSTC(161,0);        // no analysis from jam internal subr.
  jam->setMSTC(162,1);        // Output collision histroy
  jam->setMSTC(165,1);        //
  //jam->setMSTC(41,0);         // 0:no resonance decay after simulation.

  jam->setNumberOfTestParticle(ntest);
  if (sw_weakdecay) jam->setWeakDecay(); //allow weak decays

  cout << "jam event generation start" << endl;
  switch (args.jamInitType) {
  case InitialType_PHASE:
    {
      IParticleSample* psamp = new ParticleSampleFromPhasespaceDat(args.fnameInitialPhasespaceData);
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  case InitialType_C0LRF:
    hadronicCascadeC0Lrf(iparam, jam, args.fnameInitialPhasespaceData);
    break;
  case InitialType_HYDRO:
    hadronicCascadeHydrojet(iparam, jam, args.fnameInitialPhasespaceData);
    break;
  case InitialType_Debug201304:
    {
      IParticleSample* psamp = new ParticleSampleDebugJamAsymmetry;
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  case InitialType_KAWAGUCHI:
    {
      IParticleSample* psamp = new ParticleSampleFromOversampledPhasespace(args.fnameInitialPhasespaceData);
      jam->generateEvent(psamp);
      delete psamp;
    }
    break;
  default: // default dir = "test"
    jam->generateEventFromHypersurfaceFiles(
      dir+"/freezeout.dat",
      dir+"/position.dat",
      baryonfree,
      deltat, deltax, deltay, deltah);
    break;
  }

  std::cout << "Average initial particle number from hydrojet:"
            << " before decay= " << jam->getIniAverageParticleNumber1()
            << " after decay= " << jam->getIniAverageParticleNumber2()
            << std::endl;

  delete jam;
#endif
}

//-----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  int const ext = args.read(argc, argv);
  if (ext) return ext;

  // Set Random number generator.
  Random* rand = new PyRand(randomSeed);
  // Random* rand = new Random(randomSeed);
  Random::setRandom(rand);

  // If you do not have output "freezeout.dat", "position.dat"
  if (args.makeElem) getHydroSpectra();

  hadronicCascade(args.nev);
  delete rand;

  return 0;
}

void generatePhasespaceData20141020(int ibase, std::string dirJAM) {
  static const int NEvent = 1000;

  Hydro2JamInitParams iparam;
  iparam.mevent              = NEvent; // not used
  iparam.seed                = randomSeed;
  iparam.dir_reso            = dir;
  iparam.kintmp              = freezeoutTemp;
  iparam.eos_pce             = eos_pce;
  iparam.dir                 = dirJAM;
  iparam.dpd                 = dumpPhaseSpace;
  iparam.fnamePS             = fnamePS;
  iparam.fnamePS0            = fnamePS0;
  iparam.outputHistogramFlag = outputHistogramFlag;

  {
    ResonanceListPCE reso(iparam.kintmp, iparam.eos_pce, resodata);

    ParticleSampleViscous* psamp = new ParticleSampleViscous(&reso, args.fnameInitialPhasespaceData);
    psamp->setOverSamplingFactor(NEvent);
    if (kashiwa::getenvAsBool("hydro2jam_turnsOffViscousEffect", false))
      psamp->setTurnsOffViscousEffect(true);
    psamp->update();

    // 振り分け
    std::vector<HydroParticleCF*> const& plist = psamp->getParticleList();
    std::vector<HydroParticleCF*> phases[NEvent];
    for (std::vector<HydroParticleCF*>::const_iterator i = plist.begin(); i != plist.end(); ++i)
      phases[std::min(int(Random::getRand() * NEvent), NEvent - 1)].push_back(*i);

    // 保存
    for (int i = 0; i < NEvent; i++) {
      char fn[200];
      std::sprintf(fn, "%s/dens%06d_phasespace0.dat", dirJAM.c_str(), ibase + i);
      savePhasespaceData(fn, phases[i], psamp->getParticleIdType());
    }

    delete psamp;
  }
}
