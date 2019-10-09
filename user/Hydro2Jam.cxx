//...A main program to use the initial condition of hadronic cascade
//...from hydro
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef USE_JAM
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "uty/Random.h"
#include "Hydro2Jam.h"
//#include "HydroSpec.h"
//#include "jam/JAMParticle.h"

#include "ksh/util.h"

using namespace std;

static std::string elementOutputFilenames[151] = {
	"ELEMENT.A0.PC170",
	"ELEMENT.DELTA.PC170",
	"ELEMENT.DELTABAR.PC170",
	"ELEMENT.ETA.PC170",
	"ELEMENT.ETAP.PC170",
	"ELEMENT.F0.PC170",
	"ELEMENT.KBAR.PC170",//previously "ELEMENT.K0S.PC170"
	"ELEMENT.KSTAR.PC170",
	"ELEMENT.KSTARBAR.PC170",
	"ELEMENT.LAMBDA.PC170",
	"ELEMENT.LAMBDABAR.PC170",
	"ELEMENT.OMEGA.PC170",
	"ELEMENT.PHI.PC170",
	"ELEMENT.RHO.PC170",
	"ELEMENT.SIGMA.PC170",
	"ELEMENT.SIGMAB.PC170",
	"ELEMENT.SIGMABBAR.PC170",
	"ELEMENT.PI.PC170",
	"ELEMENT.K.PC170",
	"ELEMENT.PRO.PC170",
	"ELEMENT.PBAR.PC170",
	"ELEMENT.A1_1260.PC170",
	"ELEMENT.A2_1320.PC170",
	"ELEMENT.B1_1235.PC170",
	"ELEMENT.PI_1300.PC170",
	"ELEMENT.PI2_1670.PC170",
	"ELEMENT.RHO_1465.PC170",
	"ELEMENT.F2_1270.PC170",
	"ELEMENT.F2P_1525.PC170",
	"ELEMENT.H1_1170.PC170",
	"ELEMENT.F0P.PC170",
	"ELEMENT.H1P.PC170",
	"ELEMENT.ETA_1295.PC170",
	"ELEMENT.F1_1285.PC170",
	"ELEMENT.F1_1420.PC170",
	"ELEMENT.F0_1300.PC170",
	"ELEMENT.OMEGA_1420.PC170",
	"ELEMENT.F1_1510.PC170",
	"ELEMENT.OMEGA_1600.PC170",
	"ELEMENT.K1_1270.PC170",
	"ELEMENT.K1BAR_1270.PC170",
	"ELEMENT.K2S_1430.PC170",
	"ELEMENT.K2SBAR_1430.PC170",
	"ELEMENT.K0S_1430.PC170",
	"ELEMENT.K0SBAR_1430.PC170",
	"ELEMENT.K1_1400.PC170",
	"ELEMENT.K1BAR_1400.PC170",
	"ELEMENT.KSTAR_1410.PC170",
	"ELEMENT.KSTARBAR_1410.PC170",
	"ELEMENT.N_1440.PC170",
	"ELEMENT.NBAR_1440.PC170",
	"ELEMENT.N_1520.PC170",
	"ELEMENT.NBAR_1520.PC170",
	"ELEMENT.N_1535.PC170",
	"ELEMENT.NBAR_1535.PC170",
	"ELEMENT.N_1650.PC170",
	"ELEMENT.NBAR_1650.PC170",
	"ELEMENT.DELTA_1600.PC170",
	"ELEMENT.DELTABAR_1600.PC170",
	"ELEMENT.DELTA_1620.PC170",
	"ELEMENT.DELTABAR_1620.PC170",
	"ELEMENT.LAMBDA_1405.PC170",
	"ELEMENT.LAMBDABAR_1405.PC170",
	"ELEMENT.LAMBDA_1520.PC170",
	"ELEMENT.LAMBDABAR_1520.PC170",
	"ELEMENT.LAMBDA_1600.PC170",
	"ELEMENT.LAMBDABAR_1600.PC170",
	"ELEMENT.LAMBDA_1670.PC170",
	"ELEMENT.LAMBDABAR_1670.PC170",
	"ELEMENT.SIGMAB_1385.PC170",
	"ELEMENT.SIGMABBAR_1385.PC170",
	"ELEMENT.SIGMAB_1660.PC170",
	"ELEMENT.SIGMABBAR_1660.PC170",
	"ELEMENT.SIGMAB_1670.PC170",
	"ELEMENT.SIGMABBAR_1670.PC170",
	"ELEMENT.XI.PC170",
	"ELEMENT.XIBAR.PC170",
	"ELEMENT.XI_1530.PC170",
	"ELEMENT.XIBAR_1530.PC170",
	"ELEMENT.OMEGAB.PC170",
	"ELEMENT.OMEGABBAR.PC170",
	"ELEMENT.PHI_1680.PC170",
	"ELEMENT.RHO_1700.PC170",
	"ELEMENT.KSTAR_1680.PC170",
	"ELEMENT.KSTARBAR_1680.PC170",
	"ELEMENT.K_3_1780.PC170",
	"ELEMENT.K_3BAR_1780.PC170",
	"ELEMENT.K_2_1770.PC170",
	"ELEMENT.K_2BAR_1770.PC170",
	"ELEMENT.K_2_1820.PC170",
	"ELEMENT.K_2BAR_1820.PC170",
	"ELEMENT.N_1675.PC170",
	"ELEMENT.NBAR_1675.PC170",
	"ELEMENT.N_1680.PC170",
	"ELEMENT.NBAR_1680.PC170",
	"ELEMENT.N_1700.PC170",
	"ELEMENT.NBAR_1700.PC170",
	"ELEMENT.N_1710.PC170",
	"ELEMENT.NBAR_1710.PC170",
	"ELEMENT.N_1720.PC170",
	"ELEMENT.NBAR_1720.PC170",
	"ELEMENT.N_1990.PC170",
	"ELEMENT.NBAR_1990.PC170",
	"ELEMENT.DELTA_1700.PC170",
	"ELEMENT.DELTABAR_1700.PC170",
	"ELEMENT.DELTA_1900.PC170",
	"ELEMENT.DELTABAR_1900.PC170",
	"ELEMENT.DELTA_1905.PC170",
	"ELEMENT.DELTABAR_1905.PC170",
	"ELEMENT.DELTA_1910.PC170",
	"ELEMENT.DELTABAR_1910.PC170",
	"ELEMENT.DELTA_1920.PC170",
	"ELEMENT.DELTABAR_1920.PC170",
	"ELEMENT.DELTA_1930.PC170",
	"ELEMENT.DELTABAR_1930.PC170",
	"ELEMENT.DELTA_1950.PC170",
	"ELEMENT.DELTABAR_1950.PC170",
	"ELEMENT.LAMBDA_1690.PC170",
	"ELEMENT.LAMBDABAR_1690.PC170",
	"ELEMENT.LAMBDA_1800.PC170",
	"ELEMENT.LAMBDABAR_1800.PC170",
	"ELEMENT.LAMBDA_1810.PC170",
	"ELEMENT.LAMBDABAR_1810.PC170",
	"ELEMENT.LAMBDA_1820.PC170",
	"ELEMENT.LAMBDABAR_1820.PC170",
	"ELEMENT.LAMBDA_1830.PC170",
	"ELEMENT.LAMBDABAR_1830.PC170",
	"ELEMENT.LAMBDA_1890.PC170",
	"ELEMENT.LAMBDABAR_1890.PC170",
	"ELEMENT.LAMBDA_2100.PC170",
	"ELEMENT.LAMBDABAR_2100.PC170",
	"ELEMENT.LAMBDA_2110.PC170",
	"ELEMENT.LAMBDABAR_2110.PC170",
	"ELEMENT.SIGMAB_1750.PC170",
	"ELEMENT.SIGMABBAR_1750.PC170",
	"ELEMENT.SIGMAB_1775.PC170",
	"ELEMENT.SIGMABBAR_1775.PC170",
	"ELEMENT.SIGMAB_1915.PC170",
	"ELEMENT.SIGMABBAR_1915.PC170",
	"ELEMENT.SIGMAB_1940.PC170",
	"ELEMENT.SIGMABBAR_1940.PC170",
	"ELEMENT.SIGMAB_2030.PC170",
	"ELEMENT.SIGMABBAR_2030.PC170",
	"ELEMENT.XI_1690.PC170",
	"ELEMENT.XIBAR_1690.PC170",
	"ELEMENT.XI_1820.PC170",
	"ELEMENT.XIBAR_1820.PC170",
	"ELEMENT.XI_1950.PC170",
	"ELEMENT.XIBAR_1950.PC170",
	"ELEMENT.XI_2030.PC170",
	"ELEMENT.XIBAR_2030.PC170",
};

//=========Set input values and switches ========================
Hydro2Jam::Hydro2Jam(Hydro2JamInitParams const& iparam)
{
  this->initialize(iparam);
}

Hydro2Jam::Hydro2Jam(int mevent, int seed, string dir_reso,int kintmp,
	int eos_pce, string dir, int dpd, string fnamePS, string fnamePS0)
{
  Hydro2JamInitParams iparam;
  iparam.mevent  =mevent;
  iparam.seed    =seed;
  iparam.dir_reso=dir_reso;
  iparam.kintmp  =kintmp;
  iparam.eos_pce =eos_pce;
  iparam.dir     =dir;
  iparam.dpd     =dpd;
  iparam.fnamePS =fnamePS;
  iparam.fnamePS0=fnamePS0;
  this->initialize(iparam);
}

void Hydro2Jam::initialize(Hydro2JamInitParams const& iparam)
{
  dirReso           =iparam.dir_reso;
  kinTmp            =iparam.kintmp;
  eosPCE            =iparam.eos_pce;
  nEvent            =iparam.mevent;
  dumpPhaseSpaceData=iparam.dpd;

  jam=new Jam1();

  //....Initialize JAM
  jam->setMSTC(1,iparam.seed);      // random seed.
  jam->setMSTC(2,iparam.mevent);    // number of event.
  //jam->setMSTC(38,6);      // io number for jamlist.
  jam->setMSTC(8,0);         // job mode.
  jam->setMSTC(16,0);        // display on/off.
  jam->setPARC(6,5.0);       // scale of display
  jam->setMSTC(54,0);        //avoid first coll inside the same nucleus off
  //jam->setMSTC(39,0);      //no output fname(4)

  //....Switch on some analysis.
  jam->setMSTC(156,0);        // analysis of collision distribution
  jam->setMSTC(161,0);        // no analysis from jam internal subr.
  jam->setMSTC(162,0);        // Output collision histroy
  jam->setMSTC(165,0);        //
  jam->setPARC(7, 1.0);     // Output time interval (fm/c)
  jam->setMSTC(81,0);         // 1:hard scattering on/off
  jam->setMSTC(4,100);       // user defined frame.

  //jam->setMSTC(61,0);       // isotropic resonance decay option

  //....Initial setting for JAM.
  string frame="user";  // comp. frame in this case, user defined
  double dt=100.0;          // collision time(fm/c)
  int nstep=1;            // time step (i.e. no time step)

  //...dummy in this case.
  double bmin=0.0;         // minimum impact parameter (dummy)
  double bmax=0.0;         // maximum impact parameter (dummy)

  jam->setPARD(16,10.0);       // user defined frame.

  std::string fn2 = jam->getFNAME(2);
  std::string fn3 = jam->getFNAME(3);
  std::string fn4 = jam->getFNAME(4);
  if(iparam.dir.size()>0) {
    std::string dir = iparam.dir + "/";
    fn2 = dir + fn2;
    fn3 = dir + fn3;
    fn4 = dir + fn4;
    /*
      if(dir.size() > 70) {
	    cerr << "(Hydro2Jam:)Sorry file name is too long, change";
	    cerr << " the array of fname in jam2.inc" << endl;
	    cerr << dir.size() << endl;
	    exit(1);
      }
    */
    jam->setFNAME(2,fn2.c_str());
    jam->setFNAME(3,fn3.c_str());
    jam->setFNAME(4,fn4.c_str());
    jam->setFNAME(8,dir.c_str());
  }

  //...Initialize jam.
  jam->jamInit(iparam.mevent,bmin,bmax,dt,nstep, "user", "p ", "p ", "2gev");

  isFile=0;
  numberTestParticle=1;

  if(dumpPhaseSpaceData) {
    if(iparam.dir.size()>0){
      ofs.open((iparam.dir + "/" + iparam.fnamePS).c_str());
      ofs0.open((iparam.dir + "/" + iparam.fnamePS0).c_str());
    }else{
      ofs.open(iparam.fnamePS.c_str());
      ofs0.open(iparam.fnamePS0.c_str());
    }
  }

  //jam->setMDCY(jam->jamComp(111) ,1,0);   // no pi0 decay
  //jam->setMDCY(jam->jamComp(3122),1,1);   // Lambda decay
  //jam->setMDCY(jam->jamComp(3222),1,1);   // Sigma- decay
  //jam->setMDCY(jam->jamComp(3212),1,1);   // Sigma0 decay
  //jam->setMDCY(jam->jamComp(3112),1,1);   // Sigma+ decay

  if(!kashiwa::getenvAsBool("hydro2jam_phi_decays",true))
    jam->setMDCY(jam->jamComp(333),1,0); // no phi decay

  numOutputHist=10;
}

Hydro2Jam::~Hydro2Jam()
{
  if(dumpPhaseSpaceData) {
    ofs << -999 << endl;
    ofs.close();
    ofs0 << -999 << endl;
    ofs0.close();
  }
  delete jam;
}

void Hydro2Jam::generateEvent(IParticleSample* psamp){
  //int nevent=jam->getMSTC(2);
  int const nEvent=this->nEvent;
  aveNumberPart1=0.0;
  aveNumberPart2=0.0;
  jam->setMSTC(5,numberTestParticle);

  if(nEvent*numberTestParticle>0)
    psamp->setAdviceNumberOfExpectedEvents(nEvent*numberTestParticle);

  int nprint=1;

  //...Simulation start.
  for(int iev=1;iev<= nEvent;iev++) {

    // if(iev%nprint == 0)
    //   std::cout << "Hydro2Jam.cxx(Hydro2Jam::generateEvent): event=" << iev << std::endl;

    //...Set initial particle momentum and coordinate in JAM.
    nv=0;
    nbary=0;
    nmeson=0;

    // Sample particles from hydro output.
    for(int i=0;i<numberTestParticle;i++){
      psamp->update();
      if(psamp->getIsOutput())
        initJam(psamp->getFileNamePos());
      else
        initJam(psamp);
    }

    if(iev%nprint==0){
      std::cout << "Hydro2Jam.cxx(Hydro2Jam::generateEvent):"<< "iev="<<iev<<": "
                << "sampling done. The number of initial test particles is nv=" << nv << "." << std::endl;
    }

    aveNumberPart1 += (double)nv/(nEvent*numberTestParticle);

    //...C.M.correction.
    //cmCorrection(nv);

    if(dumpPhaseSpaceData)
      printPhaseSpaceData(ofs0); // output the distribution to phasespace0.dat

    if(kashiwa::getenvAsBool("hydro2jam_decay_only",false)){
      // perform resonance decay (no rescatterings will be performed.)
      jam->finalResonanceDecay();

      std::cout <<"  env(hydro2jam_decay_only): resonance decay is performed without rescattering." << std::endl;

    }else{
      // Simulate one hadronic cascade event.
      jam->jamEvt(iev);

      if(iev%nprint==0){
        std::cout << "Hydro2Jam.cxx(Hydro2Jam::generateEvent):"<< "iev="<<iev<<": "
                  << "cascade done. The expected number of collisions is "
                  << (jam->getMSTD(41)+jam->getMSTD(42))/numberTestParticle << "." << std::endl;
      }
    }

    if(dumpPhaseSpaceData)
      printPhaseSpaceData(ofs); // output the distribution to phasespace0.dat

  }  // end event simulation loop.

  //...Final output.
  jam->jamFin();
}

void Hydro2Jam::generateEventFromHypersurfaceFiles(string fn_freezeout_dat, string fn_position_dat, int baryonfree, double deltat, double deltax, double deltay, double deltah)
{
  // Now initilize the sampling of the particles from hydro simulation.
  ParticleSampleHydrojet* psamp =
    new ParticleSampleHydrojet(dirReso, elementOutputFilenames, kinTmp, eosPCE, resodata);
  psamp->setDtau(deltat);
  psamp->setDx(deltax);
  psamp->setDy(deltay);
  psamp->setDh(deltah);
  psamp->setBaryonFree(baryonfree);

  if(isFile)
    psamp->setIsOutput(1);
  else
    psamp->setIsOutput(0);

  psamp->setHypersurfaceFilenames(fn_freezeout_dat,fn_position_dat);

  this->generateEvent(psamp);
  delete psamp;
}

void Hydro2Jam::printPhaseSpaceData(ofstream& output)
{
  int nv = jam->getNV();
  output << nv << "  " << numberTestParticle << endl;
  for(int i=1;i<=nv;i++) {
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
    output << setw(5) << ks
           << setw(10) << kf
           << setw(14) << px
           << setw(14) << py
           << setw(14) << pz
           << setw(14) << m
           << setw(14) << x
           << setw(14) << y
           << setw(14) << z
           << setw(14) << t
           << endl;
  }
}

//#define DBG20140822_HYDRO2JAM
#ifdef DBG20140822_HYDRO2JAM
struct dbg_counter{
  int value;
  std::string name;
  dbg_counter(std::string const& name):value(0),name(name){}
  void operator()(){value++;}
  ~dbg_counter(){
    std::cerr<<"count:"<<name<<" = "<<value<<";"<<std::endl;
  }
};
#endif

void Hydro2Jam::initJam(IParticleSample* psamp)
{
  std::vector<Particle*> const& plist = psamp->getParticleList();

  //std::cerr<<"dbg20140822: "<<psamp->getParticleIdType()<<" (HydroParticleID = "<<ParticleIDType::HydroParticleID<<")"<<std::endl;

  std::vector<Particle*>::const_iterator mp;
  for(mp=plist.begin(); mp != plist.end(); mp++) {
    Particle const* const particle=*mp;

    int kf; // PDG particle code.
    switch(psamp->getParticleIdType()){
    case ParticleIDType::HydroParticleID:
      {
        int ir = particle->id;
        kf = getJamID(ir+1);
      }
      break;
    case ParticleIDType::PDGCode:
      kf=particle->id;
      break;
    default:
      std::cerr<<"Hydro2Jam::initJam: invalid value of psamp->getParticleIdType()."<<std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(kf == 0) continue;

    int kc=jam->jamComp(kf);      // internal particle code.
    //int ibary=isign(kchg(kc,6),kf)  // baryon number
    int ibary=jam->getBaryonNumber(kc,kf);  // baryon number
    if(ibary == 0)
	    nmeson++;
    else
	    nbary++;

    nv++;

    //...Zero the vector.
    jam->jamZero(nv);

    int ks=2;
    if(jam->getPMAS(kc,2) <= 1e-7 || jam->getMDCY(kc,1) == 0
       || jam->getMDCY(kc,2) == 0 || jam->getMDCY(kc,3) == 0) ks=1;
    jam->setK(1,nv, ks);
    jam->setK(2,nv, kf);
    jam->setK(3,nv, 0);
    jam->setK(4,nv, 0);
    jam->setK(5,nv, -1);
    jam->setK(6,nv, 0);
    jam->setK(7,nv, 1);
    jam->setK(8,nv, 1);
    jam->setK(9,nv, ibary);
    jam->setK(10,nv,0);
    jam->setK(11,nv,0);

    double const x  = particle->x;
    double const y  = particle->y;
    double const z  = particle->z;
    double const t  = particle->t;
    double const px = particle->px;
    double const py = particle->py;
    double const pz = particle->pz;
    double pe = particle->e;
    double pm;
    if (pe<0.0) {
      pm = jam->jamMass(kf);
      pe = sqrt(px*px+py*py+pz*pz+pm*pm);
    }else{
      pm = sqrt(pe*pe-(px*px+py*py+pz*pz));
    }

    jam->setP(1,nv, px);
    jam->setP(2,nv, py);
    jam->setP(3,nv, pz);
    jam->setP(4,nv, pe);
    jam->setP(5,nv, pm);

    jam->setR(1,nv, x);
    jam->setR(2,nv, y);
    jam->setR(3,nv, z);
    jam->setR(4,nv, t);
    jam->setR(5,nv, t);         // formation time

    //...Vertex
    jam->setV(1,nv, x);
    jam->setV(2,nv, y);
    jam->setV(3,nv, z);
    jam->setV(4,nv, t);

    //.....Set resonance decay time.
    jam->setV(5,nv, 1.e+35);
    if(jam->getK(1,nv) == 2)
      jam->setV(5,nv,t+jam->jamDecayTime(1,kf,kc,ks,pm,pe));

    //jam->print(nv);

  }

  jam->setNV( nv );          // set total number of particles.
  jam->setNBARY( nbary );    // set total number of baryon.
  jam->setNMESON( nmeson );  // set total number of mesons.
}

// void LoadHydroData(std::vector<Particle*>& plist,std::string const& fname){
//   std::ifstream fdata(fname.c_str());
//   if(!fdata){
//     std::cerr<<"Hydro2Jam.cxx (LoadHydroData): failed to open the file '"<<fname<<"'"<<std::endl;
//     std::exit(EXIT_FAILURE);
//   }

//   double px,py,pz,e,em,tau,rx,ry,eta;
//   int ir;
//   std::string line;
//   int lp=1;
//   while(std::getline(fdata,line)){
//     ip++;
//     if(line.find('#')>=0)continue;

//     std::istringstream is(line);
//     if(!(is >> px >> py >> pz >> e >> em >> ir >> tau >> rx >> ry >> eta)) {
//       std::cerr<<fname<<":"<<ip<<": Invalid data at line "<<lp<<std::endl;
//       std::exit(EXIT_FAILURE);
//     }

//     // type==ParticleIDType::HydroParticleID
//     Particle* particle=new Particle;
//     particle->setID(ir);

//     particle->setPx(px);
//     particle->setPy(py);
//     particle->setPz(pz);
//     particle->setPe(-1.0); // jam->jamMass で自動決定

//     particle->setX(rx);
//     particle->setY(ry);
//     particle->setZ(tau*sinh(eta));
//     particle->setT(tau*cosh(eta));

//     plist.push_back(particle);
//   }
// }

void Hydro2Jam::initJam(string fname)
{
  ifstream fdata;
  if(fdata.is_open()) {
    cerr << "funny! (Hydro2Jam:) file areadly opend"
         << fname << endl;
    exit(1);
  }

  fdata.open(fname.c_str(),ios::in);
  if(!fdata) {
    cerr << "Hydro2Jam: Error: unable to open file "
         << fname  << endl;
    exit(1);
  }
  double px,py,pz,e,em,tau,rx,ry,eta;
  int ir;
  string line;
  int lp=1;
  while(getline(fdata,line)) {
    int com = line.find('#');
    if(com >=0) continue;
    istringstream is(line);
    if(!(is >> px >> py >> pz >> e >> em >> ir >> tau >> rx >> ry >> eta)) {
	    cerr << "Invalid data at line " << lp << endl;
	    exit(1);
    }
    lp++;

    /*
      cout << "px= " << px
	    << " py=" << py
	    << " pz= "<< pz
	    << " e= " << e
	    << " em= " << em
	    << " ir= " << ir
	    << " tau= " << tau
	    << " rx= " << rx
	    << " ry= " << ry
	    << " eta= " << eta
	    << endl;
      cin.get();
    */

    int kf=getJamID(ir+1);   // PDG particle code.
    if(kf == 0) continue;
    int kc=jam->jamComp(kf);      // internal particle code.
    //int ibary=isign(kchg(kc,6),kf)  // baryon number
    int ibary=jam->getBaryonNumber(kc,kf);  // baryon number
    if(ibary == 0)
	    nmeson++;
    else
	    nbary++;

    nv++;

    //...Zero the vector.
    jam->jamZero(nv);

    //....Particle mass.
    double pm=jam->jamMass(kf);

    if(jam->getPMAS(kc,2) <= 1e-7 || jam->getMDCY(kc,1) == 0
       || jam->getMDCY(kc,2) == 0 || jam->getMDCY(kc,3) == 0)
      jam->setK(1,nv,1);
    else
      jam->setK(1,nv,2);

    jam->setK(2,nv, kf);
    jam->setK(3,nv, 0);
    jam->setK(4,nv, 0);
    jam->setK(5,nv, -1);
    jam->setK(6,nv, 0);
    jam->setK(7,nv, 1);
    jam->setK(8,nv, 1);
    jam->setK(9,nv, ibary);
    jam->setK(10,nv, 0);
    jam->setK(11,nv, 0);
    jam->setP(1,nv, px);
    jam->setP(2,nv, py);
    jam->setP(3,nv, pz);
    jam->setP(4,nv, sqrt(px*px+py*py+pz*pz+pm*pm));
    jam->setP(5,nv, pm);

    jam->setR(1,nv, rx);
    jam->setR(2,nv, ry);
    jam->setR(3,nv, tau*sinh(eta));
    jam->setR(4,nv, tau*cosh(eta));
    jam->setR(5,nv, jam->getR(4,nv)  );         // formation time

    //...Vertex
    jam->setV(1,nv, jam->getR(1,nv));
    jam->setV(2,nv, jam->getR(2,nv));
    jam->setV(3,nv, jam->getR(3,nv));
    jam->setV(4,nv, jam->getR(4,nv));

    //.....Set resonance decay time.
    jam->setV(5,nv, 1.e+35);
    if(jam->getK(1,nv) == 2)
      jam->setV(5,nv, jam->getR(4,nv)
                + jam->jamDecayTime(1,kf,kc,
                                    jam->getK(1,nv),jam->getP(5,nv),jam->getP(4,nv)));

    //jam->print(nv);

  }

  fdata.close();

  jam->setNV( nv );          // set total number of particles.
  jam->setNBARY( nbary );    // set total number of baryon.
  jam->setNMESON( nmeson );  // set total number of mesons.

}

void Hydro2Jam::cmCorrection()
{
	double cx=0.0;
	double cy=0.0;
	double cz=0.0;
	double px=0.0;
	double py=0.0;
	double pz=0.0;
	double s=0.0;
	for(int i=1;i<=nv;i++) {
    px += jam->getP(1,i);
    py += jam->getP(2,i);
    pz += jam->getP(3,i);
    cx += jam->getR(1,i)*jam->getP(5,i);
    cy += jam->getR(2,i)*jam->getP(5,i);
    cz += jam->getR(3,i)*jam->getP(5,i);
    s  += jam->getP(5,i);
	}
	cx = -cx/s;
	cy = -cy/s;
	cz = -cz/s;
	px = -px/nv;
	py = -py/nv;
	pz = -pz/nv;

	//cout << "cx= " << cx << " cy= " << cy << " cz= " << cz << endl;
	//cout << "px= " << px << " py= " << py << " pz= " << pz << endl;
	//cin.get();

	for(int i=1;i<=nv;i++) {
    jam->setR(1,i, jam->getR(1,i)+cx);
    jam->setR(2,i, jam->getR(2,i)+cy);
    jam->setR(3,i, jam->getR(3,i)+cz);
    jam->setV(1,i, jam->getR(1,i));
    jam->setV(2,i, jam->getR(2,i));
    jam->setV(3,i, jam->getR(3,i));
    double p1 = jam->getP(1,i)+px;
    double p2 = jam->getP(2,i)+py;
    double p3 = jam->getP(3,i)+pz;
    double m  = jam->getP(5,i);
    jam->setP(1,i, p1);
    jam->setP(2,i, p2);
    jam->setP(3,i, p3);
    double e=sqrt(m*m + p1*p1 + p2*p2 + p3*p3);
    jam->setP(4,i,e);
	}
}

#if 0
// int Hydro2Jam::getJamID(int irshift); alternative implementation, not yet tested

class HydroParticleCodeTable{
  static const int nresonance=151;
  std::vector<int> table[nresonance];

  void registerMapping(int hydroParticleCode,int jamCode1,int jamCode2=0,int jamCode3=0,int jamCode4=0){
    if(1<=hydroParticleCode&&hydroParticleCode<=nresonance){
      std::vector<int>& list=table[hydroParticleCode-1];
      list.clear();
      list.push_back(jamCode1);
      if(jamCode2)list.push_back(jamCode2);
      if(jamCode3)list.push_back(jamCode3);
      if(jamCode4)list.push_back(jamCode4);
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
    if(1<=hydroParticleCode&&hydroParticleCode<=nresonance){
      std::vector<int>& list=table[hydroParticleCode-1];
      if(list.size()==1)
        return list[0];
      else if(list.size()>0){
        int index=int(Random::getRand()*list.size());
        if(index==list.size())index--;
        return list[index];
      }
    }

    std::cerr << "Hydro2Jam.cxx: unexpected hydroParticleCode=" << hydroParticleCode << std::endl;
    return 0;
  }
};

int Hydro2Jam::getJamID(int irshift){
  static HydroParticleCodeTable table;
  return table.generateJamID(irshift);
}
#endif

//...Convert hydro particle code to JAM code. Isospin is also generated.
/// @param irshift = ir + 1
/// @return JAM code
int Hydro2Jam::getJamID(int irshift)
{
  double rand;
  switch(irshift) {
  case 1:  // a0
    rand = Random::getRand();
    if(rand < 0.333333)
      return 10111;
    else if(rand < 0.666666)
      return 10211;
    else
      return -10211;
  case 2:    // Delta
    rand = Random::getRand();
    if(rand <  0.25)
      return 1114;
    else if(rand < 0.5)
      return 2114;
    else if(rand < 0.75)
      return 2214;
    else
      return 2224;
  case 3:    // Delta-bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -1114;
    else if(rand < 0.5)
      return -2114;
    else if(rand < 0.75)
      return -2214;
    else
      return -2224;
  case 4:   // eta
    return 221;
  case 5:   // eta'
    return 331;
  case 6:   // f_0(980)
    return 10221;
  case 7:   // Kbar
    if(Random::getRand() < 0.5) return -311;
    return -321;
  case 8:   // K*(892)
    rand = Random::getRand();
    if(rand <  0.5)
      return 313;
    else
      return 323;
  case 9:   // K*(892)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -313;
    else
      return -323;
  case 10:   // Lambda
    return 3122;
  case 11: // Lambda-bar
    return -3122;
  case 12:   // omega
    return 223;
  case 13:  // phi(1020)
    return 333;
  case 14:  // rho(770)
    rand = Random::getRand();
    if(rand <  0.333333)
      return 113;
    else if(rand < 0.666666)
      return 213;
    else
      return -213;
  case 15:  // f_0(600) i.e. sigma meson
    return 10220;
  case 16:  // Sigma
    rand = Random::getRand();
    if(rand < 0.333333)
      return 3112;
    else if(rand <0.666666)
      return 3212;
    else
      return 3222;
  case 17:  // Sigma-bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -3112;
    else if(rand <0.666666)
      return -3212;
    else
      return -3222;
  case 18:  // pion
    rand = Random::getRand();
    if(rand < 0.333333)
      return 111;
    else if(rand < 0.666666)
      return 211;
    else
      return -211;
  case 19: // K=(K_0,K+)
    if(Random::getRand() < 0.5) return 321;
    return 311;
  case 20: // nucleon
    if(Random::getRand() < 0.5) return 2112; // n
    return 2212; // proton
  case 21: // anti-nucleon
    if(Random::getRand() < 0.5) return -2112; // n
    return -2212; // proton
  case 22://a_1(1260)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 20113;
    else if(rand < 0.666666)
      return 20213;
    else
      return -20213;
  case 23://a_2(1320)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 115;
    else if(rand < 0.666666)
      return 215;
    else
      return -215;
  case 24://b_1(1235)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 10113;
    else if(rand < 0.666666)
      return 10213;
    else
      return -10213;
  case 25://pi(1300)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 20111;
    else if(rand < 0.666666)
      return 20211;
    else
      return -20211;
  case 26://pi_2(1670)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 10115;
    else if(rand < 0.666666)
      return 10215;
    else
      return -10215;
  case 27://rho(1465)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 30113;
    else if(rand < 0.666666)
      return 30213;
    else
      return -30213;
  case 28://f_2(1270)
    return 225;
  case 29://f_2'(1525)
    return 335;
  case 30://h_1(1170)
    return 10223;
  case 31://f'_0
    return 10331;
  case 32://h'_1
    return 10333;
  case 33://eta(1295)
    return 20221;
  case 34://f_1(1285)
    return 20223;
  case 35://f_1(1420)
    return 20333;
  case 36://f_0(1300)
    return 30221;
  case 37://omega(1420)
    return 30223;
  case 38://f_1(1510)
    return 50223;
  case 39://omega(1600)
    return 60223;
  case 40://K_1(1270)
    rand = Random::getRand();
    if(rand <  0.5)
      return 10313;
    else
      return 10323;
  case 41://K_1(1270)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -10313;
    else
      return -10323;
  case 42://K_2*(1430)
    rand = Random::getRand();
    if(rand <  0.5)
      return 315;
    else
      return 325;
  case 43://K_2*(1430)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -315;
    else
      return -325;
  case 44://K_0*(1430)
    rand = Random::getRand();
    if(rand <  0.5)
      return 10311;
    else
      return 10321;
  case 45://K_0*(1430)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -10311;
    else
      return -10321;
  case 46://K_1(1400)
    rand = Random::getRand();
    if(rand <  0.5)
      return 20313;
    else
      return 20323;
  case 47://K_1(1400)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -20313;
    else
      return -20323;
  case 48://K*(1410)
    rand = Random::getRand();
    if(rand <  0.5)
      return 30313;
    else
      return 30323;
  case 49://K*(1410)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -30313;
    else
      return -30323;
  case 50://N(1440)
    if(Random::getRand() < 0.5) return 12112;
    return 12212;
  case 51://N(1440) bar
    if(Random::getRand() < 0.5) return -12112;
    return -12212;
  case 52://N(1520)
    if(Random::getRand() < 0.5) return 1214;
    return 2124;
  case 53://N(1520) bar
    if(Random::getRand() < 0.5) return -1214;
    return -2124;
  case 54://N(1535)
    if(Random::getRand() < 0.5) return 22112;
    return 22212;
  case 55://N(1535) bar
    if(Random::getRand() < 0.5) return -22112;
    return -22212;
  case 56://N(1650)
    if(Random::getRand() < 0.5) return 32112;
    return 32212;
  case 57://N(1650) bar
    if(Random::getRand() < 0.5) return -32112;
    return -32212;
  case 58://Delta(1600)
    rand = Random::getRand();
    if(rand <  0.25)
      return 31114;
    else if(rand < 0.5)
      return 32114;
    else if(rand < 0.75)
      return 32214;
    else
      return 32224;
  case 59://Delta(1600) bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -31114;
    else if(rand < 0.5)
      return -32114;
    else if(rand < 0.75)
      return -32214;
    else
      return -32224;
  case 60://Delta(1620)
    rand = Random::getRand();
    if(rand <  0.25)
      return 1112;
    else if(rand < 0.5)
      return 1212;
    else if(rand < 0.75)
      return 2122;
    else
      return 2222;
  case 61://Delta(1620) bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -1112;
    else if(rand < 0.5)
      return -1212;
    else if(rand < 0.75)
      return -2122;
    else
      return -2222;
  case 62://Lambda(1405)
      return 13122;
  case 63://Lambda(1405) bar
      return -13122;
  case 64://Lambda(1520)
      return 3124;
  case 65://Lambda(1520) bar
      return -3124;
  case 66://Lambda(1600)
      return 23122;
  case 67://Lambda(1600) bar
      return -23122;
  case 68://Lambda(1670)
      return 33122;
  case 69://Lambda(1670) bar
      return -33122;
  case 70://Sigma(1385)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 3114;
    else if(rand < 0.666666)
      return 3214;
    else
      return 3224;
  case 71://Sigma(1385) bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -3114;
    else if(rand < 0.666666)
      return -3214;
    else
      return -3224;
  case 72://Sigma(1660)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 13112;
    else if(rand < 0.666666)
      return 13212;
    else
      return 13222;
  case 73://Sigma(1660) bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -13112;
    else if(rand < 0.666666)
      return -13212;
    else
      return -13222;
  case 74://Sigma(1670)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 13114;
    else if(rand < 0.666666)
      return 13214;
    else
      return 13224;
  case 75://Sigma(1670) bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -13114;
    else if(rand < 0.666666)
      return -13214;
    else
      return -13224;
  case 76://Xi
    if(Random::getRand() < 0.5) return 3322;
    return 3312;
  case 77://Xi bar
    if(Random::getRand() < 0.5) return -3322;
    return -3312;
  case 78://Xi(1530)
    if(Random::getRand() < 0.5) return 3314;
    return 3324;
  case 79://Xi(1530) bar
    if(Random::getRand() < 0.5) return -3314;
    return -3324;
  case 80://Omega
    return 3334;
  case 81://Omega bar
    return -3334;
  case 82://phi(1680)
    return 30333;
  case 83://rho(1700)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 40113;
    else if(rand < 0.666666)
      return 40213;
    else
      return -40213;
  case 84://K*(1680)
    rand = Random::getRand();
    if(rand <  0.5)
      return 40313;
    else
      return 40323;
  case 85://K*(1680)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -40313;
    else
      return -40323;
  case 86://K_3(1780)
    rand = Random::getRand();
    if(rand <  0.5)
      return 317;
    else
      return 327;
  case 87://K_3(1780)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -317;
    else
      return -327;
  case 88://K_2(1770)
    rand = Random::getRand();
    if(rand <  0.5)
      return 10315;
    else
      return 10325;
  case 89://K_2(1770)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -10315;
    else
      return -10325;
  case 90://K_2(1820)
    rand = Random::getRand();
    if(rand <  0.5)
      return 20315;
    else
      return 20325;
  case 91://K_2(1820)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -20315;
    else
      return -20325;
  case 92://N(1675)
    rand = Random::getRand();
    if(rand <  0.5)
      return 2116;
    else
      return 2216;
  case 93://N(1675)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -2116;
    else
      return -2216;
  case 94://N(1680)
    rand = Random::getRand();
    if(rand <  0.5)
      return 12116;
    else
      return 12216;
  case 95://N(1680)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -12116;
    else
      return -12216;
  case 96://N(1700)
    rand = Random::getRand();
    if(rand <  0.5)
      return 21214;
    else
      return 22124;
  case 97://N(1700)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -21214;
    else
      return -22124;
  case 98://N(1710)
    rand = Random::getRand();
    if(rand <  0.5)
      return 42112;
    else
      return 42212;
  case 99://N(1710)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -42112;
    else
      return -42212;
  case 100://N(1720)
    rand = Random::getRand();
    if(rand <  0.5)
      return 31214;
    else
      return 32124;
  case 101://N(1720)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -31214;
    else
      return -32124;
  case 102://N(1990)
    rand = Random::getRand();
    if(rand <  0.5)
      return 11218;
    else
      return 12128;
  case 103://N(1990)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -11218;
    else
      return -12128;
  case 104://D(1700)
    rand = Random::getRand();
    if(rand <  0.25)
      return 11114;
    else if(rand < 0.5)
      return 12114;
    else if(rand < 0.75)
      return 12214;
    else
      return 12224;
  case 105://D(1700)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -11114;
    else if(rand < 0.5)
      return -12114;
    else if(rand < 0.75)
      return -12214;
    else
      return -12224;
  case 106://D(1900)
    rand = Random::getRand();
    if(rand <  0.25)
      return 11112;
    else if(rand < 0.5)
      return 11212;
    else if(rand < 0.75)
      return 12122;
    else
      return 12222;
  case 107://D(1900)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -11112;
    else if(rand < 0.5)
      return -11212;
    else if(rand < 0.75)
      return -12122;
    else
      return -12222;
  case 108://D(1905)
    rand = Random::getRand();
    if(rand <  0.25)
      return 1116;
    else if(rand < 0.5)
      return 1216;
    else if(rand < 0.75)
      return 2126;
    else
      return 2226;
  case 109://D(1905)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -1116;
    else if(rand < 0.5)
      return -1216;
    else if(rand < 0.75)
      return -2126;
    else
      return -2226;
  case 110://D(1910)
    rand = Random::getRand();
    if(rand <  0.25)
      return 21112;
    else if(rand < 0.5)
      return 21212;
    else if(rand < 0.75)
      return 22122;
    else
      return 22222;
  case 111://D(1910)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -21112;
    else if(rand < 0.5)
      return -21212;
    else if(rand < 0.75)
      return -22122;
    else
      return -22222;
  case 112://D(1920)
    rand = Random::getRand();
    if(rand <  0.25)
      return 21114;
    else if(rand < 0.5)
      return 22114;
    else if(rand < 0.75)
      return 22214;
    else
      return 22224;
  case 113://D(1920)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -21114;
    else if(rand < 0.5)
      return -22114;
    else if(rand < 0.75)
      return -22214;
    else
      return -22224;
  case 114://D(1930)
    rand = Random::getRand();
    if(rand <  0.25)
      return 11116;
    else if(rand < 0.5)
      return 11216;
    else if(rand < 0.75)
      return 12126;
    else
      return 12226;
  case 115://D(1930)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -11116;
    else if(rand < 0.5)
      return -11216;
    else if(rand < 0.75)
      return -12126;
    else
      return -12226;
  case 116://D(1950)
    rand = Random::getRand();
    if(rand <  0.25)
      return 1118;
    else if(rand < 0.5)
      return 2118;
    else if(rand < 0.75)
      return 2218;
    else
      return 2228;
  case 117://D(1950)_bar
    rand = Random::getRand();
    if(rand <  0.25)
      return -1118;
    else if(rand < 0.5)
      return -2118;
    else if(rand < 0.75)
      return -2218;
    else
      return -2228;
  case 118://L(1690)0
      return 13124;
  case 119://L(1690)0_bar
      return -13124;
  case 120://L(1800)0
      return 43122;
  case 121://L(1800)0_bar
      return -43122;
  case 122://L(1810)0
      return 53122;
  case 123://L(1810)0_bar
      return -53122;
  case 124://L(1820)0
      return 3126;
  case 125://L(1820)0_bar
      return -3126;
  case 126://L(1830)0
      return 13126;
  case 127://L(1830)0_bar
      return -13126;
  case 128://L(1890)0
      return 23124;
  case 129://L(1890)0_bar
      return -23124;
  case 130://L(2100)0
      return 3128;
  case 131://L(2100)0_bar
      return -3128;
  case 132://L(2110)0
      return 23126;
  case 133://L(2110)0_bar
      return -23126;
  case 134://S(1750)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 23112;
    else if(rand < 0.666666)
      return 23212;
    else
      return 23222;
  case 135://S(1750)_bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -23112;
    else if(rand < 0.666666)
      return -23212;
    else
      return -23222;
  case 136://S(1775)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 3116;
    else if(rand < 0.666666)
      return 3216;
    else
      return 3226;
  case 137://S(1775)_bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -3116;
    else if(rand < 0.666666)
      return -3216;
    else
      return -3226;
  case 138://S(1915)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 13116;
    else if(rand < 0.666666)
      return 13216;
    else
      return 13226;
  case 139://S(1915)_bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -13116;
    else if(rand < 0.666666)
      return -13216;
    else
      return -13226;
  case 140://S(1940)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 23114;
    else if(rand < 0.666666)
      return 23214;
    else
      return 23224;
  case 141://S(1940)_bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -23114;
    else if(rand < 0.666666)
      return -23214;
    else
      return -23224;
  case 142://S(2030)
    rand = Random::getRand();
    if(rand < 0.333333)
      return 3118;
    else if(rand < 0.666666)
      return 3218;
    else
      return 3228;
  case 143://S(2030)_bar
    rand = Random::getRand();
    if(rand < 0.333333)
      return -3118;
    else if(rand < 0.666666)
      return -3218;
    else
      return -3228;
  case 144://X(1690)
    rand = Random::getRand();
    if(rand <  0.5)
      return 13312;
    else
      return 13322;
  case 145://X(1690)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -13312;
    else
      return -13322;
  case 146://X(1820)
    rand = Random::getRand();
    if(rand <  0.5)
      return 13314;
    else
      return 13324;
  case 147://X(1820)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -13314;
    else
      return -13324;
  case 148://X(1950)
    rand = Random::getRand();
    if(rand <  0.5)
      return 23312;
    else
      return 23322;
  case 149://X(1950)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -23312;
    else
      return -23322;
  case 150://X(2030)
    rand = Random::getRand();
    if(rand <  0.5)
      return 3316;
    else
      return 3326;
  case 151://X(2030)_bar
    rand = Random::getRand();
    if(rand <  0.5)
      return -3316;
    else
      return -3326;
  default:
    cerr << "funny ir+1=" << irshift << endl;
    return 0;

  }
}
#endif
