#include <cmath>
#include <cstdlib>
#include <new>
#include <algorithm>

#include "uty/Random.h"
#include <physicsbase/Const.h>
#include <ksh/phys/Minkowski.h>
#include <spectra/IntegratedCooperFrye.h>
#include "ParticleSampleHydrojet.h"

#define CHG20110804
//#define DBG20141010_TestViscousCorrection 202
//#define DBG20141010_TestViscousCorrection 302

using namespace std;

static const double FreezeoutSkipTemperature=0.01; // unit: [/fm]

ParticleSampleHydrojet::ParticleSampleHydrojet(std::string const& dir, std::string* outf,int kin, int eos_pce, std::string const& fname)
  : ElementReso(dir,outf,kin,eos_pce,fname)
{
	outfilepos = "particlesample_pos.dat";
	if(dir.size()>0)outfilepos = dir + "/"+ "particlesample_pos.dat";
	outfileneg = "particlesample_neg.dat";
	if(dir.size()>0)outfileneg = dir + "/"+ "particlesample_neg.dat";

	//	di = 10;
	di = 20; // iteration number in bisection method
	isOutput=1;
	plist.clear();

  mode_delayed_cooperfrye=false;

  // constants
	tmpf = 0.16/sctr*1000.0;
	mubf = 1.6/sctr*1000.0;
	meanf = 0.45/sctr*1000.0;

  // 2013/04/23, KM, reverse z axis
  {
    const char* env=std::getenv("ParticleSample_ReverseParticleList");
    this->fReverseParticleList=env&&std::atoi(env);
    if(this->fReverseParticleList)
      std::cout<<"ParticleSampleHydrojet: ReverseParticleList mode enabled!"<<std::endl;
  }

  // 2013/04/30, KM, shuffle the particle list
  {
    const char* env=std::getenv("ParticleSample_ShuffleParticleList");
    this->fShuffleParticleList=env&&std::atoi(env);
    if(this->fShuffleParticleList)
      std::cout<<"ParticleSampleHydrojet: ShuffleParticleList mode enabled!"<<std::endl;
  }
}

ParticleSampleHydrojet::~ParticleSampleHydrojet()
{
  if(plist.size()>0) {
    std::vector<Particle*>::iterator cp;
    for(cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }
}

void ParticleSampleHydrojet::initialize(std::string const& fn_freezeout_dat, std::string const& fn_p)
{
  nreso_loop = ResonanceListPCE::nreso;
  //    if(baryonfree)nreso_loop = 20;

  // Open files for input.
  openDataFile2(fn_freezeout_dat, fn_p);

#ifdef CHG20110804
  {
    // Try to read ELEMENT.* cache files
    //**************************
    std::cout<<"ParticleSampleHydrojet.cxx(ParticleSampleHydrojet::initialize): checking Cooper-Frye cache files (.POS/.NEG)... "<<std::flush;
    resDataPos = new ifstream [nreso_loop];
    for(int i=0;i<nreso_loop;i++) {
      string fnpos = elemFile[i]+".POS";
      resDataPos[i].open(fnpos.c_str(), ios::in);
      if (!resDataPos[i])  {
        mode_delayed_cooperfrye=true;
        break;
      }
    }

    resDataNeg = new ifstream [nreso_loop];
    for(int i=0;i<nreso_loop;i++) {
      string fnneg = elemFile[i]+".NEG";
      resDataNeg[i].open(fnneg.c_str(), ios::in);
      if (!resDataNeg[i])  {
        mode_delayed_cooperfrye=true;
        break;
      }
    }

    if(mode_delayed_cooperfrye){
      std::cout
        <<"no(incomplete).\n"
        <<"ParticleSampleHydrojet.cxx(ParticleSampleHydrojet::initialize): entering delayed Cooper-Frye evaluation mode."<<std::endl;

      // Close files if ELEMENT.* is not a complete set.
      if(resDataPos){
        for(int i=0;i<nreso_loop;i++){
          if(resDataPos[i].is_open())
            resDataPos[i].close();
        }
      }
      if(resDataNeg){
        for(int i=0;i<nreso_loop;i++){
          if(resDataNeg[i].is_open())
            resDataNeg[i].close();
        }
      }
    }else{
      std::cout<<"yes"<<std::endl;
    }
    //**************************
  }
#else
  resDataPos = new ifstream [nreso_loop];
  for(int i=0;i<nreso_loop;i++) {
    string fnpos = elemFile[i]+".POS";
    resDataPos[i].open(fnpos.c_str(), ios::in);
    if (!resDataPos[i])  {
      cerr << "ParticleSampleHydrojet::initialize! unable to open file " << fnpos << endl;
      mode_delayed_cooperfrye=true;
      break;
    }
  }

  resDataNeg = new ifstream [nreso_loop];
  for(int i=0;i<nreso_loop;i++) {
    string fnneg = elemFile[i]+".NEG";
    resDataNeg[i].open(fnneg.c_str(), ios::in);
    if (!resDataNeg[i])  {
      cerr << "ParticleSampleHydrojet::initialize! unable to open file " << fnneg << endl;
      mode_delayed_cooperfrye=true;
      break;
    }
  }
#endif

  // open files for output.
  if(isOutput!=0){
    outdatpos.open(outfilepos.c_str());
    if (!outdatpos)  {
      cerr << "ParticleSampleHydrojet::initialize! unable to open file " << outfilepos << endl;
      exit(1);
    }
    outdatpos << "# p_x(GeV) p_y(GeV) p_z(GeV) E(GeV) m(GeV)";
    outdatpos << " ir tau(fm) x(fm) y(fm) eta"  << endl;

    outdatneg.open(outfileneg.c_str());
    if (!outdatneg)  {
      cerr << "ParticleSampleHydrojet::initialize! unable to open file " << outfileneg << endl;
      exit(1);
    }
    outdatneg << "# p_x(GeV) p_y(GeV) p_z(GeV) E(GeV) m(GeV)";
    outdatneg << " ir tau(fm) x(fm) y(fm) eta"  << endl;
  }
}

//void ParticleSampleHydrojet::analyze(string fn_freezeout_dat, string fn_p, string fn_ecc)
void ParticleSampleHydrojet::analyze(string fn_freezeout_dat, string fn_p)
{
  nreso_loop = ResonanceListPCE::nreso;

  // Open files.
  initialize(fn_freezeout_dat, fn_p);

  if(plist.size()>0) {
    vector<Particle*>::iterator cp;
    for(cp = plist.begin(); cp != plist.end();cp++) delete *cp;
    plist.clear();
  }

  double ran;
  int nsamp=1;
  int ipos = 1;
  double numResPos;
  double numResNeg;

  while(!readData()){
    if(!baryonfree){
      for(int i=0;i<nreso_loop;i++){
        rlist[i].mu=0.0;
        if(rlist[i].bf==1) rlist[i].mu = mubf*sqrt(1.0-tf*tf/tmpf/tmpf)-meanf*nbf;
        if(rlist[i].anti) rlist[i].mu = -mubf*sqrt(1.0-tf*tf/tmpf/tmpf)+meanf*nbf;
        //if(rlist[i].anti) rlist[i].mu = -mubf*sqrt(1.0-tf*tf/tmpf/tmpf)-meanf*nbf;
        //if(rlist[i].bf==1) rlist[i].mu = mub;
        //if(rlist[i].anti) rlist[i].mu = -mub;
      }
    }

    readPData();

#ifdef CHG20110804
    if(tf == 0.0) {
      std::cerr << "ParticleSampleHydrojet.cxx(ParticleSampleHydrojet::analyze)! TF=ZERO" << std::endl;
      continue;
    }

    // 2014-07-30
    if(tf<FreezeoutSkipTemperature){
      //std::cerr << "ParticleSampleHydrojet.cxx(ParticleSampleHydrojet::analyze): skipped lowT surface." << std::endl;
      continue;
    }
#endif

#if DBG20141010_TestViscousCorrection==202
    // for jam202
    {
      double gamma = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
      kashiwa::phys::vector4 u(gamma,vx*gamma,vy*gamma,vz*gamma);
      kashiwa::phys::vector4 ds(ds0,-dsx,-dsy,-dsz);
      std::fprintf(stderr,"v2:ds2:uds %g %g %g\n",gamma*gamma*(vx*vx+vy*vy+vz*vz),ds*ds,u*ds);
    }
#endif

    // Loop over all particles.
    for(int ir=0;ir<nreso_loop;ir++) {
#ifdef CHG20110804
      if(!mode_delayed_cooperfrye){
        //**************************
        resDataPos[ir] >> numResPos;
        if(numResPos > 1.0)
          cout << "Funny fluid element! ir =" << ir
               << " " << numResPos << endl;

        resDataNeg[ir] >> numResNeg;
        if(numResNeg > 1.0)
          cout << "Funny fluid element! ir =" << ir << endl;
        //**************************
      }else{
        double npos=0.0;
        double nneg=0.0;

        double gamma = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
        kashiwa::phys::vector4 u(gamma,vx*gamma,vy*gamma,vz*gamma);
        kashiwa::phys::vector4 ds(ds0,-dsx,-dsy,-dsz);
        double beta = 1.0/tf;

        if(rlist[ir].bf==-1){
          hydrojet::spectra::IntegrateBosonCooperFrye(npos,nneg,u,ds,beta,rlist[ir].mass,rlist[ir].mu);
        }else{
          hydrojet::spectra::IntegrateFermionCooperFrye(npos,nneg,u,ds,beta,rlist[ir].mass,rlist[ir].mu);
        }
#if DBG20141010_TestViscousCorrection==302
        if(ir==18){
          std::fprintf(stderr,"jam30x:npos=%g\n",npos);
          std::fflush(stderr);
        }
#endif
        double n=(nneg+npos)*rlist[ir].degeff;
        numResPos=npos*rlist[ir].deg;
        numResNeg=nneg*rlist[ir].deg;
      }
#else
      resDataPos[ir] >> numResPos;
      if(numResPos > 1.0)
        cout << "Funny fluid element! ir =" << ir
             << " " << numResPos << endl;

      resDataNeg[ir] >> numResNeg;
      if(numResNeg > 1.0)
        cout << "Funny fluid element! ir =" << ir << endl;
#endif

      if(bulk == 1){
        ds0 = dss * cosh(hh);
        dsz = -dss * sinh(hh);
      }else{
        ds0 = dss * sinh(hh);
        dsz = -dss * cosh(hh);
      }

      //	    if(!bulk){
      for(int isamp=0;isamp<nsamp;isamp++){
        if((iw == 1)  ||  (iw == 5)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,-vy,yv,ds0,dsx,-dsy,dsz,ir,ipos,tau,xx,-yy,eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx,vy,-yv,ds0,-dsx,dsy,-dsz,ir,ipos,tau,-xx,yy,-eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx,-vy,-yv,ds0,-dsx,-dsy,-dsz,ir,ipos,tau,-xx,-yy,-eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,-vy,yv,ds0,dsx,-dsy,dsz,ir,ipos,tau,xx,-yy,eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx,vy,-yv,ds0,-dsx,dsy,-dsz,ir,ipos,tau,-xx,yy,-eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx,-vy,-yv,ds0,-dsx,-dsy,-dsz,ir,ipos,tau,-xx,-yy,-eta);

        } else if((iw == 3)  ||  (iw == 7)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(-vx,vy,-yv,ds0,-dsx,dsy,-dsz,ir,ipos,tau,-xx,yy,-eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(-vx,vy,-yv,ds0,-dsx,dsy,-dsz,ir,ipos,tau,-xx,yy,-eta);

        } else if((iw == 2)  ||  (iw == 6)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,-vy,yv,ds0,dsx,-dsy,dsz,ir,ipos,tau,xx,-yy,eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,-vy,yv,ds0,dsx,-dsy,dsz,ir,ipos,tau,xx,-yy,eta);

        } else if((iw == 4)  ||  (iw == 8)) {

          //out going
          ipos = 1;
          ran = Random::getRand();
          if(ran < numResPos)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);

          //in coming
          ipos = 0;
          ran = Random::getRand();
          if(ran < numResNeg)
            getSample(vx,vy,yv,ds0,dsx,dsy,dsz,ir,ipos,tau,xx,yy,eta);

        }

      }
      //	    }//if(bulk)
    }

  }

  // 2013/04/30, KM, shuffle the particle list
  if (this->fShuffleParticleList) {
#if __cplusplus >= 201703L
    std::shuffle(this->plist.begin(), this->plist.end(), RandomURGB());
#else
    std::random_shuffle(this->plist.begin(), this->plist.end());
#endif
  }

  finish();
}

void ParticleSampleHydrojet::finish()
{
  closeDataFile();
  closePDataFile();

#ifdef CHG20110804
  if(!mode_delayed_cooperfrye)
    for(int i=0;i<nreso_loop;i++) {
      //**************************
      resDataPos[i].close();
      resDataNeg[i].close();
      //**************************
    }
#else
  for(int i=0;i<nreso_loop;i++) {
    resDataPos[i].close();
    resDataNeg[i].close();
  }
#endif

  if(isOutput) {
    outdatpos.close();
    outdatneg.close();
  }

  //**************************
  delete [] resDataPos;
  delete [] resDataNeg;
  //**************************
}

void ParticleSampleHydrojet::
getSample(double vx,double vy,double yv,
	  double ds0,double dsx,double dsy,double dsz,int ir, int ipos,
	  double tau,double x0,double y0,double eta0)
{

  double  p[58],pw[58];
  double  p1[12],pw1[12];
  //    double  p1[38],pw1[38];

  double ptmid=1e3/sctr;
  double dx = getDx();
  double dy = getDy();
  double dh = getDh();
  double dtau = getDtau();
  double vz = tanh(yv);
  double gamma =  cosh(yv)/sqrt(1.-(vx*vx+vy*vy)*cosh(yv)*cosh(yv));
  double beta= 1./tf;
  double mres = rlist[ir].mass;
  double mres2 = mres*mres;
  double prds,prx,pry,prz,er,pu;

  GauLag(0.0,ptmid,p,pw);
  double fm = 0.0;
  for(int ip=0;ip<58;ip++) {
    double eee=sqrt(p[ip]*p[ip]+mres2);
    double aaa=(eee-rlist[ir].mu)*beta;
    if(aaa < 30.0) {
      aaa=exp(aaa);
      fm += p[ip]*p[ip]*pw[ip]/(aaa + rlist[ir].bf);
    }
  }

  double ranemis;
  // double facranmax = 1.6;
  // double facranmax = 1.7;//06/28/2010, lower switching T, larger radial flow
  // double facranmax = 1.8;//08/27/2019, LHC, larger radial flow
  double facranmax = 2.0; // 2014-07-30 for RFH
  double ranmax = dx*dy*dh*tau*facranmax;

  // surface:
  //dsx = tau*dy*deta*dtau
  //dsy = tau*dx*deta*dtau
  //dss = dtau*dx*dy
  if(bulk == 0){
    if(dsx !=0.0 || dsy !=0.0){
      ranmax = dx*dh*tau*dtau*facranmax;
      // Assuming dx = dy
      //	ranmax = dx*dx*dx*2.0*tau*facranmax;
      //                ^^^Coming from MABIKI/3
      //                      in FreezeOutHyperSurface.cxx

    }else{
      ranmax = dtau*dx*dy*facranmax;
      //	ranmax = dx*dx*dx*2.0*facranmax;
      //                ^^^Coming from MABIKI/3
      //                      in FreezeOutHyperSurface.cxx
    }
  }
  do{
    do{

      //Generate momentum [0:6GeV/c]
      //according to Bose/Fermi distribution
      //in local rest frame
      //using bisection method
      double r1 = Random::getRand()*fm;
      double pmax = 6000.0/sctr;
      double pmin = 0.0;
      double ppp = (pmax+pmin)*0.5;

      for(int id=0;id<di;id++) {
        ppp = (pmax+pmin)*0.5;
        Gauss12(0.0,ppp,p1,pw1);
        //	    Gauss38(0.0,ppp,p1,pw1);
        double fp = 0.0;
        for(int ip=0;ip<12;ip++) {
          //	    for(int ip=0;ip<38;ip++) {
          double eee=sqrt(p1[ip]*p1[ip]+mres2);
          double aaa=(eee-rlist[ir].mu)*beta;
          if(aaa <100.0) {
            aaa=exp(aaa);
            fp += p1[ip]*p1[ip]*pw1[ip]/(aaa+rlist[ir].bf);
          }
        }
        double f = fp - r1;
        if(f > 0.0) pmax = ppp;
        else pmin = ppp;
      }
      // random variable on unit sphere
      double r2 = -2*Random::getRand() + 1;
      double theta = acos(r2);
      double phd = 2*M_PI*Random::getRand();

      // uniform random number on surface
      double prxd = ppp*sin(theta)*cos(phd);
      double pryd = ppp*sin(theta)*sin(phd);
      double przd = ppp*cos(theta);
      double erd = sqrt(ppp*ppp+mres2);

      //Lorentz transformation by flow velocity
      double ddd=gamma*(erd+(prxd*vx+pryd*vy+przd*vz)*gamma/(1.+gamma));

      //Momentum in lab. frame
      prx = prxd + vx*ddd;
      pry = pryd + vy*ddd;
      prz = przd + vz*ddd;

      double prt = sqrt(prx*prx+pry*pry);
      er = sqrt(prt*prt+prz*prz+mres2);

      //	double mrt = sqrt(prt*prt+mres2);
      //	double yr = log((er+prz)/(er-prz))*0.5;

      pu = gamma*(er-vx*prx-vy*pry-vz*prz); //pu = erd
      prds = er*ds0+prx*dsx+pry*dsy+prz*dsz;
      if(!ipos)prds=-(er*ds0+prx*dsx+pry*dsy+prz*dsz);

    } while (prds < 0.0);

    if(prds/pu/gamma > ranmax){
      cout << "Warning: prds/pu/gamma is greater than maximum random number. "
           << "Please increase 'facranmax'"
           << " at least "
           << prds/pu/gamma/ranmax*facranmax
           << " [at ParticleSampleHydrojet::getSample]"
           << endl;
      //continue;
    }

    ranemis = ranmax*Random::getRand();

  }while (ranemis >  prds/pu/gamma);

  //↓ は bulk emission の時だけしか正しく無い気がする by KM
  //Uniformly distributed in a fluid element in coordinate space
  double ran1 = Random::getRand();
  double ran2 = Random::getRand();
  double ran3 = Random::getRand();
  double xx = x0+dx*(ran1-0.5);
  double yy = y0+dy*(ran2-0.5);
  double eta= eta0+dh*(ran3-0.5);

  // 2013/04/23, KM, reverse z axis
  if(this->fReverseParticleList){
    prz=-prz;
    eta=-eta;
  }

  if(isOutput) {
    outputData(prx,pry,prz,er,mres,ir,tau,xx,yy,eta,ipos);
  } else {
    putParticle(prx,pry,prz,er,mres,ir,tau,xx,yy,eta,ipos);
  }
}

void ParticleSampleHydrojet::putParticle(double px,double py,double pz,
	double e,double m, int ir, double tau,double x,
	double y, double eta,int ipos)
{
  if(!ipos) return;
  Particle* jp = new(std::nothrow) Particle(ir);
  if (!jp) {
    std::cerr << "(ParticleSampleHydrojet::putPartile:) No more memory" << std::endl;
    exit(1);
  }

  jp->px = px * hbarC;
  jp->py = py * hbarC;
  jp->pz = pz * hbarC;
  //jp->setPe(std::sqrt(m*m+px*px+py*py+pz*pz)*hbarC);
  jp->e = -1.0; // 自動で計算する様にする
  jp->x = x;
  jp->y = y;
  jp->t = tau * std::cosh(eta);
  jp->z = tau * std::sinh(eta);
  this->plist.push_back(jp);
}

void ParticleSampleHydrojet::outputData(double prx,double pry,double prz,
	double er,double mres, int ir, double tau,double xx,
	double yy, double eta, int ipos)
{

	if (ipos) {
	  outdatpos << std::setw(14) << prx*hbarC
              << std::setw(14) << pry*hbarC
              << std::setw(14) << prz*hbarC
              << std::setw(14) << er*hbarC
              << std::setw(14) << mres*hbarC
              << std::setw(4)  << ir
              << std::setw(6)  << tau
              << std::setw(14) << xx
              << std::setw(14) << yy
              << std::setw(14) << eta
              << std::endl;
	} else {
	  outdatneg << std::setw(14) << prx*hbarC
              << std::setw(14) << pry*hbarC
              << std::setw(14) << prz*hbarC
              << std::setw(14) << er*hbarC
              << std::setw(14) << mres*hbarC
              << std::setw(4)  << ir
              << std::setw(6)  << tau
              << std::setw(14) << xx
              << std::setw(14) << yy
              << std::setw(14) << eta
              << std::endl;
	}

}
