#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <ksh/phys/Minkowski.hpp>
#include <spectra/IntegratedCooperFrye.hpp>

#include "ElementReso.hpp"

namespace idt {
namespace hydro2jam {

ElementReso::ElementReso(std::string dir, std::string* outf, int kint, int eos_pce, std::string fname):
  HydroSpectrum(kint, eos_pce), rlist(kint, eos_pce, fname)
{
  tmpf = 0.16 / hbarc_MeVfm * 1000.0;
  mubf = 1.6 / hbarc_MeVfm * 1000.0;
  meanf = 0.45 / hbarc_MeVfm * 1000.0;

  int const nreso_loop = this->rlist.numberOfResonances();
  this->elemFile.resize(nreso_loop);
  this->outdat.resize(nreso_loop);
  this->outdatPos.resize(nreso_loop);
  this->outdatNeg.resize(nreso_loop);

  for (int i = 0; i < nreso_loop; i++) {
    if (dir.size() >0)
      elemFile[i] = dir + "/" + outf[i];
    else
      elemFile[i] = outf[i];
  }
}

ElementReso::~ElementReso()
{
    // delete [] mass;
    // delete [] deg;
    // delete [] degeff;
    // delete [] mu;
}

void ElementReso::initialize()
{
  int const nreso_loop = this->rlist.numberOfResonances();
  //    nreso_loop = nreso;
  //    if(baryonfree)nreso_loop = 20;

  // Open output files.
  for (int i = 0; i < nreso_loop; i++) {
    if (i < 21)	outdat[i].open(elemFile[i].c_str());
    std::string elemFilepos = elemFile[i] + ".POS";
    outdatPos[i].open(elemFilepos.c_str());
    std::string elemFileneg = elemFile[i] + ".NEG";
    outdatNeg[i].open(elemFileneg.c_str());
    if (!outdat[i]) {
      // ??? outdat[i] is not opened if(i>=21).
      std::cerr << "Error: unable to open file"
                << elemFile[i] << std::endl;
      std::exit(1);
    }
  }

  //...integral region
  //    double ymin = -6.0;
  //    double ymax = 6.0;
  double ptmin = 1e3 / hbarc_MeVfm;

  // numerical integral
  Gauss12(0.0, 2. * pi, phi, phiw);
  //    Gauss38(ymin, ymax, y, yw);
  GauLag(0.0, ptmin, pt, ptw);
}

void ElementReso::analyze(std::string fnameFreezeoutDat) {
  int const nreso_loop = this->rlist.numberOfResonances();
  for(int ireso = 0; ireso < nreso_loop; ireso++)
    ElementReso::integrateForResonance(fnameFreezeoutDat, ireso);
}

void ElementReso::integrateForResonance(std::string const& fnameFreezeoutDat, int ireso) {
  //---------------------------------------------------------------------------
  // initialize files

  // open output files
  std::ofstream ostr_elm;
  if (ireso < 21) {
    ostr_elm.open((elemFile[ireso]).c_str());
    if (!ostr_elm) {
      std::cerr << "ElementReso::analyze_impl2! failed to create file " << elemFile[ireso] << std::endl;
      std::exit(1);
    }
  }
  std::ofstream ostr_pos((elemFile[ireso] + ".POS").c_str());
  if (!ostr_pos) {
    std::cerr << "ElementReso::analyze_impl2! failed to create file " << elemFile[ireso] << ".POS" << std::endl;
    std::exit(1);
  }
  std::ofstream ostr_neg((elemFile[ireso] + ".NEG").c_str());
  if (!ostr_neg) {
    std::cerr << "ElementReso::analyze_impl2! failed to create file " << elemFile[ireso] << ".NEG" << std::endl;
    std::exit(1);
  }

  ResonanceListPCE::resonance& recreso = this->rlist[ireso];

  openDataFile(fnameFreezeoutDat);

  //---------------------------------------------------------------------------
  while (!readData()) {
    if (!baryonfree) {
      recreso.mu = 0.0;
      if (recreso.bf == 1) recreso.mu = mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
      if (recreso.anti) recreso.mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) + meanf * nbf;
      //if (recreso.anti) recreso.mu = -mubf * sqrt(1.0 - tf * tf / tmpf / tmpf) - meanf * nbf;
      //if (recreso.bf==1) recreso.mu = mub;
      //if (recreso.anti) recreso.mu = -mub;
    }

    if (tf == 0.0) {
      std::cerr << "(ElementReso::analyze:) TF=ZERO" << std::endl;
      //write(24, 9000)0.0
      continue;
    }

    //-------------------------------------------------------------------------
    // integration

    double npos = 0.0;
    double nneg = 0.0;

    double const gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    kashiwa::phys::vector4 u(gamma, vx * gamma, vy * gamma, vz * gamma);
    kashiwa::phys::vector4 ds(ds0, -dsx, -dsy, -dsz);
    double const beta = 1.0 / tf;

    if (recreso.bf == -1) {
      idt::hydro2jam::IntegrateBosonCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
    } else {
      idt::hydro2jam::IntegrateFermionCooperFrye(npos, nneg, u, ds, beta, recreso.mass, recreso.mu);
    }

    double const n = (nneg + npos) * recreso.degeff;
    npos *= recreso.deg;
    nneg *= recreso.deg;
    //-------------------------------------------------------------------------

    if (ireso < 21) ostr_elm << n << "\n";
    ostr_pos << npos << "\n";
    ostr_neg << nneg << "\n";
  }

  //  close files.
  if (ireso < 21) ostr_elm.close();
  ostr_pos.close();
  ostr_neg.close();

  closeDataFile();
}

}
}
