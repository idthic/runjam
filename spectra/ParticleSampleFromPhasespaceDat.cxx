#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "ParticleSampleFromPhasespaceDat.h"

void ParticleSampleFromPhasespaceDat::update(){
  this->clearParticleList();

  int iline;
  {
    std::ifstream ifs(this->fname_phasespace_dat.c_str());
    if(!ifs)goto error_failed_to_open;

    std::string line;
    iline=1;
    if(!std::getline(ifs,line))goto error_invalid_format;

    int npart,dummy;
    if(2!=std::sscanf(line.c_str()," %d %d",&npart,&dummy))
      goto error_invalid_format;

    for(int i=0;i<npart;i++){
      iline++;
      if(!std::getline(ifs,line))goto error_invalid_format;

      int kc,kf;
      double px,py,pz,m;
      double x,y,z,t;
      if(10!=std::sscanf(
        line.c_str()," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
        &kc,&kf,&px,&py,&pz,&m,&x,&y,&z,&t
      ))goto error_invalid_format;

      this->addParticleMinkowski(kf,px,py,pz,m,x,y,z,t);
    }

    iline++;
    if(!std::getline(ifs,line))
      goto error_invalid_format;
    if(1!=std::sscanf(line.c_str()," %d",&npart)||npart!=-999)
      goto error_invalid_format;
    return;
  }

 error_failed_to_open:
  std::cerr
    <<"spectra/ParticleSampleFromPhasespaceDat.cxx(ParticleSampleFromPhasespaceDat::update): failed to open the file ("
    <<this->fname_phasespace_dat<<")"<<std::endl;
  std::exit(EXIT_FAILURE);
  return; /*NOTREACHED*/

 error_invalid_format:
  std::cerr
    <<this->fname_phasespace_dat<<":"<<iline<<": invalid format (@ParticleSampleFromPhasespaceDat::update)"
    <<std::endl;
  std::exit(EXIT_FAILURE);
  return; /*NOTREACHED*/
}
