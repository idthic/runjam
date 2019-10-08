#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include "ParticleSampleFromOversampledPhasespace.h"

void ParticleSampleFromOversampledPhasespace::readPhasespaceDat(){
  this->clearParticleList();

  bool eventCountKnown = m_overSamplingFactor >= 0;
  this->pcache.clear();
  if(eventCountKnown)
    this->pcache.reserve(this->m_overSamplingFactor);

  int iline;
  {
    std::ifstream ifs(this->fname_phasespace_dat.c_str());
    if(!ifs)goto error_failed_to_open;

    int isample;
    std::string line;
    iline=0;
    for(isample=0;;isample++){
      iline++;
      if(!std::getline(ifs,line))goto error_invalid_format;
      int npart,dummy;
      if(2!=std::sscanf(line.c_str()," %d %d",&npart,&dummy)) break;

      if(isample==this->m_overSamplingFactor) goto error_invalid_format;

      if(npart < 0) goto error_invalid_format;

      this->pcache.resize(isample+1);
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
        this->pcache[isample].push_back(this->plist.back());
      }
    }

    // 最終行の中身が -999 かどうかチェックする。
    //   条件: line に最終行が入っていること。
    int lastNumber;
    if(1!=std::sscanf(line.c_str()," %d",&lastNumber)||lastNumber!=-999)
      goto error_invalid_format;
    if(eventCountKnown && isample != m_overSamplingFactor)
      goto error_invalid_format;
    return;
  }

 error_failed_to_open:
  std::cerr
    <<"spectra/ParticleSampleFromOversampledPhasespace.cxx(ParticleSampleFromOversampledPhasespace::update): failed to open the file ("
    <<this->fname_phasespace_dat<<")"<<std::endl;
  std::exit(EXIT_FAILURE);
  return; /*NOTREACHED*/

 error_invalid_format:
  std::cerr
    <<this->fname_phasespace_dat<<":"<<iline<<": invalid format (@ParticleSampleFromOversampledPhasespace::update)"
    <<std::endl;
  std::exit(EXIT_FAILURE);
  return; /*NOTREACHED*/
}

void ParticleSampleFromOversampledPhasespace::update(){
  if(this->m_currentSampleIndex<0){
    this->readPhasespaceDat();
    this->m_currentSampleIndex=0;
  }else{
    this->m_currentSampleIndex++;
  }
}
