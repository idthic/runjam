#include <cmath>
#include <algorithm>
#include "ParticleSampleDebugJamAsymmetry.h"

void ParticleSampleDebugJamAsymmetry::update(){
  const double hbarC=0.197327053;
  this->clearParticleList();

  // // 002
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);

  // 003
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.101,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-0.099,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+0.099,1);
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.101,1);

  // // 003a
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-5.001,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-4.999,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+4.999,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+5.001,1);

  // // 003b
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-5.001,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-4.999,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+4.999,1);
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+5.001,1);

  // // 003c
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.003,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-0.001,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+0.001,1);
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.003,1);

  // // 004b
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+5.001,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+4.999,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-4.999,1);
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-5.001,1);

  // // 005
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  // }

  // // 005a 005b
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+3.0,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-1.0,0.13957, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(17, 0,0,-3.0,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,+1.0,0.13957, 0,0,-0.001,1);
  // }

  // // 005c
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.001,1);
  // }

  // // 005d
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+5.0,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-0.0,0.13957, 0,0, 0.000,1);
  // }else{
  //   this->addParticle(17, 0,0,-5.0,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,+0.0,0.13957, 0,0, 0.000,1);
  // }

  // // 005e
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+1.0,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(17, 0,0,-1.0,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  // }

  // // 005f
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+2.0,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(17, 0,0,-2.0,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  // }

  // // 006 -> crash?
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.002,1);
  // this->addParticle(17, 0,0, 0.0,0.13957, 0,0, 0.000,1);
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.002,1);

  // // 006a
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.002,1);
  //   this->addParticle(17, 0,0, 0.0,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.002,1);
  // }else{
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.002,1);
  //   this->addParticle(17, 0,0, 0.0,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.002,1);
  // }

  // // 006b
  // if(this->count++&1){
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.002,1);
  //   this->addParticle(17, 0,0,-0.5,0.13957, 0,0,+0.001,1);
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.002,1);
  // }else{
  //   this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.002,1);
  //   this->addParticle(17, 0,0,+0.5,0.13957, 0,0,-0.001,1);
  //   this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.002,1);
  // }

  // // 206
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.101,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+0.099,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-0.099,1);
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.101,1);

  // // 207
  // this->addParticle(17, 0,0,+1.5,0.13957, 0,0,-0.101,1);
  // this->addParticle(17, 0,0,-0.5,0.13957, 0,0,-0.099,1);
  // this->addParticle(17, 0,0,+0.5,0.13957, 0,0,+0.099,1);
  // this->addParticle(17, 0,0,-1.5,0.13957, 0,0,+0.101,1);
  // std::random_shuffle(this->plist.begin(),this->plist.end());

  //---------------------------------------------------------------------------
  // nucleons

  // // 301
  // // double efactor=1.0; // 301
  // double efactor=2.0; // 301a
  // // double efactor=0.5; // 301b
  // // double efactor=4.0; // 301c
  // if(this->count++&1){
  //   this->addParticle(19, 0,0,+1.0/hbarC*efactor,0.94, 0,0,-0.001,1);
  //   this->addParticle(19, 0,0,-0.5/hbarC*efactor,0.94, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(19, 0,0,-1.0/hbarC*efactor,0.94, 0,0,+0.001,1);
  //   this->addParticle(19, 0,0,+0.5/hbarC*efactor,0.94, 0,0,-0.001,1);
  // }

  // // 302
  // double efactor=1.0; // 302
  // // double efactor=2.0; // 302a
  // // double efactor=0.5; // 302b
  // // double efactor=4.0; // 302c
  // if(this->count++&1){
  //   this->addParticle(19, +1.0/hbarC*efactor,+1.0/hbarC*efactor,+1.0/hbarC*efactor,0.94, 0,0,-0.001,1);
  //   this->addParticle(19, -0.5/hbarC*efactor,-0.5/hbarC*efactor,-0.5/hbarC*efactor,0.94, 0,0,+0.001,1);
  // }else{
  //   this->addParticle(19, -1.0/hbarC*efactor,-1.0/hbarC*efactor,-1.0/hbarC*efactor,0.94, 0,0,+0.001,1);
  //   this->addParticle(19, +0.5/hbarC*efactor,+0.5/hbarC*efactor,+0.5/hbarC*efactor,0.94, 0,0,-0.001,1);
  // }

  // // 303
  // double efactor=1.0; // 302
  // // double efactor=2.0; // 302a
  // // double efactor=0.5; // 302b
  // // double efactor=4.0; // 302c
  // if(this->count++&1){
  //   this->addParticle(19, +1.0/hbarC*efactor,+1.0/hbarC*efactor,+1.0/hbarC*efactor,0.94, -0.002,-0.002,-0.002,1);
  //   this->addParticle(19, +0.5/hbarC*efactor,-0.5/hbarC*efactor,-0.5/hbarC*efactor,0.94, -0.001,+0.001,+0.001,1);
  // }else{
  //   this->addParticle(19, -1.0/hbarC*efactor,-1.0/hbarC*efactor,-1.0/hbarC*efactor,0.94, +0.002,+0.002,+0.002,1);
  //   this->addParticle(19, -0.5/hbarC*efactor,+0.5/hbarC*efactor,+0.5/hbarC*efactor,0.94, +0.001,-0.001,-0.001,1);
  // }
}
