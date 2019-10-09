// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_UTIL_H
#define KASHIWA_UTIL_H
#include <cstdlib>
#include <cstring>
namespace kashiwa{

  static inline bool getenvAsBool(const char* envname,bool defaultValue){
    const char* str=std::getenv(envname);
    if(!str||!*str)
      return defaultValue;
    else if(isdigit(*str))
      return std::atoi(str)!=0;
    else
      return std::strcmp(str,"false")!=0;
  }

  static inline int getenvAsInt(const char* envname,int defaultValue){
    const char* str=std::getenv(envname);
    if(str&&isdigit(*str))
      return std::atoi(str);
    else
      return defaultValue;
  }

}
#endif
