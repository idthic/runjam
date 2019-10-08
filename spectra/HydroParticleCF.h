// -*- mode:c++;indent-tabs-mode:nil -*-
#ifndef HydroParticleCF_h
#define HydroParticleCF_h

#include <uty/Vector4.h>

// This is a class for particles from Cooper-Frye
class HydroParticleCF
{
private:
  int id;
  double mass;
  Vector4 r;
  Vector4 p;
public:
  HydroParticleCF(int i) {
    id=i;
    r.setTXYZ(0.0, 0.0, 0.0, 0.0);
    p.setTXYZ(0.0, 0.0, 0.0, 0.0);
  }
  ~HydroParticleCF() { }

  int    getID() const {return id;}
  void   setID(int i) {id=i;}

  double getPx() const  {return p.getX();}
  double getPy() const  {return p.getY();}
  double getPz() const  {return p.getZ();}
  double getPe() const  {return p.getT();}
  double getX()  const  {return r.getX();}
  double getY()  const  {return r.getY();}
  double getZ()  const  {return r.getZ();}
  double getT()  const  {return r.getT();}

  void   setPx(double px)    {p.setX(px);}
  void   setPy(double py)    {p.setY(py);}
  void   setPz(double pz)    {p.setZ(pz);}
  void   setPe(double e)     {p.setT(e);}
  void   setX(double x)      {r.setX(x);}
  void   setY(double y)      {r.setY(y);}
  void   setZ(double z)      {r.setZ(z);}
  void   setT(double t)      {r.setT(t);}
  void   setCoordinate(Vector4& rr) {r=rr;}
  void   setMomentum(Vector4& pp)   {p=pp;}
};
#endif
