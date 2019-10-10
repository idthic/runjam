// -*- C++ -*-
#ifndef kashiwa_phys_Minkowski_hpp
#define kashiwa_phys_Minkowski_hpp
namespace kashiwa{
namespace phys{
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  class vector4{
  public:
    double t;
    double x;
    double y;
    double z;
  public:
    vector4():t(0),x(0),y(0),z(0){}
    vector4(double t,double x,double y,double z):t(t),x(x),y(y),z(z){}
  public:
    // internal product
    double operator*(const vector4& r) const{
      return this->t*r.t-this->x*r.x-this->y*r.y-this->z*r.z;
    }
    // dual-representation
    vector4 operator~() const{
      return vector4(t,-x,-y,-z);
    }
  };
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
}
}
#endif
