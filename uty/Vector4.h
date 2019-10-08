#ifndef Vector4_h
#define Vector4_h

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>

using namespace std;

class Vector4
{

protected:
    double t;   // time or energy
    double x;   // px
    double y;   // py
    double z;   // pz

    static const double eps1;

public:
    //Vector4() {x=0; y=0; z=0; t=0;}
    Vector4(double t1=0.0, double x1=0.0,double y1=0.0,double z1=0.0);
    Vector4(const double* v);
    Vector4(const Vector4& p) : t(p.t),x(p.x),y(p.y),z(p.z) {}
    ~Vector4() { }

    void setXYZT(double x1, double y1,double z1,double t1) {
	x = x1; y = y1; z = z1; t = t1;
    }

    void   setTXYZ(double a,double b,double c,double d) {
	t=a; x=b; y=c; z=d;
    }
    void   set(double a,double b,double c,double d) {
	t=a; x=b; y=c; z=d;
    }
    void   set(Vector4& a) {t=a(0);x=a(1);y=a(2);z=a(3);}
    void   setT(double a) {t=a;}
    void   setX(double a) {x=a;}
    void   setY(double a) {y=a;}
    void   setZ(double a) {z=a;}

    double getT() const {return t;}
    double getX() const {return x;}
    double getY() const {return y;}
    double getZ() const {return z;}

    double T() const {return t;}
    double X() const {return x;}
    double Y() const {return y;}
    double Z() const {return z;}

    double  operator () (int i) const;
    double  operator [] (int i) const {return (*this)(i);}

    double& operator () (int i);
    double& operator [] (int i) {return (*this)(i);}

    Vector4& operator = (const Vector4& p) {
	t = p.t; x = p.x; y = p.y; z = p.z;
	return *this;
    }
    bool operator == (const Vector4& v) const {
	return (v.t == t && v.x == x && v.y == y && v.z == z)? true:false;
    }
    bool operator != (const Vector4& v) const {
	return (v.t != t || v.x != x || v.y != y || v.z != z)? true:false;
    }
    Vector4 operator - () const {
	return Vector4(-t,-x,-y,-z);
    }
    Vector4& operator -= (const Vector4& p) {
	t -= p.t; x -=  p.x; y -= p.y; z -= p.z;
	return *this;
    }

    Vector4& operator *= (double a) {
	t *=a; x *= a; y *= a; z *= a;
	return *this;
    }
    Vector4& operator /= (double a) {
	t /=a; x /= a; y /= a; z /= a;
	return *this;
    }
    Vector4& operator /= (int a) {
	t /=a; x /= a; y /= a; z /= a;
	return *this;
    }
    Vector4& operator += (Vector4& a) {
	t += a.t; x += a.x; y += a.y; z += a.z;
	return *this;
    }
    Vector4& operator += (double a) {
	t +=a; x += a; y += a; z += a;
	return *this;
    }
    Vector4& operator -= (double a) {
	t -=a; x -= a; y -= a; z -= a;
	return *this;
    }

    /*
    Vector4 operator + (const Vector4& a) {
	return Vector4(a.T()+t, a.X()+x, a.Y()+y, a.Z()+z);
    }
    Vector4 operator - (const Vector4& a) {
	return Vector4(t-a.T(), x-a.X(), a.y-Y(), z-a.Z());
    }
    Vector4 operator * (const double a) {
	return Vector4(t*a, x*a, y*a, z*a);
    }
    */
    Vector4 operator / (const double a) {
	return Vector4(t/a, x/a, y/a, z/a);
    }
    Vector4 operator / (const int a) {
	return Vector4(t/a, x/a, y/a, z/a);
    }


    double pt2()    const {return x*x + y*y;}
    double pt()     const {return sqrt(pt2());}
    double dot3()   const {return x*x + y*y + z*z;}
    double dot4a()  const {return t*t + x*x + y*y + z*z;}
    double pp()     const {return sqrt(dot3());}
    double sqdot3() const {return sqrt(dot3());}
    double dot4(Vector4& p2) { return t*p2.t - x*p2.x -y*p2.y - z*p2.z;}
    double dot4()   const {return t*t - x*x - y*y - z*z;}
    double dot4sq() const {return sqrt( std::max(0.0,dot4()));}
    void   setE(const double m)  { t = sqrt(m*m + dot3());}

    double rapidity() const {
	return 0.5*log( std::max(t + z,1e-8)/std::max(t - z,1e-8) );
    }
    double pseudoRapidity() const {
	return 0.5*log( std::max(sqdot3() + z,1e-8)/std::max(sqdot3() - z,1e-8) );
    }
    double et() const { return t*pt()/std::max(pp(),1e-8);}
    double v2() const { return (x*x - y*y)/std::max(pt2(),1e-10);}
    double v4() const { 
      return (x*x*x*x+y*y*y*y-6.0*x*x*y*y)/std::max(pt2()*pt2(),1e-10);
    }

    double Angle();
    void   LorentzBoost(const double bex,const double bey,const double bez);
    void   Rotate(const double the, const double phi);
    void   RotateBoost(const double the, const double phi,
		const double bex,const double bey,const double bez) {
	Rotate(the,phi); LorentzBoost(bex,bey,bez);
    }

    friend ostream& operator<<(ostream& os, const Vector4& p) {
	char str[1024];
	sprintf(str,"(%g, %g, %g, %g)",p.t,p.x,p.y,p.z);
	os << str;
	return os;
    }

};

Vector4 operator + (const Vector4& a, const Vector4& b);
Vector4 operator - (const Vector4& a, const Vector4& b);
Vector4 operator * (const Vector4& a, const double b);
Vector4 operator * (const double a, const Vector4& b);


#endif // Vector4_h
