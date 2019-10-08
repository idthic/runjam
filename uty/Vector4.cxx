#include "Vector4.h"
#include <cstdlib>

using namespace std;

const double Vector4::eps1=1-1e-12;

Vector4::Vector4(double t1, double x1, double y1,double z1)
	:t(t1), x(x1),y(y1),z(z1) { }

Vector4::Vector4(const double* v)
    : t(v[0]), x(v[1]), y(v[2]), z(v[3]) { }


double Vector4::operator () (int i) const
{
    switch(i) {
	case 0:
	    return t;
	case 1:
	    return x;
	case 2:
	    return y;
	case 3:
	    return z;
	default:
	    cerr << "Vector4::operator() bad index i= " << i << endl;
	    abort();
    }
    return 0;
}

double& Vector4::operator () (int i)
{
    switch(i) {
	case 0:
	    return t;
	case 1:
	    return x;
	case 2:
	    return y;
	case 3:
	    return z;
	default:
	    cerr << "Vector4::operator() bad index i= " << i << endl;
	    abort();
    }
    return t;
}

Vector4 operator + (const Vector4& a, const Vector4& b)
{
    return Vector4(a.T()+b.T(), a.X()+b.X(), a.Y()+b.Y(),a.Z()+b.Z());
}
Vector4 operator - (const Vector4& a, const Vector4& b)
{
    return Vector4(a.T()-b.T(), a.X()-b.X(), a.Y()-b.Y(),a.Z()-b.Z());
}
Vector4 operator * (const Vector4& a, const double b)
{
    return Vector4(a.T()*b, a.X()*b, a.Y()*b,a.Z()*b);
}

Vector4 operator * (const double a, const Vector4& b)
{
    return Vector4(b.T()*a, b.X()*a, b.Y()*a,b.Z()*a);
}

double Vector4::Angle()
{
    double angl=0.0;
    double r=sqrt(x*x+y*y);
    if(r < 1e-20) return angl;

    if(fabs(x)/r < 0.8) {
        //angl=sign(acos(x/r),y)
	angl = y>0 ? fabs(acos(x/r)): -fabs(acos(x/r));
    }else {
        angl=asin(y/r);
        if(x < 0.0 && angl >= 0.0)
          angl=M_PI-angl;
        else if(x < 0.0)
          angl=-M_PI-angl;
    }
 
    return angl;
}

void Vector4::LorentzBoost(const double bex,const double bey,const double bez)
{
    if(bex*bex+bey*bey+bez*bez <= 1e-20) return;

    double dbx=bex;
    double dby=bey;
    double dbz=bez;
    double db=sqrt(dbx*dbx+dby*dby+dbz*dbz);

//...Rescale boost vector if too close to unity.
    if(db > eps1) {
        cerr << "(Boost:) boost vector too large " << db << endl;
        dbx=dbx*(eps1/db);
        dby=dby*(eps1/db);
        dbz=dbz*(eps1/db);
        db=eps1;
    }
    double dga=1.0/sqrt(1.0-db*db);
    Vector4 dp(t,x,y,z);

    double dbp = dbx*x + dby*y + dbz*z;
    double dgabp=dga*(dga*dbp/(1.0+dga) + t);
    x = dp.x + dgabp*dbx;
    y = dp.y + dgabp*dby;
    z = dp.z + dgabp*dbz;
    t = dga*(dp.t + dbp);
} 

void Vector4::Rotate(const double the, const double phi)
{
//...Local arrays.
      double rot[3][3],pr[3];
 
//...Rotate, typically from z axis to direction (theta,phi).
    if(the*the+phi*phi > 1e-20) {
        rot[0][0]=cos(the)*cos(phi);
        rot[0][1]=-sin(phi);
        rot[0][2]=sin(the)*cos(phi);
        rot[1][0]=cos(the)*sin(phi);
        rot[1][1]=cos(phi);
        rot[1][2]=sin(the)*sin(phi);
        rot[2][0]=-sin(the);
        rot[2][1]=0.0;
        rot[2][2]=cos(the);
	pr[0]=x;
	pr[1]=y;
	pr[2]=z;
	x=rot[0][0]*pr[0]+rot[0][1]*pr[1]+rot[0][2]*pr[2];
	y=rot[1][0]*pr[0]+rot[1][1]*pr[1]+rot[1][2]*pr[2];
	z=rot[2][0]*pr[0]+rot[2][1]*pr[1]+rot[2][2]*pr[2];
    }

}
