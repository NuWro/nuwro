#ifndef _vect_h_
#define _vect_h_
#include <iostream>
#include <cassert>
#include <cmath>
#include "vec.h"
#include "vecrand.h"
#include "jednostki.h"

using namespace std;


/// Minkowski 4-vector
class vect
{
public:
  double t, x, y, z;
//  friend class spinor;

public:
  inline vect (double x0 = 0, double x1 = 0, double x2 = 0, double x3 = 0);
  inline vect (vec v0, double t0=0);
  inline vect (double t0, vec v0);
  inline vec operator= (vec v);
  inline vect& operator+= (vect v);
  inline vect& operator-= (vect v);
  inline vect& operator+= (vec v);
  inline vect& operator-= (vec v);
  inline operator vec ();
  inline double &operator[] (int i);
  inline friend double operator* (vect a, vect b);    ///< scalar product   
  inline friend vect operator* (double a, vect b);    ///< multiply 4-vector by a number
  inline friend vect operator* (vect b, double a);    ///< multiply 4-vector by a number
  inline friend vect operator/ (vect b, double a);    ///< divide 4-vector by a number
  inline friend vect operator+ (vect a, vect b);      ///< add 4-vectors
  inline friend vect operator- (vect a, vect b);      ///< subtract 4-vectors
  inline vect operator- ();                           ///< negate 
//  inline vec speed ();                                ///< momentum / energy
  inline vec v ()    ;                                ///< momentum / energy
  inline vect & boost (vec v);                        ///< boost to a frame moving with velocity v
  inline vect & boost1 (vec v);                       ///< boost to a frame moving with velocity v
  inline vect & boost2 (vec v);                       ///< boost from a frame moving with velocity v
  inline vect & boost3 (vec v);                       ///< boost to a frame moving with velocity v
  inline vect & boost4 (vec v);                       ///< boost to a frame moving with velocity v
  inline friend ostream & operator<< (ostream & o, vect v); ///< print a 4 vector
  inline double length();                             ///< length of the space component
};

inline vec speed (vect a, vect b);            ///< CMS speed of a and b
inline vec rand_from_ball (double r);         ///< random vector from ball of radius r
inline vec rand_dir ();                       ///< random unit vector in 3D space  
inline vec bodek_rand_from_ball (double r);   ///< punkt ze sfery z ogonem Bodka
inline vec spectral_choice (int p, int n);
inline vec deuterium ();
inline double q2 (vect k, vect kprim);
inline double q0 (vect k, vect kprim);


////////////////////////////////////////////////////////////////
///            I M P L E M E N T A C J A                     ///  
////////////////////////////////////////////////////////////////

vect::vect (double x0, double x1, double x2, double x3)
           :t (x0), x (x1), y (x2), z (x3)
  {
  }
vect::vect (vec v0, double t0):t (t0), x (v0.x), y (v0.y), z (v0.z)
  {
  }
vect::vect (double t0, vec v0):t (t0), x (v0.x), y (v0.y), z (v0.z)
  {
  }
double & vect::operator[] (int i)
  {
    return *(&t + i);
  }
double operator* (vect a, vect b)
  {
    return a.t * b.t - a.x * b.x - a.y * b.y - a.z * b.z;
  }
vect operator+ (vect a, vect b)
  { 
    return vect (a.t + b.t, a.x + b.x, a.y + b.y, a.z + b.z);
  }
vect operator- (vect a, vect b)
  {
    return vect (a.t - b.t, a.x - b.x, a.y - b.y, a.z - b.z);
  }
vect vect::operator- ()
  {
    return vect (-t, -x, -y, -z);
  }
vec vect::operator= (vec v)
  {
    x = v.x;
    y = v.y;
    z = v.z;
    return v;
  }
vect& vect::operator+= (vect v)
  {
    x+=v.x;
    y+=v.y;
    z+=v.z;
    t+=v.t;
    return *this;
  }
vect& vect::operator-= (vect v)
  {
    x-=v.x;
    y-=v.y;
    z-=v.z;
    t-=v.t;
    return *this;
  }
vect& vect::operator+= (vec v)
  {
    x+=v.x;
    y+=v.y;
    z+=v.z;
    return *this;
  }
vect& vect::operator-= (vec v)
  {
    x-=v.x;
    y-=v.y;
    z-=v.z;
    return *this;
  }
vect operator* (double a, vect b)
  {
    return vect (a * b.t, a * b.x, a * b.y, a * b.z);
  }
vect operator* (vect b, double a)
  {
    return vect (a * b.t, a * b.x, a * b.y, a * b.z);
  }
vect operator/ (vect b, double a)
  {
    return vect (b.t / a, b.x / a, b.y / a, b.z / a);
  }
ostream & operator<< (ostream & o, vect v)
  {
    return o << '(' << v.t << ',' << v.x << ',' << v.y << ',' << v.z << ')';
  }
vect::operator vec ()
  {
    return vec (x, y, z);
  }  
/*vec vect::speed ()
  { 
    return vec (x / t, y / t, z / t);
  }
*/
vec vect::v ()
  { 
    return vec (x / t, y / t, z / t);
  }
// Alternative boost implementations 
vect & vect::boost1 (vec v)  // orig boost
  {//cout<<v<<endl;				
    double vv = v * v;
    //cout<<vv<<endl;
    if(vv==0) return *this;
    if(vv>=1) 
       cout<<v<<vv<<flush;
    assert (vv < 1);
//     double old_norm=*this * *this;
    vec p = *this;
    //double dv=sqrt(vv);
    vec zmienne = (p * v) * v / vv;
    vect stale = vect (p - zmienne);
    *this = *this - stale;
    double gamma = 1 / sqrt (1 - vv);
    *this =
      vect (t + x * v.x + y * v.y + z * v.z, x + t * v.x, y + t * v.y,
	    z + t * v.z) * gamma + stale;
//     double new_norm=*this * *this;
//     cout<<"Niedokadność boostowania="<<old_norm-new_norm<<"/"<<old_norm<<endl;
    return *this;
  }

vect & vect::boost2(vec v)   //inverse boost
  {				
    return boost(-v);
  }

double vect::length()
  { 
    return sqrt(x*x+y*y+z*z);
  }

vect & vect::boost3 (vec v)    //boost3 (robi błędy)
  {
    double vv = v * v;
    if(vv==0) return *this;
    if(vv>=1) 
       cout<<v<<vv<<flush;
    assert (vv < 1);
    double pv = x*v.x+y*v.y+z*v.z;
    double gamma = 1 / sqrt (1 - vv);
    double A=(gamma-1)*pv/vv+gamma*t;
    t=gamma*(t+pv);
    x+= A*v.x;
    y+= A*v.y;
    z+= A*v.z;
    return *this;
  }
  
/// Boost as two Lorentz reflections
vect & vect::boost4 (vec v) // boost 4
  {
    double vv = v * v;
    if(vv==0) return *this; // not necessary in this version of boost
    if(vv>=1) 
       cout<<v<<vv<<flush;
    assert (vv < 1);
    t=-t;
    vect s(v,sqrt(1-vv)+1);
    return *this-=2*(*this*s)/(s.t*s.t-vv)*s;
  }

vect & vect::boost (vec v)    
{ return boost1(v);
}


///////////////////////////////////////////////////////////////////////////
///                  FUNKCJE POMOCNICZE
///////////////////////////////////////////////////////////////////////////

vec
speed (vect a, vect b)
{
  vect c = a + b;
  return c.v ();
}


    
double
q2 (vect k, vect kprim)
{
  vect q = k - kprim;
  return q * q;
}

double
q0 (vect k, vect kprim)
{
  return k.t - kprim.t;
}



#endif
