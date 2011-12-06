#ifndef _quat_h_
#define _quat_h_
#include <iostream>
#include <cassert>
#include <cmath>
#include "vec.h"

using namespace std;

class quat
{
public:
  double w, x, y, z;

public:
  inline quat (double x0 = 0, double x1 = 0, double x2 = 0, double x3 = 0);
  inline quat (double t, vec a):w(t),x(a.x),y(a.y),z(a.z){}
  inline double &operator[] (int i);
  inline friend quat operator* (quat a, quat b);      // quaternion product   
  inline friend quat operator* (double a, quat b);    // multiply quat by a number
  inline friend quat operator* (quat b, double a);    // multiply quat by a number
  inline friend quat operator/ (quat b, double a);    // divide quat by a number
  inline friend quat operator/ (quat a, quat b);      // divide quat by a quat 
  inline friend quat operator/ (double a, quat b);    // divide number by a quat
  inline friend quat operator+ (quat a, quat b);      // add quats
  inline friend quat operator- (quat a, quat b);      // subtract quats
  inline quat operator- ();                           // negate 
  inline quat operator~ ();                           // conjugate
  inline quat conj ();                                // conjugate
  inline double norm2 ();
  inline double norm ();
  inline quat normalize ();
  inline quat inv ();
  inline friend ostream & operator<< (ostream & o, quat v); // print a quat
  inline vec operator*(vec v)
  {quat a=*this*quat(0,v)*~*this/norm2();
   return vec(a.x,a.y,a.z);
  }
};


inline quat rot(vec os,double alfa)
{os*=sin(alfa/2)/os.length();
 return quat(cos(alfa/2),os.x,os.y,os.z);
} 

/// quat transforming a into b (unit vectors)
inline quat rot(vec a,vec b)
{ a=a.dir();
  b=b.dir();
  vec c=(a+b).dir();
  if(c.norm2()==0)
   
 return quat(a*c,vecprod(a,c));
} 



////////////////////////////////////////////////////////////////
///            I M P L E M E N T A C J A                     ///  
////////////////////////////////////////////////////////////////

quat::quat (double x0, double x1, double x2, double x3)
           :w (x0), x (x1), y (x2), z (x3)
  {
  }
double & quat::operator[] (int i)
  {
    return *(&w + i);
  }
quat operator* (quat a, quat b)
  { return quat(
       (a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z),
       (a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y),
       (a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x),
       (a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w));
  }
quat operator+ (quat a, quat b)
  { 
    return quat (a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
  }
quat operator- (quat a, quat b)
  {
    return quat (a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z);
  }
quat quat::operator- ()
  {
    return quat (-w, -x, -y, -z);
  }
quat quat::operator~ ()
  {
    return quat (w, -x, -y, -z);
  }
quat quat::conj ()
  {
    return quat (w, -x, -y, -z);
  }
double quat::norm2 ()
  {
    return w*w+x*x+y*y+z*z;
  }
double quat::norm ()
  {
    return sqrt(w*w+x*x+y*y+z*z);
  }
quat quat::normalize ()
  { double n=norm();
    w/=n;
    x/=n;
    y/=n;
    z/=n;
    return *this;
  }
quat quat::inv ()
  { double n2=w*w+x*x+y*y+z*z;
    return quat(w/n2,x/n2,y/n2,z/n2);
  }
quat operator* (double a, quat b)
  {
    return quat (a * b.w, a * b.x, a * b.y, a * b.z);
  }
quat operator* (quat b, double a)
  {
    return quat (a * b.w, a * b.x, a * b.y, a * b.z);
  }
quat operator/ (quat b, double a)
  {
    return quat (b.w / a, b.x / a, b.y / a, b.z / a);
  }
quat operator/ (quat a, quat b)
  {
    return a*b.inv();
  }
quat operator/ (double a, quat b)
  {
    return a*b.inv();
  }
ostream & operator<< (ostream & o, quat v)
  {
    return o << '(' << v.w << ',' << v.x << ',' << v.y << ',' << v.z << ')';
  }
    
#endif
