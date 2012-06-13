#ifndef _vec_h_
#define _vec_h_
#include <iostream>
#include <sstream>
#include <cmath>
//#include <stdlib.h>
//#include <string>
//#include <string.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////
///   vec is a 3D vector                                                 ///
////////////////////////////////////////////////////////////////////////////

class vec
{
public:
  double x, y, z;                                 ///< Cartesian coordinates 

public:
   inline double &operator[] (int i);              ///< treat vec as double[3]
   inline vec();                                   ///< vec=0 
   inline vec (double x1, double x2, double x3);   ///< vec from coordinates
   inline vec (string v);                    ///< vec from string ex. "1 2 3"
   inline double length ();                        ///< vector length
   inline double norm () const;                    ///< vector length
   inline double norm2 () const;                   ///< vector length squared
   inline void normalize ();                       ///< vector length
   inline vec dir();                               ///< normalized vector
   inline friend double operator* (vec a, vec b);  ///< scalar product
   inline friend double angle (vec a, vec b);      ///< angle between a and b
   inline friend vec vecprod(vec a, vec b);        ///< vector product
   inline friend vec operator+ (vec a, vec b);     ///<  a + b
   inline        vec operator+= (vec b);           ///<  this += b
   inline friend vec operator- (vec a, vec b);     ///< a - b
   inline        vec operator-= (vec b);           ///<  this -= b
   inline vec operator- ();                        ///< -a
   inline friend vec operator* (double a, vec b);  ///< mult by number
   inline friend vec operator* (vec b, double a);  ///< mult by number
   inline vec operator*= ( double a);              ///< mult by number
   inline friend vec operator/ (vec b, double a);  ///< div by number
   inline vec operator/= (double a);               ///< div by number
   inline vec fromZto(vec new_axis);               ///< rotate self so that 0,0,1 goes to new_axis
   inline friend vec min (vec a, vec b);           ///< min(a.x,b.x) min(a.y,b.y)  min(a.z,b.z) 
   inline friend vec max (vec a, vec b);           ///< max(a.x,b.x) max(a.y,b.y)  max(a.z,b.z) 
   inline friend vec abs (vec a);                  ///< vec(abs(x),abs(y),abs(z))
   inline friend ostream & operator<< (ostream & o, vec v); ///< print vector
   inline friend double costhetalab (vec k, vec kprim); ///< 
   inline friend double cos (vec k, vec kprim);
   inline friend double thetalab (vec k, vec kprim);
};

////////////////////////////////////////////////////////////////////////////
//             I M P L E M E N T A T I O N                                //  
////////////////////////////////////////////////////////////////////////////

double & vec::operator[] (int i)
  {
    return *(&x + i - 1);
  }
  vec::vec():x(0),y(0),z(0)
  {
  }
  vec::vec (double x1, double x2, double x3):x (x1), y (x2), z (x3)
  {
  }
  vec::vec(string v) 
  {
	stringstream(v)>>x>>y>>z;
  }
/*  vec::vec(const string v) {
      char * str = strdup(v.c_str());
      char * pch;
      pch = strtok (str,", ");
      x=atof(pch);
      pch = strtok(NULL,", ");
      y=atof(pch);
      pch = strtok(NULL,", ");
      z=atof(pch);
}
*/
  double vec::length ()
  {
    return sqrt (x * x + y * y + z * z);
  }
  double vec::norm () const
  {
    return sqrt (x * x + y * y + z * z);
  }
  void vec::normalize () 
  { double d=sqrt (x * x + y * y + z * z);
    if(d==0) return;
    x/=d;
    y/=d;
    z/=d;
  }
  double vec::norm2 () const
  {
    return  (x * x + y * y + z * z);
  }
  vec vec::dir()
  {
     return *this/length();
  }
double operator* (vec a, vec b)
  {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
double angle (vec a, vec b)
  {
    return acos(a*b/sqrt((a*a)*(b*b)));
  }
vec vecprod(vec a, vec b)
  {
    return vec( a.y * b.z - b.y * a.z,
                a.z * b.x - b.z * a.x,
	            a.x * b.y - b.x * a.y);
  }
vec operator+ (vec a, vec b)
  {
    return vec (a.x + b.x, a.y + b.y, a.z + b.z);
  }
vec vec::operator+= (vec a)
  { 
	x+=a.x;
	y+=a.y;
	z+=a.z;
    return *this;
  }
vec operator- (vec a, vec b)
  {
    return vec (a.x - b.x, a.y - b.y, a.z - b.z);
  }
vec vec::operator-= (vec a)
  { 
	x-=a.x;
	y-=a.y;
	z-=a.z;
    return *this;
  }
vec vec::operator- ()
  {
    return vec (-x, -y, -z);
  }
vec operator* (double a, vec b)
  {
    return vec (a * b.x, a * b.y, a * b.z);
  }
vec operator* (vec b, double a)
  {
    return vec (a * b.x, a * b.y, a * b.z);
  }
vec vec::operator*= (double a)
  { 
	x*=a;
    y*=a;
    z*=a;
    return *this;
  }
vec operator/ (vec b, double a)
  {
    return vec (b.x / a, b.y / a, b.z / a);
  }
vec vec::operator/= (double a)
  {
	x/=a;
    y/=a;
    z/=a;
    return *this;
  }
vec vec::fromZto(vec new_axis) ///< rotate self so that 0,0,1 goes to new_axis
  {
 	vec zvec=vec(0,0,1);
    if(new_axis.z<0)
		   zvec.z=-1;
	   else 
		   z*=-1; 
    vec s=(zvec+new_axis.dir()).dir();
    return *this-=2*(*this*s)*s;
  }
vec min (vec a, vec b)
  {
  	return vec(min(a.x,b.x),min(a.y,b.y),min(a.z,b.z));
  }
vec max(vec a, vec b)
  {
  	return vec(max(a.x,b.x),max(a.y,b.y),max(a.z,b.z));
  }
vec abs (vec a)
  {
  	return vec(abs(a.x),abs(a.y),abs(a.z));
  }
ostream & operator<< (ostream & o, vec v)
  {
    return o << '(' << v.x << ',' << v.y << ',' << v.z << ')';
  }
double costhetalab (vec k, vec kprim)
  {
    return k * kprim / sqrt ((k * k) * (kprim * kprim));
  }

double cos (vec k, vec kprim)
  {
    return k * kprim / sqrt ((k * k) * (kprim * kprim));
  }

double thetalab (vec k, vec kprim)
  {
    return acos (k * kprim / sqrt ((k * k) * (kprim * kprim))); 
  }

////////////////////////////////////////////////////////////////////////////
//              I M P L E M E N T A T I O N  -  E N D                     //  
////////////////////////////////////////////////////////////////////////////
#endif
