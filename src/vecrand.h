#ifndef _vecrand_h_
#define _vecrand_h_
#include "vec.h"
#include "jednostki.h"

using namespace std;

double frandom();

inline vec rand_from_ball (double r);         ///< random vector from ball of radius r
inline vec rand_dir ();                       ///< random unit vector in 3D space  
inline vec bodek_rand_from_ball (double r);   ///< punkt ze sfery z ogonem Bodka
inline vec spectral_choice (int p, int n);
inline vec deuterium ();
inline vec rand3(vec A,double b,double c);
inline vec rand_ort (vec a);
//inline vec rand_ort2 (vec a);


////////////////////////////////////////////////////////////////
///            I M P L E M E N T A C J A                     ///  
////////////////////////////////////////////////////////////////

/// Random R such that |R| < r
vec rand_from_ball (double r)
{				//losujemy kierunek punkt z kuli o promieniu r
  double xx, yy, zz, d;
  do
    {
      xx = 2 * frandom () - 1;
      yy = 2 * frandom () - 1;
      zz = 2 * frandom () - 1;
    }
  while ((d = xx * xx + yy * yy + zz * zz) > 1);
  // i mnoymy przez r
  return vec (r * xx, r * yy, r * zz);
}

////////////////////////////////////////////////////////////////
/// momentum of nucleon in Carbon or Oxygen
vec spectral_choice (int p, int n)   //dependence on nucleus (p,n)
{ 
	double pp=0;

	if (p==8 && n==8)//oxygen
	{
		do{
			double x=frandom();
			if (x>=0 && x<0.1)
			pp = 188.06*pow(x,0.32279);

			//1062.09*sqrt(x+304.04*x*x-833.65*x*x*x)-1249.86*log(1+13.2414*x+105.915*x*x);

			if (x>=0.1 && x<0.92)
			pp = -5.141772e+01 
			     +2.854000e+03*x 
			     -2.261280e+04*x*x 
			     +9.965746e+04*x*x*x 
			     -2.469936e+05*x*x*x*x 
			     +3.453132e+05*x*x*x*x*x 
			     -2.535928e+05*x*x*x*x*x*x 
			     +7.605974e+04*x*x*x*x*x*x*x;
			//61.9663+x*326.887-x*x*412.57+x*x*x*360.456;

			if (x>=0.92 && x<=1)
			pp = 0.0019118/(1.00394212-x)/(1.00394212-x) + 25.294/(1.057584-x) - 672.12 + 907.51*x;

			//256.238/(1.0-x+0.17818) + 0.00419236/(1-x+0.0056208)/(1-x+0.0056208) + 952.842 -1726.59*x;

		}
		while (pp>728);
	}

	if (p==6 && n==6)//carbon
	{
		do{
			double x=frandom();
			if (x>=0 && x<0.1)
			pp = 187.47*pow(x,0.33837);

			if (x>=0.1 && x<0.92)
			pp = -5.124758e+01 
			     +2.782654e+03*x 
			     -2.192732e+04*x*x  
			     +9.627319e+04*x*x*x  
			     -2.372355e+05*x*x*x*x  
			     +3.291887e+05*x*x*x*x*x  
			     -2.397192e+05*x*x*x*x*x*x  
			     +7.126416e+04*x*x*x*x*x*x*x;

			if (x>=0.92 && x<=1)
			pp = 0.014243/(1.010629-x)/(1.010629-x) + 212.06/(1.16607-x) + 307.58 - 913.36*x;
		}
		while (pp>688);
	}

  return pp*rand_dir();
}

////////////////////////////////////////////////////////////////
/// momentum of nucleon in deuterium
vec deuterium ( )
{
	double x=frandom();
	double pp = 113.049*sqrt(x) - 129.565*x +166.956*x*x + 3.21967/(1.01-x);
    return pp*rand_dir();
}
////////////////////////////////////////////////////////////////
/// random vec R such that |R|=1
vec rand_dir ()
{
  double xx, yy, zz, d;
  do
    {
      xx = 2 * frandom () - 1;
      yy = 2 * frandom () - 1;
      zz = 2 * frandom () - 1;
    }
  while ((d = xx * xx + yy * yy + zz * zz) > 1 || (d<0.25));
  d=sqrt(d);
  return vec (xx/d, yy/d, zz/d);
}
////////////////////////////////////////////////////////////////
/// random unit vec B orthogonal to A ( |B|=1 and B*A=0 )
vec rand_ort (vec A) 
{				
  double phi=2*Pi*frandom();
  return vec(cos(phi),sin(phi),0).fromZto(A);
}

////////////////////////////////////////////////////////////////
/// momentum of nucleon in nucleus with Bodek "tail" 
/// r = fermi momentum
vec bodek_rand_from_ball (double r)
{
  double xx, yy, zz, d;
  double x=frandom();
  double a=2/GeV;
  if(x<1-6*(r*r*a*a/M_PI/M_PI))
    return rand_from_ball(r);    
  do
    {
      xx = 2 * frandom () - 1;
      yy = 2 * frandom () - 1;
      zz = 2 * frandom () - 1;
      d = xx * xx + yy * yy + zz * zz;
    }
  while (d > 1 || d <0.25 );
  x=frandom()*6*(r*r*a*a/M_PI/M_PI);
  double R=1/(1-r/(4*GeV));
  double S=6*R*a*a/Pi/Pi*r*r*r;
  double p=4*GeV*S/(S+4*GeV*x);
//  cout<<p<<endl;
  d=sqrt(d);   
  return vec (p/d * xx, p/d * yy, p/d * zz);
}

////////////////////////////////////////////////////////////////
/// random vec B such that |B|=b and |B-A|=c
vec rand3(vec A,double b,double c)
{  
	double a=A.length();
	double cgamma=(a*a+b*b-c*c)/2/a/b;
	double sgamma=sqrt(1-cgamma*cgamma);
	return b/a*cgamma*A+b*sgamma*rand_ort(A);
}


#endif
