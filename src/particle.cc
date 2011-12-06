#include "particle.h"
#define echo(x)

/// particle decay;
/// p4 is totalfourmomentum 
/// p1, p2 particles to decay into (their masses must be set beforehand)
bool
decay (vect p4, particle & p1, particle & p2)
{
  double MM = p4 * p4;
  double Msum=p1.mass () + p2.mass ();
  if (Msum*Msum >= MM)
    return false;
  // random vector form the unit ball 
  double xx, yy, zz, rr;
  do
    {
      xx = 2 * frandom () - 1;
      yy = 2 * frandom () - 1;
      zz = 2 * frandom () - 1;
    }
  while ((rr = xx * xx + yy * yy + zz * zz) > 1);
  // value of momentum of decay products in cms frame
//  double m1m1 = p1.mass () * p1.mass ();
//  double m2m2 = p2.mass () * p2.mass ();
  double pp = cms_momentum2(MM,p1.mass2(),p2.mass2());
//    (MM * MM + m1m1 * m1m1 + m2m2 * m2m2 -
//     2 * (MM * m1m1 + MM * m2m2 + m1m1 * m2m2)) / (4 * MM);
  double k = sqrt (pp / rr);
  // adjust the length of random vector to the calculated value
  vec p (k * xx, k * yy, k * zz);
//     cout<<"p="<<p<<endl;
  p1.set_momentum (p);
  p2.set_momentum (-p);
  vec cmsspeed = p4.v ();
  // go back to the lab frame - in which initial particle has velocity p4.v();
  p1.boost (cmsspeed);
  p2.boost (cmsspeed);
  return true;
}

///  speed of a in the frame of b
double
relative_speed (particle & a, particle & b)
{
  vect plab = a.p4 ();
  plab.boost (-b.v ());
  return plab.v ().length ();
}

/// energy of a in the frame of b
double Ek_in_frame (particle & p, vec v)
{
  vect plab = p.p4 ();
  return plab.boost (-v).t-p.mass();
}

///  
double
get_cos (double a, double b, double c, double d, double e, double f, double g, double h)
{		
	double max = fabs(a) + fabs(b) + fabs(c) + fabs(d) + fabs(e) + fabs(f) + fabs(g) + fabs(h);
	double x, x2, x3, x4, x5, x6, x7;
		
	do
	{
		x = -1.0 + 2.0*frandom();
		x2 = x*x;
		x3 = x*x2;
		x4 = x*x3;
		x5 = x*x4;
		x6 = x*x5;
		x7 = x*x6;
	
	}while (a*x7 + b*x6 + c*x5 + d*x4 + e*x3 + f*x2 + g*x + h < max*frandom());
		
	return x;
}

bool decay2 (vect k,vect p4, particle & p1, particle & p2,double &coef)
{
  double pp=cms_momentum2(p4*p4,p1.mass2(),p2.mass2());
  if(pp<0) return false;     
  double x=frandom();
  double z;
  if(true)
  {
	z=1-2*x*x;
	coef=4*x;
  }
/*  if(false)
  { 
	  z=1-2*x;
	  coef=2;
  }  
  if(false)
  {
    z=cos(x*Pi);
    coef=sin(x*Pi)*Pi;
  }
*/
  double phi=frandom()*2*Pi;
  double st=sqrt(1-z*z);

  vec p(st*cos(phi),st*sin(phi),z);
  vec cmsspeed = p4.v();
  k.boost(-cmsspeed);
  p.fromZto(k);  
  p*=sqrt (pp);
//  cout<<p<<' '<<theta<<endl;
  p1.set_momentum (p);
  p2.set_momentum (-p);
  // go back to the lab frame 
  p1.boost (cmsspeed);
  p2.boost (cmsspeed);
  return true;
}

