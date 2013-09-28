#include "scatter.h"
#include "particle.h"
#include <TGenPhaseSpace.h>

/// transate root TLorentzVector to vect  
static inline vect makevect (TLorentzVector * v)
{
  return vect (v->T (), v->X (), v->Y (), v->Z ());
}

/// Scatter p1,p2 into p3,p4  isotropically in cms 
// if f is given returns f(s,q2) for chosen kinematics
double scatter_2 (particle p1, particle p2,
		         particle & p3, particle & p4, 
	             double (*f) (double, double))
{
  vect suma = vect (p1) + vect (p2);
  if(suma*suma<=0) return 0;
  if (!::decay (suma, p3, p4))
	{
	//cerr << "The process is kinematically impossible because"<<endl;
    //cerr << "suma ="<<vect (suma)<<"  p3 ="<< vect (p3)<<"  p4 ="<<vect (p4)<<endl;
    return 0;
    }
  vect p13 = p1 - p3;
  double q2 = p13 * p13;
  if (f)
    return (*f) (suma * suma, q2);
  else
    return q2;
}

/// Scatter p1,p2 into p3,p4 according to distribution 
/// given by f=Ax^3+Bx where x=\cos \theta
bool scatterAB (particle p1, particle p2, 
		 	      particle & p3, particle & p4, 
			      double A, double B, double C, double D, double E, double F, double G, double H)
{
  vect sum = vect (p1) + vect (p2);
  double MM =sum*sum;
  double pp =cms_momentum2(MM, p3.mass2(), p4.mass2());
  if (pp<0 )
  {
//   Invariant mass to small for the reaction
     return 0;// 0 means failure
  }
  vec v = sum.v ();
  vect P1 = p1;
  P1.boost (-v);  // go to CMS
  if (vect(P1).length () == 0)
  {  
    cerr<<"scatterAB: p1.Ek()==0"<<endl;
    return 0; // 0 means failure
  }

  
  double z = get_cos (A, B, C, D, E, F, G, H);
  double r = sqrt (1 - z*z);
  
  double phi=frandom()*2*Pi;
  vec P3=vec(r*cos(phi), r*sin(phi), z)*sqrt (pp);
  
  P3.fromZto(P1);

  p3.set_momentum (P3);
  p4.set_momentum (-P3);

  p3.boost (v);
  p4.boost (v);

  return 1;  // 1 means success
}


/// Scatter p1,p2 into p[0],..,p[n-1] uniformly in phase space
int scatter_n (int n, particle p1, particle p2, 
		    particle p[])
  { 
    vect in = p1 + p2;
    
    double masses[4];
    double sm = 0;
    for (int i = 0; i < n; i++)
      {
	   sm += masses[i] = p[i].mass ();
      }

    if (sm * sm > in * in )
      {
 	  // cerr<< "scatter: to small invariant mass for the reaction"<<endl;
      // cerr<<"in*in="<<in*in<<endl;
      // cerr<<"out*out="<<sm*sm<<endl;
	   return 0;
      }

    TGenPhaseSpace event;
    TLorentzVector v1 (in.x, in.y, in.z, in.t);
    event.SetDecay (v1, n, masses); // give the total momentum, number of particles, 
                                    // and masses of particles in the final state

    double x;
    do	// Probability of each kinematics should be the same
      {
		x = event.Generate ();
      }
    while (x < frandom () * event.GetWtMax ());

    for (int i = 0; i < n; i++)
      {
		p[i].p4 () = makevect (event.GetDecay (i));
		p[i].r = p1.r;
     if(!p[i].is_valid())
	{
	  cerr <<"Scatter: invalid particle:"<< p[i]<<endl;
	  exit(26);
	}
      }
    return 1;
  }


