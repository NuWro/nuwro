
#include "scatter.h"
#include "particle.h"
#include <TGenPhaseSpace.h>
#include <iomanip>

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

    double masses[5];
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


// NOT CURRENTLY USED
// Solve kinematics with hyperon BE using self consistency method:

// Boost initial particles into CMS
// Generate a scattering angle for the hyperon
// Guess an effective mass for the hyperon
// Boost back to lab
// Check if eff mass is consistent with on shell hyperon modified for binding energy
// If not, guess a new value and repeat until they match

double scatter2_with_BE_SC(particle p1,particle p2, particle &p3, particle &p4, double Y_Eb){

  double range = 0.01; // Maximum allowed disagreement between two calculations of effective mass

  // Compute the 4 momentum of cms
  vect suma = vect (p1) + vect (p2);
  if(suma*suma<=0) return 0;

  double s = suma*suma; // Squared cms energy
  double pcms2=0;

  particle hyp = p4;
  particle lep = p3;

  vec vcms = suma.v();

  // Make initial guess of hyperon effective mass as nucleon mass minus lepton mass
  // this gaurentees this initial attempt will always succeed

  hyp.set_mass(sqrt(s) - lep.mass()-3);

  // Perform the decay IN THE CMS using this initial guess

  ::decay(vect(p1)+vect(p2), hyp, lep);

  // Get the direction of the 3 momentum of the hyperon in the cms
  vec unit_vec = vec(hyp) / vec(hyp).length();

  // Boost to lab

  hyp.boost(vcms);
  lep.boost(vcms);

  double diff=0;

  // Cap at 50 iterations (usually finds eff mass in <5) 
  for(int i=0;i<50;i++){

     // Check is hyperon's 4 momentum consistent with BE adjustment

     double calc_energy = hyp.t;
     double target_energy =  sqrt(p4.mass()*p4.mass() + hyp.length()*hyp.length()) - Y_Eb;

     // Break the loop when the values are close enough or if step size falls below target size
     if(calc_energy < target_energy + range/2 && calc_energy > target_energy - range/2) break; 

     diff = calc_energy - target_energy; // Difference between two energy calculations - when this is zero you have correct effective mass

     // Boost back to cms
     hyp.boost(-vcms);
     lep.boost(-vcms);

     double hyp_mass = hyp.mass();

     label:

     // Get new estimate of hyperon effective mass

     hyp.set_mass(hyp_mass - diff);

     // Recalculate the cms momentum
     // If this fails then pcms2=-1

     pcms2 = cms_momentum2(s,lep.mass2(),hyp.mass2());

     // If unable to generate kinematics for this 
     if(pcms2 < 0 && abs(diff) > range/2){ 

        diff /= 2;
        goto label;

     }
     else if(abs(diff)<range/2) return 0;
 
     hyp.set_momentum(unit_vec * sqrt(pcms2) );
     lep.set_momentum(-unit_vec * sqrt(pcms2) );

     // Boost back to lab
     hyp.boost(vcms);
     lep.boost(vcms);


  }


  // Having found the correct effective mass, set
  // the outgoing particle 4 momenta and masses

  double hyp_mass = p4.mass();

  p3 = lep;
  p4 = hyp;

  double q2 = (lep.t - p1.t)*(lep.t - p1.t) - (vec(p1)-vec(lep))*(vec(p1)-vec(lep));

  return q2;
}


// NOT CURRENTLY USED
// Modified cms energy trick

double scatter2_with_BE(particle p1,particle p2, particle &p3, particle &p4, double Y_Eb, vec &cms_dir){

  particle hyp = p4;
  particle lep = p3;

  // Compute the 4 momentum of cms blob
  vect suma = vect (p1) + vect (p2);
  if(suma*suma<=0) return 0;

  // To boost to cms
  vec vcms = suma.v();

  // Generate random unit vector in true cms

  double xx,yy,zz,rr;

 do
    {
      xx = 2 * frandom () - 1;
      yy = 2 * frandom () - 1;
      zz = 2 * frandom () - 1;
      
      
    }
 while ((rr = sqrt(xx * xx + yy * yy + zz * zz)) > 1);

 rr=sqrt(xx*xx+yy*yy+zz*zz);

 vec dir(xx/rr,yy/rr,zz/rr);

 vect true_cms_dir(1,dir);

 // Boost this vector to the lab frame
 true_cms_dir.boost(vcms);

 // Take the total initial 4 momentum and add the hyperon BE
 vect suma_be = suma;
 suma_be.t += Y_Eb;

 // To boost to modified cms
 vec vcms_be = suma_be.v();

 // Boost direction vector into modified cms
 true_cms_dir.boost(-vcms_be);

 // Try to decay the CMS

 double s = suma_be*suma_be; // Squared cms energy
 double pp = cms_momentum2(s,p3.mass2(),p4.mass2());


 if(pp <= 0 ) return 0; // Process not allowed by kinematics

 lep.set_momentum( vec(true_cms_dir)/vec(true_cms_dir).length() * sqrt(pp)  );
 hyp.set_momentum( - vec(true_cms_dir)/vec(true_cms_dir).length() * sqrt(pp)  );

 lep.t = sqrt(lep.mass()*lep.mass() + pp );

 hyp.t = sqrt(hyp.mass()*hyp.mass() + pp );


 // Record direction of outgoing lepton in cms for use in rescale_momenta
 cms_dir = dir;

 // Boost back to lab
 lep.boost(vcms_be);
 hyp.boost(vcms_be);

 // Subtract BE off hyperon energy 
 hyp.t -= Y_Eb;


 // Calculate the q2

 double q2 = (lep.t - p1.t)*(lep.t - p1.t) - (vec(p1)-vec(lep))*(vec(p1)-vec(lep));

 p3 = lep;
 p4 = hyp;

 return q2;


}

// Takes total 4 momentum pcms, generates kinematics for particles
// p3 and p4 with p3 travelling along cms_dir in the CMS

//bool rescale_momenta(vect pcms, vec cms_dir, particle &p3, particle &p4,double Y_Eb){
bool rescale_momenta(vect pcms,vec cms_dir,particle &p3, particle &p4){

   // Just to make sure
   cms_dir *= 1.0/cms_dir.length();

   vec vcms_be = pcms.v();

   double s = pcms*pcms;

   double pp = cms_momentum2(s,p3.mass2(),p4.mass2());

   if(pp < 0) return false; // Process not allowed by kinematics

   p3.set_momentum( cms_dir * sqrt(pp) );
   p4.set_momentum( -cms_dir * sqrt(pp) );

   p3.t = sqrt(p3.mass()*p3.mass() + pp );
   p4.t = sqrt(p4.mass()*p4.mass() + pp );

   p3.boost(vcms_be);
   p4.boost(vcms_be);

   /*
   // Check 4 momentum
      std::cout << std::endl;
      std::cout << pcms << std::endl;
      std::cout << vect(p3) << std::endl;
      std::cout << vect(p4) << std::endl;
   */

   return true;

}
