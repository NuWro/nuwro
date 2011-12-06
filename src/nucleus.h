#ifndef _nucleus_h_
#define _nucleus_h_

#include "particle.h"
#include "pdg.h"
#include "params.h"
#include "generatormt.h"
#include <stack>

//#include <algorithm>  

using namespace std;

/// this class implements the target nucleus
/// it should take into account of the its radius, charge ,
/// proton and neutron density as well as local fermi momentum and Pauli blocking
/// the number of nucleons can change as a result of scattering and the density 
/// and Fermi momentum should be consistently updated - this is still NOT implemented

class nucleus  
{
public:	
  int p;                ///< initial number of protons 
  int n;                ///< initial number of neutrons
  particle *spectator;   ///< nucleon on which absorbtion took place
protected:
  int pr;			    ///< real number of protons 
  int nr;			    ///< real number of neutrons
  double _r;			///< radius
  double Eb;			///< nucleon bounding energy
  double local_kf;		///< local Fermi momentum (if used)
  double kf;			///< global Fermi momentum (if used)
  bool bodek;			///< "bodek tail" i moementum distribution (yes or no)
  double _V;			///< Volume 
  double kf_deuter;     ///< ???????????? 
  int target_model;      ///< Possible values:
						 /// 0 - free target; 
						 /// 1 - Fermi gas; 
						 /// 2 - local Fermi gas; 
						 /// 3 - Bodek; 
						 /// 4 - spectral function; 
						 /// 5 - deuterium; 
						 /// 6 - only proton in deuterium (for tests only);
		 
  vect p4;               ///< minus total fourmomentum of lost nucleons  
  
public:

	nucleus(params &par):
		p   (par.nucleus_p),
		n   (par.nucleus_n),
		Eb  (par.nucleus_E_b*MeV),
		kf  (par.nucleus_kf*MeV),
		target_model(par.nucleus_target),
		spectator(NULL)
    {
		reset1(par);
    }
        
  virtual void reset1(params &par);///< Reset nuclues to a state described by par
  virtual int  A(){return n+p;}                    ///< A initial
  virtual int  Z(){return p;}                      ///< Z initial
  virtual int  N(){return n;}                      ///< N initial
  virtual int  Zr(){return pr;}                    ///< Z real
  virtual int  Nr(){return nr;}                    ///< N real
  virtual int  Ar(){return pr+nr;}                    ///< N real
  virtual int  add_charge(int);                    /// convert neutron to proton or back
  virtual double radius()=0;			           ///< end of nucleus
  virtual double V(){return _V;};	               ///< bounding energy
  virtual double frac_proton ();                   ///< percentage of protons 
  virtual double frac_neutron ();                  ///< percentage of neutrons
  virtual double density (double r)=0;             ///< nucleon density at dist r from center
  virtual double get_random_r ()=0;                ///< random distance from the center 
  virtual void remove_nucleon (particle P);	       ///< remove nucleon from the nucleous
  virtual void insert_nucleon (particle P);	       ///< insert nucleon back to the nucleous
  virtual bool pauli_blocking (particle & pa);     ///< true = the particle is blocked
  virtual bool pauli_blocking (int pdg, double p, double r);	///< true = particle is blocked
  virtual bool pauli_blocking (particle p[], int n);	///< true = any of the particles is blocked
  virtual bool pauli_blocking_old (particle &pa, double ped);  ///< TRUE = PARTICLE IS BLOCKED
  virtual bool pauli_blocking_old (int pdg, double p, double r, double ped);  ///< TRUE = PARTICLE IS BLOCKED
  virtual double mstar (double r);	               ///< nucleon effective mass inside nucleus  
  virtual particle get_nucleon (vec r);            ///< random nucleon located at r
  virtual particle shoot (vec r);                  ///< get nucleon with (random) momentum not exceeding localkf(r)
  virtual particle shoot ();                       ///< get nucleon with random position and momentum not exceeding localkf(r)
  virtual double localkf (particle & pa);          ///< local Fermi momentum for particle (pdg and position dependent)
  virtual double localkf (int pdg, double r);      ///< local Fermi momentum from pdg code and  dist r from nucleus center 
  virtual double kf_from_density (double dens);    ///< local Fermi momentum from density 
};

////////////////////////////////////////////////////////////////////////
///                I m p l e m e n t a t i on
////////////////////////////////////////////////////////////////////////

inline int nucleus::add_charge (int x)   
  { 
	pr+=x;
	nr-=x;
	assert(pr>0);
	assert(nr>0);   
	return x;
  }

////////////////////////////////////////////////////////////////////////
inline double nucleus::localkf (particle & pa)   ///< local Fermi momentum for particle (pdg and position dependent)
  {
    return localkf (pa.pdg, pa.r.length ());
  }

///////////////////////////////////////////////////////////////////////
inline double nucleus::kf_from_density (double dens)
  { static const double C= 0.5 * 3 * Pi * Pi * fermi3;
    assert(dens==dens);  // exclude not a number
    assert(dens>=0);     // exclude negative 
    if(dens<=0) 
     {
      return 0;}
    return pow ( C*dens , 1./3.) * 197.4 * MeV;
  }

///////////////////////////////////////////////////////////////////////
inline  double nucleus::localkf (int pdg, double r)
  {
    if(pr==1 && nr==0) return 0;
    double dens=density(r);
    if(dens<0)
      cout<<"p="<<p<<" n="<<n<<"pr="<<pr<<" nr="<<nr<<" r="<<r<<" dens="<<dens*fermi3<<endl;
    assert(dens==dens);
    assert(p+n>0);
    if(dens==0) return 0;	
	return kf_from_density (dens*(pdg==pdg_proton?pr:nr)/(p+n));
  }

///////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking_old (particle &pa, double ped)    // TRUE = PARTICLE IS BLOCKED
    {
      return pauli_blocking_old(pa.pdg, pa.momentum(), pa.r.length(), ped);
    }

///////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking_old (int pdg, double p, double r, double ped)    // TRUE = PARTICLE IS BLOCKED
    { if(p+n==1 or not nucleon(pdg))
      return false;
	if (target_model==2)
	  return p<localkf(pdg,r); 
      	if (target_model==5)
	 {return p<ped;}
	else
	return p<kf;
	}
///////////////////////////////////////////////////////////////////////
inline particle nucleus::shoot ()
  { 
	  if(p+n==1)	//free nucleon
         {int pdg= ( p ? PDG::pdg_proton : PDG::pdg_neutron);  
	  return particle(pdg,PDG::mass(pdg));
         }
        
 		//if more nucleons in the target
      int pdg=(frandom()*(pr+nr)<pr ? PDG::pdg_proton : PDG::pdg_neutron);
      particle N(pdg,PDG::mass(pdg));
      
      if(p==1 && n==1) //deuterium    
	     {
	     	N.set_momentum(deuterium ()); 
	     	N.r = get_random_r(); 
	     	return N; 
	     }

/* model_target=
 * 0 for free target; 
 * 1 for Fermi gas; 
 * 2 for local Fermi gas; 
 * 3 for Bodek; 
 * 4 for spectral function; 
 * 5 for deuterium; 
 * 6 for only proton in deuterium (for tests only);
 */
	
//	double lkf=localkf(N.pdg, N.r.length());

	switch(target_model)
	{
//	case 0: if(p+n==1) {N.set_momentum(vec(0,0,0));	N.r = vec(0,0,0); return N;}//free target 
//	case 1: N.set_momentum(rand_from_ball(kf));	break; // fermi gas
//	case 2: N.set_momentum(rand_from_ball(localkf(N.pdg, N.r.length())));	break; // local fermi gas
//	case 3: N.set_momentum(bodek_rand_from_ball(kf));	break; //Bodek
//	case 4: N.set_momentum(spectral_choice(p,n));	break; //spectral function 
//	case 5: N.set_momentum(deuterium ()); N.r = get_random_r(); return N; //deuterium 
//	case 6: N.set_momentum(deuterium ()); N.r = get_random_r(); return N; //proton in deuterium
	
//in the future it is possible to add SF for other nuclei as well
	default: return get_nucleon (get_random_r () * rand_dir ());
	}

  }
  
////////////////////////////////////////////////////////////////////////
inline particle nucleus::shoot (vec r)
  {	 
    return get_nucleon (r);	
  }

////////////////////////////////////////////////////////////////////////
inline particle nucleus::get_nucleon (vec r) 
  {
    particle p0;
    p0.r=r;
    if (pr + nr == 1)
      {
		if (pr == 1)
		  p0.set_proton ();
		else
		  p0.set_neutron ();
		p0.set_momentum(vec(p4));
		return p0;
      }

    double x = frandom ();
    if (x * (pr + nr) < pr)
      p0.set_proton ();
    else
      p0.set_neutron ();

    double lkf;

    if (local_kf == 0)
      lkf = kf;
    else
      lkf = localkf (p0.pdg, r.length ());
/* model_target=
 * 0 for free target; 
 * 1 for Fermi gas; 
 * 2 for local Fermi gas; 
 * 3 for Bodek; 
 * 4 for spectral function; 
 * 5 for deuterium; 
 * 6 for only proton in deuterium (for tests only);
 */
//	cout<<"target_model="<<target_model<<endl;
	switch(target_model)
	{
		case 1: p0.set_momentum(rand_from_ball(kf));	break; // fermi gas
		case 2: {
//			     double lkf= localkf(p0.pdg, p0.r.length());
			     p0.set_momentum(rand_from_ball(lkf));	break; // local fermi gas
			     } 
		case 3: p0.set_momentum(bodek_rand_from_ball(kf));	break; //Bodek
		case 4: p0.set_momentum(spectral_choice(p,n));	break; //spectral function 
		default: p0.set_momentum(rand_from_ball(localkf(p0.pdg, p0.r.length())));   //local fermi gas
	}
	//account for the fact that the nucleus is in motion
	p0.set_momentum(p0.p()+vec(p4)/Ar()); ////???????
	if(not (p0.v2()<1))
	cerr<<target_model<<" p="<<this->p<<" n="<<this->n<<" lkf="<<lkf<<" v2="<< p0.v2()<<" "<<p0.v()<<endl;
	//assert(p0.v2()<1);
    return p0;
  }

////////////////////////////////////////////////////////////////////////
inline void nucleus::reset1(params &par)
    {
       pr=p=par.nucleus_p;
       nr=n=par.nucleus_n;	 
       Eb=par.nucleus_E_b*MeV;
       kf=par.nucleus_kf*MeV;
  	   target_model=par.nucleus_target;
  	   spectator=NULL;
  	   
       if(p+n==1) 
          Eb=kf=_r=0; /// ?????? 	    
       ///Potential -- 
       using namespace PDG;
       _V = sqrt(mass_proton*mass_proton + kf*kf) - mass_proton + 7*MeV;
    }
 
////////////////////////////////////////////////////////////////////////
inline double nucleus::frac_proton ()  ///< percentage of protons 
  { 
	  return double (pr) / (pr + nr); 
  }

////////////////////////////////////////////////////////////////////////  
inline double nucleus::frac_neutron () ///< percentage of neutrons
  { 
	  return double (nr) / (pr + nr); 
  }

///////////////////////////////////////////////////////////////////////////////
// check if particle is Pauli blocked - depends on localkf at particle position
///////////////////////////////////////////////////////////////////////////////
inline  bool nucleus::pauli_blocking (particle & pa) 
  {
    return pauli_blocking (pa.pdg, pa.momentum (), pa.r.length ());
  }

///////////////////////////////////////////////////////////////////////////////
// check if any of p[0]...p[n-1] is Pauli blocked 
///////////////////////////////////////////////////////////////////////////////
inline  bool nucleus::pauli_blocking (particle p[], int n)
  {
    for (int i = 0; i < n; i++)
      if (pauli_blocking (p[i]))
        return true;
    return false;
  }

///////////////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking (int pdg, double p, double r)    // TRUE = PARTICLE IS BLOCKED
  {
	if(pr+nr==1 or not nucleon(pdg))   return false;
	if(local_kf==0)                    return p<kf;
	else                               return p<localkf(pdg,r);
  }

///////////////////////////////////////////////////////////////////////////////
inline double nucleus::mstar (double r)
  {				// na razie nie jet to zbyt mądre
    return mass_proton;
  }

///////////////////////////////////////////////////////////////////////////////
/// remove nucleon P from the nucleus 
///////////////////////////////////////////////////////////////////////////////
inline void nucleus::remove_nucleon(particle P)
  {
   switch(P.pdg)
    {case pdg_proton:  //if(pr==0) cout<<proc<<endl; 
             assert(pr>0); pr--;break;
     case pdg_neutron: //if(nr==0) cout<<proc<<endl; 
			 assert(nr>0);  nr--;break;
     default: return;
    }
    p4=p4-P.p4();
  }

///////////////////////////////////////////////////////////////////////////////
/// insert nucleon P into the nucleus 
///////////////////////////////////////////////////////////////////////////////
inline void nucleus::insert_nucleon(particle P)
  { switch(P.pdg)
    {case pdg_proton: pr++;break;
     case pdg_neutron: nr++;break;
     default: return;
    }
    p4=p4+P.p4();
  }

#endif
