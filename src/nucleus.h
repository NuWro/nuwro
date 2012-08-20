#ifndef _nucleus_h_
#define _nucleus_h_

#include "particle.h"
#include "pdg.h"
#include "params.h"
#include "generatormt.h"
#include "nucleus_data.h"
#include "isotopes.h"
#include <stack>

using namespace std;

/// this class implements the target nucleus
/// it should take into account of the its radius, charge ,
/// proton and neutron density as well as local fermi momentum and Pauli blocking
/// the number of nucleons can change as a result of scattering and the density 
/// and Fermi momentum should are consistently updated - this is still NOT implemented

class nucleus  
{
	public:	

	int p;                ///< initial number of protons 
	int n;                ///< initial number of neutrons
	particle *spectator;  ///< nucleon on which absorbtion took place

	private:

	nucleus_data *d;      ///< density profiles from experimental data
	public:
	isotope      *i;      ///< isotope data from: "The Ame2003 atomic mass evaluation (II)"  by G.Audi, A.H.Wapstra and C.Thibault
					      ///< Nuclear Physics A729 p. 337-676, December 22, 2003.
	
	int pr;			      ///< real number of protons 
	int nr;			      ///< real number of neutrons
	vect _p4;             ///< minus total fourmomentum of lost nucleons  
	double _r;			  ///< nucleus radius
	double _Eb;			  ///< binding energy per nucleon (from exp data)
	double _kf;			  ///< global Fermi momentum
	double _V;			  ///< nucleon potential (Fermi energy) 
	int kMomDist;   	  ///< Type of nucleon momentum distribution:
						  /// 0 - free nucleon (forced for H1); 
						  /// 1 - Fermi gas (mean Fermi momentum); 
						  /// 2 - local Fermis gas; 
						  /// 3 - "bodek tail" in momentum distribution
						  /// 4 - spectral function (for carbon and oxygen, else lFG); 
						  /// 5 - deuterium (forced for H2); 
						  /// 6 - proton in deuterium (for tests only);
		 
	public:

	nucleus(params &par);                    ///< construct nucleus from params
	int  A(){return n+p;}                    ///< A initial
	int  Z(){return p;}                      ///< Z initial
	int  N(){return n;}                      ///< N initial
	int  Zr(){return pr;}                    ///< Z real
	int  Nr(){return nr;}                    ///< N real
	int  Ar(){return pr+nr;}                 ///< N real
	double radius(){return _r;}	             ///< radius of nucleus
	double r(){return _r;}	                 ///< radius of nucleus
	double V(){return _V;};	                 ///< potential of nucleons
	double V(particle &p);         			 ///< potential of nucleon p
	double frac_proton ();                   ///< percentage of protons 
	double frac_neutron ();                  ///< percentage of neutrons
	double density (double r);               ///< nucleon density at dist r from center
	double get_random_r ();                  ///< random distance from the center 
	void remove_nucleon (particle P);	     ///< remove nucleon from the nucleous
	void insert_nucleon (particle P);	     ///< insert nucleon back to the nucleous
	double localkf (particle & pa);          ///< local Fermi momentum for particle (pdg and position dependent)
	double localkf (int pdg, double r);      ///< local Fermi momentum from pdg code and  dist r from nucleus center 
	double kF(){return _kf;}	             ///< global Fermi momentum
	double Mf();			                 ///< nucleon effective mass inside nucleus  
	double Ef();                             ///< nucleon Fermi energy
	double Eb(){return _Eb;}                 ///< nucleon binding energy (from experimantal data)
	particle get_nucleon (vec r);            ///< random nucleon located at r (used in Interaction.cc)
	particle get_nucleon ();                 ///< random nucleon (used in nuwro.cc)
	bool pauli_blocking (particle & pa);     ///< true = the particle is blocked
	bool pauli_blocking (int pdg, double p, double r);	                ///< true = particle is blocked
	bool pauli_blocking (particle p[], int n);	                        ///< true = any of the particles is blocked
	bool pauli_blocking_old (particle &pa, double ped);                 ///< TRUE = PARTICLE IS BLOCKED
	bool pauli_blocking_old (int pdg, double p, double r, double ped);  ///< TRUE = PARTICLE IS BLOCKED
};

////////////////////////////////////////////////////////////////////////
///                I m p l e m e n t a t i on
////////////////////////////////////////////////////////////////////////


inline double nucleus::Ef()
{
	double const M=0.5*(PDG::mass_proton+PDG::mass_neutron);
	return sqrt(_kf*_kf+M*M)-M+_Eb;
}


inline double kf_from_density (double dens)
{ 
	static const double C= 1.5 * Pi * Pi * fermi3;
	if(dens>0)
		return cbrt( C*dens) * 197.4 * MeV;
	else
		return 0;
}


////////////////////////////////////////////////////////////////////////
inline double nucleus::localkf (particle & pa)   ///< local Fermi momentum for particle (pdg and position dependent)
{
	return localkf (pa.pdg, pa.r.length ());
}


///////////////////////////////////////////////////////////////////////
inline  double nucleus::localkf (int pdg, double r)
{
	if(pr+nr==1) 
		return 0;
	double dens=density(r);
	assert(dens==dens);
	assert(p+n>0);
	if(dens==0) 
		return 0;	
	return kf_from_density (dens*(pdg==pdg_proton?pr:nr)/(p+n));
//	return kf_from_density (2*dens*(pdg==pdg_proton?pr:nr)/(p+n));
//	return FermiMomentum (2*dens*(pdg==pdg_proton?pr:nr)/(p+n));
}

///////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking_old (particle &pa, double ped)    // TRUE = PARTICLE IS BLOCKED
{
	return pauli_blocking_old(pa.pdg, pa.momentum(), pa.r.length(), ped);
}

///////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking_old (int pdg, double p, double r, double ped)    // TRUE = PARTICLE IS BLOCKED
{
	if(p+n==1 or not nucleon(pdg))
		return false;
	switch(kMomDist)
	{
		case 2:	 return p<localkf(pdg,r); 
		case 5:	 return p<ped;
		default: return p<_kf;
	}
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
	if(pr+nr==1 or not nucleon(pdg))   
		return false;
	switch(kMomDist)
	{
		case 2:
			return p<localkf(pdg,r);
		default:
			return p<_kf;
		
	}
}

///////////////////////////////////////////////////////////////////////////////
inline double nucleus::Mf ()
{				
	if(d)
		return d->Mf();
	else
		return Meff(_kf);
}

///////////////////////////////////////////////////////////////////////////////
/// remove nucleon P from the nucleus 
///////////////////////////////////////////////////////////////////////////////
inline void nucleus::remove_nucleon(particle P)
{
	switch(P.pdg)
    {
		case pdg_proton:  //if(pr==0) cout<<proc<<endl; 
			assert(pr>0); pr--;break;
		case pdg_neutron: //if(nr==0) cout<<proc<<endl; 
			assert(nr>0);  nr--;break;
		default: return;
    }
    _p4=_p4-P.p4();
}

///////////////////////////////////////////////////////////////////////////////
/// insert nucleon P into the nucleus 
///////////////////////////////////////////////////////////////////////////////
inline void nucleus::insert_nucleon(particle P)
{ 
	switch(P.pdg)
    {
		case pdg_proton:  pr++; break;
		case pdg_neutron: nr++; break;
		default: return;
    }
    _p4=_p4+P.p4();
}

///////////////////////////////////////////////////////////////////////////////
/// potential of a nucleon
///////////////////////////////////////////////////////////////////////////////
inline double nucleus::V(particle &p)
{
	double lkf = localkf(p);
	return (sqrt(p.mass2() + lkf*lkf) - p.mass());// + 7*MeV);// + 5*MeV;
}

#endif
