#ifndef _nucleus_h_
#define _nucleus_h_

#include "particle.h"
#include "pdg.h"
#include "params.h"
#include "generatormt.h"
#include "nucleus_data.h"
#include "isotopes.h"
#include <stack>


#include "MEC_data/C12/C12_pp.h"
#include "MEC_data/C12/C12_np.h"
#include "MEC_data/C12/C12_pn.h"
#include "MEC_data/C12/C12_3p3h.h"
#include "MEC_data/O16/O16_pp.h"
#include "MEC_data/O16/O16_np.h"
#include "MEC_data/O16/O16_pn.h"
#include "MEC_data/O16/O16_3p3h.h"
#include "MEC_data/Ca40/Ca40_pp.h"
#include "MEC_data/Ca40/Ca40_np.h"
#include "MEC_data/Ca40/Ca40_pn.h"
#include "MEC_data/Ca40/Ca40_3p3h.h"

using namespace std;



// A class to load Valencia 2020 model hadronic tables Phys. Rev. C 102 (2020) 024601 J.E. Sobczyk
// Class Implemented by Hemant Prasad on 23 / 11 / 2024
class Response_function {
  public:
    double *_channel_pp;
    double *_channel_np;
    double *_channel_pn;
    double *_channel_3p3h;
    Response_function() {};
    ~Response_function() {};

    void Attach_ptr(double* ptr1, double* ptr2, double* ptr3, double* ptr4)   {
      _channel_pp = ptr1;
      _channel_np = ptr2;
      _channel_pn = ptr3;
      _channel_3p3h = ptr4;
    }
    double channel(int l, int k){
      double value=0;
      //int k = (int)( 0.5*(240 - n + 1) + (m-1) );
      switch(l) {
        case 0: value=_channel_pp[k];
                break;
        case 1: value=_channel_np[k];
                break;
        case 2: value=_channel_pn[k];
                break;
        case 3: value=_channel_3p3h[k];
                break;
      }
      return value;
    }

};

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

        int MAX_EVENT_REDRAWS = 10;
        double effective_mfp_scale = 0.25;
        double cosine_threshold = 0.001;

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
	int kMomDist;   	  ///< Type of nucleon momentum distribution:
						  /// 0 - free nucleon (forced for H1);
						  /// 1 - Fermi gas (mean Fermi momentum);
						  /// 2 - local Fermis gas;
						  /// 3 - "bodek tail" in momentum distribution
						  /// 4 - effective spectral function (for carbon and oxygen, else lFG);
						  /// 5 - deuterium (forced for H2);
              /// 6 - effective potential

	// Hyperon potential at r=0
	double Lambda_Eb;
	double Sigma_Eb;
  Response_function W00;          /// Component of hadronic tensor in (w,q)
  Response_function W03;          /// Component of hadronic tensor in (w,q)
  Response_function W11;          /// Component of hadronic tensor in (w,q)
  Response_function W12;          /// Component of hadronic tensor in (w,q)
  Response_function W33;          /// Component of hadronic tensor in (w,q)

	public:

	nucleus(params &par);                    ///< construct nucleus from params
	void reset(){pr=p;nr=n;_p4=vect();}      ///< recover  initial state
	int  A(){return n+p;}                    ///< A initial
	int  Z(){return p;}                      ///< Z initial
	int  N(){return n;}                      ///< N initial
	int  Zr(){return pr;}                    ///< Z real
	int  Nr(){return nr;}                    ///< N real
	int  Ar(){return pr+nr;}                 ///< N real
	double radius(){return _r;}	             ///< radius of nucleus
	double r(){return _r;}	                 ///< radius of nucleus
	double V(particle &p);         			 ///< potential of nucleon p
	double frac_proton ();                   ///< percentage of protons
	double frac_neutron ();                  ///< percentage of neutrons
	double density (double r);               ///< nucleon density at dist r from center
	double get_random_r ();                  ///< random distance from the center
	bool remove_nucleon (particle P);	     ///< remove nucleon from the nucleus
	void insert_nucleon (particle P);	     ///< insert nucleon back to the nucleus
	double localkf (particle & pa);          ///< local Fermi momentum for particle (pdg and position dependent)
	double localkf_ (int pdg, double r);      ///< local Fermi momentum from pdg code and  dist r from nucleus center
	double kF(){return _kf;}	             ///< global Fermi momentum
	double Ef (particle & pa);		 		 ///< nucleon Fermi energy dependent on kMomDist
	double Mf();			                 ///< nucleon effective mass inside nucleus
	double Ef();                             ///< nucleon Fermi energy
	double Eb(){return _Eb;}                 ///< nucleon binding energy (from experimantal data)
	particle get_nucleon (vec r);            ///< random nucleon located at r (used in Interaction.cc)
	particle get_nucleon ();                 ///< random nucleon (used in nuwro.cc)
	bool pauli_blocking (particle & pa);     ///< true if the particle is blocked
	bool pauli_blocking (particle p[], int n); ///< true if any of p[0]..p[n-1] is blocked (used in cascade)
	bool pauli_blocking_old (particle &pa, double p0); ///< TRUE = PARTICLE IS BLOCKED (p0 - momentum of initial nucleon)

	double hyp_BE(double r,int pdg);

  double GetValue_W00(int l, int k) {return W00.channel(l,k);}
  double GetValue_W03(int l, int k) {return W03.channel(l,k);}
  double GetValue_W11(int l, int k) {return W11.channel(l,k);}
  double GetValue_W12(int l, int k) {return W12.channel(l,k);}
  double GetValue_W33(int l, int k) {return W33.channel(l,k);}
};

////////////////////////////////////////////////////////////////////////
///                I m p l e m e n t a t i on
////////////////////////////////////////////////////////////////////////

inline double nucleus::Ef()
{
	double const M=0.5*(PDG::mass_proton+PDG::mass_neutron);

	return sqrt(_kf*_kf+M*M)-M;
}

inline double kf_from_density (double dens)
{
	static const double C= 3 * Pi * Pi * fermi3;
	if(dens>0)
		return cbrt( C*dens) * 197.4 * MeV;
	else
		return 0;
}

////////////////////////////////////////////////////////////////////////
inline double nucleus::localkf (particle & pa)   ///< local Fermi momentum for particle (pdg and position dependent)
{
	return localkf_ (pa.pdg, pa.r.length ());
}

///////////////////////////////////////////////////////////////////////
inline  double nucleus::localkf_ (int pdg, double r)
{
	if(pr+nr==1)
		return 0;
	double dens=density(r);
	assert(dens==dens);
	assert(p+n>0);
	if(dens==0)
		return 0;
	return kf_from_density (dens*(pdg==pdg_proton?pr:nr)/(pr+nr));
}

///////////////////////////////////////////////////////////////////////
inline bool nucleus::pauli_blocking_old (particle &pa, double p0)    // TRUE = PARTICLE IS BLOCKED
{
	if(p+n==1 or not pa.nucleon())
		return false;
	switch(kMomDist)
	{
		case 0:	 return false;
		case 2:	 return pa.momentum()<localkf(pa); // Local Fermi Gas
		case 4:	 return pa.momentum()<localkf(pa); // Local Fermi Gas
		case 5:	 return pa.momentum()<p0;                            // deuteron
		default: return pa.momentum()<_kf;                           // Global Fermi Gas
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
// check if particle is Pauli blocked - depends on localkf at particle position
///////////////////////////////////////////////////////////////////////////////
inline  bool nucleus::pauli_blocking (particle & pa)
{
	if(pr+nr==1 or not pa.nucleon())
		return false;
	switch(kMomDist)
	{
		case 0:
			return false;
		case 2:
		case 4:
			return pa.momentum()<localkf(pa);// Local FG
		default:
			return pa.momentum()<_kf;                           // Global FG
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
inline bool nucleus::remove_nucleon(particle P)
{
	switch(P.pdg)
    {
		case pdg_proton:  //if(pr==0) cout<<proc<<endl;
			if(pr<=0)
				return false;
			pr--;
			break;
		case pdg_neutron: //if(nr==0) cout<<proc<<endl;
			if(nr<=0)
				return false;
			nr--;break;
		default: return false;
    }
    _p4=_p4-P.p4();
    return true;
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

///////////////////////////////////////////////////////////////////////////////
/// potential of a hyperon
///////////////////////////////////////////////////////////////////////////////
inline double nucleus::hyp_BE(double r,int pdg)
{
	switch(pdg)
	{
		case 3122: return Lambda_Eb*density(r)/density(0);
		case 3212: case 3112: case 3222: return Sigma_Eb*density(r)/density(0);
		default: std::cout << "This is not a hyperon! Returning 0 for hyperon potential" << std::endl; return 0;
	}

	std::cout << "Should never reach here!" << std::endl;
	return 0;
}

#endif
