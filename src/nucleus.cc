#include "nucleus.h"
#include "nucleus_data.h"
#include "particle.h"

#include <cassert>
   
using namespace std;

///////////////////////////////////////////////////////////////////////////////
/// Construct nucleus from params given in the input file 
///////////////////////////////////////////////////////////////////////////////
nucleus::nucleus(params &par):
	p   (par.nucleus_p),
	n   (par.nucleus_n),
	spectator(NULL),
	d(NULL),
	i(NULL),
	_Eb  (par.nucleus_E_b*MeV),
	_kf  (par.nucleus_kf*MeV),
	kMomDist(par.nucleus_target)
{
	using namespace PDG;
	pr=p;
	nr=n;	 
	if(p+n==1) 
		_Eb=_kf=_r=0; /// ?????? 	    
    i=isotope_find(p,n);     
    if(i->Z==p && i->N==n)
		_Eb=i->binding_energy*keV;
	///Potential -- 
	using namespace PDG;
	
	if(par.nucleus_model==1)
	{	
		d=best_data(p,n);
		_r=d->r() * cbrt( A() / d->A() );
		_kf=d->kF();
	}
	else
	{
		_r=cbrt(p+n)*1.2*fermi;
	}
	if((par.target_type==0 || par.target_type==1)  and par.nucleus_kf>0)
	{
		_kf=par.nucleus_kf;
	}
	//set hyperon binding energy
	Lambda_Eb = par.hyp_Lambda_Eb;
	Sigma_Eb = par.hyp_Sigma_Eb;


  //Attaching appropriate hadronic grids according to the nucleus configuration 
  if (p > 12) {  // Calcium grid 
    W00.Attach_ptr(Ca40_W00pp, Ca40_W00np, Ca40_W00pn, Ca40_W00_3p3h);
    W03.Attach_ptr(Ca40_W03pp, Ca40_W03np, Ca40_W03pn, Ca40_W03_3p3h);
    W11.Attach_ptr(Ca40_W11pp, Ca40_W11np, Ca40_W11pn, Ca40_W11_3p3h);
    W12.Attach_ptr(Ca40_W12pp, Ca40_W12np, Ca40_W12pn, Ca40_W12_3p3h);
    W33.Attach_ptr(Ca40_W33pp, Ca40_W33np, Ca40_W33pn, Ca40_W33_3p3h);
  } else if ((p > 6) and (p<=12)) { // Oxygen grid
      W00.Attach_ptr(O16_W00pp, O16_W00np, O16_W00pn, O16_W00_3p3h);
      W03.Attach_ptr(O16_W03pp, O16_W03np, O16_W03pn, O16_W03_3p3h);
      W11.Attach_ptr(O16_W11pp, O16_W11np, O16_W11pn, O16_W11_3p3h);
      W12.Attach_ptr(O16_W12pp, O16_W12np, O16_W12pn, O16_W12_3p3h);
      W33.Attach_ptr(O16_W33pp, O16_W33np, O16_W33pn, O16_W33_3p3h);
    } else {    // Carbon Grid
      W00.Attach_ptr(C12_W00pp, C12_W00np, C12_W00pn, C12_W00_3p3h);
      W03.Attach_ptr(C12_W03pp, C12_W03np, C12_W03pn, C12_W03_3p3h);
      W11.Attach_ptr(C12_W11pp, C12_W11np, C12_W11pn, C12_W11_3p3h);
      W12.Attach_ptr(C12_W12pp, C12_W12np, C12_W12pn, C12_W12_3p3h);
      W33.Attach_ptr(C12_W33pp, C12_W33np, C12_W33pn, C12_W33_3p3h);
    }
}

double nucleus::density (double r)
{	
	static double const Vol0=4.0/3*Pi*pow(1.2*fermi,3);		

	if(d)
		return max(d->dens(r/_r*d->r())*Ar()/A(),0.);
	else
		return r < _r ? 1/Vol0*Ar()/A() : 0;
}
  
///////////////////////////////////////////////////////////////////////////////
/// random distance from the center weighted with nuclear density
///////////////////////////////////////////////////////////////////////////////
double nucleus::get_random_r()
{
	if(d)
		return _r*(d->random_r()/d->r());
	else
		return _r*frandom_sqr();
} 
///////////////////////////////////////////////////////////////////////
particle nucleus::get_nucleon ()
{ 
	if(p+n==1)	//free nucleon
	{
		int pdg= ( p==1 ? PDG::pdg_proton : PDG::pdg_neutron);  
		return particle(pdg,PDG::mass(pdg)); //r=(0,0,0) //p=(0,0,0)
	}

	if(p==1 && n==1) //deuterium    
	{
		int pdg=(frandom()<frac_proton() ? PDG::pdg_proton : PDG::pdg_neutron);
		particle N(pdg,PDG::mass(pdg));
		N.set_momentum(deuterium ()); 
		N.r = get_random_r(); 
		return N; 
	}

	//if more nucleons in the target
	return get_nucleon (get_random_r () * rand_dir ());
}
 

////////////////////////////////////////////////////////////////////////
particle nucleus::get_nucleon (vec r) 
{
	particle p0;
	p0.r=r;
	if(nr==0 || (pr!=0 && frandom () < frac_proton()))
		p0.set_proton ();
	else
		p0.set_neutron ();
	//~ if(pr+nr==1)
	//~ {
		//~ p0.set_momentum(vec(_p4));
		//~ return p0;
	//~ }
	
	if(p==1 and n==1)
	{
		p0.set_momentum(deuterium());
		return p0;
	}

	switch(kMomDist)
	{
		case  1: p0.set_momentum(rand_from_ball(_kf)); break; // global fermi gas
		case  2: p0.set_momentum(rand_from_ball(localkf(p0)));break; // local fermi gas		
		case  3: p0.set_momentum(bodek_rand_from_ball(_kf)); break; //Bodek Ritchie
		case  4: assert ("effective SF can be run only for carbon and oxygen" && p==n && (p==6 || p==8) ); //effective SF is defined only for two targets
					p0.set_momentum(spectral_choice(p,n));	break;
		case  0: // free nucleon
		case  5: p0.set_momentum(deuterium()); break; //deuterium
		case  6: p0.set_momentum(rand_from_ball(localkf(p0)));//effective potential

		default: p0.set_momentum(rand_from_ball(localkf(p0)));  //local Fermi Gas
	}

	// account for the fact that the nucleus is in motion
	//p0.set_momentum(p0.p()+vec(_p4)/Ar()); ////???????
	if(not (p0.v2()<1))
		cerr<<kMomDist
			<<" p="<<this->p
			<<" n="<<this->n
			<<" lkf="<<localkf(p0)
			<<" v2="<< p0.v2()
			<<" "<<p0.v()<<endl;
	// assert(p0.v2()<1);
	
	return p0;
}

double nucleus :: Ef (particle &pa)
{
	double const M = 0.5*(PDG::mass_proton + PDG::mass_neutron);
	double kmom = 0;
	
	switch (kMomDist)
	{
		case 0: case 5: case 6: return 0; break;
		case 2: kmom = localkf (pa); break;
		default: kmom = kF(); break;
	};

	
	return sqrt(kmom*kmom + M*M) - M;
}
