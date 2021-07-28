#include "nucleus.h"
#include "nucleus_data.h"
#include "particle.h"
   
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
	

	/* model_target=
	* 0 for free target; 
	* 1 for Fermi gas; 
	* 2 for local Fermi gas; 
	* 3 for Bodek; 
	* 4 for spectral function; 
	* 5 for deuterium; 
	* 6 for only proton in deuterium (for tests only);
	*/

	switch(kMomDist)
	{
		case  1: p0.set_momentum(rand_from_ball(_kf)); break; // fermi gas
		case  2: p0.set_momentum(rand_from_ball(localkf(p0)));break; // local fermi gas		
		case  3: p0.set_momentum(bodek_rand_from_ball(_kf)); break; //Bodek
		case  4: if(p==n && (p==6 || p==8))
				{
					p0.set_momentum(spectral_choice(p,n));	//spectral function for carbon and oxygen
					break;
				}
		case 0: // proton
		case 5: // deuterium
		case 6: // deuterium
				// use local fermi gas if H1 or H2 options used for other nucleus
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
