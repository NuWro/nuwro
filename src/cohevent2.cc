#include <cassert>
#include "particle.h"
#include "jednostki.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "pdg.h"
#include "nucleus.h"
#include "cohdynamics2.h"
#include "event1.h"
#include "quat.h"


//      Coherent single pion production
//      Implementation of the Rein-Sehgal model
//      Nucl. Phys. B 



////////////////////////////////////////////////////////////////////////
void
cohevent2 (params & p, event & e, nucleus & t, bool cc)
////////////////////////////////////////////////////////////////////////
{
  vec lepton3mom;
  vec pion3mom;
  e.flag.coh = true;
  e.flag.cc = cc;
  e.flag.nc = !cc;
  e.flag.dis = false;
  e.flag.qel = false;

  int cohA = t.A();


  particle nu = e.in[0];
  particle lepton;
  particle pion;

//      Identification of final states
//      In NC reaction always pi0 is produced
//      In CC neutrino reaction pi+ is produced and in antineutrino reaction pi- is produced

  if (!cc)
    {
      lepton = nu;
      pion.pdg = 111;
    }
  else if (nu.pdg > 0)
    {
      lepton.pdg = nu.pdg - 1;
      pion.pdg = 211;
    }
  else
    {
      lepton.pdg = nu.pdg + 1;
      pion.pdg = -211;
    }

	// Setting final particles masses
	lepton.set_mass (PDG::mass (lepton.pdg));
	pion.set_mass (PDG::mass (pion.pdg));

	double cohm = lepton.mass ();
	double wynik = cohwaga2 (nu.t, cc, cohA, p.coh_mass_correction, cohm,
							 lepton3mom,pion3mom);	//weight
	//cout<<cohnu.t<<" "<<cc<<" "<<cohA<<" "<<p.coh_mass_correction<<"  "<<cohm<<" "<<wynik<<endl;
	e.weight = wynik;
	{
	/// apply T to lepton3mom and pion3mom 
	/// T = transformation taking vec(0,0,1) to nu.speed();
	/// it is one reflection ( if nu.z<0 ) or two reflections (if nu.z>=0)
       vec zvec=vec(0,0,1);
       if(nu.z<0)
		   zvec.z=-1;
	   else 
	   {
		   lepton3mom.z*=-1; 
		   pion3mom.z*=-1;
	   }
       vec s=(zvec+nu.p().dir()).dir();
       lepton3mom-=2*(lepton3mom*s)*s;
       pion3mom-=2*(pion3mom*s)*s;
	}
	lepton.set_momentum (lepton3mom);
	pion.set_momentum (pion3mom);

	e.out.push_back (lepton);
	e.out.push_back (pion);
//	cout<<"--------------"<<endl;
//	cout<<vect(e.in[0])<<endl;
//	cout<<e.out[0]+e.out[1]<<endl;
	//return wynik;
}
