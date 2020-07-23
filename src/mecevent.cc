#include <cassert>
#include "particle.h"
#include "jednostki.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "beam.h"
#include "pdg.h"
#include "nucleus.h"
#include "mecdynamics.h"
#include "event1.h"
using namespace NUWRO;

//      MEC
//      Implementation of the TEM
// only muon neutrino
//  flux direction is (0,0,1)

////////////////////////////////////////////////////////////////////////
void
mecevent (params & p, event & e, nucleus & t, bool cc)
////////////////////////////////////////////////////////////////////////
{
	bool nu;					 //neutrino or antineutrino
	e.par = p;
	e.flag.cc = cc;
	e.flag.nc = !cc;
	e.flag.dis = false;
	e.flag.qel = false;
	e.flag.coh = false;
	e.flag.mec = true;
	bool fsi = p.kaskada_on;
//	double fermimom = p.nucleus_kf;							 
//	double potwell = sqrt(fermimom*fermimom + 939*939) - 939;//Fermi energy
//  double ebinding = 8*MeV;
	double potwell = t.Ef();	//Fermi energy
	double ebinding= t.Eb();	//Binding energy

	if(t.A()<4)
	{
		e.weight=0;
		return;
	}

	//      Initial neutrino
	particle mecnu = e.in[0];

	//      Final lepton; neutrino for nc and charged lepton for cc
	particle meclepton;

	//      Final nucleons
	particle mecnucleon1, mecnucleon2;

	//      Identification of final states
	//      In NC reaction always pi0 is produced
	//      In CC neutrino reaction pi+ is produced and in antineutrino reaction pi- is produced

	if (!cc)					 //not really necessary as will have CC reaction only
	{
		meclepton = mecnu;
	}
	else if (mecnu.pdg > 0)
	{
		meclepton.pdg = mecnu.pdg - 1;
	}
	else
	{
		meclepton.pdg = mecnu.pdg + 1;
	}

	if (meclepton.pdg >0)
		nu=true;
	else
		nu=false;

	//      Setting final lepton masses
	meclepton.set_mass (PDG::mass (meclepton.pdg));

	double mecm = meclepton.mass ();

	//cout<<mecnu.t<<endl;

								 //weight
	double wynik = mecweight (mecnu.t, nu, t, p, meclepton, mecnucleon1, mecnucleon2, fsi, potwell);
	//cout<<"sleep4"<<endl;
	e.weight = wynik;
    
	e.out.push_back (meclepton);
	//cout<<"sleep5"<<"  "<<meclepton<<endl;

	if (!fsi)
	{
		if (mecnucleon1.Ek() >0.1)
			e.out.push_back (mecnucleon1);
		//cout<<"sleep6"<<"  "<<mecnucleon1<<endl;
		if (mecnucleon2.Ek() >0.1)
			e.out.push_back (mecnucleon2);
		//cout<<"sleep7"<<"  "<<mecnucleon2<<endl;
	}
	else
	{
		if (mecnucleon1.Ek() >potwell+ebinding)
			e.out.push_back (mecnucleon1);
	
		if (mecnucleon2.Ek() >potwell+ebinding)
			e.out.push_back (mecnucleon2);
	
	}
}
