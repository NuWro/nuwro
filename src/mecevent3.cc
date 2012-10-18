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
#include "mecdynamics3.h"
#include "event1.h"

//      MEC
//      Implementation of the TEM-FG
// only muon neutrino
//  flux direction is (0,0,1)

////////////////////////////////////////////////////////////////////////
void
mecevent3 (params & p, event & e, nucleus & t, bool cc)
////////////////////////////////////////////////////////////////////////
{
	e.par = p;
	e.flag.mec = true;
	e.flag.cc = cc;
	e.flag.nc = !cc;
	e.flag.dis = false;
	e.flag.qel = false;
	e.flag.coh = false;

	nucleus jadro (p);
	int mecA = jadro.p + jadro.n;// number of nucleons

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

	//      Setting final lepton masses
	meclepton.set_mass (PDG::mass (meclepton.pdg));

	double mecm = meclepton.mass ();

	//cout<<mecnu.t<<endl;

								 //weight
	double wynik = mecweight3 (mecnu.t, mecnu.pdg > 0, cc, mecA, meclepton, mecnucleon1, mecnucleon2);
	//cout<<"sleep4"<<endl;
	e.weight = wynik;

	e.out.push_back (meclepton);
	//cout<<"sleep5"<<"  "<<meclepton<<endl;
	if (mecnucleon1.momentum() >0.1)
		e.out.push_back (mecnucleon1);
	//cout<<"sleep6"<<"  "<<mecnucleon1<<endl;
	if (mecnucleon2.momentum() >0.1)
		e.out.push_back (mecnucleon2);
	//cout<<"sleep7"<<"  "<<mecnucleon2<<endl;
}
