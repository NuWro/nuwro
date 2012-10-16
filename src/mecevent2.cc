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
#include "mecdynamics2.h"
#include "event1.h"


//      MEC 
//      Implementation of the np-nh
// only muon neutrino
//  flux direction is (0,0,1)

////////////////////////////////////////////////////////////////////////
void
mecevent2 (params & p, event & e, nucleus & t, bool cc)
////////////////////////////////////////////////////////////////////////
{
  e.par = p;
  e.flag.mec = true;
  e.flag.cc = cc;
  e.flag.nc = !cc;
  e.flag.dis = false;
  e.flag.qel = false;
  e.flag.coh = false;
  bool fsi = p.kaskada_on;
  double fermimom = p.nucleus_kf;
  double potwell = sqrt(fermimom*fermimom + 939*939) - 939;//Fermi energy

nucleus jadro (p);
int mecA = jadro.p + jadro.n;// number of nucleons

//      Initial neutrino
particle mecnu = e.in[0];

//      Final lepton; neutrino for nc and charged lepton for cc
particle meclepton;

//      Final nucleons
particle mecnucleon1, mecnucleon2, mecnucleon3;

//      Identification of final states
//      In NC reaction always pi0 is produced
//      In CC neutrino reaction pi+ is produced and in antineutrino reaction pi- is produced

  if (!cc)//not really necessary as will have CC reaction only 
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

if (mecnu.t>150)// for lower energies MM does not work
{//150 loop

double wynik = mecweight2 (mecnu.t, cc, mecA, meclepton, mecnucleon1, mecnucleon2, mecnucleon3, fsi, potwell);	//weight
//cout<<"sleep4"<<endl;
e.weight = wynik;

//cout<<mecnucleon1<<"  "<<mecnucleon2<<"  "<<mecnucleon3<<endl;

e.out.push_back (meclepton);
//cout<<"sleep5"<<"  "<<meclepton<<endl;

if (!fsi)
{
if (mecnucleon1.Ek() >0.1)
e.out.push_back (mecnucleon1);

if (mecnucleon2.Ek() >0.1)
e.out.push_back (mecnucleon2);

if (mecnucleon3.Ek() >0.1)
e.out.push_back (mecnucleon3);
}
else
{
if (mecnucleon1.Ek() >potwell+8)
e.out.push_back (mecnucleon1);

if (mecnucleon2.Ek() >potwell+8)
e.out.push_back (mecnucleon2);

if (mecnucleon3.Ek() >potwell+8)
e.out.push_back (mecnucleon3);
}
}//end 150 loop

else
  e.weight = 0.0;

}
