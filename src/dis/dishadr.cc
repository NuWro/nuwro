#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "params.h"
#include "pdg_name.h"
#include "parameters.h"
#include "LeptonMass.h"
//#include "jednostki.h"
#include "grv94_bodek.h"
#include "dis_cr_sec.h"
#include "fragmentation.h"
#include "vect.h"
#include "charge.h"
#include "event1.h"
#include <TMCParticle.h>
#include <TPythia6.h>

TPythia6 *pythia2 = new TPythia6 ();
extern "C" int pycomp_ (const int *);

void
dishadr (event & e, bool current, double hama, double entra)
{
////////////////////////////////////////////////////////////////////////////
//      Setting Pythia parameters
//      Done by Jaroslaw Nowak
//////////////////////////////////////////////

//stabilne pi0
  pythia2->SetMDCY (pycomp_ (&pizero), 1, 0);
//C Thorpe: Adding Hyperons as stable dis particles
  pythia2->SetMDCY (pycomp_ (&Lambda), 1, 0);
  pythia2->SetMDCY (pycomp_ (&Sigma), 1, 0);
  pythia2->SetMDCY (pycomp_ (&SigmaP), 1, 0);
  pythia2->SetMDCY (pycomp_ (&SigmaM), 1, 0);

  // C Thorpe: Stablize kaons
  pythia22->SetMDCY ( pycomp_ (&DIS_PDG::Kplus) , 1, 0);
  pythia22->SetMDCY ( pycomp_ (&DIS_PDG::Kzero) , 1, 0);
  pythia22->SetMDCY ( pycomp_ (&DIS_PDG::Kminus) , 1, 0);


  pythia2->SetMSTU (20, 1);	//advirsory warning for unphysical flavour switch off
  pythia2->SetMSTU (23, 1);	//It sets counter of errors at 0
  pythia2->SetMSTU (26, 0);	//no warnings printed 

  // PARJ(32)(D=1GeV) is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet parton system
  pythia2->SetPARJ (33, 0.1);

  // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below which 
  //the fragmentation of a parton system is stopped and two final hadrons formed.
  pythia2->SetPARJ (34, 0.5);
  pythia2->SetPARJ (35, 1.0);

  //PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point of the 
  //fragmentation. Strongly corlated with PARJ(33-35)

  pythia2->SetPARJ (37, 1.);	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  //MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and thus allow 
  // a cluster to decay rather than collaps
  pythia2->SetMSTJ (18, 3);	//do not change

/////////////////////////////////////////////////
//              End of setting Pythia parameters
////////////////////////////////////////////////

  double Mtrue = (PDG::mass_proton+PDG::mass_neutron)/2;
  double Mtrue2 = Mtrue * Mtrue;
  double W2 = hama * hama;
  double nu = entra;

  vect nuc0 = e.in[1];
  vect nu0 = e.in[0];
  nu0.boost (-nuc0.speed ());	//neutrino 4-momentum in the target rest frame

  vec nulab = vec (nu0.x, nu0.y, nu0.z);	//can be made simpler ???
  int lepton = e.in[0].pdg;
  int nukleon = e.in[1].pdg;
  double m = lepton_mass (abs (lepton), current);	//mass of the given lepton (see pgd header file)
  double m2 = m * m;

  double E = nu0.t;
  double E2 = E * E;

  int nParticle = 0;
  int nCharged = 0;
  int NPar = 0;
  Pyjets_t *pythiaParticle;	//deklaracja event recordu
  double W1 = hama / GeV;	//W1 w GeV-ach potrzebne do Pythii

  while (NPar < 5)
    {
      hadronization (E, hama, entra, m, lepton, nukleon, current);
      pythiaParticle = pythia2->GetPyjets ();
      NPar = pythia2->GetN ();
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////Kinematics
///////////////////////////////////////////////
//////////// The input data is: neutrinoo 4-momentum, invariant hadronic mass and energy transfer
////////////  With this data the kinematics is resolved
////////////////////////////////////////////////////
/////////// We know that nu^2-q^2= (k-k')^2=m^2-2*k.k'= m^2-2*E*(E-nu)+ 2*E*sqrt((E-nu)^2-m^2)*cos theta
///////////
///////////              (M+nu)^2-q^2=W^2
////////////
////////////             (k-q)^2= m^2 = nu^2-q^2 -2*E*nu + 2*E*q*cos beta
/////////////
/////////////             theta is an angle between leptons and beta is an angle between neutrino and momentum transfer
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////              Of course it is not necessary to calculate vectors k' and q separately because of momentum conservation
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double q = sqrt (kwad (Mtrue + nu) - W2);
  double kprim = sqrt (kwad (E - nu) - m2);
  double cth = (E2 + kprim * kprim - q * q) / 2 / E / kprim;


  vec kkprim;			//the unit vector in the direction of scattered lepton
  kinfinder (nulab, kkprim, cth);	//hopefully should produce kkprim 

  kkprim = kprim * kkprim;	//multiplied by its length

  vect lepton_out = vect (E - nu, kkprim.x, kkprim.y, kkprim.z);

  vec momtran = nulab - kkprim;

  vec hadrspeed = momtran / sqrt (W2 + q * q);
  nParticle = pythia2->GetN ();

  if (nParticle == 0)
    {
      cout << "nie ma czastek" << endl;
      cin.get ();
    }

  vect par[100];
  double ks[100];		//int czy double ???

  par[0] = lepton_out;

//powrot do ukladu spoczywajacej tarczy
  par[0] = par[0].boost (nuc0.speed ());	//ok

  particle lept (par[0]);

  if (current == true && lepton > 0)
    {
      lept.pdg = lepton - 1;
    }
  if (current == true && lepton < 0)
    {
      lept.pdg = lepton + 1;
    }
  if (current == false)
    {
      lept.pdg = lepton;
    }

  e.out.push_back (lept);	//final lepton; ok

  for (int i = 0; i < nParticle; i++)
    {
      par[i].t = pythiaParticle->P[3][i] * GeV;
      par[i].x = pythiaParticle->P[0][i] * GeV;
      par[i].y = pythiaParticle->P[1][i] * GeV;
      par[i].z = pythiaParticle->P[2][i] * GeV;
      rotation (par[i], momtran);
      ks[i] = pythiaParticle->K[0][i];
      par[i] = par[i].boost (hadrspeed);	//correct direction ???
      par[i] = par[i].boost (nuc0.speed ());
      particle part (par[i]);



      part.ks = pythiaParticle->K[0][i];
      part.pdg = pythiaParticle->K[1][i];
      part.orgin = pythiaParticle->K[2][i];

      e.temp.push_back (part);
      if (ks[i] == 1)		//condition for a real particle in the final state
	{
	  e.out.push_back (part);
	}
    }
}
