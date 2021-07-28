#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "params.h"
#include "pdg_name.h"
#include "parameters.h"
#include "LeptonMass.h"
#include "jednostki.h"
#include "grv94_bodek.h"
#include "dis_cr_sec.h"
#include "fragmentation.h"
#include "vect.h"
#include "charge.h"
#include "dis2res.h"
#include "event1.h"
#include <TMCParticle.h>
#include <TPythia6.h>

TPythia6 *pythia7 = new TPythia6 ();
extern "C" int pycomp_ (const int *);

void
reshadr (event & e, bool current, double hama, double entra, double rel0,
	 double reldelta0, double reldelta1, double reldelta2, double reldis0,
	 double reldis1, double reldis2)
{
////////////////////////////////////////////////////////////////////////////
//      Setting Pythia parameters
//      Done by Jaroslaw Nowak
//////////////////////////////////////////////

std::cout << "Setting RES PYTHIA parameters" << std::endl;

//stabilne pi0
  pythia7->SetMDCY (pycomp_ (&pizero), 1, 0);
//C Thorpe: Adding Hyperons as stable dis particles
  pythia7->SetMDCY (pycomp_ (&Lambda), 1, 0);
  pythia7->SetMDCY (pycomp_ (&Sigma), 1, 0);
  pythia7->SetMDCY (pycomp_ (&SigmaP), 1, 0);
  pythia7->SetMDCY (pycomp_ (&SigmaM), 1, 0);

  pythia7->SetMSTU (20, 1);	//advirsory warning for unphysical flavour switch off
  pythia7->SetMSTU (23, 1);	//It sets counter of errors at 0
  pythia7->SetMSTU (26, 0);	//no warnings printed 

  // PARJ(32)(D=1GeV) is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet parton system
  pythia7->SetPARJ (33, 0.1);

  // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below which 
  //the fragmentation of a parton system is stopped and two final hadrons formed.
  pythia7->SetPARJ (34, 0.5);
  pythia7->SetPARJ (35, 1.0);

  //PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point of the 
  //fragmentation. Strongly corlated with PARJ(33-35)

  pythia7->SetPARJ (37, 1.);	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  //MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and thus allow 
  // a cluster to decay rather than collaps
  pythia7->SetMSTJ (18, 3);	//do not change

/////////////////////////////////////////////////
//              End of setting Pythia parameters
////////////////////////////////////////////////

  double Mtrue = M12;
  double Mtrue2 = Mtrue * Mtrue;
  double W2 = hama * hama;
  double nu = entra;

  vect nuc0 = e.in[1];
  vect nu0 = e.in[0];
  nu0.boost (-nuc0.speed ());	//neutrino 4-momentum in the target rest frame

  vec numom = vec (nu0.x, nu0.y, nu0.z);	//can be made simpler ???
  int lepton = e.in[0].pdg;
  int nukleon = e.in[1].pdg;
  double m = lepton_mass (abs (lepton), current);	//mass of the given lepton (see pgd header file)
  double m2 = m * m;
  double E = nu0.t;

  double E2 = E * E;

//in the MC method we determine which process actually happened
  int channel;
  double los = frandom ();
  if (reldelta0 > los)
    channel = 1;
  else
    {
      if (reldelta0 + reldelta1 > los)
	channel = 2;
      else
	{
	  if (reldelta0 + reldelta1 + reldelta2 > los)
	    channel = 3;
	  else
	    {
	      if (reldelta0 + reldelta1 + reldelta2 + rel0 > los)
		channel = 0;
	      else
		{
		  if (reldelta0 + reldelta1 + reldelta2 + rel0 + reldis0 >
		      los)
		    channel = -1;
		  else
		    {
		      if (reldelta0 + reldelta1 + reldelta2 + rel0 + reldis0 +
			  reldis1 > los)
			channel = -2;
		      else
			channel = -3;
		    }
		}
	    }
	}
    }

  if (hama < 1210 && channel < 0)
    channel *= -1;

  int meson = 0;

  if (channel == 1 || channel == -1)
    meson = 211;
  if (channel == 2 || channel == -2)
    meson = 111;
  if (channel == 3 || channel == -3)
    meson = -211;

  double q = sqrt (kwad (Mtrue + nu) - W2);
  double kprim = sqrt (kwad (E - nu) - m2);
  double cth = (E2 + kprim * kprim - q * q) / 2 / E / kprim;

  vec kkprim;			//the unit vector in the direction of scattered lepton
  kinfinder (numom, kkprim, cth);	//produces kkprim 

  kkprim = kprim * kkprim;	//multiplied by its length

  vect lepton_out = vect (E - nu, kkprim.x, kkprim.y, kkprim.z);

  vec momtran = numom - kkprim;

  if (E - nu < 0)
    {
      cout << "ojej2" << endl;
      cin.get ();
    }

  vec hadrspeed = momtran / sqrt (W2 + q * q);

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

  int nParticle = 0;
  int nCharged = 0;
  int NPar = 0;
  Pyjets_t *pythiaParticle;	//deklaracja event recordu
  double W1 = hama / GeV;	//W1 w GeV-ach potrzebne do Pythii

  if (channel <= 0)
    {
      if (channel == 0)		//more inelastic event
	{
	  while (NPar < 6)
	    {
	      hadronization (E, hama, entra, m, lepton, nukleon, current);
	      pythiaParticle = pythia7->GetPyjets ();
	      NPar = pythia7->GetN ();
	    }
	}

      if (channel < 0)		//dis SPP event
	{
	  while (!
		 (NPar == 5
		  && (pythiaParticle->K[1][3] == meson
		      || pythiaParticle->K[1][4] == meson)))
	    {
	      hadronization (E, hama, entra, m, lepton, nukleon, current);
	      pythiaParticle = pythia7->GetPyjets ();
	      NPar = pythia7->GetN ();
	    }
	}

      nParticle = pythia7->GetN ();
      if (nParticle == 0)
	{
	  cout << "nie ma czastek" << endl;
	  cin.get ();
	}

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
	  if (ks[i] == 1)	//condition for a real particle in the final state
	    {
	      e.out.push_back (part);
	    }
	}
    }
  else				//Delta decay 
    {
      int nukleon2 = nukleon_out_ (hama, lepton, nukleon, meson, current);	//which nucleon in the final state

      vect finnuk, finpion;

      kin2part (hama, nukleon2, meson, finnuk, finpion);	//produces 4-momenta of final pair: nucleon + pion

      finnuk = finnuk.boost (hadrspeed);
      finnuk = finnuk.boost (nuc0.speed ());

      finpion = finpion.boost (hadrspeed);
      finpion = finpion.boost (nuc0.speed ());

      particle ppion (finpion);
      particle nnukleon (finnuk);

      ppion.pdg = meson;
      nnukleon.pdg = nukleon2;

      e.out.push_back (ppion);
      e.out.push_back (nnukleon);
    }

}
