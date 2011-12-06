#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "jednostki.h"
#include "params.h"
//#include "hist.h"
#include "pdg_name.h"
#include "parameters.h"
#include "LeptonMass.h"
//#include "bodek_old.h"
#include "grv94_bodek.h"
//#include "wiazka.h"
#include "dis_cr_sec.h"
//#include "delta.h"
#include "fragmentation.h"
#include "vect.h"
#include "charge.h"
//#include "lorentz.h"
//#include "dis2res.h"
#include "event1.h"


void
disweight (params & p, event & e, bool cc)
{

/////////////////initial parameters/////////////////////
  bool current = cc;		//cc==true for charge current , nc == false
  int nukleon = e.in[1].pdg;
  int lepton = e.in[0].pdg;

  double m = lepton_mass (abs (lepton), current);	//mass of the given lepton (see pgd header file)
  double m2 = m * m;

  vect nu0 = e.in[0];
  vect nuc0 = e.in[1];

  nu0.boost (-nuc0.speed ());
  double Mtrue;

  double cut = p.res_dis_cut;


  Mtrue = M12;
  double Mtrue2 = Mtrue * Mtrue;
  double E = nu0.t;

  if (E < ((cut + m) * (cut + m) - Mtrue2) / 2 / Mtrue)
    {
      e.weight = 0;
      return;
    }

  else
    {
/////////////////////////////////////////////////////////////
//      Selection of points in W, nu plane
////////////////////////////////////////////////////////////
      double W = (sqrt (Mtrue2 + 2 * Mtrue * E) - m - cut) * frandom () + cut;

      double W2 = W * W;
      double E2 = E * E;

      double wminus =
	((Mtrue + E) * (W2 - Mtrue2 - m2) + 2 * Mtrue * E2 -
	 E * sqrt (kwad (W2 - Mtrue2 - m2 - 2 * Mtrue * E) -
		   4 * m2 * Mtrue * (Mtrue + 2 * E))) / 2 / Mtrue / (Mtrue +
								     2 * E);

      double wplus =
	((Mtrue + E) * (W2 - Mtrue2 - m2) + 2 * Mtrue * E2 +
	 E * sqrt (kwad (W2 - Mtrue2 - m2 - 2 * Mtrue * E) -
		   4 * m2 * Mtrue * (Mtrue + 2 * E))) / 2 / Mtrue / (Mtrue +
								     2 * E);

      double numin = max (wminus, m);
      double numax = min (wplus, E - m);

      double nu = numin + (numax - numin) * frandom ();	//losujemy jednorodnie nu 

      double waga = cr_sec_dis (E, W, nu, lepton, nukleon, current);

      if (waga < 0)
	{
	  e.weight = 0;
	  return;
	}

      double przedzial =
	(numax - numin) * (sqrt (Mtrue2 + 2 * Mtrue * E) - m - cut);



      e.weight = waga * 1e-38 * przedzial;
      e.HadrMass = W;
      e.EnTran = nu;
      return;

    }
}
