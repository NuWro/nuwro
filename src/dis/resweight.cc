////////////////This function calculates RES #Sect
////////////////RES region is defined in the hadronic mass from 1080 to the parameter res_dis_cut introduced in params.txt
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////The following model of the #Sect is postulated
///////// XSect (nu + N -> mu + N' + pi)/dW d omega  = [XSect (nu + N -> mu + N' + pi)/dW d omega ]_Delta * alfadelta (W) +
//////////                                             [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_SPP* SPP(W) * alfados(W) +
//////////                                             [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_nonSPP
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// The main idea behind is that DIS contribution simulates non-resonant part
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// The functions alfa and beta are a priori arbitrary but they have to satisfy the conditions
/////////// alfadis(res_dis_cut)=1, ....alfadelta(res_dis_cut=0,............alfadis (1080)=0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////alfadis and alfadelta are defined in alfa.cc
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////In principle for each channel the functions should be different. They should be fitted to experimental data
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////Channels are labelled by 4 integers corresponding to: nu/anu   cc/nc   proton/neutron    pi+/pi0/pi-
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "params.h"
//#include "hist.h"
#include "pdg_name.h"
#include "LeptonMass.h"
#include "parameters.h"
#include "jednostki.h"
//#include "bodek_old.h"
#include "grv94_bodek.h"
//#include "wiazka.h"
#include "dis_cr_sec.h"
#include "delta.h"
#include "fragmentation.h"
#include "vect.h"
#include "charge.h"
//#include "lorentz.h"
#include "event1.h"
#include "alfa.h"
#include "charge.h"


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void
resweight (params & p, event & e, bool cc)
///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
{
  double M2 = M12 * M12;
  bool current = cc;		//cc==true for charge current 
  int nukleon = e.in[1].pdg;
  int lepton = e.in[0].pdg;

  double m = lepton_mass (abs (lepton), current);	//mass of the produced lepton (see pgd header file)
  double m2 = m * m;

  vect nu0 = e.in[0];
  vect nuc0 = e.in[1];
  nu0.boost (-nuc0.speed ());

  double cut = p.res_dis_cut;

  double E = nu0.t;

  if (E < ((1080 + m) * (1080 + m) - M2) / 2 / M12)
    {
      e.weight = 0;
      return;
    }

  else
    {
/////////////////////////////////////////////////////////////
//      Selection of points in W, nu plane
/////////////////////////////////////////////////////////////

      double Wmax = min (cut, sqrt (M2 + 2 * M12 * E) - m);

      double W = 1080 + (Wmax - 1080) * frandom ();

      double W2 = W * W;
      double E2 = E * E;

      double wminus =
	((M12 + E) * (W2 - M2 - m2) + 2 * M12 * E2 -
	 E * sqrt (kwad (W2 - M2 - m2 - 2 * M12 * E) -
		   4 * m2 * M12 * (M12 + 2 * E))) / 2 / M12 / (M12 + 2 * E);

      double wplus =
	((M12 + E) * (W2 - M2 - m2) + 2 * M12 * E2 +
	 E * sqrt (kwad (W2 - M2 - m2 - 2 * M12 * E) -
		   4 * m2 * M12 * (M12 + 2 * E))) / 2 / M12 / (M12 + 2 * E);

      double numin = max (wminus, m);
      double numax = min (wplus, E - m);

      double nu = numin + (numax - numin) * frandom ();	//losujemy jednorodnie nu 

      double przedzial = (numax - numin) * (Wmax - 1080);

///////////////////////////////////////////////////////////////
//      End of selection of points in W, nu plane
/////////////////////////////////////////////////////////////

      double fromdis = cr_sec_dis (E, W, nu, lepton, nukleon, current);

      if (fromdis < 0)
	fromdis = 0;
      int j, k, l, t;

      if (lepton > 0)		//translation into "spp language"
	j = 0;
      else
	{
	  j = 1;
	}

      if (current)
	k = 0;
      else
	{
	  k = 1;
	}

      if (nukleon == 2212)
	l = 0;
      else
	{
	  l = 1;
	}

      double spp0 = SPPF (j, k, l, 0, W);
      double spp1 = SPPF (j, k, l, 1, W);
      double spp2 = SPPF (j, k, l, 2, W);

      double nonspp = fromdis * (1 - spp0 - spp1 - spp2);

      double dis0 = fromdis * spp0 * alfadis (j, k, l, 0, W);
      double dis1 = fromdis * spp1 * alfadis (j, k, l, 1, W);
      double dis2 = fromdis * spp2 * alfadis (j, k, l, 2, W);	//can be made simpler !!!

//below we have to determine first which nukleon is in the final state (we want to the function cr_sec_delta )
//it is not very convenient but I do not want to change everything at once
//we calculate the total electric charge of the final nucleon and pion system

      double delta0, delta1, delta2;
      int finalcharge = charge (nukleon) + (1 - k) * (1 - 2 * j);

      double adel0 = alfadelta (j, k, l, 0, W);
      double adel1 = alfadelta (j, k, l, 1, W);
      double adel2 = alfadelta (j, k, l, 2, W);

      if (finalcharge == 2)
	{
	  delta0 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2212, 211,
			  current) * adel0;
	  delta1 = delta2 = 0;
	}

      if (finalcharge == 1)
	{
	  delta0 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2112, 211,
			  current) * adel0;
	  delta1 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2212, 111,
			  current) * adel1;
	  delta2 = 0;
	}

      if (finalcharge == 0)
	{
	  delta1 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2112, 111,
			  current) * adel1;
	  delta2 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2212, -211,
			  current) * adel2;
	  delta0 = 0;
	}

      if (finalcharge == -1)
	{
	  delta2 =
	    cr_sec_delta (E, W, nu, lepton, nukleon, 2112, -211,
			  current) * adel2;
	  delta0 = delta1 = 0;
	}

//we arrived at the overall strength !!!
      double resstrength =
	nonspp + dis0 + dis1 + dis2 + delta0 + delta1 + delta2;

      double rel0 = nonspp / resstrength;
      double reldis0 = dis0 / resstrength;
      double reldis1 = dis1 / resstrength;
      double reldis2 = dis2 / resstrength;
      double reldelta0 = delta0 / resstrength;
      double reldelta1 = delta1 / resstrength;
      double reldelta2 = delta2 / resstrength;

      e.weight = resstrength * 1e-38 * przedzial;
      e.HadrMass = W;
      e.EnTran = nu;
      e.rel0 = rel0;
      e.reldis0 = reldis0;
      e.reldis1 = reldis1;
      e.reldis2 = reldis2;
      e.reldelta0 = reldelta0;
      e.reldelta1 = reldelta1;
      e.reldelta2 = reldelta2;

      return;
    }

}
