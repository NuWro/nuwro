////////////////This function calculates RES event
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
#include "masses.h"
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
#include "dis2res.h"
#include<TMCParticle.h>
#include<TPythia6.h>
#include "pauli.h"

TPythia6 *pythia77 = new TPythia6 ();
extern "C" int pycomp_ (const int *);

///////////////////////////////////////////////////////////////////
void
resevent (params & p, event & e, bool cc)
///////////////////////////////////////////////////////////////////
{				//cout<<"hello1"<<endl;

  int FFset = p.delta_FF_set;
  double M2 = M12 * M12;
  bool current = cc;		//cc==true for charge current 
  int nukleon = e.in[1].pdg;
  int lepton = e.in[0].pdg;

  double m = lepton_mass (abs (lepton), current);	//mass of the produced lepton (see pgd header file)
  double m2 = m * m;

  vect nu0 = e.in[0];
  vect nuc0 = e.in[1];
//cout<<nuc0.t<<" "<<nuc0.x<<" "<<nuc0.y<<" "<<nuc0.z<<" "<<nuc0.speed().length()<<endl;
  nu0.boost (-nuc0.speed ());

  vec numom = vec (nu0.x, nu0.y, nu0.z);	//can be made simpler ???

  double cut = p.res_dis_cut;
  target jadro (p);

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

      double Mtrue = M12;

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
//cout<<"hello2"<<endl;
///////////////////////////////////////////////////////////////
//      End of selection of points in W, nu plane
/////////////////////////////////////////////////////////////

      double fromdis = cr_sec_dis (E, W, nu, lepton, nukleon, current);

      if (fromdis < 0)
	fromdis = 0;
//cout<<"hello3"<<endl;
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

//const int j=jj;
//const int k=kk;
//const int l=ll;

//cout<<"hello4"<<endl;
//cout<<"a1"<<j<<" "<<k<<" "<<l<<endl;

      double spp0 = SPPF (j, k, l, 0, W);
      double spp1 = SPPF (j, k, l, 1, W);
      double spp2 = SPPF (j, k, l, 2, W);

      double nonspp = fromdis * (1 - spp0 - spp1 - spp2);

      double dis0 = fromdis * spp0 * alfadis (j, k, l, 0, W);
//cout<<"a2"<<j<<" "<<k<<" "<<l<<endl;
      double dis1 = fromdis * spp1 * alfadis (j, k, l, 1, W);
//cout<<"a3"<<j<<" "<<k<<" "<<l<<endl;
      double dis2 = fromdis * spp2 * alfadis (j, k, l, 2, W);	//can be made simpler !!!
//cout<<"a4"<<j<<" "<<k<<" "<<l<<endl;

//cout<<W<<"  "<<nukleon<<"  "<<fromdis<<"  "<<alfadis(j,k,l,0,W)<<"  "<<dis0<<"  "<<alfadis(j,k,l,1,W)<<"  "<<dis1<<"  "<<alfadis(j,k,l,2,W)<<"  "<<dis2<<"  "<<endl;
//cin.get();

//below we have to determine first which nukleon is in the final state (we want to the function cr_sec_delta )
//it is not very convenient but I do not want to change everything at once
//we calculate the total electric charge of the final nucleon and pion system

      double delta0, delta1, delta2;
      int finalcharge = charge1 (nukleon) + (1 - k) * (1 - 2 * j);

      double adel0 = alfadelta (j, k, l, 0, W);
      double adel1 = alfadelta (j, k, l, 1, W);
      double adel2 = alfadelta (j, k, l, 2, W);

      if (finalcharge == 2)
	{
	  delta0 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2212, 211,
			  current) * adel0;
	  delta1 = delta2 = 0;
	}

      if (finalcharge == 1)
	{
	  delta0 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2112, 211,
			  current) * adel0;
	  delta1 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2212, 111,
			  current) * adel1;
	  delta2 = 0;
	}

      if (finalcharge == 0)
	{
	  delta1 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2112, 111,
			  current) * adel1;
	  delta2 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2212, -211,
			  current) * adel2;
	  delta0 = 0;
	}

      if (finalcharge == -1)
	{
	  delta2 =
	    cr_sec_delta (FFset, E, W, nu, lepton, nukleon, 2112, -211,
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

////////////////////////////////////////////////////////////////////////////
//      Setting Pythia parameters
//      Done by Jaroslaw Nowak
//////////////////////////////////////////////

//stabilne pi0
      pythia77->SetMDCY (pycomp_ (&pizero), 1, 0);

      pythia77->SetMSTU (20, 1);	//advirsory warning for unphysical flavour switch off
      pythia77->SetMSTU (23, 1);	//It sets counter of errors at 0
      pythia77->SetMSTU (26, 0);	//no warnings printed 

      // PARJ(32)(D=1GeV) is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet parton system
      pythia77->SetPARJ (33, 0.1);

      // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below which 
      //the fragmentation of a parton system is stopped and two final hadrons formed.
      pythia77->SetPARJ (34, 0.5);
      pythia77->SetPARJ (35, 1.0);

      //PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point of the 
      //fragmentation. Strongly corlated with PARJ(33-35)

      pythia77->SetPARJ (37, 1.);	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      //MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and thus allow 
      // a cluster to decay rather than collaps
      pythia77->SetMSTJ (18, 3);	//do not change

/////////////////////////////////////////////////
//              End of setting Pythia parameters
////////////////////////////////////////////////

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
			  if (reldelta0 + reldelta1 + reldelta2 + rel0 +
			      reldis0 + reldis1 > los)
			    channel = -2;
			  else
			    channel = -3;
			}
		    }
		}
	    }
	}

      if (W < 1210 && channel < 0)
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

      vec kkprim;		//the unit vector in the direction of scattered lepton
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
      double W1 = W / GeV;	//W1 w GeV-ach potrzebne do Pythii




      if (channel <= 0)
	{
	  if (channel == 0)	//more inelastic event
	    {			//cout<<"MoreInelastic with"<<hama<<"  "<<current<<"  "<<lepton<<"  "<<nukleon<<endl;
	      while (NPar < 6)
		{
		  hadronization (E, W, nu, m, lepton, nukleon, current);
		  pythiaParticle = pythia77->GetPyjets ();
		  NPar = pythia77->GetN ();
		}
	    }

	  if (channel < 0)	//dis SPP event
	    {			//cout<<"DIS_SPP with"<<hama<<"  "<<current<<"  "<<lepton<<"  "<<nukleon<<"  "<<channel<<"  "<<meson<<endl;
	      while (!
		     (NPar == 5
		      && (pythiaParticle->K[1][3] == meson
			  || pythiaParticle->K[1][4] == meson)))
		{
		  hadronization (E, W, nu, m, lepton, nukleon, current);
		  pythiaParticle = pythia77->GetPyjets ();
		  NPar = pythia77->GetN ();
		}
	    }

	  nParticle = pythia77->GetN ();
//cout<<nParticle<<endl;
	  if (nParticle == 0)
	    {
	      cout << "nie ma czastek" << endl;
	      cin.get ();
	    }
//vect test = vect(0,0,0,0);
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
		{		//test=test+par[i];
		  e.out.push_back (part);
		}
	    }
//cout<<"resdis  "<<W<<"  "<<sqrt(test*test)<<endl;
	}
      else			//Delta decay 
	{			//cout<<"DeltaDecay with"<<hama<<"  "<<current<<"  "<<lepton<<"  "<<nukleon<<"  "<<channel<<"  "<<meson<<endl;
	  int nukleon2 = nukleon_out_ (W, lepton, nukleon, meson, current);	//which nucleon in the final state

	  vect finnuk, finpion;

	  kin2part (W, nukleon2, meson, finnuk, finpion);	//produces 4-momenta of final pair: nucleon + pion

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

//cout<<"delta  "<<W<<"  "<<sqrt((finpion+finnuk)*(finpion+finnuk))<<endl;
	}

     /* if (p.pauli_blocking)
	{
	  for (int i = 0; i < e.out.size (); i++)
	    if (((e.out[i].pdg == 2112) || (e.out[i].pdg == 2212))
		&& jadro.pauli_blocking (e.out[i]))
	      {
		e.weight = 0;
	      }
	}*/

    }

}
