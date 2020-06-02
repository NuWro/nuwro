#include <iostream>
#include <stdio.h>
#include <vector>
#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
//#include "qelevent.h"
#include "pdg.h"
//#include "nucleus.h"
#include "generatormt.h"

using namespace NSNWRO;



class ND280_event
{
public:
  Int_t Nev;
  Float_t Pos[3];
  Int_t Mode;
  Int_t Numnu;
  Int_t Ipnu[100];		//[Numnu
  Float_t Abspnu[100];		//[Numnu]
  Float_t Pnu[100][3];		//[Numnu]
  Int_t Npar;
  Int_t Ipv[100];		//[Npar]
  Int_t Icrnv[100];		//[Npar]
  Float_t Pmomv[100][3];	//[Npar]
  Float_t xpi[3];
  Float_t npi[3];
  Float_t cospibm;
  Float_t ppi;
  Int_t ppid;
  Float_t xpi0[3];
  Float_t npi0[3];
  Float_t cospi0bm;
  Float_t ppi0;
  Float_t Abspv[100];		//Npar
  Float_t Enu;
  Int_t Iflgv[100];		// Npar
  Int_t Iorgv[100];		// Npar
  Int_t idfd;
  Int_t modef;
  Float_t nnu[3];
  Float_t norm;
  Float_t rnu;
  Float_t xnu;
  Float_t ynu;

  void Create_Branches (TTree *);
  void Translate_From_NuWro (event * e);
  int ND280_Mode (event * e);
};


void
ND280_event::Create_Branches (TTree * tt2)
{

  tt2->Branch ("Nev", &Nev, "Nev/I");
  tt2->Branch ("Pos", Pos, "Pos[3]/F");
  tt2->Branch ("Mode", &Mode, "Mode/I");
  tt2->Branch ("Numnu", &Numnu, "Numnu/I");
  tt2->Branch ("Ipnu", Ipnu, "Ipnu[Numnu]/I");
  tt2->Branch ("Abspnu", Abspnu, "Abspnu[Numnu]/F");
  tt2->Branch ("Pnu", Pnu, "Pnu[Numnu][3]/F");
  tt2->Branch ("Npar", &Npar, "Npar/I");
  tt2->Branch ("Ipv", Ipv, "Ipv[Npar]/I");
  tt2->Branch ("Icrnv", Icrnv, "Icrnv[Npar]/I");
  tt2->Branch ("Pmomv", Pmomv, "Pmomv[Npar][3]/F");
  tt2->Branch ("xpi", xpi, "xpi[3]/F");
  tt2->Branch ("npi", npi, "npi[3]/F");
  tt2->Branch ("cospibm", &cospibm, "cospibm/F");
  tt2->Branch ("ppi", &ppi, "ppi/F");
  tt2->Branch ("ppid", &ppid, "ppid/I");
  tt2->Branch ("xpi0", xpi0, "xpi0[3]/F");
  tt2->Branch ("npi0", npi0, "npi0[3]/F");
  tt2->Branch ("cospi0bm", &cospi0bm, "cospi0bm/F");
  tt2->Branch ("ppi0", &ppi0, "ppi0/F");
  ////
  tt2->Branch ("Abspv", Abspv, "Abspv[Npar]/F");
  tt2->Branch ("Enu", &Enu, "Enu/F");
  tt2->Branch ("Iflgv", Iflgv, "Iflgv[Npar]/I");
  tt2->Branch ("Iorgv", Iorgv, "Iorgv[Npar]/I");
  tt2->Branch ("idfd", &idfd, "idfd/I");
  tt2->Branch ("modef", &modef, "modef/I");
  tt2->Branch ("nnu", nnu, "nnu[3]/F");
  tt2->Branch ("norm", &norm, "norm/F");
  tt2->Branch ("rnu", &rnu, "rnu/F");
  tt2->Branch ("xnu", &xnu, "xnu/F");
  tt2->Branch ("ynu", &ynu, "ynu/F");

}

int
ND280_event::ND280_Mode (event * e)
{


  Int_t proton_pdg, neutron_pdg, pion_pdg, pion_plus_pdg, pion_minus_pdg,
    lambda_pdg, eta_pdg, kaon_pdg, kaon_plus_pdg;
  proton_pdg = 2212;
  eta_pdg = 221;
  neutron_pdg = 2112;
  pion_pdg = 111;
  pion_plus_pdg = 211;
  pion_minus_pdg = -211;
  //O_16_pdg = 100069;   // oznacznie z Neuta
  lambda_pdg = 3122;
  kaon_pdg = 311;
  kaon_plus_pdg = 321;


  if (e->flag.qel)		// kwiazielastyczne oddziaływanie
    {
      if (e->flag.anty)		// jeśli jest to oddziaływanie z antyneutrinem
	{
	  if (e->flag.cc)
	    return -1;
	  else
	    {
	      if (e->nof (proton_pdg))
		return -51;
	      else if (e->nof (neutron_pdg))
		return -52;	// sprawdzam dodatkowo ?
	    }
	}
      else			// oddziaływanie z neutrinem
	{
	  if (e->flag.cc)
	    return 1;
	  else
	    {
	      if (e->nof (proton_pdg))
		return 51;
	      else if (e->nof (neutron_pdg))
		return 52;
	    }
	}
    }


  if (e->flag.res)		//rezonansowa produkcja: pojedynczy pion, pojed.eta, kaon, multipiony  
    {

      Int_t liczba_pionow, liczba_kaonow;

      liczba_pionow =
	e->nof (pion_pdg) + e->nof (pion_plus_pdg) + e->nof (pion_minus_pdg);
      liczba_kaonow = e->nof (kaon_pdg) + e->nof (kaon_pdg);

      if (liczba_pionow > 1 || liczba_pionow == 0)	// multipiony
	{
	  if (e->flag.anty)
	    {
	      if (e->flag.cc)
		return -21;
	      else
		return -41;
	    }
	  else
	    {
	      if (e->flag.cc)
		return 21;
	      else
		return 41;
	    }
	}

      if (liczba_pionow == 1)
	{
	  if (e->flag.anty)	// jeśli jest to oddziaływanie z antyneutrinem
	    {
	      if (e->flag.cc)
		{
		  if (e->nof (neutron_pdg) && e->nof (pion_minus_pdg))
		    return -11;
		  if (e->nof (neutron_pdg) && e->nof (pion_pdg))
		    return -12;
		  if (e->nof (proton_pdg) && e->nof (pion_minus_pdg))
		    return -13;
		}
	      else
		{
		  if (e->nof (proton_pdg))
		    {
		      if (e->nof (pion_minus_pdg))
			return -33;
		      else if (e->nof (pion_pdg))
			return -32;
		    }
		  else if (e->nof (neutron_pdg))
		    {
		      if (e->nof (pion_plus_pdg))
			return -34;
		      else if (e->nof (pion_pdg))
			return -31;
		    }
		}
	    }
	  else			// oddziaływanie z neutrinem
	    {
	      if (e->flag.cc)
		{
		  if (e->nof (proton_pdg) && e->nof (pion_plus_pdg))
		    return 11;
		  if (e->nof (proton_pdg) && e->nof (pion_pdg))
		    return 12;
		  if (e->nof (neutron_pdg) && e->nof (pion_plus_pdg))
		    return 13;
		}
	      else
		{
		  if (e->nof (proton_pdg))
		    {
		      if (e->nof (pion_minus_pdg))
			return 33;
		      else if (e->nof (pion_pdg))
			return 32;
		    }
		  else if (e->nof (neutron_pdg))
		    {
		      if (e->nof (pion_plus_pdg))
			return 34;
		      else if (e->nof (pion_pdg))
			return 31;
		    }
		}
	    }
	}

      if (e->nof (eta_pdg))	// produkcja rezonansowa ety
	{
	  if (e->flag.anty)	// jeśli jest to oddziaływanie z antyneutrinem
	    {
	      if (e->flag.cc)
		return -22;
	      else
		{
		  if (e->nof (neutron_pdg))
		    return -42;
		  else if (e->nof (proton_pdg))
		    return -43;	// sprawdzam dodatkowo ?
		}
	    }
	  else			// oddziaływanie z neutrinem
	    {
	      if (e->flag.cc)
		return 22;
	      else
		{
		  if (e->nof (neutron_pdg))
		    return 42;
		  else if (e->nof (proton_pdg))
		    return 43;
		}
	    }
	}

      if (e->nof (lambda_pdg) == 1 && liczba_kaonow == 1)	// produkcja rezonansowa kaonu
	{
	  if (e->flag.anty)	// jeśli jest to oddziaływanie z antyneutrinem
	    {
	      if (e->flag.cc && e->nof (kaon_pdg))
		return -23;
	      else
		{
		  if (e->nof (kaon_pdg))
		    return -44;
		  else if (e->nof (kaon_plus_pdg))
		    return -45;
		}
	    }
	  else			// oddziaływanie z neutrinem
	    {
	      if (e->flag.cc && e->nof (kaon_plus_pdg))
		return 23;
	      else
		{
		  if (e->nof (kaon_pdg))
		    return 44;
		  else if (e->nof (kaon_plus_pdg))
		    return 45;
		}
	    }


	}

    }

  if (e->flag.coh)		// koherentne  oddziaływanie tylko na O(16) 
    {
      Int_t _target;
      _target = e->par.nucleus_p + e->par.nucleus_n;	// liczba masowa  O(16) 

      if (_target == 16)
	{
	  if (e->flag.anty)	// jeśli jest to oddziaływanie z antyneutrinem
	    {
	      if (e->flag.cc && e->nof (pion_minus_pdg))
		return -16;
	      else if (e->nof (pion_pdg))
		return -36;
	    }
	  else			// oddziaływanie z neutrinem
	    {
	      if (e->flag.cc && e->nof (pion_plus_pdg))
		return 16;
	      else if (e->nof (pion_pdg))
		return 36;
	    }
	}
    }

  // gleboko nieelastyczne rozpraszanie               
  if (e->flag.dis)
    {
      if (e->flag.anty)
	{
	  if (e->flag.cc)
	    return -26;
	  else
	    return -46;
	}
      else
	{
	  if (e->flag.cc)
	    return 26;
	  else
	    return 46;
	}
    }

  return 9999;
}




void
ND280_event::Translate_From_NuWro (event * e)
{
  Int_t nie_dotyczyI = -3;
  Float_t nie_dotyczyF = -3.3;
  Int_t nie_wiadomoI = -13;
  Float_t nie_wiadomoF = -13.13;

  Pos[0] = e->in[1].r.x;
  Pos[1] = e->in[1].r.y;
  Pos[2] = e->in[1].r.z;	// wartość 0 w wynikach

  Mode = ND280_Mode (e);
  Int_t Numnu_in = e->in.size ();
  Int_t Numnu_out = e->out.size ();
  Numnu = Numnu_in + Numnu_out;

  for (Int_t j = 0; j < Numnu_in; j++)
    {

      Ipnu[j] = e->in[j].pdg;
      Abspnu[j] = 0.001 * vec (e->in[j].x, e->in[j].y, e->in[j].z).length ();
      Pnu[j][0] = 0.001 * e->in[j].x;
      Pnu[j][1] = 0.001 * e->in[j].y;
      Pnu[j][2] = 0.001 * e->in[j].z;
    }

  for (Int_t j = 0; j < Numnu_out; j++)
    {

      Ipnu[j + Numnu_in] = e->out[j].pdg;
      Abspnu[j + Numnu_in] =
	0.001 * vec (e->out[j].x, e->out[j].y, e->out[j].z).length ();
      Pnu[j + Numnu_in][0] = 0.001 * e->out[j].x;
      Pnu[j + Numnu_in][1] = 0.001 * e->out[j].y;
      Pnu[j + Numnu_in][2] = 0.001 * e->out[j].z;
    }

  Int_t Npar_in = e->in.size ();
  Int_t Npar_out = e->out.size ();
  Npar = Npar_in + Npar_out;

  for (Int_t k = 0; k < Npar_in; k++)
    {
      Ipv[k] = e->out[k].pdg;
      Icrnv[k] = 0;
      Pmomv[k][0] = e->out[k].x;
      Pmomv[k][1] = e->out[k].y;
      Pmomv[k][2] = e->out[k].z;
      Abspv[k] = vec (e->out[k].x, e->out[k].y, e->out[k].z).length ();
      Iflgv[k] = nie_wiadomoI;	// przyjmowane wartości -1,0,2,3,4,7
      Iorgv[k] = nie_wiadomoI;	// przyjmowane wartości 0,1,2,4,7    
    }

  for (Int_t k = 0; k < Npar_out; k++)
    {
      Ipv[k + Npar_in] = e->out[k].pdg;
      if (e->out[k].pdg == 14)
	Icrnv[k + Npar_in] = 0;
      else
	Icrnv[k + Npar_in] = 1;
      Pmomv[k + Npar_in][0] = e->out[k].x;
      Pmomv[k + Npar_in][1] = e->out[k].y;
      Pmomv[k + Npar_in][2] = e->out[k].z;
      Abspv[k + Npar_in] =
	vec (e->out[k].x, e->out[k].y, e->out[k].z).length ();
      Iflgv[k + Npar_in] = nie_wiadomoI;	// przyjmowane wartości -1,0,2,3,4,7
      Iorgv[k + Npar_in] = nie_wiadomoI;	// przyjmowane wartości 0,1,2,4,7    
    }

  xpi[0] = nie_dotyczyF;
  xpi[1] = nie_dotyczyF;
  xpi[2] = nie_dotyczyF;
  npi[0] = nie_dotyczyF;	// nie dotyczy Nuwro
  npi[1] = nie_dotyczyF;	// nie dotyczy NuWro
  npi[2] = nie_dotyczyF;	// nie dotyczy NuWro                                                               
  cospibm = nie_dotyczyF;
  ppi = nie_dotyczyF;
  ppid = nie_dotyczyI;
  xpi0[0] = nie_dotyczyF;
  xpi0[1] = nie_dotyczyF;
  xpi0[2] = nie_dotyczyF;
  npi0[0] = nie_dotyczyF;	// nie dotyczy NuWro
  npi0[1] = nie_dotyczyF;	// nie dotyczy NuWro
  npi0[2] = nie_dotyczyF;	// nie dotyczy NuWro           
  cospi0bm = nie_dotyczyF;
  ppi0 = nie_dotyczyF;
  Enu = nie_wiadomoF;		// nie wiadomo 
  idfd = nie_wiadomoI;		// przypisana wartosc 5
  modef = nie_wiadomoI;		// przypisana wartosc 11 ,12                                                                      
  nnu[0] = nie_wiadomoF;
  nnu[1] = nie_wiadomoF;
  nnu[2] = nie_wiadomoF;
  norm = nie_wiadomoF;
  rnu = nie_wiadomoF;
  xnu = nie_wiadomoF;
  ynu = nie_wiadomoF;

}


void
NuWro_To_ND280 (const char *plik1, const char *plik2)
{

  event *e = new event;
  TFile *ff1 = new TFile (plik1);
  TTree *tt1 = (TTree *) ff1->Get ("treeout");
  tt1->SetBranchAddress ("e", &e);
  int n = tt1->GetEntries ();
  TFile *ff2 = new TFile (plik2, "recreate");
  TTree *tt2 = new TTree ("h10", "NEUT like event tree");
  ND280_event e2;
  e2.Create_Branches (tt2);

  for (int i = 0; i < n; i++)
    {
      tt1->GetEntry (i);
      e2.Translate_From_NuWro (e);
      e2.Nev = i;
      tt2->Fill ();
    }

  ff2->Write ();
  ff2->Close ();
  ff1->Close ();
  //delete tt2;
  delete ff2;
  //delete tt1;
  delete ff1;
}




int
main (int argc, char *argv[])
{

  if (argc == 1)
    {
      NuWro_To_ND280 ("eventsout.root", "ND280.root");
      cout << "Plik wynikowy: \"ND280.root\"" << endl;
    }
  else if (argc == 2)
    {
      NuWro_To_ND280 (argv[1], "ND280.root");
      cout << "Plik wynikowy: \"ND280.root\"" << endl;
    }

  else if (argc == 3)
    {
      NuWro_To_ND280 (argv[1], argv[2]);
      cout << "Plik wynikowy: " << argv[2] << endl;
    }



}
