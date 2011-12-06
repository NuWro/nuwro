#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
//#include "sobek.h"
#include<TMCParticle.h>
#include<TPythia6.h>

#include "generatormt.h"
//routines called from pythia6

#include "pdg_name.h"
#include "parameters.h"
#include "dis_nc.h"
#include "dis_cr_sec.h"
#include <cstdlib>
#include "masses.h"


extern "C" void py3ent_(int *, const int *, const int *, const int*, double *, double *,double *);
extern "C" void py2ent_(int *, const int *, const int *, double *);
extern "C" void py1ent_(int *, const int *, double *, double *, double *);
extern "C" void pylist_(int *);
extern "C" void pyedit_(int *);
extern "C" int  pycomp_(const int *);
extern "C" void pydecy_(int *);
extern "C" void pyerrm_(int *,int *);



//choosing the hit parton and spectator diquark
//written by J.Nowak


/////////////////////Hadronization for neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////

int hit_parton_nc_nu_p(double E, double W, double nu, double m)
{

  double l = frandom();

  double kwark_d    = cr_sec_nc_nu_p_d(E,W,nu,m);
  double kwark_u    = cr_sec_nc_nu_p_u(E,W,nu,m);
  double kwark_dbar = cr_sec_nc_nu_p_dbar(E,W,nu,m);
  double kwark_ubar = cr_sec_nc_nu_p_ubar(E,W,nu,m);
  double kwark_s    = cr_sec_nc_nu_p_s(E,W,nu,m);
  double kwark_sbar = cr_sec_nc_nu_p_sbar(E,W,nu,m);

  double suma = kwark_d + kwark_u + kwark_dbar + kwark_ubar + kwark_s + kwark_sbar;
  
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_u)/suma){return quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s)/suma){return quark_s;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s+kwark_sbar)/suma){return anti_quark_s;}

 else
    {
      cerr<<"Imposible choice of quark 1"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_nc_nu_p(int hit_parton)
{

    double r = frandom();

    if(hit_parton == quark_d){return quark_d;}
    if(hit_parton == quark_u){return quark_u;}
    if(hit_parton == anti_quark_u) {return anti_quark_u;}
    if(hit_parton == anti_quark_d) {return anti_quark_d;}
    if(hit_parton == quark_s){return quark_s;}
    if(hit_parton == anti_quark_s){return anti_quark_s;}

    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_nc_nu_p(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_d && transfer == quark_d){return quark_d;}
  if(hit_parton == quark_u && transfer == quark_u){return quark_u;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_d){return quark_d;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u){return quark_u;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  if(hit_parton == quark_s && transfer == quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  else {
    cerr<<hit_parton<<endl;
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_nc_nu_p(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();
//quark d & anti-d
    if(abs(hit_parton) == quark_d && (quark_mass  + diquark_mass < W)){return diquark_uu_1;}
//quark u & anti-u
    if(abs(hit_parton) == quark_u && (quark_mass  + diquark_mass_0 < W)){return diquark_ud_0;}

// quark s
  
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}

// anti-quark s
  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}
  
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_nc_nu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u){return Kzero;}

    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u){return Kzero;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

double hadronization_nc_nu_p(double E, double W, double nu, double m)
{

//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 =0;

  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  hit_parton = hit_parton_nc_nu_p(E,W,nu,m);
  transfer0 = transfer_nc_nu_p(hit_parton);
  frag_parton = frag_parton_nc_nu_p(hit_parton, transfer0);
  spectator1 = spectator_diquark1_nc_nu_p(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_nc_nu_p(hit_parton,transfer0,frag_parton,spectator1,W1);
/*cout<<"hit_parton="<<hit_parton<<endl;
cout<<"transfer  ="<<transfer0<<endl;
cout<<"frag_parto="<<frag_parton<<endl;
cout<<"diquark   ="<<spectator1<<endl;

cout<<"mezon     ="<<spectator2<<endl;
//cin.get();
*/

  }
  

  if(frag_parton == quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  
  
  
/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////////////QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }

    }	

//////////////////////anti_QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      pylist_(&lista);
}

///antineurtinos proton

int hit_parton_nc_anu_p(double E, double W, double nu, double m)
{

  double l = frandom();

  double kwark_d    = cr_sec_nc_anu_p_d(E,W,nu,m);
  double kwark_u    = cr_sec_nc_anu_p_u(E,W,nu,m);
  double kwark_dbar = cr_sec_nc_anu_p_dbar(E,W,nu,m);
  double kwark_ubar = cr_sec_nc_anu_p_ubar(E,W,nu,m);
  double kwark_s    = cr_sec_nc_anu_p_s(E,W,nu,m);
  double kwark_sbar = cr_sec_nc_anu_p_sbar(E,W,nu,m);

  double suma = kwark_d + kwark_u + kwark_dbar + kwark_ubar + kwark_s + kwark_sbar;
  
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_u)/suma){return quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s)/suma){return quark_s;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s+kwark_sbar)/suma){return anti_quark_s;}

 else
    {
      cerr<<"Imposible choice of quark 2"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_nc_anu_p(int hit_parton)
{

    double r = frandom();

    if(hit_parton == quark_d){return quark_d;}
    if(hit_parton == quark_u){return quark_u;}
    if(hit_parton == anti_quark_u) {return anti_quark_u;}
    if(hit_parton == anti_quark_d) {return anti_quark_d;}
    if(hit_parton == quark_s){return quark_s;}
    if(hit_parton == anti_quark_s){return anti_quark_s;}

    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_nc_anu_p(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_d && transfer == quark_d){return quark_d;}
  if(hit_parton == quark_u && transfer == quark_u){return quark_u;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_d){return quark_d;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u){return quark_u;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  if(hit_parton == quark_s && transfer == quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_nc_anu_p(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();
 

//quark d & anti-d
    if(abs(hit_parton) == quark_d &&  (1.08< W) ){return diquark_uu_1;}
//quark u & anti-u
    if(abs(hit_parton) == quark_u &&  (1.08 < W)){return diquark_ud_0;}

// quark s
  
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}

  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}
  
  else{
    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_nc_anu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u){return Kzero;}

    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u){return Kzero;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

double hadronization_nc_anu_p(double E, double W, double nu, double m)
{

//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 =0;

  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  hit_parton = hit_parton_nc_nu_p(E,W,nu,m);
  transfer0 = transfer_nc_nu_p(hit_parton);
  frag_parton = frag_parton_nc_nu_p(hit_parton, transfer0);
  spectator1 = spectator_diquark1_nc_nu_p(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_nc_nu_p(hit_parton,transfer0,frag_parton,spectator1,W1);

/*cout<<hit_parton<<endl;
cout<<  transfer0 <<endl;//= transfer_cc_anu_p(hit_parton,W, nu);
cout<<  frag_parton<<endl;// = frag_parton_cc_anu_p(hit_parton, transfer0);
cout<<  spectator1 <<endl;//= spectator_diquark1_cc_anu_p(hit_parton, transfer0, frag_parton, W1);
cout<<  spectator2 <<endl;//= spectator_meson_cc_anu_p(hit_parton,transfer0,frag_parton,spectator1,W1);
  */
  
  }
  

  if(frag_parton == quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  
  
  
/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////////////QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass_0 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass_0)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
//	      cout<<"totutututu"<<endl;
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

//////////////////////anti_QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass_0 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass_0)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      pylist_(&lista);
}






///////////NEUTRON///////////////////////



int hit_parton_nc_nu_n(double E, double W, double nu, double m)
{

  double l = frandom();

  double kwark_d    = cr_sec_nc_nu_n_d(E,W,nu,m);
  double kwark_u    = cr_sec_nc_nu_n_u(E,W,nu,m);
  double kwark_dbar = cr_sec_nc_nu_n_dbar(E,W,nu,m);
  double kwark_ubar = cr_sec_nc_nu_n_ubar(E,W,nu,m);
  double kwark_s    = cr_sec_nc_nu_n_s(E,W,nu,m);
  double kwark_sbar = cr_sec_nc_nu_n_sbar(E,W,nu,m);

  double suma = kwark_d + kwark_u + kwark_dbar + kwark_ubar + kwark_s + kwark_sbar;
  
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_u)/suma){return quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s)/suma){return quark_s;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s+kwark_sbar)/suma){return anti_quark_s;}

 else
    {//cout<<"W="<<W<<"   nu="<<nu<<"   m="<<m<<"   E="<<E<<"  suma"<<suma<<endl;
//cout<<kwark_d<<"  "<<kwark_u<<"  "<<kwark_dbar<<"  "<<kwark_ubar<<"  "<<kwark_s<<"  "<<kwark_sbar<<endl;
//cout<< (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu)<<endl;
//cout<< nu / E<<endl;
  

      cerr<<"Impossible choice of quark 3"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_nc_nu_n(int hit_parton)
{

    double r = frandom();

    if(hit_parton == quark_d){return quark_d;}
    if(hit_parton == quark_u){return quark_u;}
    if(hit_parton == anti_quark_u) {return anti_quark_u;}
    if(hit_parton == anti_quark_d) {return anti_quark_d;}
    if(hit_parton == quark_s){return quark_s;}
    if(hit_parton == anti_quark_s){return anti_quark_s;}

    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_nc_nu_n(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_d && transfer == quark_d){return quark_d;}
  if(hit_parton == quark_u && transfer == quark_u){return quark_u;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_d){return quark_d;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u){return quark_u;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s)
	{
	  if(frandom()<1/3.){return quark_s;}
	  else{return quark_u;}
	}

  if(hit_parton == quark_s && transfer == quark_s)
	{
	  if(frandom()<1/3.){return quark_s;}
	  else{return quark_u;}
	}

  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_nc_nu_n(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();


  
  if(hit_parton == quark_d && transfer == quark_d && frag_parton == quark_d && (quark_mass  + diquark_mass < W)){return diquark_ud_0;}
  if(hit_parton == quark_u && transfer == quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass < W)){return diquark_dd_1;}

  if(hit_parton == anti_quark_d && transfer == anti_quark_d && frag_parton == quark_d && (quark_mass  + diquark_mass < W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass < W)){return diquark_dd_1;}

// quark s
  
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_dd_1;}
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_sd_0;}

  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_sd_0;}
  
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_nc_nu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u){return Kzero;}

    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u){return Kzero;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

double hadronization_nc_nu_n(double E, double W, double nu, double m)
{

//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 =0;

  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  hit_parton = hit_parton_nc_nu_n(E,W,nu,m);
  transfer0 = transfer_nc_nu_n(hit_parton);
  frag_parton = frag_parton_nc_nu_n(hit_parton, transfer0);
  spectator1 = spectator_diquark1_nc_nu_n(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_nc_nu_n(hit_parton,transfer0,frag_parton,spectator1,W1);
  }
  

  if(frag_parton == quark_d && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == quark_u && spectator1 == diquark_dd_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_d && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_u && spectator1 == diquark_dd_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  
  
  
/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////////////QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_dd_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_sd_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

//////////////////////anti_QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_dd_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_sd_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
           }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      cout<<"scenariusz="<<scenariusz<<endl;
//      pylist_(&lista);
}

///antineurtinos proton

int hit_parton_nc_anu_n(double E, double W, double nu, double m)
{

  double l = frandom();

  double kwark_d    = cr_sec_nc_anu_n_d(E,W,nu,m);
  double kwark_u    = cr_sec_nc_anu_n_u(E,W,nu,m);
  double kwark_dbar = cr_sec_nc_anu_n_dbar(E,W,nu,m);
  double kwark_ubar = cr_sec_nc_anu_n_ubar(E,W,nu,m);
  double kwark_s    = cr_sec_nc_anu_n_s(E,W,nu,m);
  double kwark_sbar = cr_sec_nc_anu_n_sbar(E,W,nu,m);

  double suma = kwark_d + kwark_u + kwark_dbar + kwark_ubar + kwark_s + kwark_sbar;
  
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_u)/suma){return quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s)/suma){return quark_s;}
  if(l<(kwark_d+kwark_u+kwark_dbar+kwark_ubar+ kwark_s+kwark_sbar)/suma){return anti_quark_s;}

 else
    {
      cerr<<"Imposible choice of quark 4"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_nc_anu_n(int hit_parton)
{

    double r = frandom();

    if(hit_parton == quark_d){return quark_d;}
    if(hit_parton == quark_u){return quark_u;}
    if(hit_parton == anti_quark_u) {return anti_quark_u;}
    if(hit_parton == anti_quark_d) {return anti_quark_d;}
    if(hit_parton == quark_s){return quark_s;}
    if(hit_parton == anti_quark_s){return anti_quark_s;}

    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_nc_anu_n(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_d && transfer == quark_d){return quark_d;}
  if(hit_parton == quark_u && transfer == quark_u){return quark_u;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_d){return quark_d;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u){return quark_u;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  if(hit_parton == quark_s && transfer == quark_s)
	{
	  if(frandom()<2/3.){return quark_s;}
	  else{return quark_u;}
	}

  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}

//    cout<<"hit_parton="<<hit_parton<<endl;

}


int spectator_diquark1_nc_anu_n(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();

  if(hit_parton == quark_d && transfer == quark_d && frag_parton == quark_d && (quark_mass  + diquark_mass < W)){return diquark_ud_0;}
  if(hit_parton == quark_u && transfer == quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass < W)){return diquark_dd_1;}

  if(hit_parton == anti_quark_d && transfer == anti_quark_d && frag_parton == quark_d && (quark_mass  + diquark_mass < W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass < W)){return diquark_dd_1;}

// quark s
  
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}

  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s && (quark_s_mass  + diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass  + diquark_s_mass + Kzero_mass + 0.1*l < W)){return diquark_su_0;}
  
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_nc_anu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == quark_s && transfer == quark_s && frag_parton == quark_u){return Kzero;}

    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_s){return Kplus;}
    if(hit_parton == anti_quark_s && transfer == anti_quark_s && frag_parton == quark_u){return Kzero;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

double hadronization_nc_anu_n(double E, double W, double nu, double m)
{

//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 =0;

  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  hit_parton = hit_parton_nc_nu_n(E,W,nu,m);
  transfer0 = transfer_nc_nu_n(hit_parton);
  frag_parton = frag_parton_nc_nu_n(hit_parton, transfer0);
  spectator1 = spectator_diquark1_nc_nu_n(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_nc_nu_n(hit_parton,transfer0,frag_parton,spectator1,W1);

/*cout<<"hit_parton ="<<hit_parton<<endl;
cout<<"frag_parton="<<frag_parton<<endl;
cout<<"diquark    ="<<spectator1<<endl;
cout<<"meson      ="<<spectator2<<endl;
*/
  }
//cout<<"W="<<W<<endl;
//cin.get();  

if(spectator2 == 0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}

/*  if(frag_parton == quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_d && spectator1 == diquark_uu_1){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
  if(frag_parton == anti_quark_u && spectator1 == diquark_ud_0){py2ent_(&ip, &frag_parton, &spectator1, &W1);}
*/
  
  
/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////////////QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass_0)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

//////////////////////anti_QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_s_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass_0)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	}
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_s 
    && frag_parton == quark_u && spectator1 == diquark_su_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_s_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          
    }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      cout<<"cenariusz="<<scenariusz<<endl;
//      pylist_(&lista);
}



