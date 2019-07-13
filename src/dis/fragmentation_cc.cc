#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
//#include "sobek.h"
#include <TMCParticle.h>
#include <TPythia6.h>

#include "generatormt.h"
//#include "hit_parton.nc.h"

#include "pdg_name.h"
#include "parameters.h"
#include "dis_cc_neutron.h"
#include "dis_cc_proton.h"
#include "masses.h"
#include "dis_cr_sec.h"


//routines called from pythia6

extern "C" void py3ent_(int *, const int *, const int *, const int*, double *, double *,double *);
extern "C" void py2ent_(int *, const int *, const int *, double *);
extern "C" void py1ent_(int *, const int *, double *, double *, double *);
extern "C" void pylist_(int *);
extern "C" void pyedit_(int *);
extern "C" int  pycomp_(const int *);
extern "C" void pydecy_(int *);
extern "C" void pyerrm_(int *,int *);

int scenariusz=-1;

//choosing the hit parton and spectator diquark
//written by J.Nowak


/////////////////////Hadronization for neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////

int hit_parton_cc_nu_p(double E, double W, double nu, double m)
{

  double l = frandom();
  double kwark_d = cr_sec_cc_nu_p_d_BY(E,W,nu,m);
  double kwark_sb = cr_sec_cc_nu_p_sb_BY(E,W,nu,m);
  double kwark_ubar = cr_sec_cc_nu_p_ubar_BY(E,W,nu,m);

if(kwark_d<0)kwark_d=0;
if(kwark_sb<0)kwark_sb=0;
if(kwark_ubar<0)kwark_ubar=0;

  double suma = kwark_d + kwark_sb + kwark_ubar;


if(kwark_d +kwark_sb + kwark_ubar<0)
{  
cout<<"W="<<W<<" nu="<<nu<<" kwark_d="<<kwark_d<<" kwark_s="<<kwark_sb<<" kwark_ubar="<<kwark_ubar<<endl;
//cin.get();
} 
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_ubar+kwark_sb)/suma){return quark_s;}

 else
    {
      cerr<<"Impossible choice of quark 5"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_cc_nu_p(int hit_parton, double W, double nu)
{
    double Q2 = M2 + 2*M12*nu - W*W;
    double xx = Q2/(2*M12*nu);
        
    double r = frandom();

    if(hit_parton == anti_quark_u && r < sin_2_theta_C) {return anti_quark_d;}
    if(hit_parton == anti_quark_u && r > sin_2_theta_C) {return anti_quark_s;}

if(hit_parton==quark_d && x_d2c(W,nu)< xx) {return quark_u;}
if(hit_parton==quark_d && x_d2c(W,nu)> xx && r > sin_2_theta_C){return quark_u;}
if(hit_parton==quark_d && x_d2c(W,nu)> xx && r < sin_2_theta_C){return quark_c;}


if(hit_parton == quark_s && x_s2c(W,nu)<xx){return quark_u;}
if(hit_parton == quark_s && x_s2c(W,nu)>xx && r < sin_2_theta_C){return quark_u;}
if(hit_parton == quark_s && x_s2c(W,nu)>xx && r > sin_2_theta_C){return quark_c;}


    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_cc_nu_p(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  //nu + d (uu) --> mu + u (uu) 
  //interactiong quark is quark d (valence or from sea) and outgoing quark is quark u or c
  
  if(hit_parton == quark_d && transfer == quark_c){return quark_c;}
  if(hit_parton == quark_d && transfer == quark_u){return quark_u;}

  //nu + u_bar (uuud) --> mu + dbar (uuud) --> mu + u (uu)
  //nu + u_bar (uuud) --> mu + sbar (uuud) --> mu + (sbar u) + u (ud) || mu + (sbar u) + d (uu) || mu + (sbar d) + u (uu)
  
  if(hit_parton == anti_quark_u && transfer == anti_quark_d){return quark_u;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s)
	{
	  if(frandom()<0.25){return quark_d;}
	  else{return quark_u;}
	}

  //nu + s (sbar uud) --> mu + u + (sbar u) + (ud) || mu +u + (sbar d) + (uu)
  
  if(hit_parton == quark_s && transfer == quark_c)
	{
	    double l=frandom();
	  if(l<0.5){return quark_u;}
	  if(l<0.75){return quark_c;}    
	  else{return quark_d;}
	}
  if(hit_parton == quark_s && transfer == quark_u)
	{
	  if(frandom()<0.25){return quark_d;}
	  else{return quark_u;}
	}
  
  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_cc_nu_p(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();
//quark d
//cout<<W<<" "<<quark_d_mass+diquark_mass<<endl;
  if(hit_parton == quark_d && transfer == quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass <= W)){return diquark_uu_1;}
  if(hit_parton == quark_d && transfer == quark_c && frag_parton == quark_c && (W>3.+ 0.1*l )){return diquark_uu_1;}
  
//antiquark u  

  if(hit_parton == anti_quark_u && transfer == anti_quark_d && frag_parton == quark_u && (quark_mass + diquark_mass<W)){return diquark_uu_1;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_d && (quark_mass + diquark_mass + Kplus_mass + 0.1*l<W)){return diquark_uu_1;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass + diquark_mass + Kzero_mass +0.1*l<W) ){return diquark_uu_1;}

// quark s
  
  if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_u && (quark_mass  +diquark_mass + Kzero_mass + 0.1*l < W)){return diquark_uu_1;}
  if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_d && (quark_mass  +diquark_mass + Kplus_mass + 0.1*l < W)){return diquark_uu_1;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_u && (W>3.+0.1*l)){return diquark_cd_0;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_c && (W>3.+0.1*l)){return diquark_uu_1;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_d && (quark_mass  +diquark_mass +Dsplus_mass +0.1*l<W)){return diquark_uu_1;}
  
    
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_cc_nu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_d) {return Kplus;}
	
    if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_u) {return Kzero;}
	
//quark s

    if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_d){return Kplus;}
	
    if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_u){return Kzero;}
    
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_u){return Kplus;}
	
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_c){return Kzero;}
    
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_d
    	&& (W> 3.5 )){return Dsplus;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

void hadronization_cc_nu_p(double E, double W, double nu, double m)
{

//int scenariusz=-1;


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
  hit_parton = hit_parton_cc_nu_p(E,W,nu,m);
  transfer0 = transfer_cc_nu_p(hit_parton,W,nu);
  frag_parton = frag_parton_cc_nu_p(hit_parton, transfer0);
  spectator1 = spectator_diquark1_cc_nu_p(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_cc_nu_p(hit_parton,transfer0,frag_parton,spectator1,W1);
  }
  

  if(frag_parton == quark_u && spectator1 == diquark_uu_1)
    {
//      cout<<"Interaction type u + (uu)"<<endl;
      py2ent_(&ip, &frag_parton, &spectator1, &W1);

if(hit_parton == 1){scenariusz =1;}
if(hit_parton == 2){scenariusz =0;}

    }

  if(frag_parton == quark_c && spectator1 == diquark_uu_1 && W1>3.0)
    {
//      cout<<"Interaction type c + (uu)"<<endl;
      py2ent_(&ip, &frag_parton, &spectator1, &W1);
          scenariusz = 2;
    }

 

/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////ANTIQUARK U///////////////////////////////////////  
  if(hit_parton == anti_quark_u && transfer0 == anti_quark_s
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

/*cout<<"c2="<<cthe2<<endl;
cout<<"c3="<<cthe3<<endl;
cout<<"W="<<W1<<endl;
*/
	      
	    }
          scenariusz = 3;
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_u && transfer0 == anti_quark_s 
    && frag_parton == quark_u && spectator1 == diquark_uu_1 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	      
//	      cout<<"TUTAJJJJJJJJJJJJJJJJJJJ"<<endl;
	      
	    }
          scenariusz = 4;
    }	

//////////////////////QUARK S///////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_u
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
          scenariusz =5;
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == quark_s && transfer0 == quark_u 
    && frag_parton == quark_u && spectator1 == diquark_uu_1 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
          scenariusz = 6;
    }	

///////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_u && spectator1 == diquark_cd_0 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_c_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_c_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
          scenariusz = 7;
    }	

///////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_c && spectator1 == diquark_uu_1 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_c_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_c_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
          scenariusz = 8;
    }	

///////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Dsplus )
    {
      while(x2*W1/2 < 2.1 || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dsplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 


}

    scenariusz = 9;
}
      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      cout<<"scenariusz="<<scenariusz<<endl;
//      pylist_(&lista);

}


/////////////////////Hadronization for anti-neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////
///////////////do sprawdzenia kat cabibo

int hit_parton_cc_anu_p(double E, double W, double nu, double m)
{

  double l = frandom();

  double kwark_u = cr_sec_cc_anu_p_u_BY(E,W,nu,m);
  double kwark_sb = cr_sec_cc_anu_p_sb_BY(E,W,nu,m);
  double kwark_dbar = cr_sec_cc_anu_p_dbar_BY(E,W,nu,m);

if(kwark_u<0)kwark_u=0;
if(kwark_sb<0)kwark_sb=0;
if(kwark_dbar<0)kwark_dbar=0;

  double suma = kwark_u + kwark_sb + kwark_dbar;

if(kwark_u<0 || kwark_sb<0 || kwark_dbar<0)
{
cout<<"W="<<W<<" nu="<<nu<<" kwark_u="<<kwark_u<<" kwark_s="<<kwark_sb<<" kwark_dbar="<<kwark_dbar<<endl;
} 

  if(l<kwark_u/suma){return quark_u;}
  if(l<(kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_u+kwark_dbar+kwark_sb)/suma){return anti_quark_s;}

 else
    {
      cerr<<"Impossible choice of quark 6"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_cc_anu_p(int hit_parton, double W, double nu)
{
    double Q2 = M2 + 2*M12*nu - W*W;
    double xx = Q2/(2*M12*nu);

    double r = frandom();

if(hit_parton == anti_quark_d && x_d2c(W,nu)<xx){return anti_quark_u;}
if(hit_parton == anti_quark_d && x_d2c(W,nu)>xx && r > sin_2_theta_C) {return anti_quark_u;}
if(hit_parton == anti_quark_d && x_d2c(W,nu)>xx && r < sin_2_theta_C) {return anti_quark_c;}

if(hit_parton == quark_u && r > sin_2_theta_C){return quark_d;}
if(hit_parton == quark_u && r < sin_2_theta_C){return quark_s;}

if(hit_parton == anti_quark_s && x_s2c(W,nu)<xx){return anti_quark_u;}
if(hit_parton == anti_quark_s && x_s2c(W,nu)>xx && r < sin_2_theta_C){return anti_quark_u;}
if(hit_parton == anti_quark_s && x_s2c(W,nu)>xx && r > sin_2_theta_C){return anti_quark_c;}

    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_cc_anu_p(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_u && transfer == quark_s){return quark_s;}
  if(hit_parton == quark_u && transfer == quark_d){return quark_d;}


  if(hit_parton == anti_quark_s && transfer == anti_quark_u){return quark_s;}

  
  if(hit_parton == anti_quark_d && transfer == anti_quark_u){return quark_u;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_c)
	{
	  if(frandom()<0.5){return quark_u;}
	  else{return quark_d;}
	}

  
  if(hit_parton == anti_quark_s && transfer == anti_quark_c)
	{
	    double l=frandom();
	  if(l<0.5){return quark_u;}
	  if(l<0.75){return quark_s;}    
	  else{return quark_d;}
	}
  
  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_cc_anu_p(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();
//quark u
//cout<<W<<" "<<quark_d_mass+diquark_mass<<endl;

  if(hit_parton == quark_u && transfer == quark_s && frag_parton == quark_s && (quark_mass + diquark_mass + Kplus_mass  +0.1*l<= W)){return diquark_ud_1;}
  if(hit_parton == quark_u && transfer == quark_d && frag_parton == quark_d && (1.08 <= W)){return diquark_ud_0;}
  
//antiquark d

  if(hit_parton == anti_quark_d && transfer == anti_quark_u && frag_parton == quark_u && (1.08 <= W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_u && (3.0 + 0.1*l<W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_d && (3.0 + 0.1*l<W) ){return diquark_uu_1;}

// antiquark s
  
  if(hit_parton == anti_quark_s && transfer == anti_quark_u && frag_parton == quark_s && (quark_mass + diquark_mass + Kplus_mass +0.1*l <= W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_u && (3.0 + 0.1*l < W)){return diquark_sd_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_s && (3.0 +0.1*l <=W)){return diquark_uu_1;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_d && (3.0 +0.1*l<= W)){return diquark_uu_1;}

    
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_cc_anu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_u) {return Dzero_bar;}
	
    if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_d) {return Dminus;}
	
//quark s

    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_u){return Dzero_bar;}
	
    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_s){return Dminus;}
    
    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_d){return Dsminus;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}



void hadronization_cc_anu_p(double E, double W, double nu, double m)
{
//int scenariusz=-1;
//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 =0;

//cout<<"tutaj"<<endl;

  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  hit_parton = hit_parton_cc_anu_p(E,W,nu,m);
  transfer0 = transfer_cc_anu_p(hit_parton,W, nu);
  frag_parton = frag_parton_cc_anu_p(hit_parton, transfer0);
  spectator1 = spectator_diquark1_cc_anu_p(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_cc_anu_p(hit_parton,transfer0,frag_parton,spectator1,W1);
/*
cout<<hit_parton<<endl;
cout<<  transfer0 <<endl;//= transfer_cc_anu_p(hit_parton,W, nu);
cout<<  frag_parton<<endl;// = frag_parton_cc_anu_p(hit_parton, transfer0);
cout<<  spectator1 <<endl;//= spectator_diquark1_cc_anu_p(hit_parton, transfer0, frag_parton, W1);
cout<<  spectator2 <<endl;//= spectator_meson_cc_anu_p(hit_parton,transfer0,frag_parton,spectator1,W1);
*/
  }




  

if(spectator2 == 0)
{
  py2ent_(&ip, &frag_parton, &spectator1, &W1);

}

/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////ANTIQUARK d///////////////////////////////////////  
  if(hit_parton == anti_quark_d && transfer0 == anti_quark_c
    && frag_parton == quark_u && spectator1 == diquark_dd_1 &&  spectator2 == Dzero_bar )
    {
      while(x2*W1/2 < Dzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	      
	    }
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_d && transfer0 == anti_quark_c 
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Dminus )
    {
      while(x2*W1/2 < Dplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	      
	      
	    }
    }	

//////////////////////AntiQUARK S///////////////////////////////  

  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_u && spectator1 == diquark_sd_0 &&  spectator2 == Dzero_bar )
    {
      while(x2*W1/2 < Dzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_s && spectator1 == diquark_uu_1 &&  spectator2 == Dminus )
    {
      while(x2*W1/2 < Dplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
    }	

///////////////////////////////////////////
  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Dsminus )
    {
      while(x2*W1/2 < Dsplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dsplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
    }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 } 
//      pylist_(&lista);
}






/////////////////////////NEUTRON/////////////////////////////////////////

int hit_parton_cc_nu_n(double E, double W, double nu, double m)
{

  double l = frandom();
  double kwark_d = 0;
  double kwark_sb =0;
  double kwark_ubar =0;
  double suma = kwark_d + kwark_sb + kwark_ubar;

  kwark_d = cr_sec_cc_nu_n_d_BY(E,W,nu,m);
  kwark_sb = cr_sec_cc_nu_n_sb_BY(E,W,nu,m);
  kwark_ubar = cr_sec_cc_nu_n_ubar_BY(E,W,nu,m);

if(kwark_d<0)kwark_d=0;
if(kwark_sb<0)kwark_sb=0;
if(kwark_ubar<0)kwark_ubar=0;

  suma = kwark_d + kwark_sb + kwark_ubar;

if(kwark_d<0 || kwark_sb<0 || kwark_ubar<0)
{
cout<<"W="<<W<<" nu="<<nu<<" kwark_d="<<kwark_d<<" kwark_s="<<kwark_sb<<" kwark_ubar="<<kwark_ubar<<endl;
//cin.get();
}
 
  if(l<kwark_d/suma){return quark_d;}
  if(l<(kwark_d+kwark_ubar)/suma){return anti_quark_u;}
  if(l<(kwark_d+kwark_ubar+kwark_sb)/suma){return quark_s;}

 else
    {
      cerr<<"Impossible choice of quark 7"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_cc_nu_n(int hit_parton, double W, double nu)
{
    double Q2 = M2 + 2*M12*nu - W*W;
    double xx = Q2/(2*M12*nu);

    double r = frandom();


if(hit_parton == anti_quark_u && r < sin_2_theta_C) {return anti_quark_d;}
if(hit_parton == anti_quark_u && r > sin_2_theta_C) {return anti_quark_s;}
    
    
if(hit_parton == quark_d && x_d2c(W,nu)<xx) {return quark_u;}
if(hit_parton == quark_d && x_d2c(W,nu)>xx && r > sin_2_theta_C){return quark_u;}
if(hit_parton == quark_d && x_d2c(W,nu)>xx && r < sin_2_theta_C){ return quark_c;}


if(hit_parton == quark_s && x_s2c(W,nu)<xx) {return quark_u;}
if(hit_parton == quark_s && x_s2c(W,nu)>xx && r < sin_2_theta_C){return quark_u;}
if(hit_parton == quark_s && x_s2c(W,nu)>xx && r > sin_2_theta_C){return quark_c;}


    else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;}
    
    
}

int frag_parton_cc_nu_n(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  //nu + d (du) --> mu + u (du) or   mu + c (du)
  //interactiong quark is quark d (valence or from sea) and outgoing quark is quark u or c
  
  if(hit_parton == quark_d && transfer == quark_c){return quark_c;}
  if(hit_parton == quark_d && transfer == quark_u){return quark_u;}

  //nu + u_bar (uuud) --> mu + dbar (uuud) --> mu + u (uu)
  //nu + u_bar (uuud) --> mu + sbar (uuud) --> mu + (sbar u) + u (ud) || mu + (sbar u) + d (uu) || mu + (sbar d) + u (uu)
  
  if(hit_parton == anti_quark_u && transfer == anti_quark_d){return quark_u;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s)
	{
	  if(frandom()<0.5){return quark_d;}
	  else{return quark_u;}
	}

  //nu + s (sbar uud) --> mu + u + (sbar u) + (ud) || mu +u + (sbar d) + (uu)
  
  if(hit_parton == quark_s && transfer == quark_c)
	{
	    double l=frandom();
	  if(l<0.5){return quark_d;}
	  if(l<0.75){return quark_c;}    
	  else{return quark_u;}
	}
  if(hit_parton == quark_s && transfer == quark_u)
	{
	  if(frandom()<0.5){return quark_d;}
	  else{return quark_u;}
	}
  
  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_cc_nu_n(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}

//quark d
//cout<<W<<" "<<quark_d_mass+diquark_mass<<endl;
  if(hit_parton == quark_d && transfer == quark_u && frag_parton == quark_u && (quark_mass  + diquark_mass_0<=W)){return diquark_ud_0;}
  if(hit_parton == quark_d && transfer == quark_c && frag_parton == quark_c && (W>3.)){return diquark_ud_0;}
  
//antiquark u  

  if(hit_parton == anti_quark_u && transfer == anti_quark_d && frag_parton == quark_u && (quark_mass + diquark_mass_0<=W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_d && (quark_mass + diquark_mass + Kzero_mass + 0.005 <=W)){return diquark_uu_1;}
  if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_u && (quark_mass + diquark_mass + Kplus_mass + 0.005 <=W)){return diquark_dd_1;}
    
// quark s
  
  if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_u && (quark_mass  +diquark_mass + Kplus_mass + 0.005 <W)){return diquark_dd_1;}
  if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_d && (quark_mass  +diquark_mass + Kzero_mass + 0.005 <W)){return diquark_uu_1;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_u && (W>3.)){return diquark_cd_0;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_c && (W>3.)){return diquark_dd_1;}
  if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_d && (W>3.25)){return diquark_ud_0;}
  
    
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_cc_nu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_d 
	&& (quark_mass + diquark_mass + Kzero_mass < W )) {return Kzero;}
	
    if(hit_parton == anti_quark_u && transfer == anti_quark_s && frag_parton == quark_u
    	&& (quark_mass + diquark_mass + Kplus_mass < W )) {return Kplus;}
	
//quark s

    if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_d
    	&& (quark_mass + diquark_mass + Kzero_mass < W )){return Kzero;}
	
    if(hit_parton == quark_s && transfer == quark_u && frag_parton == quark_u
    	&& (quark_mass + diquark_mass + Kplus_mass < W )){return Kplus;}
    
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_u
    	&& ( W>3. )){return Kzero;}
	
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_c
    	&& ( W>3. )){return Kplus;}
    
    if(hit_parton == quark_s && transfer == quark_c && frag_parton == quark_d
    	&& (W> 3.5 )){return Dsplus;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}

void hadronization_cc_nu_n(double E, double W, double nu, double m)
{

//int scenariusz=-1;
//number of line of first parton in py2ent_ routine

    int ip = 0;
    int lista = 1;

  double W1 = W/1000;
  int hit_parton = 0;
  int transfer0 = 0;
  int frag_parton = 0;
  int spectator1 = 0;
  int spectator2 = 0;
  int dmeson = 0;
        
  while(hit_parton == 0 || transfer0 == 0 || frag_parton == 0 || spectator1 == 0)
  {
  
  hit_parton = hit_parton_cc_nu_n(E,W,nu,m);
  transfer0 = transfer_cc_nu_n(hit_parton, W,nu);
  frag_parton = frag_parton_cc_nu_n(hit_parton, transfer0);
  spectator1 = spectator_diquark1_cc_nu_n(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_cc_nu_n(hit_parton,transfer0,frag_parton,spectator1,W1);


}    



  if(frag_parton == quark_u && spectator1 == diquark_ud_0)
    {
      py2ent_(&ip, &frag_parton, &spectator1, &W1);

if(hit_parton == 1){scenariusz =1;}
if(hit_parton == 2){scenariusz =0;}

    }

  if(frag_parton == quark_c && spectator1 == diquark_ud_0 && W1>3.0)
    {
      py2ent_(&ip, &frag_parton, &spectator1, &W1);

scenariusz = 2;

    }

 

/////////////Additional meson///////////////////////////////////////////////////////

if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////ANTIQUARK U////////////////////////////////////////////////////////////

  if(hit_parton == anti_quark_u && transfer0 == anti_quark_s
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
scenariusz = 3;
      
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_u && transfer0 == anti_quark_s 
    && frag_parton == quark_u && spectator1 == diquark_dd_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
      scenariusz = 4;

    }	

//////////////////////QUARK S///////////////////////////////////////////////////////////


  if(hit_parton == quark_s && transfer0 == quark_u
    && frag_parton == quark_d && spectator1 == diquark_uu_1 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
      scenariusz = 5;

    }	

////////////////////////////////////////////////////////////////////////////////////////  

  if(hit_parton == quark_s && transfer0 == quark_u 
    && frag_parton == quark_u && spectator1 == diquark_dd_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
      scenariusz = 6;

    }	

/////////////////////////////////////////////////////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_u && spectator1 == diquark_cd_0 &&  spectator2 == Kzero )
    {
      while(x2*W1/2 < Kzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_c_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_c_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
      scenariusz = 7;

    }	

////////////////////////////////////////////////////////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_c && spectator1 == diquark_dd_1 &&  spectator2 == Kplus )
    {
      while(x2*W1/2 < Kplus_mass || x1*W1/2 < quark_c_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_c_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Kplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
      scenariusz = 8;

    }	

/////////////////////////////////////////////////////////////////////////////////////////////

  if(hit_parton == quark_s && transfer0 == quark_c 
    && frag_parton == quark_d && spectator1 == diquark_ud_0 &&  spectator2 == Dsplus )
    {
      while(x2*W1/2 < 2.1 || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass_0 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(kwad(x1*W1/2)-kwad(quark_mass));
	      pa2 = sqrt(kwad(x2*W1/2)-kwad(Dsplus_mass));
	      pa3 = sqrt(kwad(x3*W1/2)-kwad(diquark_mass_0));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 


}    

scenariusz = 9;

}	

py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
} 

}


/////////////////////Hadronization for anti-neutrinos////////////////////////
/////////////////////////Neutron/////////////////////////////////////////
///////////////do sprawdzenia kat cabibo

int hit_parton_cc_anu_n(double E, double W, double nu, double m)
{

  double l = frandom();
  double kwark_u = 0;
  double kwark_sb =0;
  double kwark_dbar =0;
  double suma = kwark_u + kwark_sb + kwark_dbar;

  kwark_u = cr_sec_cc_anu_n_u_BY(E,W,nu,m);
  kwark_sb = cr_sec_cc_anu_n_sb_BY(E,W,nu,m);
  kwark_dbar = cr_sec_cc_anu_n_dbar_BY(E,W,nu,m);

if(kwark_u<0)kwark_u=0;
if(kwark_sb<0)kwark_sb=0;
if(kwark_dbar<0)kwark_dbar=0;

  
  suma = kwark_u + kwark_sb + kwark_dbar;

if(kwark_u<0 || kwark_sb<0 || kwark_dbar<0){cout<<"W="<<W<<" nu="<<nu<<" kwark_d="<<kwark_u<<" kwark_s="<<kwark_sb<<" kwark_ubar="<<kwark_dbar<<endl;}

  
  if(l<kwark_u/suma){return quark_u;}
  if(l<(kwark_u+kwark_dbar)/suma){return anti_quark_d;}
  if(l<(kwark_u+kwark_dbar+kwark_sb)/suma){return anti_quark_s;}

 else
    {
      cerr<<"Imposible choice of quark 8"<<endl;
      cin.get();
      return 0;
    }
}

int transfer_cc_anu_n(int hit_parton, double W, double nu)
{



    double Q2 = M2 + 2*M12*nu - W*W;
    double xx = Q2/(2*M12*nu);

    double r = frandom();

//cout<<"tutajd "<<x_d2c(W,nu)<<" "<<xx<<endl;
//cout<<"tutajs "<<x_s2c(W,nu)<<" "<<xx<<endl;

if(hit_parton == anti_quark_d && x_d2c(W,nu)<xx) {return anti_quark_u;}
if(hit_parton == anti_quark_d && x_d2c(W,nu)>xx && r > sin_2_theta_C) {return anti_quark_u;}
if(hit_parton == anti_quark_d && x_d2c(W,nu)>xx && r < sin_2_theta_C) {return anti_quark_c;}


if(hit_parton == quark_u && r > sin_2_theta_C){return quark_d;}
if(hit_parton == quark_u && r < sin_2_theta_C){return quark_s;}


if(hit_parton == anti_quark_s && x_s2c(W,nu)<xx ) {return anti_quark_u;}
if(hit_parton == anti_quark_s && x_s2c(W,nu)>xx && r < sin_2_theta_C){return anti_quark_u;}
if(hit_parton == anti_quark_s && x_s2c(W,nu)>xx && r > sin_2_theta_C){return anti_quark_c;}

//else{
    cerr<<"Implsible flavour transition"<<endl;
    return 0;
//    }
    

}

int frag_parton_cc_anu_n(int hit_parton, int transfer)
{
  
  if (hit_parton == 0) {return 0;}
  
  
  if(hit_parton == quark_u && transfer == quark_s){return quark_s;}
  if(hit_parton == quark_u && transfer == quark_d){return quark_d;}


  if(hit_parton == anti_quark_s && transfer == anti_quark_u){return quark_s;}

  
  if(hit_parton == anti_quark_d && transfer == anti_quark_u){return quark_d;}
  
  if(hit_parton == anti_quark_d && transfer == anti_quark_c)
	{
	  if(frandom()<0.75){return quark_d;}
	  else{return quark_u;}
	}

  
  if(hit_parton == anti_quark_s && transfer == anti_quark_c)
	{
	  double l=frandom();
	  if(l<0.5){return quark_d;}
	  if(l<0.75){return quark_s;}    
	  else{return quark_u;}
	}
  
  else {
    cerr<<"Non parton for fragmentation"<<endl;
    return 0;}
}


int spectator_diquark1_cc_anu_n(int hit_parton, int transfer, int frag_parton, double W)
{
  
  if(hit_parton == 0) {return 0;}
double l=frandom();
//quark u
//cout<<W<<" "<<quark_d_mass+diquark_mass<<endl;

  if(hit_parton == quark_u && transfer == quark_s && frag_parton == quark_s && (quark_mass + diquark_mass + Kplus_mass + 0.1*l <= W)){return diquark_dd_1;}
  if(hit_parton == quark_u && transfer == quark_d && frag_parton == quark_d && (quark_mass  + diquark_mass <= W)){return diquark_dd_1;}
  
//antiquark d

  if(hit_parton == anti_quark_d && transfer == anti_quark_u && frag_parton == quark_d && (quark_mass + diquark_mass <= W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_d && (quark_mass + diquark_mass + Dzero_mass + 0.1*l<W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_u && (quark_mass + diquark_mass + Dplus_mass + 0.1*l<W) ){return diquark_dd_1;}

// antiquark s
  
  if(hit_parton == anti_quark_s && transfer == anti_quark_u && frag_parton == quark_s && (quark_mass + diquark_mass + Kplus_mass + 0.1*l <= W)){return diquark_dd_1;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_d && (quark_mass + diquark_s_mass + Dzero_mass + 0.1*l <=W)){return diquark_sd_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_s && (quark_s_mass + diquark_mass + Dplus_mass + 0.1*l <=W)){return diquark_ud_0;}
  if(hit_parton == anti_quark_s && transfer == anti_quark_c && frag_parton == quark_u && (quark_mass + diquark_mass + Dsplus_mass +0.1*l <=W)){return diquark_ud_0;}

    
  else{
//    cerr<<"error in first spectator choice"<<endl;
    return 0;
   }
}


int spectator_meson_cc_anu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W)
{
  
  
    if(hit_parton == 0) {return 0;}
    
    if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_d) {return Dzero_bar;}
	
    if(hit_parton == anti_quark_d && transfer == anti_quark_c && frag_parton == quark_u) {return Dminus;}
	
//quark s

    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_d){return Dzero_bar;}
	
    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_s){return Dminus;}
    
    if(hit_parton == anti_quark_s && transfer ==  anti_quark_c && frag_parton == quark_d){return Dsminus;}
    
else
{
//  cerr<<"error in second spectator choice"<<endl;
  return 0;
}
}



void hadronization_cc_anu_n(double E, double W, double nu, double m)
{

//int scenariusz=-1;
//number of line of first parton in py2ent_ routine
//    int scenariusz=-1;

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
  hit_parton = hit_parton_cc_anu_n(E,W,nu,m);
  transfer0 = transfer_cc_anu_n(hit_parton,W,nu);
  frag_parton = frag_parton_cc_anu_n(hit_parton, transfer0);
  spectator1 = spectator_diquark1_cc_anu_n(hit_parton, transfer0, frag_parton, W1);
  spectator2 = spectator_meson_cc_anu_n(hit_parton,transfer0,frag_parton,spectator1,W1);
  }
  

if(spectator2 == 0)
{

if(hit_parton == 2){scenariusz =1;}
if(hit_parton == -1){scenariusz =0;}
//cout<<"TUTAJ 2"<<endl;
//cout<<"scen ="<<scenariusz<<endl;

  py2ent_(&ip, &frag_parton, &spectator1, &W1);
}

/////////////Additional meson////////////////////// 
if(spectator2 != 0)
{  
      double x1 = 0, x2 = 0, x3 = 0;
      double pa1 = 0, pa2 =0, pa3 = 0;
      double cthe2 = 10, cthe3 = 10;

//////////////ANTIQUARK d///////////////////////////////////////  
  if(hit_parton == anti_quark_d && transfer0 == anti_quark_c
    && frag_parton == quark_d && spectator1 == diquark_dd_1 &&  spectator2 == Dzero_bar )
    {
      while(x2*W1/2 < Dzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	      
	    }
scenariusz = 2;
    }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_d && transfer0 == anti_quark_c 
    && frag_parton == quark_u && spectator1 == diquark_dd_1 &&  spectator2 == Dminus )
    {
      while(x2*W1/2 < Dplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	      
	      
	    }
scenariusz = 3;
    }	

//////////////////////AntiQUARK S///////////////////////////////  

  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_d && spectator1 == diquark_sd_0 &&  spectator2 == Dzero_bar )
    {
      while(x2*W1/2 < Dzero_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dzero_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_s_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
scenariusz = 4; 
   }	

/////////////////////////////////////////////////////  
  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_s && spectator1 == diquark_ud_0 &&  spectator2 == Dminus )
    {
      while(x2*W1/2 < Dplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_s_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass_0)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 

	    }
scenariusz = 5;
    }	

///////////////////////////////////////////
  if(hit_parton == anti_quark_s && transfer0 == anti_quark_c
    && frag_parton == quark_u && spectator1 == diquark_ud_0 &&  spectator2 == Dsminus )
    {
      while(x2*W1/2 < Dsplus_mass || x1*W1/2 < quark_mass || x3*W1/2 < diquark_mass 
		|| cthe2 > 1 || cthe2 < -1 || cthe3 > 1 || cthe3 < -1)
	    {
	      x2 = 2*hadr();
	      x1 = 2*hadr();
	      x3 = 2-x1-x2;
	      
	      pa1 = sqrt(max2(1e-10,kwad(x1*W1/2)-kwad(quark_mass)));
	      pa2 = sqrt(max2(1e-10,kwad(x2*W1/2)-kwad(Dsplus_mass)));
	      pa3 = sqrt(max2(1e-10,kwad(x3*W1/2)-kwad(diquark_mass)));
	      
	      cthe2 = (pa3*pa3-pa1*pa1-pa2*pa2)/2/pa1/pa2;
	      cthe3 = (pa2*pa2-pa1*pa1-pa3*pa3)/2/pa1/pa3; 
	    }
scenariusz = 6; 
   }	

      py3ent_(&ip, &frag_parton, &spectator1, &spectator2, &W1, &x1,&x2);
  
 }
// pylist_(&lista);
// cout<<"scenariusz="<<scenariusz<<endl; 

}


