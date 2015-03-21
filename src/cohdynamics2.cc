#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "vec.h"
#include "params.h"
#include "nucleus.h"
//#include "sampler.h"

//      Parametrization of pion-nucleon Xsections 
//      taken from the original Rein-Sehgal paper
static const int SN=70;
/////////////////////////////////////////////////////////////////////
double cohincl (double en) 
{
  if (en < 0.2)
    return en * 25 / 0.2;
  if (en < 0.3)
    return en * (56.5 - 25.0) / 0.1 + 25.0 - (56.5 - 25.0) / 0.1 * 0.2;
  if (en < 0.4)
    return en * (17.0 - 56.5) / 0.1 + 56.5 - (17.0 - 56.5) / 0.1 * 0.3;
  if (en < 0.6)
    return en * (8.0 - 17.0) / 0.2 + 17.0 - (8.0 - 17.0) / 0.2 * 0.4;
  if (en < 1.0)
    return en * (25.5 - 8.0) / 0.4 + 8.0 - (25.5 - 8.0) / 0.4 * 0.6;
  if (en < 1.3)
    return en * (20.5 - 25.5) / 0.3 + 25.5 - (20.5 - 25.5) / 0.3 * 1.0;
  if (en < 1.5)
    return en * (26.0 - 20.5) / 0.2 + 20.5 - (26.0 - 20.5) / 0.2 * 1.3;
  if (en < 1.8)
    return en * (23.0 - 26.0) / 0.3 + 26.0 - (23.0 - 26.0) / 0.3 * 1.5;
  if (en < 2.0)
    return en * (24.0 - 23.0) / 0.2 + 23.0 - (24.0 - 23.0) / 0.2 * 1.8;
  if (en < 5.0)
    return en * (22.0 - 24.0) / 3.0 + 24.0 - (22.0 - 24.0) / 3.0 * 2.0;
  
  else
    return 22;
}


/////////////////////////////////////////////////////////////////////
double cohtot (double en) 
{
  if (en < 0.25)
    return en * 109 / 0.25;
  if (en < 0.6)
    return en * (26.0 - 109.0) / 0.35 + 109 - (26.0 - 109) / 0.35 * 0.25;
  if (en < 0.7)
    return en * (29.0 - 26.0) / 0.1 + 26.0 - (29.0 - 26.0) / 0.1 * 0.6;
  if (en < 0.8)
    return en * (27.0 - 29.0) / 0.1 + 29.0 - (27.0 - 29.0) / 0.1 * 0.7;
  if (en < 1.0)
    return en * (38.0 - 27.0) / 0.2 + 27.0 - (38.0 - 27.0) / 0.2 * 0.8;
  if (en < 1.2)
    return en * (31.5 - 38.0) / 0.2 + 38.0 - (31.5 - 38.0) / 0.2 * 1.0;
  if (en < 1.5)
    return en * (36.0 - 31.5) / 0.3 + 31.5 - (36.0 - 31.5) / 0.3 * 1.2;
  if (en < 2.0)
    return en * (30.5 - 36.0) / 0.5 + 36.0 - (30.5 - 36.0) / 0.5 * 1.5;
  if (en < 5.0)
    return en * (26.0 - 30.5) / 3.0 + 30.5 - (26.0 - 30.5) / 3.0 * 2.0;
  
  else
    return 26;
}


double pi = M_PI;
double pi2 = pi * pi;
double pi3 = pi2 * pi;

//      Conversion constants
//      All the quantities are expressed in MeV
//      Few exceptions are mentioned explicitly
double przelcm = 197.326968e-13;	// MeV*cm
double przelcm2 = przelcm * przelcm;
double przelfermi = 197.326968;	//MeV fm
double przelfermi2 = przelfermi * przelfermi;

//      Parameters used in the Rein-Sehgal model
//      This part should be simplified within the entire Monte Carlo
double cohG = 1.16637e-11;	// MeV-2
double cohG2 = cohG * cohG;
double cohM = 938.0;
double cohM2 = cohM * cohM;
double cohmpi = 138.0;
double cohmpi2 = cohmpi * cohmpi;

//      This has to be done more general
//double cohm=105.0;
//double cohm2=cohm*cohm;
double cohMA = 1000;		// RS choice
double cohMA2 = cohMA * cohMA;
double cohfpi = 0.93 * cohmpi;
double cohfpi2 = cohfpi * cohfpi;
double cohr0 = 1.0;		// RS choice;  in fermi
double cohr02 = cohr0 * cohr0;
double cohweight;
double cohx, cohy, cohz, cohphi;
double cohQ2;

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//      cohwaga returns Xsec per nucleon
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
  
//      cohprocess = 1 for CC ; cohprocess = 0 for NC
//      coh_mass_correction = 1 for Rein-Sehgal m/neq 0 correction
double cohwaga2 (double E, bool cohprocess, int cohA, 
                 bool coh_mass_correction, double cohm,
                 vec& cohlepton3mom,vec& cohpion3mom) 
{ //cohm=0;
  double cohm2 = cohm * cohm;
  int A2 = cohA * cohA;
  if (E < 140)
    return 0;
  double cohb = cohr02 * pow (cohA, 2.0 / 3.0) / 3.0;	//in fermi^2
  double E2 = E * E;
  double E3 = E2 * E;
//  static sampler<SN> sY("losy",0);
//	double losy=sY.random();
	double losy=frandom();

  cohy = cohmpi / E + losy * (1 - cohmpi / E);
  double cohywaga = 1 - cohmpi / E;
  
    
//      Integration range in x variable
  double cohxmax = 2 * E * (1 - cohy) / cohM / cohy;

 
	double px=0.5;
//	static sampler<SN> sX("losx",0);
//	double losx=sX.random();
    double losx=frandom();
	double cohxwaga=cohxmax;
	double cohxmid;
    cohx=cohxmax*losx;
    
	if (cohxmax>1.0)
	{
		cohxmid=0.3;
		if (losx<px)
		{
			cohx= cohxmid*losx/px;
			cohxwaga=cohxmid/px;
		}
		else
		{
			cohx= ( cohxmid - cohxmax*px + losx*(cohxmax-cohxmid) )/(1-px);
			cohxwaga= (cohxmax-cohxmid) /(1-px);}
		}
	else
	{
		cohxmid=cohxmax/5.0;
		if (losx<px)
		{
			cohx= cohxmid*losx/px;
			cohxwaga=cohxmid/px;
		}
		else
		{
			cohx= ( cohxmid - cohxmax*px + losx*(cohxmax-cohxmid) )/(1-px);
			cohxwaga= (cohxmax-cohxmid) /(1-px);
		}
	}


//  double cohxwaga = cohxmax;
  double cohQ2 = 2.0 * cohM * cohx * E * cohy;
  double cohy2 = cohy * cohy;
  double cohx2 = cohx * cohx;
  double cohsigma_piN = cohtot ((cohy * E)/ 1000) * 1e-27;
  double cohFabs =
    exp (-9.0 * pow (cohA, 1.0 / 3) * cohincl ((cohy * E) / 1000) / 10.0 /
	 16.0 / pi / cohr02);
//      This parameter is a fit to the plot shown in the original RS paper
//      At very small values of neutrino energy (~ 300-500) xsec is still underestimated 
  double cohr2 = 0.2;
  
/*Czarek1
//static sampler<SN> sZ("losz",1);
//cohz=-1.0+2.0*sZ.random();
cohz=-1.0+2.0*frandom();
double cohzwaga=2.0;
Czarek1*/
  double cohzwaga;
  double zpar = 1.01;
//  static sampler<SN> sampz("losz",0);
//  double losz = sampz.random ();
  double losz = frandom ();
  cohz = zpar - (zpar + 1) * pow ((zpar - 1) / (zpar + 1), losz);
  cohzwaga =
    (1 + zpar) * log ((zpar + 1) / (zpar - 1)) * pow ((zpar - 1) / (zpar + 1),
						      losz);
/*Czarek1*/

  cohphi = 2.0 * pi * frandom ();
  double cohphiwaga = 2.0 * pi;
  double cohz2 = cohz * cohz;
  double cohwyrminus = sqrt (1.0 - cohmpi2 / E2 / cohy2);
  double cohwyrplus = sqrt (1.0 + 2.0 * cohM * cohx / E / cohy);
  double cohwyrkin =
    sqrt (2.0 * cohM * cohx * (1.0 - cohy) / E / cohy -
	  cohM2 * cohx2 / E2);
  double cohprop = cohMA2 / (cohMA2 + cohQ2);
  double cohprop2 = cohprop * cohprop;
  double cohczek1 = (-2.0 * cohb * E2 * cohy2) / przelfermi2 * 
    (1 + cohM * cohx / E / cohy - cohmpi2 / 2.0 / E2 / cohy2 -
     cohwyrminus * cohz * (1.0 + cohM * cohx / E) );
  double cohczek2 =
    cos (cohphi) * (-2.0 * cohb * E2 * cohy2) / przelfermi2 *
    cohwyrminus  *sqrt (1.0 - cohz2) * cohwyrkin;
  double cohprzekroj = 1.0 / przelcm2 * cohG2 * cohM * cohfpi2 * A2 * E3 * cohy2 * (1 - cohy) 
        / 16.0 / pi3 * cohsigma_piN * cohsigma_piN * (1 + cohr2) 
        * cohprop2 * cohwyrplus * cohwyrminus * cohFabs  /2.0 / pi	//miara w calkowaniu po phi
        * exp (cohczek1 + cohczek2);
  double cohQ2min = cohm2 * cohy / (1 - cohy);
  if (!cohprocess)
  {
      cohweight = cohprzekroj * cohxwaga * cohywaga * cohzwaga * cohphiwaga;
  }
  else
  {
      double cohCorr =
			(1.0 - 0.5 * cohQ2min / (cohQ2 + cohmpi2)) 
			* (1.0 -  0.5 * cohQ2min / (cohQ2 + cohmpi2))  
			+cohy *	0.25 * cohQ2min * (cohQ2 - cohQ2min) 
			 / (cohQ2 + cohmpi2) / (cohQ2 + cohmpi2);
      cohweight =
			2.0 * cohprzekroj  *cohxwaga * cohywaga * cohzwaga * cohphiwaga * cohCorr ;
  }// if(false)
   //    cout<<"Old "<<cohm<<' '<<E<<' '<<cohQ2<<' '<<cohQ2min<<endl;
  if ((coh_mass_correction) && 
	  ((cohy > 1.0 - cohm / E) || (cohQ2 < cohQ2min)))
     cohweight = 0.0;
  
//      From the values : x, y, cos(theta) momenta of final states are reconstructed
    
//      Energy and momentum transfer
  double cohomega = cohy * E;
  double cohq = sqrt (cohomega * cohomega + cohQ2);
  
//      Leton momentum; separately CC and NC; can be made shorter
  if (cohprocess)
  {
      if ((E - cohomega) * (E - cohomega) - cohm2 < 0)
		return 0;
     double cohscattcos = (-cohQ2 - cohm2 + 2.0 * E * (E - cohomega))  
		/2.0 / E / sqrt ((E -	cohomega) *  (E - cohomega) - cohm2);
      double cohleptonmom =
	sqrt ((E - cohomega) * (E - cohomega) - cohm2);
      double cohleptonphi = 2 * pi * frandom ();
      cohlepton3mom =cohleptonmom *
	vec (sqrt (1.0 - cohscattcos * cohscattcos) * sin (cohleptonphi),
	     sqrt (1.0 - cohscattcos * cohscattcos) * cos (cohleptonphi), 
	     cohscattcos);
    }
  
  else
    
    {
      double cohscattcos = (-cohQ2 + 2.0 * E * (E - cohomega)) 
							/2.0 / E / (E - cohomega);
      double cohleptonmom = E - cohomega;
      double cohleptonphi = 2.0 * pi * frandom ();
      cohlepton3mom =cohleptonmom * 
		vec (sqrt (1.0 - cohscattcos * cohscattcos) * sin (cohleptonphi),
			 sqrt (1.0 - cohscattcos * cohscattcos) * cos (cohleptonphi), 
			 cohscattcos);
    } 
//      Pion momentum
  double cohpionmom = sqrt (E2 * cohy2 - cohmpi2);
  double cohpionphi = 2 * pi * frandom ();
  cohpion3mom =cohpionmom * 
		vec (sqrt (1.0 - cohz2) * sin (cohpionphi),
			 sqrt (1.0 - cohz2) * cos (cohpionphi),
			 cohz);
  
//if (cohweight/cohA > 1e-36)
//cout<<"waga="<<cohweight/cohA<<"  E="<<E<<"  y="<<cohy<<"  x="<<cohx<<"  z="<<cohz<<"  fi="<<cohphi<<endl;
    double res=cohweight / cohA;// *2.0/3.0;
//    sY.report(res,res/cohA);  
//    sY.report(res);  
//    sX.report(res);
//    sZ.report(res);
    return res;
    return cohweight / cohA *2.0/3.0;//the last multiplicative factor is taken from Nuance MC
}

