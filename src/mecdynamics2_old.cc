#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "vec.h"
#include "params.h"
#include "nucleus.h"
#include "jednostki.h"

//double pi = M_PI;;
//double Pi2 = Pi * Pi;
double Pi3 = Pi2 * Pi;

//double przelcm = 197.326968e-13;	// MeV*cm
//double przelcm2 = przelcm * przelcm;
//double przelfermi = 197.326968;	//MeV fm
//double przelfermi2 = przelfermi * przelfermi;

double mecM = ( 938.272013 + 939.565346 )/2.0;
double mecM2 = mecM * mecM;

double mecm = 105.6583668;
double mecm2 = mecm*mecm;

double meckf = 220.0;
double mecEf = sqrt(meckf*meckf + mecM2);

double mecMv2=710000;
double mecMA=1016;
double mecMA2 = mecMA*mecMA;
double mecGA = -1.267;			 // negative!!!

double stala=5.07*1e-6;

static double mecmax2(double a, double b)
{
	if (a<b)
		return b;
	else
		return a;
}


static double mecGe(double Q2)
{
	if (opcja == 1)
		return 1/(1+Q2/mecMv2)/(1+Q2/mecMv2);

	if ( opcja > 1 )
	{
		double tau = Q2/4.0/939/939;
		return ( 1.0-0.0578*tau )/( 1.0 + (11.1+(13.6+33.0*tau)*tau)*tau ) -
			( 1.25+1.30*tau )*tau/( 1.0 +(-9.86+(305.0+(-758.0+802.0*tau)*tau)*tau)*tau );
	}
}


static double mecGm(double Q2)
{
	if (  opcja == 1  )
		return 4.706*mecGe(Q2,opcja);

	if ( (opcja == 2 ) || ( opcja == 4 ) )
	{
		double tau = Q2/4.0/939/939;
		return 2.792847351*( 1.0+0.15*tau )/( 1.0 + (11.1+(19.6+7.54*tau)*tau)*tau ) +
			1.91304273*( 1.0+1.81*tau )/( 1.0 + (14.1+(20.7+68.7*tau)*tau)*tau );
	}
	if (opcja == 3 )
	{
		double tau = Q2/4.0/939/939;
		return ( 2.792847351*( 1.0+0.15*tau )/( 1.0 + (11.1+(19.6+7.54*tau)*tau)*tau ) +
			1.91304273*( 1.0+1.81*tau )/( 1.0 + (14.1+(20.7+68.7*tau)*tau)*tau ) )
			*sqrt(1 + 6.0*Q2/1e6*exp(-Q2/0.35/1e6)) ;
	}
}


static double mecGa(double Q2, int opcja)
{
	if ( opcja < 4 )
		return mecGA/(1+Q2/mecMA2)/(1+Q2/mecMA2);

	if (opcja == 4)
	{
		double Maxial = 1350;
		double Maxial2 = Maxial*Maxial;
		return mecGA/(1+Q2/Maxial2)/(1+Q2/Maxial2);
	}
}


static double mecGp(double Q2, int opcja)
{return 2.0*939*939*mecGa(Q2,opcja)/(140*140 + Q2);}

static double mecF2(double Q2, int opcja)
{return (mecGm(Q2,opcja)-mecGe(Q2,opcja))/(1+Q2/4/mecM2);}

static double mecF1(double Q2, int opcja)
{return ( mecGe(Q2,opcja) + mecGm(Q2,opcja)*Q2/4/mecM2 )/(1+Q2/4/mecM2);}

///////////////////////////////////////////////////////////////////////////
								 //subtract energy by pot and adjusts momentum
static void spowalniacz (double pot, vect &pp)
{
	if ( pp.t-pot < sqrt(pp*pp) )
	{
		vec zero = vec(0,0,0);
		pp=vect (zero,sqrt(pp*pp));
	}

	else
	{
		double en = pp.t;
		vec kom = vec (pp);
		double scale = sqrt( kom.normSq() - pot*(2*en-pot) )/kom.length();
		pp=vect(kom*scale,en-pot);
	}
}


////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void model_2body2_old (double E, double w, double q, double Bin, particle &meclep, particle &nuc1, particle &nuc2)
{
	//it is assumed that neutrino direction is (0,0,1); perhaps should be relaxed...

	double muonmom = sqrt ((E-w)*(E-w) - mecm2);
	double cosmuon = (2.0*E*(E-w) - mecm2 - q*q + w*w)/2.0/E/muonmom;
	double phi = 2.0*Pi*frandom();

								 //momentum transfer
	vec qq(cos(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), sin(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon) , (E - muonmom*cosmuon));

								 //muon momentum
	vec kprim(-cos(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), -sin(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), muonmom*cosmuon);
	vect qqq(qq,w);				 //muon 4-momentum

	vec N1, N2;
	vect N11, N22, sumadd;
	vec maksio_vec;

	double length1, length2;
	double length;
	double Winv;

	vec trans;

	do
	{
		N1=rand_from_ball (meckf);
		N2=rand_from_ball (meckf);
		//cout<<E<<"  "<<w<<"  "<<q<<"  "<<N1<<"  "<<N2<<endl;
		length1= sqrt(mecM2 + N1.normSq());
		length2= sqrt(mecM2 + N2.normSq());

		N11=vect(N1, length1-mecEf+mecM);
								 //Fermi energy is subtracted
		N22=vect(N2, length2-mecEf+mecM);

		vect sum=N11+N22;

		sumadd = sum+qqq;
	}

								 //to be able to make Lorentz boost and to "decay"
	while ( ( sumadd.t < sumadd.length() )   ||   ( sumadd*sumadd <4.0*mecM2 ) );
	//cout<<"po"<<endl;
	trans= sumadd.speed();
	sumadd.boost2(trans);		 //boost to the CM frame

	Winv= sumadd.t;

	vec F1=rand_dir ()*sqrt(Winv*Winv-4.0*mecM2)/2.0;
	vec F2=-F1;

	vect nuc11=(Winv/2.0,F1);	 //hadron four momenta
	nuc11.t=Winv/2.0;
	vect nuc22=(Winv/2.0,F2);
	nuc22.t=Winv/2.0;

	nuc11.boost2(-trans);
	nuc22.boost2(-trans);		 //nucleons in the LAB frame

	spowalniacz(Bin,nuc11);
	spowalniacz(Bin,nuc22);

	vec nucmom1 = vec (nuc11);
	vec nucmom2 = vec (nuc22);
	//cout<<"ped_x="<<cohleptonmom_x<<"   wynik="<<endl;
	meclep.set_momentum(kprim);
	nuc1.set_momentum(nucmom1);
	nuc2.set_momentum(nucmom2);
	//cout<<"dormo"<<endl;

}


///////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//      mecweight returns Xsec per nucleon for TEM model
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double mecweight2_old (double E, bool cc, int mecA, particle &meclepton, particle &mecnucleon1, particle &mecnucleon2, particle &mecnucleon3)
{

	double w = (E-mecm)*frandom();
	double q = qmin + (qmax-qmin)*frandom();
	vect nuc1, nuc2, nuc3;
	// here is the isospin model; I assume that 3/5 times a pair is p-p and 2/5 times it is p-n
	mecnucleon1.pdg = 2212;		 //there must be at least one proton

	double losso = frandom();
	if (losso<0.6)
	{
		mecnucleon2.pdg = 2212;
	}
	else
	{
		mecnucleon2.pdg = 2112;
	}
	losso = frandom();
	//czarek  if (losso<??)
	if (losso<0.5)				 //czarek
	{
		mecnucleon3.pdg = 2212;
	}
	else
	{
		mecnucleon3.pdg = 2112;
	}

	mecnucleon1.set_mass (PDG::mass (mecnucleon1.pdg));
	mecnucleon2.set_mass (PDG::mass (mecnucleon2.pdg));
	mecnucleon3.set_mass (PDG::mass (mecnucleon3.pdg));

	double weight2 = (qmax-qmin*(E-mecm)*marteau_dwa (E,??)/12.0/1e38;
		if (weight2<0)
		{weight2 = 0;}
		double weight3 = (qmax-qmin*(E-mecm)*marteau_trzy (E,??)/12.0/1e38;
		if (weight3<0)
		{weight3 = 0;}

		if ( weight2>(weight2+weight3)*frandom() )
		{model_2body2_old (E, w, q, 8, meclepton, mecnucleon1, mecnucleon2);}
			else
				{model_3body (E, w, q, 8, meclepton, mecnucleon1, mecnucleon2, mecnucleon3);}

				return weight2+weight3;

}
