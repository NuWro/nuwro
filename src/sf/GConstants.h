#ifndef G_CONSTANTS_H
#define G_CONSTANTS_H

#include <cmath>
#include <cstdlib>
#include <pdg.h>

enum TargetElement
{	
	targC_Ben,
	targO_Ben,
	targO_GSF,
	targCa_GSF,
	targAr_GSF,
	targAr_Ben,
	targFe_Ben,
	targEnd
};

enum MomDistribs
{	
	md_C12_Ben,
	md_O16_Ben,
	md_O16_GCo,
	md_O16_CdA,	
	md_Ca40_CdA,
	md_Ca40_GCo,
	md_Ca48_GCo,
	md_Fe56_Ben,
	mdEnd
};

enum IsospinOfSF
{
	proton,
	neutron,
	isoEnd
	
};

enum FormFactors
{
	dipoleFF=1,
	BBBA05FF=2,
	BBA03FF=3,
	JLabFF=4
};


const double pi(4*atan(1.0));
//const double GFcosTheta(1.14959);//Fermi constant including the Cabibbo angle, *10^11
const double GF(1.16637);//Fermi constant *10^11
const double cosThetaC(0.97377);//Cabibbo angle, PDG 2006, p. 138
const double GFcosTheta(GF*cosThetaC);//Fermi constant including the Cabibbo angle, *10^11
const double alpha(7.297352568e+4);//fine structure constant, *10^7

const double M( 0.50*(PDG::mass_proton+PDG::mass_neutron) );

//~ const double piMass( 139.57018 );//charged pion mass

const double M2(M*M);
const double piMass2(PDG::mass_piP*PDG::mass_piP);

const double mu_p(2.792847351);
const double mu_n(-1.91304273);
//const double xi(mu_p - mu_n - 1.0);
const double gA(-1.2673);
//const double MA(1070.0);
//const double MA(1030.0);

const double MV2(710000.0);
//const double MA2(MA*MA);
////////////////////////////////////////////////////////////////////////

inline static double beta(double A){return (A-2.0)/(A-1.0);}

inline static double f34(double x, double y)
{
	return 0.75/(x*y*y);
}

struct TargetData
{const int A,N,Z;
 const double PBlock,EBlock,Beta;
 const double E1,E2;
 
 TargetData(int a, int z, double pblock, double e1=0, double e2=0):
 	A(a),N(a-z),Z(z),
	PBlock(pblock),EBlock(M2+pblock*pblock),
	Beta(beta(a)),
	E1(e1),E2(e2)
	{}
};


//                        A   Z  PBlock    E1        E2
const TargetData carbon ( 12,  6, 209);
const TargetData oxygen ( 16,  8, 209, 19.5855, 26.3294 );
const TargetData calcium( 40, 20, 217, 29.2565, 16.4631 );
const TargetData argon  ( 40, 18, 217, 29.2565, 16.4631 );
const TargetData iron   ( 56, 26, 225);

// Coefficients for the Gausian Approximation of the Spectral Function
const double OxygenN[] = 
{  // coef,mean,width
	2, 47.00,70.0,
	4, 21.80,4.0,
	2, 15.65,4.0,
	 0, 4.15,2.0,
	 0, 3.27,2.0,
	 0,-0.93,2.0,
	 0,
};


const double OxygenP[] = 
{  // coef,mean,width
	 2,45.00,70.0,
	 4,18.44,4.0,
	 2,12.11,4.0,
	 0,0.59,2.0,
	 0,0.08,2.0,
	0,-4.65,2.0,
	0
};



const double CalciumN[] = 
{ // coef,mean,width
	2,66.12,25.0,//"exact"
	4,43.80,15.0,//"exact"
	2,39.12,15.0,//"exact"
	6,22.48,6.0,//"exact"
	2,17.53,4.0,
	4,15.79,4.0,
	0
};


const double CalciumP[] = 
{ // coef,mean,width
	
	2,57.38,25.0,//"exact"
	4,36.52,15.0,//"exact"
	2,31.62,15.0,//"exact"
	6,14.95,4.0,//"exact"
	2,10.67,2.0,
	4,8.88,2.0,
	0
};


const double ArgonN[] = 
{ // coef,mean,width
	2,62,25,
	4,40,15,
	2,35,15,
	6,18,5,
	2,13.15,4,
	4,11.45,3,
	2,5.56,3,
	0
};


const double ArgonP[] = 
{ // coef,mean,width
	2,52,25,//"exact"
	4,32,15,//"exact"
	2,28,15,//"exact"
	6,11,4,//"exact"
	2,8,2,
	2,6,2,
	0,1,2,
	0
};


#endif
