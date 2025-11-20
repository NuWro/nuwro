#ifndef G_CONSTANTS_H
#define G_CONSTANTS_H

#include <cmath>
#include <cstdlib>
#include <array>

#include <pdg.h>
#include <jednostki.h>
#include <isotopes.h>

enum TargetNucleus 
{
        C12_SF,
        C12_newSF,
        O16_SF,
        Ar40_SF,
        Fe56_SF,
        Unsupported
};

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


// Atomic masses
inline double getIsotopeMass( int Z, int N ) {
    isotope *iso = isotope_find( Z, N );
    return iso->atomic_mass * Unit; // Convert micro-u to MeV/c^2
}

const double carbon12Mass = getIsotopeMass(6, 6);
const double carbon11Mass = getIsotopeMass(6, 5);
const double boron11Mass = getIsotopeMass(5, 6);
const double oxygen16Mass = getIsotopeMass(8, 8);
const double oxygen15Mass = getIsotopeMass(8, 7);
const double nitrogen15Mass = getIsotopeMass(7, 8);
const double argon40Mass = getIsotopeMass(18, 22);
const double argon39Mass = getIsotopeMass(18, 21);
const double scandium47Mass = getIsotopeMass(21, 26);
const double titanium48Mass = getIsotopeMass(22, 26);
const double manganese55Mass = getIsotopeMass(25, 30);
const double iron56Mass = getIsotopeMass(26, 30);
const double iron55Mass = getIsotopeMass(26, 29);

// Atomic numbers
const int carbonZ(6);
const int carbonN(6);
const int oxygenZ(8);
const int oxygenN(8);
const int argonZ(18);
const int argonN(22);
const int calciumZ(20);
const int calciumN(20);
const int ironZ(26);
const int ironN(30);
const int argonA(40);
const int titaniumA(48);

constexpr int CARBON = 6006;
constexpr int OXYGEN = 8008;
constexpr int ARGON = 18022;
constexpr int CALCIUM = 20020;
constexpr int IRON = 26030;

const double pi(4*atan(1.0));
// Global constants
const double GF( 1.16637e-5/GeV2 ); //Fermi constant
const double cosThetaC( 0.97418 );  //Cabibbo angle, PDG 2008, p. 145
const double cos2ThetaC( cosThetaC*cosThetaC );
const double sinWeinbergAngleSq( 0.23116 );
const double GFcosTheta(GF*cosThetaC); //Fermi constant including the Cabibbo angle, *10^11
const double alpha( 7.2973525376e+4 ); //fine structure constant, *10^7
const double reciprocalAlpha( 137.03599908 );
const double M( 0.50 * ( PDG::mass_proton + PDG::mass_neutron ) );
const double Mass2( M * M );
const double fm_to_mev = 1.0 / 197.327;
const double fm2 = fm_to_mev * fm_to_mev;
const double fm3 = fm_to_mev * fm_to_mev * fm_to_mev;
//~ const double piMass( 139.57018 ); //charged pion mass

const double M2(M*M);
const double piMass2( PDG::mass_piP * PDG::mass_piP );
const double pMass( 938.27203 * MeV );
const double nMass( 939.56536 * MeV );
const double eMass( 0.510998910*MeV );
const double mu_p(2.792847351);
const double mu_n(-1.91304273);
//const double xi(mu_p - mu_n - 1.0);
// Constants for ff.cc
const double gA( -1.2701 );
const double MV2( 710000.0 );
const double MA_cc_mec = 1014 * MeV;
const double MA_nc_mec = 1014 * MeV;
const double Dipole_Lambda = 5.6; //C Thorpe: Lambda parameter used in dipole form factors // default 1030 MeV // March 2019
//Constants associated with SU(3) representation of axial currents
const double Axial_F = 0.463;
const double Axial_D = 0.804;
static const double Axial_x = Axial_F/(Axial_F+Axial_D);
static double p_AEp[7] = {1., 0.9927, 0.9898, 0.9975, 0.9812, 0.9340, 1.};
static double p_AMp[7] = {1., 1.0011, 0.9992, 0.9974, 1.0010, 1.0003, 1.};
static double p_AEn[7] = {1., 1.1011, 1.1392, 1.0203, 1.1093, 1.5429, 0.9706};
static double p_AMn[7] = {1., 0.9958, 0.9877, 1.0193, 1.0350, 0.9164, 0.7300};
static double p_AAx[7] = {1., 0.9958, 0.9877, 1.0193, 1.0350, 0.9164, 0.7300};
static double axial_ff_beta = 0.0;
static double axial_ff_theta = 0.0;
static double axial_ff_gamma = 0.0;
static double axial_ff_alpha = 0.0;
static double MA_cc = 1030 * MeV;
static double MA_nc = 1030 * MeV;
static double MA_s = 1030 * MeV;
static double MA_hyp = 1030 * MeV; // C Thorpe Added Hyperon axial mass parameter
//static const double Dipole_Lambda = 5.6; //C Thorpe: Lambda parameter used in dipole form factors // default 1030 MeV // March 2019
static double delta_s = 0;
static int axialFFset = 0;
static const int strangeFFset = 0;
static int strange = 0;
//Second class current setup
//real and imaginary components at Q2=0
static double Rg20 = 0;
static double Ig20 = 0;
//symmetry breaking setup
static bool sym_break = false;
static int strangeEM = 0;

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
