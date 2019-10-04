#ifndef CONSTANTS_H
#define CONSTANTS_H


/*
The "Constants.h" header file, when included, provides some useful constant objects that will be used throughout the code.
*/


#include <complex>
#include "Matrix.h"
using namespace std;

// The complex number unit
const complex<double> I = complex<double> (0, 1);

// Some real constants
// const double Mp = 938.27203;  // proton mass (MeV)
// const double Mn = 939.56536;  // neutron mass (MeV)
const double MN = 938.918695;  // average
const double Mp = 938.918695;
const double MN2 = 881568.315821;//MN^2

const double Mpi_chrgd = 139.57018;  // charged-pion mass (MeV)
const double Mpi_ntrl = 134.9766;  // neutral pion mass (MeV)
//const double Mpi = Mpi_chrgd;
const double Mpi = 138.0389867; // average (2*Mpi_chrgd+Mpi_ntrl)/3
const double Mpi2 = 19054.761849; //Mpi^2

//const double muonmass = 105.658369;  // muon mass (MeV)
const double muonmass = 0.00;

const double Pi = 3.141592654;  // \pi
const double hbc = 197.3270;  // hbar*c (MeV*fm)

// // // // // // // // // // // // // // // 
const double x2Pi = 2*Pi;
const double rad_conv = Pi/180;
const double deg_conv = 180/Pi;
// // // // // // // // // // // // // // // 

// Electroweak constants
const double M_W = 80385.;  // W-boson mass (MeV)
const double M_Z = 91187.6; // Z-boson mass (MeV)
const double Cabibbo = 0.974;  // cosine Cabibbo angle
const double G_Fermi = 1.16637e-11;  // Fermi constant (MeV^{-2})
const double Alpha = 0.007297353;  // fine structure constant
const double Sin2W = 0.23122; 
const double QWeak = 0.07512; //weak charge of the proton = 1-4*Sin2W

// Delta-related constants
const double Wdth = 120.;  	// Delta-particle width (MeV)
const double MDelta = 1232.;  	// Delta-particle mass (MeV)
const double fDelta = 2.18; 	// 

// D13-related constants 
const double WdthD13 = 115.; 
const double MD13 = 1515.; 
const double fD13 = 1.62; // br = 0.60

// // S11-related constants 
const double WdthS11 = 150.; 
const double MS11 = 1535.; 
const double fS11 = 0.16; // br = 0.45

// // P11-related constants 
const double WdthP11 = 350.;
const double MP11 = 1430.; 
const double fP11 = 0.489; // br = 0.60 

// // Non-resonant constants 
const double gA = 1.26;
const double fpi = 93.;

// // // // // // // // // // 
const double pm_cutoff = 500.; // cutoff in the missing momentum: when pm>pm_cutoff the cross section is set to zero
const double k_N_PWcutoff = 600.; //if the momentum of the outgoing nucleon is larger than this quantity the projection coefficient will be zero and we recover a pure plane wave (only when DW=0)



///////////////////////////////////////////////////////////
// //    PARAMETERS RELATED WITH THE OPERATOR         // //
///////////////////////////////////////////////////////////
const int iModel = 0; // iModel=0 -> LEM without form factors for the resonances, iModel=1 -> LEM with resonance form factors, iModel=2 -> Hybrid model

//Parameters used for the operator which are not kinematic;
const int Pascalutsa = 0;// 1 to use Pascalutsa Delta decay vertex, '0' for the normal one (Delta pole)
//int medmod = 0; //Medmod is in reaction parameters

const double W_upper = 2500; //Above 2100 MeV, the LEM (chpt + resonances) is not computed, only the regge operator is used.
const double W_lower = 1200; // Below 1200 MeV the Regge is not computed

// // // Parameters of the transition function
const double W0 = 1.7E3; // MeV, center of the transition
const double leffe = 0.1E3; // MeV, width of the transition (divided by 2)
// // // // // // // // // // 
// // // // // // // // // // // // // // // // // // // // // // // // // 


// // // // // // // // // // // // // // // 
const int nr=310;		// big enough integer, actually bigger or equal than the number of lines in input files
const int np=510;		// big enough integer, actually bigger or equal than the number of lines in input files
const int degeneracy_max = 4;	// big enough integer, actually bigger or equal than (2J+1) for the maximum J involved
const int shells = 8;		// big enough integer, actually bigger or equal than number of bound shells (occupied or not) + 1
// const int shells = 16;		// big enough integer, actually bigger or equal than number of bound shells (occupied or not) + 1
const double pFermi = 230.;     // Fermi momentum, it is only used whne PB=2 => hard cutoff 
// // // // // // // // // // // // // // // 


// // // // // // // // // // // // // 
const int npunt1 = 1100;    // maximum number of radial points for the imported DW spinor
const int ntheta = 1100;      // number of cos(theta_r) points for the imported DW spinor
const double r_fin = 10.;      // nuclear size for the dw spinor
// // // // // // // // // 

// // // Integration limits for the 3d integral of the hadronic current
const double r_i = 0.0001;  //fm
const double r_f = 8.;      //fm

const double phi_i = 0.000000001;  
const double phi_f = 6.283185307;

const double costheta_i =-0.99999999;
const double costheta_f = 0.99999999;
// // // // // // // // // // // // // // // // // // // // // // 


// // // Integration limits for the 3d integral over the angles in Code.cpp
const double phiPi_i = 0.000000001;  
const double phiPi_f = 6.283185307;

const double phiN_i = 0.000000001;  
const double phiN_f = 6.283185307;

const double costhetaN_i =-0.99999999;
const double costhetaN_f = 0.99999999;

const double costhetaPi_i =-0.99999999;
const double costhetaPi_f = 0.99999999;
// // // // // // // // // // // // // // // // // // // // // // 

const int down = 0; //spin -1/2
const int up = 1;   //spin 1/2

const int proton = 1;
const int neutron = 2;

struct BWF_parameters
{
    int A, level, nucleon;
    double Em, Ma, Mb;
    int kappa, two_J, nprincipal;
    double (*r_in)[nr], (*G)[nr], (*F)[nr];
    double (*p_in)[np], (*Gp)[np], (*Fp)[np];
    double *Emiss; 
    int *L2, *J2, *nprin, *nuc, *okup;
    int level_min, level_max;
};

struct Reaction_parameters
{
    int process, decay, Helicity, medmod, nucleon;
    int Pauli_blocking;
    int irun;
    double Minf_A, lambda_rho;
    double *ABCDE, *SL_factors, *e_responses;
};

struct Lepton_kin   //4vectors and lepton tensor in lab system
{
    double k_l_inc, k_l, w, q, Qsq, leptonmass, El;
    double costhetal, sinthetal;
    
    double *kl, *kl_inc, *Q, *Pa;

    double FF_CC;
    
    double (*Lepton_S)[4], (*Lepton_A)[4];
    complex<double> (*Lepton)[4];    
};

struct Hadron_kin   //4vectors and angles in the lab system
{
    double sinthetaPi, costhetaPi;
    double sinphiPi,   cosphiPi;
    double sinthetaN,  costhetaN;
    double sinphiN,    cosphiN;
    
    double k_pi, k_N, pm, TB;
    double *kpi, *kN, *ki, *pB;         
    
};

struct CS_tabulated
{
	double W_fin, QQ_fin;
	double W_step, QQ_step;
	double *W_vec, *QQ_vec;

	double (*table)[500][5];
};

struct Hadron_prime
{
    double *Q, *minusQ;
    double *kpi, *kN, *ki;         //4vectors in the kN//z system
    
    Matrix kpiSlash, QSlash, kiSlash, kNSlash; 
    
    double *sMan, *tMan, *uMan;
    Matrix sSlash, tSlash, uSlash;
    
    double W, s, t, u, Qsq;                     //invariants

    double k_pi;  //Scalars
};


#endif
