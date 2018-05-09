#ifndef constants_h_
#define constants_h_
#include "../jednostki.h"
#include "DM.h"
#include <cmath>
//typical PDG staticants in HEP units, put all the particle masses, couplings etc. here
//Masses:
const double Mp=938.272013*MeV; //Proton mass in MeV
const double Mn=939.565346*MeV; //Neutron mass in MeV (optionally precise, different masses)
//const double M_N=0.938*GeV;
const double M_N=0.5*(Mp+Mn);//=938.9186795
const double M=M_N;
const double M2=M*M;
const double Mp2=Mp*Mp;
const double M4=M2*M2;
const double Mp4=Mp2*Mp2;
const double mpip=139.57*MeV;  //pi+- mass
const double mpi0=134.98*MeV;  //pi0 mass
//const double mpi=0.138*GeV;
const double mpi=(2*mpip+mpi0)/3;//=138.04
const double mpi2=mpi*mpi;
const double mpi4=mpi2*mpi2;
const double mrho= 775.8*MeV; //rho meson mass
const double mrho2= mrho*mrho;
const double mrho4= mrho2*mrho2;
const double Md=1232.0*MeV; // Delta mass in Mev
const double Md2=Md*Md;
const double Md4=Md2*Md2;
const double MdI=1.0/Md;
const double Md2I=1.0/Md2;
const double qc=780*MeV;
const double qc2=qc*qc;
const double mmu=105.65*MeV;
const double me=0.511*MeV;
const double mtau=1776.82*MeV;

const double W2pi=Mp+2*mpi;
//Delta stuff

const double fd33=2.16;

const double qcmdelpin0=0.5*sqrt(Md4+mpi2*mpi2+Mp2*Mp2-2*(Md2*mpi2+mpi2*Mp2+Md2*Mp2))/Md;
const double rhodelpin0=qcmdelpin0*qcmdelpin0*qcmdelpin0/(0.57*0.57*GeV2+qcmdelpin0*qcmdelpin0)/Md;

const double MA=1.015*GeV;//nucleon axial mass Juan
//const double MA=0.999*GeV;//nucleon axial mass Olga
const double MA2=MA*MA;
const double Mv=0.84*GeV;
const double lambda2_n=0.71*GeV;
const double Mn24=4*M2;
const double Mv2=Mv*Mv;
const double Mvd=0.84*GeV;
const double Mvd2=Mvd*Mvd;
const double lambda_pi=1.25*GeV;
const double lambdapi2=lambda_pi*lambda_pi;
const double lambda_rho=2.5*GeV;
const double gprime=0.65;
const double rho_0=0.17;

//couplings
const double alpha= 1.0/137.0;//3599976;//do I have to explain? ;p

//const double fstar=2.1125;//pi N Delta coupling
const double fstar=fd33;
const double fpi=93*MeV;

const double F2pnd=0.25*fstar*fstar/Pi; // fpnd^2/4pi
const double fpinn=sqrt(0.32*Pi);
const double fumpi=fpinn/mpi;
const double Crho=2.2;

const double gA=1.267;//axial nucleon coupling

//miscelanous stuff

const double cabbibo=0.974;
const double mup= 2.792847351;
const double mun= -1.91304273;
const double lambdan=5.6;

const double sqrt3=1.7320508075688772935274463415058723669428052538103806280558;
const double sqrt2=1.4142135623730950488016887242096980785696718753769480731766;

const double CGDP[5]={0,sqrt2/sqrt3,1/sqrt3,-1/sqrt3,sqrt2/sqrt3};
const double CGDPC[5]={0,sqrt2/sqrt3,-1/sqrt3,1/sqrt3,sqrt2/sqrt3};
const double CGNP[5]={0,1/sqrt2,1,1,-1/sqrt2};
const double CGNPC[5]={0,1/sqrt2,1,1,-1/sqrt2};
const double CGOTH[5]={0,0,1,-1,0};

//const double minus=-1;
const double zero=0;
const std::complex<double> delusive(0,1);
const double wmin=(M+mpi);
const double w2min=wmin*wmin;

const double gpinn=sqrt(14.8*4*Pi);
const double gpinna=Mp*gA/fpi;


//MAID 2007 stuff

const double kr=(Md*Md-M*M)/(2*Md);
const double mpw2=(M+Md)*(M+Md);
const double mmw2=(M-Md)*(M-Md);
const double kgcm0 = (Md*Md-M*M)/(2*Md);
const double normgev=sqrt(GeV);
const double helpref0=sqrt(Pi*alpha*(Md-M)*(Md-M)/(Md2-M2)/M)*normgev;
const double helpref0s=helpref0*sqrt((Md2-M2)*(Md2-M2)/6)/Md2;

// extra constants for MAID2007 output handling

const double W0=1070.0*MeV;
const double Wstep=2.5*MeV;
const double Q20=0.225*GeV2;
const double Q2step=0.05*GeV2; 

//experimental helicity amplitudes
/*
static double a12[48]={
0.000 ,-0.135 ,0.006,
0.127 ,-0.155991 ,0.006,
0.400 ,-0.126524 ,0.006,
0.530 ,-0.11414 ,0.006,
0.650 ,-0.1009 ,0.006,
0.750 ,-0.0932155 ,0.006,
0.900 ,-0.0821148 ,0.006,
1.150 ,-0.0667485 ,0.006,
1.450 ,-0.0526668 ,0.006,
2.800 ,-0.022848 ,0.006,
3.000 ,-0.00698719 ,0.006,
3.500 ,-0.0065634 ,0.006,
4.000 ,-0.00626806 ,0.006,
4.200 ,-0.00609828 ,0.006,
5.000 ,-0.00804907 ,0.006,
6.000 ,-0.00512883 ,0.006
};

static double a32[48]={
0.000 ,-0.250 ,0.008,
0.127 ,-0.296767 ,0.01,
0.400 ,-0.243465 ,0.01,
0.530 ,-0.215445 ,0.01,
0.650 ,-0.192211 ,0.01,
0.750 ,-0.171711 ,0.01,
0.900 ,-0.154633 ,0.01,
1.150 ,-0.125944 ,0.01,
1.450 ,-0.0999952 ,0.01,
2.800 ,-0.0420266 ,0.01,
3.000 ,-0.0365807 ,0.01,
3.500 ,-0.0291167 ,0.01,
4.000 ,-0.0230233 ,0.01,
4.200 ,-0.0216749 ,0.01,
5.000 ,-0.0149289 ,0.01,
6.000 ,-0.0109412 ,0.01
};

static double s12[48]={
0.000 ,0.0000000 ,0.01,
0.127 ,0.0261880 ,0.01,
0.400 ,0.0227009 ,0.0012,
0.530 ,0.0211282 ,0.0012,
0.650 ,0.0201709 ,0.0011,
0.750 ,0.0166838 ,0.0013,
0.900 ,0.0176410 ,0.0011,
1.150 ,0.0151795 ,0.0011,
1.450 ,0.0121709 ,0.0014,
2.800 ,0.0056752 ,0.0013,
3.000 ,0.0049231 ,0.01,
3.500 ,0.0040342 ,0.01,
4.000 ,0.0045812 ,0.0012,
4.200 ,0.0032821 ,0.01,
5.000 ,0.0024615 ,0.01,
6.000 ,0.0013675 ,0.01
};
*/
//Dirac algebra
//these will be constant elements of algebra:

//Dirac matrices and 1_{4x4}:
const DM gamma0=Gamma(0);
const DM gamma1=Gamma(1);
const DM gamma2=Gamma(2);
const DM gamma3=Gamma(3);
const DM gamma5=Gamma(5);

const D4V<std::complex<double > > zeros(0,0,0,0);
//default 4-vector of Dirac matrices
const D4V<DM> dirac(gamma0,gamma1,gamma2,gamma3);
//default 4-vector of identities
const D4V<DM> identity(1,1,1,1);
//default contravariant vector of identities
const D4V<DM> contridentity(1,-1,-1,-1);
//default tensor \gamma^\mu\gamma\nu
const D4T<DM> Dirac2(dirac,dirac);
//default metric tensor in 4x4 matrix space
const D4T<DM> gmunu(1,0,0,0,
                    0,-1,0,0,
                    0,0,-1,0,
                    0,0,0,-1);
//default metric tensor in 4x4 matrix space
const D4T<double> gmunud(1,0,0,0,
			0,-1,0,0,
			0,0,-1,0,
			0,0,0,-1);
// \sigma^{\mu\nu}=\frac{i}{2}\left\{\gamma^\mu,\gamma^\nu\right\}
const D4T<DM> sigma(0, 0.5*(gamma0*gamma1-gamma1*gamma0)*delusive, 0.5*(gamma0*gamma2-gamma2*gamma0)*delusive, 0.5*(gamma0*gamma3-gamma3*gamma0)*delusive, -0.5*(gamma0*gamma1-gamma1*gamma0)*delusive, 0, 0.5*(gamma1*gamma2-gamma2*gamma1)*delusive, 0.5*(gamma1*gamma3-gamma3*gamma1)*delusive, -0.5*(gamma0*gamma2-gamma2*gamma0)*delusive, -0.5*(gamma1*gamma2-gamma2*gamma1)*delusive, 0, 0.5*(gamma2*gamma3-gamma3*gamma2)*delusive, -0.5*(gamma0*gamma3-gamma3*gamma0)*delusive, -0.5*(gamma1*gamma3-gamma3*gamma1)*delusive, -0.5*(gamma2*gamma3-gamma3*gamma2)*delusive,0);


//Oset's results from NPA 468. And LDA parameters from De Vries/De Jager Do not touch!
// or we do not guarantee the results will have any physical meaning :p
//well, let's say we have nonzero density dependence in 0 (alpha,beta)
static double dataset[6][6]={{0, 0, 0, 0, 1.0, 0.31},
			   {100, 0   , 12.9,  0  , 1.0 , 0.31},
			   {200, 5.5 , 19.0,  3.7, 0.93, 0.66},
			   {300, 11.7, 16.6, 16.5, 0.47, 0.79},
			   {400, 14.5, 15.1, 21.2, 0.40, 0.85},
			   {500,  5.4, 12.0, 12.5, 0.47, 0.89}};
//pion self-energy parametrization
static double aa[5]={-5.190,  1.060, -13.46,  0.382, -0.038};
static double be[5]={ 15.35, -6.640,  46.17, -1.322,  0.204};
static double ce[5]={ 2.060,  22.66, -20.34,  1.466,  0.613};

#endif
