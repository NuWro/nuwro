#include <math.h>
#include <iostream>
#include <fstream>
#include "util2.h"
#include "calga5.h"
//#include "myfermi.h"
//#include "model.h"
//#include "oset2.h"
//#include "formfactors.h"
#include "vecrand.h"
#include "vec.h"
#include "vect.h"
#include "Dpion_a3.h"
#include "di_a3.h"
#include "dr_a3.h"
#include "nn_a3.h"
#include "nnn_a3.h"
#include <cstdlib>
#include "params.h"
#include "nucleus.h"
//#include "jednostki.h"
#include "los.h"

using namespace std;

double f1=0.6;
double g1=0.7;
double g2=0.5;
double g3=0.5;
double fkwpiNN=0.08;
double Lampi = 1000;
double Lampi2 = Lampi*Lampi;
double Lamro = 1500;
double Lamro2 = Lamro*Lamro;
double fperf=4.78;
double C2ro=2;
//double Gamma = 118;
//double GD=118;
double fstkw=4.65;
double fkw=1;
double eps=1;
double mecmrho = 775.49;
double mecmrho2 = mecmrho*mecmrho;

const double Pi3 = Pi2 * Pi;

const double mecM = ( 938.272013 + 939.565346 )/2.0;
const double mecM2 = mecM * mecM;

const double mecMD = 1232;
const double mecMD2 = mecMD*mecMD;

const double mecMDeff = 1232+20;
const double mecMDeff2 = mecMDeff * mecMDeff;

const double mecmpjon = 139;
const double mecmpjon2 = mecmpjon*mecmpjon;

static double mecm = 105.6583668;
static double mecm2 = mecm*mecm;

const double meckf = 220.0;
const double meckf2 = meckf*meckf;
const double meckf3 = meckf2*meckf;
const double mecEf = sqrt(meckf*meckf + mecM2);

const double mecMv2=710000;
const double mecMA=1016;
const double mecMA2 = mecMA*mecMA;
const double mecGA = -1.267;	 // negative!!!

const double mecstala=5.07*1e-6;

const double marteau[20][3]=
{
	{ 0,  0,  0},
	{ 0,  0,  0},
	{ 0.0096412,  0.000831472,  0.000278494},
	{ 0.0128777,  0.00267156,  0.000727129},
	{ 0.00944625,  0.00368255,  0.000831445},
	{ 0.00596565,  0.00385111,  0.000725262},
	{ 0.00357525,  0.00354324,  0.000557264},
	{ 0.00209501,  0.00302457,  0.000398252},
	{ 0.00121268,  0.00243011,  0.000271336},
	{ 0.00069551,  0.00179672,  0.000178407},
	{ 0.000395186,  0.00109438,  0.000113891},
	{ 0.000222054,  0.000324378,  7.07688e-05},
	{ 0.000123027,  0,  4.28107e-05},
	{ 6.69435e-05,  0,  2.51723e-05},
	{ 3.55937e-05,  0,  1.43389e-05},
	{ 1.83711e-05,  0,  7.87188e-06},
	{ 9.12489e-06,  0,  4.1335e-06},
	{ 4.31038e-06,  0,  2.05345e-06},
	{ 1.90418e-06,  0,  9.49729e-07},
	{ 7.67159e-07,  0,  3.99022e-07}
};

template <class F>
double supercalka1 (F f, double x1, double x2, double precision, double renorm)
{
	double results[2] = { calg20(f,x1,x2,1), calg20(f,x1,x1+(x2-x1)/2.0,1) + calg20(f,x1+(x2-x1)/2.0,x2,1)};
	if (results[0]==0 && results [1]==0)
		return 0;
	else
	{

		if ( 2.0*fabs(results[0]-results[1])< precision*fabs(results[0]+ results[1])  )
			return results[1];
		else
		{
			return supercalka1(f,x1,x1+(x2-x1)/2.0   ,precision*renorm, renorm) + supercalka1(f,x1+(x2-x1)/2.0,x2,precision*renorm, renorm);
		}

	}
}


template <class F>
double supercalka2 (F f, double x1, double x2, double precision, double renorm)
{
	double results[2] = { calg20(f,x1,x2,1), calg20(f,x1,x1+(x2-x1)/2.0,1) + calg20(f,x1+(x2-x1)/2.0,x2,1)};

	if ( 2.0*fabs(results[0]-results[1])< precision  )
		return results[1];
	else
	{
		return supercalka2(f,x1,x1+(x2-x1)/2.0   ,precision*renorm, renorm) + supercalka2(f,x1+(x2-x1)/2.0,x2,precision*renorm, renorm);
	}
}


double mecA=12;
double mecVol=3*Pi2*mecA/2.0/meckf3;

static double mecGe(double Q2)
{return 1/(1+Q2/mecMv2)/(1+Q2/mecMv2);}

static double mecGm(double Q2)
{return 4.706*mecGe(Q2);}

static double mecGa(double Q2)
{return mecGA/(1+Q2/mecMA2)/(1+Q2/mecMA2);}

static double mecGp(double Q2)
{return 2.0*940*940*mecGa(Q2)/(140*140 + Q2);}

static double mecF2(double Q2)
{return (mecGm(Q2)-mecGe(Q2))/(1+Q2/4/mecM2);}

static double mecF1(double Q2)
{return ( mecGe(Q2) + mecGm(Q2)*Q2/4/mecM2 )/(1+Q2/4/mecM2);}

								 //a parameter used in the oset_new function; for simplicity I took kf=225 as fixed
double kf_eff=220*pow(0.75, 1.0/3.0);

static double dataset[6][6]=
{
	{
		0, 0, 0, 0, 1.0, 0.31
	},
	{100, 0   , 12.9,  0  , 1.0 , 0.31},
	{200, 5.5 , 19.0,  3.7, 0.93, 0.66},
	{300, 11.7, 16.6, 16.5, 0.47, 0.79},
	{400, 14.5, 15.1, 21.2, 0.40, 0.85},
	{500,  5.4, 12.0, 12.5, 0.47, 0.89}
};

static double piondataset[9][6]=
{
	{
		5,  0,  0,  0,  0.79,  0.72
	},
	{45, 4.85, 9.45, 1.875, 0.79, 0.72},
	{85,  9.7  , 18.9,  3.7  , 0.79 , 0.72},
	{125, 11.9 , 17.7,  8.6, 0.62, 0.77},
	{165, 12.0, 16.3, 15.8, 0.42, 0.80},
	{205, 13.0, 15.2, 18.0, 0.31, 0.83},
	{245, 14.3, 14.1, 20.2, 0.36, 0.85},
	{285, 11.73, 13.53, 17.06, 0.394, 0.867},
	{325,  9.16, 12.96, 13.91, 0.429, 0.884},

};

double omega_eff_Fermi(double omega, double q)
{
								 // /sqrt(M*M+0.6*kf*kf)
	double wynik=omega+ 0.5*(omega*omega-q*q)/mecM;
	;
	if (wynik<0) wynik=0;
	return wynik;
}


double enkin_eff_Fermi(double omega, double q)
{
	double wynik=omega -mecmpjon + 0.5*(omega*omega-q*q)/mecM - mecmpjon2/2.0/mecM;
	if (wynik<0) wynik=0;
	return wynik;
}


///////Parametrization of different Delta decay channels from NPA 468:
///////Uses real photon tables
double Cgammax(int i, double omega)
{
	int j=int(omega/100);
	double wynik=0;
	if  (j<5)
	{
		double x1=j*100.0;
		double y1=dataset[j][i+1];
		double y2=dataset[j+1][i+1];
		wynik=y1+(omega-x1)*(y2-y1)/100.0;
		if (wynik<0) wynik=0;
	}
	else if(j>4)
	{
		double x1=325.0;
		double y1=dataset[4][i+1];
		double y2=dataset[5][i+1];
		wynik=y1+(omega-x1)*(y2-y1)/100.0;
		if (wynik<0) wynik=0;
	}
	return wynik;
}


double Cpiongammax(int i, double pionkin)
{
	int j=int((pionkin-5.0)/40);
	double wynik=0;
	if
		(j<8)
	{
		double x1=j*40.0+5;
		double y1=piondataset[j][i+1];
		double y2=piondataset[j+1][i+1];
		wynik=y1+(pionkin-x1)*(y2-y1)/40.0;
		if (wynik<0) wynik=0;
	}
	else if(j>7)
	{
		double x1=325.0;
		double y1=piondataset[7][i+1];
		double y2=piondataset[8][i+1];
		wynik=y1+(pionkin-x1)*(y2-y1)/40.0;
		if (wynik<0) wynik=0;
	}
	return wynik;
}


double mecc=2.355;
double mecz=0.5224;
double mecw=-0.149;

double rel_rho (double r)
{
	return
		(1.0 + mecw*r*r/mecc/mecc)/(1.0 + exp((r-mecc)/mecz))*12.0/11.1967;
}


double ratio_mec (double x)
{
	double radius = 2.2 - 0.49*tan(4.47 - 2.75*x) - 1.67*x +5.29*x*x -3.67*x*x*x;
	return rel_rho (radius);
}


//all below evaluated at 0.75 of the rho_0 !!!!! In general they depend on the nuclear density !!!!!
double oset_gamma_all (double w, double q, double stos)
{
	double weff = omega_eff_Fermi(w, q);
	return
		Cgammax(0, weff)*pow(stos, Cgammax(3,weff)) +
		Cgammax(1, weff)*pow(stos, Cgammax(4,weff)) +
		Cgammax(2, weff)*pow(stos, 2.0*Cgammax(4,weff));
}


double oset_gamma_pion (double w, double q, double stos)
{
	double weff = omega_eff_Fermi(w, q);
	return Cgammax(0, weff)*pow(stos, Cgammax(3,weff));
}


double oset_gamma_NN (double w, double q, double stos)
{
	double weff = omega_eff_Fermi(w, q);
	return Cgammax(1, weff)*pow(stos, Cgammax(4,weff));
}


double oset_gamma_NNN (double w, double q, double stos)
{
	double weff = omega_eff_Fermi(w, q);
	return Cgammax(2, weff)*pow(stos, 2.0*Cgammax(4,weff));
}


double oset_pion_all (double w, double q, double stos)
{
	double enkin_eff = enkin_eff_Fermi(w, q);
	return
		Cpiongammax(0, enkin_eff)*pow(stos, Cpiongammax(3,enkin_eff)) +
		Cpiongammax(1, enkin_eff)*pow(stos, Cpiongammax(4,enkin_eff)) +
		Cpiongammax(2, enkin_eff)*pow(stos, 2.0*Cpiongammax(4,enkin_eff));
}


double oset_pion_pion (double w, double q, double stos)
{
	double enkin_eff = enkin_eff_Fermi(w, q);
	return Cpiongammax(0, enkin_eff)*pow(stos, Cpiongammax(3,enkin_eff));
}


double oset_pion_NN (double w, double q, double stos)
{
	double enkin_eff = enkin_eff_Fermi(w, q);
	return Cpiongammax(1, enkin_eff)*pow(stos, Cpiongammax(4,enkin_eff));
}


double oset_pion_NNN (double w, double q, double stos)
{
	double enkin_eff = enkin_eff_Fermi(w, q);
	return Cpiongammax(2, enkin_eff)*pow(stos, 2.0*Cpiongammax(4,enkin_eff));
}


double oset_new (double w, double q)
{
	double x=los();
	double stos = ratio_mec(x);
	double y= 2*oset_gamma_all(w,q,stos) - oset_pion_all(w,q,stos);
	if (y>0)
		return y;
	else
		return 0;
}


double oset_pion_new (double w, double q)
{
	double x=los();
	double stos = ratio_mec(x);
	double y= 2*oset_gamma_pion(w,q,stos) - oset_pion_pion(w,q,stos);
	if (y>0)
		return y;
	else
		return 0;
}


double oset_NN_new (double w, double q)
{
	double x=los();
	double stos = ratio_mec(x);
	double y= 2*oset_gamma_NN(w,q,stos) - oset_pion_NN(w,q,stos);
	if (y>0)
		return y;
	else
		return 0;
}


double oset_NNN_new (double w, double q)
{
	double x=los();
	double stos = ratio_mec(x);
	double y= 2*oset_gamma_NNN(w,q,stos) - oset_pion_NNN(w,q,stos);
	if (y>0)
		return y;
	else
		return 0;
}


void threebody (double W, double m1, double m2, double m3, vect &p1, vect &p2, vect &p3)
{
	if (W < m1+m2+m3)
	{
		cout<<"impossible kinematics in 3-body decay"<<endl;
	}

	double W_2 = W*W;
	double m3_2 = m3*m3;
	double m1_2 = m1*m1;
	double m2_2 = m2*m2;

	double m23_2_los, m23_2_min, m23_2_max, E3;

	do
	{
		double m12_2_los = (m1+m2)*(m1+m2) + los()* ( (W-m3)*(W-m3) - (m1+m2)*(m1+m2) );
		double m23_test = (m2+m3)*(m2+m3) + los()* ( (W-m1)*(W-m1) - (m2+m3)*(m2+m3) );
		m23_2_los = m23_test;

		E3 = (W_2 + m3_2 - m12_2_los)/2.0/W;

		double E2star = (m12_2_los - m1_2 + m2_2)/2.0/sqrt(m12_2_los);
		double E3star = (W_2 - m12_2_los - m3_2)/2.0/sqrt(m12_2_los);

		m23_2_min = (E2star+E3star)*(E2star+E3star)
			- ( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) )*
			( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) );

		m23_2_max = (E2star+E3star)*(E2star+E3star)
			- ( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) )*
			( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) );

	}
	while ( (m23_2_los < m23_2_min) || (m23_2_los > m23_2_max) );

	//double m23_2_los = m23_2_min + los()* ( m23_2_max - m23_2_min);

	double E1 = (W_2 + m1_2 - m23_2_los)/2.0/W;
	double E2 = W - E1 - E3;

	double mom1 = sqrt (E1*E1 - m1_2);
	double mom2 = sqrt (E2*E2 - m2_2);
	double mom3 = sqrt (E3*E3 - m3_2);

	vec axis = rand_dir ();		 // new -- Y axis
	vec proba = rand_dir ();

	vec kier1 = vecprod (axis, proba);
	kier1 = kier1/kier1.length();//unit vector -- new Z axis!

								 // unit vetor -- new X axis
	vec kier2 = vecprod (axis, kier1);

	double kosalpha = (mom1*mom1 + mom2*mom2 - mom3*mom3)/2.0/mom1/mom2;
	double sinalpha = sqrt(1-kosalpha*kosalpha);
	double kosbeta = ( mom1 - mom2*kosalpha )/mom3;
	double sinbeta = sqrt(1 - kosbeta*kosbeta);

	vec pp1 = kier1*mom1;
	vec pp2 = mom2*(-kosalpha*kier1 + sinalpha*kier2);
	vec pp3 = mom3*(-kosbeta*kier1 - sinbeta*kier2);

	p1 = vect(pp1,E1);
	p2 = vect(pp2,E2);
	p3 = vect(pp3,E3);

}


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
		double scale = sqrt( kom.norm2() - pot*(2*en-pot) )/kom.length();
		pp=vect(kom*scale,en-pot);
	}
}


static double mecmax2(double a, double b)
{
	if (a<b)
		return b;
	else
		return a;
}


static double mecmax3(double a, double b, double c)
{
	return mecmax2(mecmax2(a,b),c)
		;
}


static double mecmin2(double a, double b)
{
	if (a<b)
		return a;
	else
		return b;
}


static double NN_elem (double w, double q, bool nowy)
{
  if (!nowy)
  {
	if (w>475)
		return 0;
	else
	{
		int n=int(w/25.0);
		double rest = w/25.0 - n;
		return rest*marteau[n+1][0]  + (1-rest)*marteau[n][0];
	}
  }
  else
  {
    double x=(q*q-w*w)/2.0/939/w;
  return  7.904*(0.0048883*x -4.11136e-6*x*x)*exp(-1.75589*x );
  }
}


static double ND_elem (double w, double q, bool nowy)
{
  if (!nowy)
  {
	if (w>475)
		return 0;
	else
	{
		int n=int(w/25.0);
		double rest = w/25.0 - n;
		return rest*marteau[n+1][1]  + (1-rest)*marteau[n][1];
	}
  }
  else
  {
    double x=(q*q-w*w)/2.0/939/w;
  return 7.904*(0.000431802 +0.000948343*x)*exp(-2.40084*x - 0.00135351*0.00135351*x*x);
  }
}

/*
double NN_newJS (double w, double q)
{
  double x=(q*q-w*w)/2.0/939/w;
  return  7.904*(0.0048883*x -4.11136e-6*x*x)*exp(-1.75589*x );
}

double ND_newJS (double w , double q)
{
  double x=(q*q-w*w)/2.0/939/w;
  return 7.904*(0.000431802 +0.000948343*x)*exp(-2.40084*x - 0.00135351*0.00135351*x*x);
}

double ND_marco (double w, double q)
{
  double x=(q*q-w*w)/2.0/939/w;
  double y=0.5*0.001*( 6.283 + 18.64*x - 16.93*x*x + 7.804*x*x*x )*exp (-x/0.472501);
  if (y<0)
    return 0;
  else
    return y;
}

double NN_marco (double w, double q)
{
  double x=(q*q-w*w)/2.0/939/w;
  double y=0.5*0.001*( -9.78 + 161.2*x - 149.8*x*x + 96.56*x*x*x )*exp (-x/0.46586);
  if (y<0)
    return 0;
  else
    return y;
}
*/
static double DD_old (double w, double q)
{
	if (w>475)
		return 0;
	else
	{
		int n=int(w/25.0);
		double rest = w/25.0 - n;
		return rest*marteau[n+1][2]  + (1-rest)*marteau[n][2];
	}
}


double NI (double w, double q)	 //wysumowana po izospinie!!!
{
	if (q==0)
		return 0;
	else
	if (w*w-q*q>=0)
		return 0;
	else
	if (-w/2+q/2*sqrt(1-4*mecM2/(w*w-q*q))>=mecEf)
		return 0;
	else
		return -mecM2/Pi/q*(
			mecEf-
			mecmax3(mecM,mecEf-w,-w/2+q/2*sqrt(1-4*mecM2/(w*w-q*q) ) ) );
}


double u2 (double q)
{
	return sqrt(
		(mecM2-meckf2-q*meckf+ mecEf*sqrt((meckf+q)*(meckf+q)+mecM2) )/
		(-mecM2-meckf2-q*meckf+ mecEf*sqrt((meckf+q)*(meckf+q)+mecM2) )
		);
}


double u1 (double q)
{
	return sqrt(
		(mecM2-meckf2+q*meckf+ mecEf*sqrt((meckf-q)*(meckf-q)+mecM2))
		/(-mecM2-meckf2+q*meckf+ mecEf*sqrt((meckf-q)*(meckf-q)+mecM2))
		);
}


double LinRelRe1 (double w, double q)
{
	return
		-mecM2/2/Pi2/q*(sqrt((q+meckf)*(q+meckf)+mecM2)-sqrt((q-meckf)*(q-meckf)+mecM2) )

		-mecM2*(2*mecEf+w)/4/Pi2/q*
		log( fabs(
		(w+mecEf-sqrt((q+meckf)*(q+meckf)+mecM2))/
		(w+mecEf-sqrt((q-meckf)*(q-meckf)+mecM2)) ) )
		-mecM2*(2*mecEf-w)/4/Pi2/q*
		log( fabs(
		(w-mecEf+sqrt((q+meckf)*(q+meckf)+mecM2))/
		(w-mecEf+sqrt((q-meckf)*(q-meckf)+mecM2)) ) )
		-mecM2/4/Pi2*(
		+log( fabs( (u2(q)-1)/(u1(q)-1) ) )
		-log( fabs( (u2(q)+1)/(u1(q)+1) ) )
		);
}


double LinRelRe2 (double w, double q)
{
	double Q2=q*q-w*w;
	if (Q2==0)
		return -mecM2/4/Pi2*2*(u2(q)-u1(q));
	else
	if (1+4*mecM2/Q2==0)
		return 0;
	else
	if ((1+4*mecM2/Q2)<0)
		return -mecM2/4/Pi2*
				2*sqrt( fabs( 1+4*mecM2/Q2 ) )*
				( atan( u2(q)/sqrt(fabs(1+4*mecM2/Q2) ))-
				atan( u1(q)/sqrt(fabs(1+4*mecM2/Q2) ))
				);
	else
		return -mecM2/4/Pi2*sqrt(1+4*mecM2/Q2 )*
			(-log(
			fabs( (u2(q)-sqrt(1+4*mecM2/Q2) )/
			(u1(q)-sqrt(1+4*mecM2/Q2) ) ) )
			+log(
			fabs( (u2(q)+sqrt(1+4*mecM2/Q2))/
			(u1(q)+sqrt(1+4*mecM2/Q2)) ) )
			);
}


double NR (double w, double q)
{
	if (q==0)
		return 0;
	else
								 //wysumowana po izospinie
		return 2*(LinRelRe1(w,q)+LinRelRe2(w,q));
}


double sw (double kat, double k, double w, double q)
{return mecM2 +w*w-q*q+2*(w*sqrt(mecM2+k*k)-k*q*kat);}

double u (double kat, double k, double w, double q)
{
	return mecM2 +w*w-q*q-
		2*(sw(kat,k,w,q)+mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2+w*w-q*q)/4/sw(kat,k,w,q);
}


double qcmf (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<=(mecM+mecmpjon)*(mecM+mecmpjon))
		return 0;
	else
		return sqrt(
			(sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2
			)/2/sqrt(sw(kat,k,w,q));
}


double blok (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)-(mecM+mecmpjon)*(mecM+mecmpjon)<=1)
		return 1.0/2;
	else
	if (sqrt((w+sqrt(mecM2+k*k))*
		(w+sqrt(mecM2+k*k))-sw(kat,k,w,q))
		*qcmf(kat,k,w,q)
		+(w+sqrt(mecM2+k*k))*sqrt(mecM2+qcmf(kat,k,w,q)*qcmf(kat,k,w,q))-
		mecEf*sqrt(sw(kat,k,w,q))
		<=0)
		return 0;
	else
	if (
		sqrt((w+sqrt(mecM2+k*k))*(w+sqrt(mecM2+k*k))
		-sw(kat,k,w,q))
		*qcmf(kat,k,w,q)
		<=(w+sqrt(mecM2+k*k))*sqrt(mecM2+qcmf(kat,k,w,q)*qcmf(kat,k,w,q))-
		mecEf*sqrt(sw(kat,k,w,q))
		)
		return 1;
	else
		return (
			sqrt((w+sqrt(mecM2+k*k))*(w+sqrt(mecM2+k*k))-sw(kat,k,w,q))
			*qcmf(kat,k,w,q)
			+(w+sqrt(mecM2+k*k))*sqrt(mecM2+qcmf(kat,k,w,q)*qcmf(kat,k,w,q))
			-mecEf*sqrt(sw(kat,k,w,q))
			)/
			2/qcmf(kat,k,w,q)/sqrt((w+sqrt(mecM2+k*k))*(w+sqrt(mecM2+k*k))-sw(kat,k,w,q));
}


double GD_new (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<0)
		return 0;
	else
		return 2*oset_new(w,q)+blok(kat,k,w,q)*
			qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*
			fstkw*mecM/6/Pi/mecmpjon2/sqrt(sw(kat,k,w,q));
}


double GD_new_noPB (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<0)
		return 0;
	else
		return 2*oset_new(w,q)+
			qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*
			fstkw*mecM/6/Pi/mecmpjon2/sqrt(sw(kat,k,w,q));
}


double GD_p_new (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<0)
		return 0;
	else
		return 2*oset_pion_new(w,q)+blok(kat,k,w,q)*
			qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*qcmf(kat,k,w,q)*
			fstkw*mecM/6/Pi/mecmpjon2/sqrt(sw(kat,k,w,q));
}


double GD_NN_new (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<0)
		return 0;
	else
		return 2*oset_NN_new(w,q);
}


double GD_NNN_new (double kat, double k, double w, double q)
{
	if (sw(kat,k,w,q)<0)
		return 0;
	else
		return 2*oset_NNN_new(w,q);
}


double dp_new (double kat, double k, double w, double q)
{
	if((sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2<0)
		return 0;
	else
	{
		double mecszer= GD_new_noPB(kat,k,w,q);
		return GD_new (kat,k,w,q)*k*k/
			((sw(kat,k,w,q)-mecMDeff2)*(sw(kat,k,w,q)-mecMDeff2)+mecMDeff2*
			mecszer*mecszer);
	}
}


double dp_p_new (double kat, double k, double w, double q)
{
	if((sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2<0)
		return 0;
	else
	{
		double mecszer= GD_new_noPB(kat,k,w,q);
		return GD_p_new(kat,k,w,q)*k*k/
			((sw(kat,k,w,q)-mecMDeff2)*(sw(kat,k,w,q)-mecMDeff2)+mecMDeff2*
			mecszer*mecszer);
	}
}


double dp_NN_new (double kat, double k, double w, double q)
{
	if((sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2<0)
		return 0;
	else
	{
		double mecszer= GD_new_noPB(kat,k,w,q);
		return GD_NN_new (kat,k,w,q)*k*k/
			((sw(kat,k,w,q)-mecMDeff2)*(sw(kat,k,w,q)-mecMDeff2)+mecMDeff2*
			mecszer*mecszer);
	}
}


double dp_NNN_new (double kat, double k, double w, double q)
{
	if((sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2<0)
		return 0;
	else
	{
		double mecszer= GD_new_noPB(kat,k,w,q);
		return GD_NNN_new (kat,k,w,q)*k*k/
			((sw(kat,k,w,q)-mecMDeff2)*(sw(kat,k,w,q)-mecMDeff2)+mecMDeff2*
			mecszer*mecszer);
	}
}


double dp_rz_new (double kat, double k, double w, double q)
{
	if((sw(kat,k,w,q)-mecM2-mecmpjon2)*(sw(kat,k,w,q)-mecM2-mecmpjon2)-4*mecM2*mecmpjon2<0)
		return 0;
	else
		return k*k*(
			(sw(kat,k,w,q)-mecMDeff2)/
			((sw(kat,k,w,q)-mecMDeff2)*(sw(kat,k,w,q)-mecMDeff2)+mecMDeff2*
			GD_new (kat,k,w,q)*GD_new (kat,k,w,q))-
			(u(kat,k,w,q)-mecMDeff2)/ ( (u(kat,k,w,q)-mecMDeff2)*(u(kat,k,w,q)-mecMDeff2)   +eps*eps)
			);
}


double ddp_new (double k, double w, double q)
{
	return calg20(fix234(dp_new,k,w,q), -1,1,1)
	//return supercalka1(fix234(dp_new,k,w,q), -1, 1, 0.025, 0.9)
		;
}


double ddp_p_new (double k, double w, double q)
{
	return calg20(fix234(dp_p_new,k,w,q), -1,1,1)
	//return supercalka1(fix234(dp_p_new,k,w,q), -1, 1, 0.025, 0.9)
		;
}


double ddp_NN_new (double k, double w, double q)
{
	return calg20(fix234(dp_NN_new,k,w,q), -1,1,1)
	//return supercalka1(fix234(dp_NN_new,k,w,q), -1, 1, 0.025, 0.9)
		;
}


double ddp_NNN_new (double k, double w, double q)
{
	return calg20(fix234(dp_NN_new,k,w,q), -1,1,1)
	//return supercalka1(fix234(dp_NNN_new,k,w,q), -1, 1, 0.025, 0.9)
		;
}


double ddp_rz_new (double k, double w, double q)
{								 //return calg20(fix234(dp_rz,k,w,q), -1,1,20)
	return supercalka2(fix234(dp_rz_new,k,w,q), -1, 1, 0.14, 1.25)
		;
}


double DI_full_new (double w, double q)
{
	return -32.0/9*mecMDeff2*calg20(fix23(ddp_new,w,q),0,meckf,1)/2/2/Pi2
	//return -32.0/9*mecMD2*supercalka1(fix23(ddp_new,w,q), 0, meckf, 0.025, 0.9)/2/2/Pi2
		;
}


double DI_pion_full_new (double w, double q)
{
	return -32.0/9*mecMDeff2*calg20(fix23(ddp_p_new,w,q),0,meckf,1)/2/2/Pi2
	//return -32.0/9*mecMD2*supercalka1(fix23(ddp_p_new,w,q), 0, meckf, 0.025, 0.9)/2/2/Pi2
		;
}


double DI_NN_full_new (double w, double q)
{
	return -32.0/9*mecMDeff2*calg20(fix23(ddp_NN_new,w,q),0,meckf,1)/2/2/Pi2
	//return -32.0/9*mecMD2*supercalka1(fix23(ddp_NN_new,w,q), 0, meckf, 0.025, 0.9)/2/2/Pi2
		;
}


double DI_NNN_full_new (double w, double q)
{
	return -32.0/9*mecMDeff2*calg20(fix23(ddp_NN_new,w,q),0,meckf,1)/2/2/Pi2
	//return -32.0/9*mecMD2*supercalka1(fix23(ddp_NNN_new,w,q), 0, meckf, 0.025, 0.9)/2/2/Pi2
		;
}


double DR_full_new (double w, double q)
{								 //return 32.0/9*mecMD/2/2/Pi2*calg20(fix23(ddp_rz,w,q),0,meckf,20)
	return 32.0/9*mecMDeff/2/2/Pi2*supercalka2(fix23(ddp_rz_new,w,q), 0, meckf, 0.14, 1.25)
		;
}


double DI (double w, double q)
{
	if ((w>3000)||(q>3600))
		return
			DI_full_new(w,q);
	else
	{
		int n=int(w/20);
		double x=w/20-n;
		int m=int(q/20);
		double y=q/20-m;
		return
			(1-x)*(1-y)*DI_a3[n][m]+
			(1-x)*y*DI_a3[n][m+1]+
			x*(1-y)*DI_a3[n+1][m]+
			x*y*DI_a3[n+1][m+1];
	}
	;
}


double DI_pion (double w, double q)
{
	if ((w>3000)||(q>3600))
		return
			DI_pion_full_new(w,q);
	else
	{
		int n=int(w/20.0);
		double x=w/20.0-n;
		int m=int(q/20);
		double y=q/20.0 -m;
		return
			(1-x)*(1-y)*Dpion_a3[n][m]+
			(1-x)*y*Dpion_a3[n][m+1]+
			x*(1-y)*Dpion_a3[n+1][m]+
			x*y*Dpion_a3[n+1][m+1];
	}
	;
}


double DI_NN (double w, double q)
{
	if ((w>3000)||(q>3600))
		return
			DI_NN_full_new(w,q);
	else
	{
		int n=int(w/20.0);
		double x=w/20.0-n;
		int m=int(q/20.0);
		double y=q/20.0-m;
		return
			(1-x)*(1-y)*NN_a3[n][m]+
			(1-x)*y*NN_a3[n][m+1]+
			x*(1-y)*NN_a3[n+1][m]+
			x*y*NN_a3[n+1][m+1];
	}
	;
}


double DI_NNN (double w, double q)
{
	if ((w>3000)||(q>3600))
		return DI_NNN_full_new(w,q);
	else
	{
		int n=int(w/20.0);
		double x=w/20.0-n;
		int m=int(q/20.0);
		double y=q/20.0-m;
		return
			(1-x)*(1-y)*NNN_a3[n][m]+
			(1-x)*y*NNN_a3[n][m+1]+
			x*(1-y)*NNN_a3[n+1][m]+
			x*y*NNN_a3[n+1][m+1];
	}
	;
}


double DR (double w, double q)
{
	if ((w>3000)||(q>3600))
		return
			DR_full_new(w,q);
	else
	{
		int n=int(w/20.0);
		double x=w/20.0-n;
		int m=int(q/20.0);
		double y=q/20.0 - m;
		return
			(1-x)*(1-y)*DR_a3[n][m]+
			(1-x)*y*DR_a3[n][m+1]+
			x*(1-y)*DR_a3[n+1][m]+
			x*y*DR_a3[n+1][m+1];
	}
	;
}


double V_NN_c(double w, double q, bool RPA)
{
	if (RPA == true)
		return f1/mecmpjon2/4.0;
	else return 0
			;
}


double V_NN_l(double w, double q, bool RPA)
{
	if (RPA == true)
		return (g1+q*q/(w*w-q*q-mecmpjon2)) *4*Pi*fkwpiNN/mecmpjon2*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)/4.0;
	else return 0
			;
}


double V_ND_l(double w, double q, bool RPA)
{
	if (RPA == true)
		return (g2+q*q/(w*w-q*q-mecmpjon2)) *4*Pi*fkwpiNN/mecmpjon2*sqrt(fperf)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)/4.0;
	else return 0
			;
}


double V_DD_l(double w, double q, bool RPA)
{
	if (RPA == true)
		return (g3+q*q/(w*w-q*q-mecmpjon2))*4*Pi*fkwpiNN/mecmpjon2*fperf*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)/4.0;
	else return 0
			;
}


double V_NN_t(double w, double q, bool RPA)
{
	if (RPA == true)
		return 4*Pi*fkwpiNN/mecmpjon2*( g1*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q) + C2ro*q*q/(w*w-q*q-mecmrho2)*
			(Lamro2-mecmrho2)/(Lamro2-w*w+q*q)*(Lamro2-mecmrho2)/(Lamro2-w*w+q*q) )/4.0;
	else return 0
			;
}


double V_ND_t(double w, double q, bool RPA)
{
	if (RPA == true)
		return 4*Pi*fkwpiNN/mecmpjon2*sqrt(fperf)*( g2*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q) +C2ro*q*q/(w*w-q*q-mecmrho2)*
			(Lamro2-mecmrho2)/(Lamro2-w*w+q*q)*(Lamro2-mecmrho2)/(Lamro2-w*w+q*q) )/4.0;
	else return 0
			;
}


double V_DD_t(double w, double q, bool RPA)
{
	if (RPA == true)
		return 4*Pi*fkwpiNN/mecmpjon2*fperf*( g3*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q)*(Lampi2-mecmpjon2)/(Lampi2-w*w+q*q) + C2ro*q*q/(w*w-q*q-mecmrho2)*
			(Lamro2-mecmrho2)/(Lamro2-w*w+q*q)*(Lamro2-mecmrho2)/(Lamro2-w*w+q*q) )/4.0;
	else return 0
			;
}


double den_c (double w, double q, bool RPA)
{
	return
		1-2*V_NN_c(w,q,RPA)*NR(w,q)+V_NN_c(w,q,RPA)*V_NN_c(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))
		;
}


double Pi_NN_qel_c(double w, double q, bool RPA)
{
	return
		NI(w,q)/den_c(w,q,RPA)
		;
}


double den_t (double w, double q, bool RPA)
{
	return
		pow( 1-V_NN_t(w,q,RPA)*NR(w,q)-V_DD_t(w,q,RPA)*DR(w,q)+
		(V_NN_t(w,q,RPA)*V_DD_t(w,q,RPA)-V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA))*
		(NR(w,q)*DR(w,q)-NI(w,q)*DI(w,q)),2.0)+
		pow( -V_NN_t(w,q,RPA)*NI(w,q)-V_DD_t(w,q,RPA)*DI(w,q)+
		(V_NN_t(w,q,RPA)*V_DD_t(w,q,RPA)-V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA))*
		(NR(w,q)*DI(w,q)+NI(w,q)*DR(w,q)),2.0);
}


double Pi_NN_qel_t(double w, double q, bool RPA)
{
	return
		NI(w,q)
		*(1+V_DD_t(w,q,RPA)*V_DD_t(w,q,RPA)*(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q))-2*V_DD_t(w,q,RPA)*DR(w,q))/den_t(w,q,RPA)
		;
}


double Pi_NN_D_t(double w, double q, bool RPA)
{
	return
		DI(w,q)* V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_t(w,q,RPA)
		;
}


double Pi_NN_pion_t(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*  V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_t(w,q,RPA)
		;
}


double Pi_NN_NN_t(double w, double q, bool RPA, bool nowy)
{
	return
		DI_NN(w,q)*  V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_t(w,q,RPA)
		+ NN_elem(w,q,nowy)*(-Pi/mecVol)*mecM/sqrt(mecM2+q*q)
		;
}


double Pi_NN_NNN_t(double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*  V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_t(w,q,RPA)
		;
}


double Pi_DD_qel_t(double w, double q, bool RPA)
{
	return
		NI(w,q)*  V_ND_t(w,q,RPA)*V_ND_t(w,q,RPA)  *(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q))/den_t(w,q,RPA)
		;
}


double Pi_DD_D_t(double w, double q, bool RPA)
{
	return
		DI(w,q)*(1+V_NN_t(w,q,RPA)*V_NN_t(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-2*V_NN_t(w,q,RPA)*NR(w,q))/den_t(w,q,RPA)
		;
}


double Pi_DD_pion_t(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*(1+V_NN_t(w,q,RPA)*V_NN_t(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-
		2*V_NN_t(w,q,RPA)*NR(w,q))/den_t(w,q,RPA)
		;
}


double Pi_DD_NN_t(double w, double q, bool RPA)
{
	return
		DI_NN(w,q)*(1+V_NN_t(w,q,RPA)*V_NN_t(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-
		2*V_NN_t(w,q,RPA)*NR(w,q))/den_t(w,q,RPA)
		+ DD_old(w,q)*(-Pi/mecVol) /fperf
		;
}


double Pi_DD_NNN_t(double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*(1+V_NN_t(w,q,RPA)*V_NN_t(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-
		2*V_NN_t(w,q,RPA)*NR(w,q))/den_t(w,q,RPA)
		;
}


double Pi_ND_qel_t(double w, double q, bool RPA)
{
	return
		NI(w,q)*V_ND_t(w,q,RPA)*(
		DR(w,q)-V_DD_t(w,q,RPA)*(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q)))/den_t(w,q,RPA)
		;
}


double Pi_ND_D_t(double w, double q, bool RPA)
{
	return
		DI(w,q)*V_ND_t(w,q,RPA)*(
		NR(w,q)-V_NN_t(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_t(w,q,RPA)
		;
}


double Pi_ND_pion_t(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*V_ND_t(w,q,RPA)*(
		NR(w,q)-V_NN_t(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_t(w,q,RPA)
		;
}


double Pi_ND_NN_t (double w, double q, bool RPA, bool nowy)
{
	return
		DI_NN(w,q)*V_ND_t(w,q,RPA)*(
		NR(w,q)-V_NN_t(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_t(w,q,RPA)
								 //*sqrt(mecM2+q*q)/mecM
		+ ND_elem(w,q,nowy)*(-Pi/mecVol) / sqrt( fperf
		)						 //* sqrt( 2.0*mecM * 2.0*mecMD * (mecM+w)/ (2.0*mecM+w) / (mecM+mecMD+w) / mecMD  )
		;
}


double Pi_ND_NNN_t (double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*V_ND_t(w,q,RPA)*(
		NR(w,q)-V_NN_t(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_t(w,q,RPA)
		;
}


double den_l (double w, double q, bool RPA)
{
	return
		pow( 1-V_NN_l(w,q,RPA)*NR(w,q)-V_DD_l(w,q,RPA)*DR(w,q)+
		(V_NN_l(w,q,RPA)*V_DD_l(w,q,RPA)-V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA))*
		(NR(w,q)*DR(w,q)-NI(w,q)*DI(w,q)),2.0)+
		pow( -V_NN_l(w,q,RPA)*NI(w,q)-V_DD_l(w,q,RPA)*DI(w,q)+
		(V_NN_l(w,q,RPA)*V_DD_l(w,q,RPA)-V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA))*
		(NR(w,q)*DI(w,q)+NI(w,q)*DR(w,q)),2.0);
}


double Pi_NN_qel_l(double w, double q, bool RPA)
{
	return
		NI(w,q)*(1+V_DD_l(w,q,RPA)*V_DD_l(w,q,RPA)*(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q))-2*V_DD_l(w,q,RPA)*DR(w,q))/den_l(w,q,RPA)
		;
}


double Pi_NN_D_l(double w, double q, bool RPA)
{
	return
		DI(w,q)*  V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_l(w,q,RPA)
		;
}


double Pi_NN_pion_l(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*  V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_l(w,q,RPA)
		;
}


double Pi_NN_NN_l(double w, double q, bool RPA, bool nowy)
{
	return
		DI_NN(w,q)*  V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_l(w,q,RPA)
		+ NN_elem(w,q,nowy)*(-Pi/mecVol)*mecM/sqrt(mecM2+q*q)
		;
}


double Pi_NN_NNN_l(double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*  V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA)  *(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))/den_l(w,q,RPA)
		;
}


double Pi_DD_qel_l(double w, double q, bool RPA)
{
	return
		NI(w,q)*  V_ND_l(w,q,RPA)*V_ND_l(w,q,RPA)  *(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q))/den_l(w,q,RPA)
		;
}


double Pi_DD_D_l(double w, double q, bool RPA)
{
	return
		DI(w,q)*(1+V_NN_l(w,q,RPA)*V_NN_l(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-2*V_NN_l(w,q,RPA)*NR(w,q))/den_l(w,q,RPA)
		;
}


double Pi_DD_pion_l(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*(1+V_NN_l(w,q,RPA)*V_NN_l(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))- 2*V_NN_l(w,q,RPA)*NR(w,q))/den_l(w,q,RPA)
		;
}


double Pi_DD_NN_l(double w, double q, bool RPA)
{
	return
		DI_NN(w,q)*(1+V_NN_l(w,q,RPA)*V_NN_l(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-2*V_NN_l(w,q,RPA)*NR(w,q))/den_l(w,q,RPA)
		+ DD_old(w,q)*(-Pi/mecVol) /fperf
		;
}


double Pi_DD_NNN_l(double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*(1+V_NN_l(w,q,RPA)*V_NN_l(w,q,RPA)*(NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q))-2*V_NN_l(w,q,RPA)*NR(w,q))/den_l(w,q,RPA)
		;
}


double Pi_ND_qel_l(double w, double q, bool RPA)
{
	return
		NI(w,q)*V_ND_l(w,q,RPA)*(DR(w,q)-V_DD_l(w,q,RPA)*(DR(w,q)*DR(w,q)+DI(w,q)*DI(w,q)))/den_l(w,q,RPA)
		;
}


double Pi_ND_D_l(double w, double q, bool RPA)
{
	return
		DI(w,q)*V_ND_l(w,q,RPA)*(NR(w,q)-V_NN_l(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_l(w,q,RPA)
		;
}


double Pi_ND_pion_l(double w, double q, bool RPA)
{
	return
		DI_pion(w,q)*V_ND_l(w,q,RPA)*(NR(w,q)-V_NN_l(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_l(w,q,RPA)
		;
}


double Pi_ND_NN_l(double w, double q, bool RPA, bool nowy)
{
	return
		DI_NN(w,q)*V_ND_l(w,q,RPA)*(NR(w,q)-V_NN_l(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_l(w,q,RPA)
								 //*sqrt(mecM2+q*q)/mecM
		+ ND_elem(w,q,nowy)*(-Pi/mecVol)/ sqrt( fperf
		)						 //* sqrt( 2.0*mecM * 2.0*mecMD * (mecM+w)/ (2.0*mecM+w) / (mecM+mecMD+w) / mecMD  )
		;
}


double Pi_ND_NNN_l(double w, double q, bool RPA)
{
	return
		DI_NNN(w,q)*V_ND_l(w,q,RPA)*(NR(w,q)-V_NN_l(w,q,RPA)*( NR(w,q)*NR(w,q)+NI(w,q)*NI(w,q) ))/den_l(w,q,RPA)
		;
}


double a00N(double w, double q)
{
	double Q2=q*q-w*w;
	return mecF1(Q2)-mecF2(Q2)*q*q/2/mecM/(mecM+sqrt(mecM2+q*q));
}


double b00N(double w, double q)
{
	double Q2=q*q-w*w;
	return q*(mecGa(Q2)*(1.0+0.5)/(mecM+sqrt(mecM2+q*q))-w/2/mecM*mecGp(Q2)/(mecM+sqrt(mecM2+q*q)));
}


double a03N(double w, double q)
{
	double Q2=q*q-w*w;
	return q*(mecF1(Q2)/(mecM+sqrt(mecM2+q*q))-w/2/mecM*mecF2(Q2)/(mecM+sqrt(mecM2+q*q)));
}


double b03N (double w, double q)
{
	double Q2=q*q-w*w;
	return mecGa(Q2) - q*q*mecGp(Q2)/2/mecM/(mecM+sqrt(mecM2+q*q));
}


double g0N(double w, double q)
{
	double Q2=q*q-w*w;
	return q*( mecF1(Q2)- mecF2(Q2)*w/2.0/mecM + mecF2(Q2)/2.0/mecM * (mecM+sqrt(mecM2+q*q)) )/( mecM+sqrt(mecM2+q*q) );
}


double d0(double w, double q)
{
	double Q2=q*q-w*w;
	return -mecGa(Q2);
}


double a00D(double w, double q)
{
	double Q2=q*q-w*w;
	return mecF1(Q2)-mecF2(Q2)*q*q/2/mecM/(mecMD+sqrt(mecMD2+q*q));
}


double b00D(double w, double q)
{
	double Q2=q*q-w*w;
	return q*(mecGa(Q2)*(1.0+0.5)/(mecMD+sqrt(mecMD2+q*q))-w/2/mecM*mecGp(Q2)/(mecMD+sqrt(mecMD2+q*q)));
}


double a03D(double w, double q)
{
	double Q2=q*q-w*w;
	return q*(mecF1(Q2)/(mecMD+sqrt(mecMD2+q*q))-w/2/mecM*mecF2(Q2)/(mecMD+sqrt(mecMD2+q*q)));
}


double b03D(double w, double q)
{
	double Q2=q*q-w*w;
	return mecGa(Q2)- q*q*mecGp(Q2)/2/mecM/(mecMD+sqrt(mecMD2+q*q));
}


double g0D(double w, double q)
{
	double Q2=q*q-w*w;
	return q*( mecF1(Q2) - mecF2(Q2)*w/2.0/mecM + mecF2(Q2)/2.0/mecM * (mecMD+sqrt(mecMD2+q*q)) )/( mecMD+sqrt(mecMD2+q*q) );
}


double d0D(double w, double q)
{
	double Q2=q*q-w*w;
	return -mecGa(Q2);
}


double SSS2(double E, double w, double q)
{
	return (mecm2+q*q-w*w)/4/E/(E-w);
}


double CSS2(double E, double w, double q)
{
	return 1-SSS2(E,w,q);
}


double L00(double E, double w, double q)
{
	return 16*E*(E-w)*CSS2(E,w,q);
}


double H00Nc(double w, double q)
{
	return pow(a00N(w,q),2.0);
}


double H00Nl(double w, double q)
{
	return pow(b00N(w,q),2.0);
}


double H00D(double w, double q)
{
	return b00D(w,q)*b00D(w,q);
}


double H00ND(double w, double q)
{
	return b00D(w,q)*b00N(w,q);
}


double L33(double E, double w, double q)
{
	return 16/q/q*(mecm2*(E*E-E*(E-w)*CSS2(E,w,q))+E*(E-w)*CSS2(E,w,q)*w*w);
}


double H33Nc(double w, double q)
{
	return pow(a03N(w,q),2.0);
}


double H33Nl(double w, double q)
{
	return pow(b03N(w,q),2.0);
}


double H33D(double w, double q)
{
	return b03D(w,q)*b03D(w,q);
}


double H33ND(double w, double q)
{
	return b03D(w,q)*b03N(w,q);
}


double L03(double E, double w, double q)
{
	return 8/q*(-mecm2*E-2*E*(E-w)*CSS2(E,w,q)*w);
}


double H03Nc(double w, double q)
{
	return a00N(w,q)*a03N(w,q);
}


double H03Nl(double w, double q)
{
	return b00N(w,q)*b03N(w,q);
}


double H03D(double w, double q)
{
	return b00D(w,q)*b03D(w,q);
}


double H03ND(double w, double q)
{
	return (b00D(w,q)*b03N(w,q)+b00N(w,q)*b03D(w,q))/2;
}


double L12(double E, double w, double q, bool nu)
{
  if (nu)
	return 16/q*(-mecm2*E+2*E*(E-w)*(2*E-w)*SSS2(E,w,q));
  else
    return -16/q*(-mecm2*E+2*E*(E-w)*(2*E-w)*SSS2(E,w,q));
}


double H12Nt(double w, double q)
{
	return 2*g0N(w,q)*d0(w,q);
}


double H12D(double w, double q)
{
	return 2*g0D(w,q)*d0(w,q);
}


double H12ND(double w, double q)
{
	return (g0N(w,q)+g0D(w,q))*d0(w,q);
}


double L11(double E, double w, double q)
{
	return 16*(2*E*(E-w)*SSS2(E,w,q)+(q*q-w*w)/q/q*E*(E-w)*CSS2(E,w,q)-
		mecm2/q/q*(E*E-E*(E-w)*CSS2(E,w,q)));
}


double H11Nt(double w, double q)
{
	return pow(g0N(w,q),2.0)+pow(d0(w,q),2.0);
}


double H11D(double w, double q)
{
	return pow(g0D(w,q),2.0)+pow(d0(w,q),2.0);
}


double H11ND(double w, double q)
{
	return g0N(w,q)*g0D(w,q)+pow(d0(w,q),2.0);
}


double trejsy_NN_qel_c (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nc(w,q)+
		L33(E,w,q)*H33Nc(w,q)+
		2*L03(E,w,q)*H03Nc(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_qel_c(w,q,RPA);
}


double trejsy_NN_qel_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nl(w,q)+
		L33(E,w,q)*H33Nl(w,q)+
		2*L03(E,w,q)*H03Nl(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_qel_l(w,q,RPA);
}


double trejsy_NN_D_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nl(w,q)+
		L33(E,w,q)*H33Nl(w,q)+
		2*L03(E,w,q)*H03Nl(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_D_l(w,q,RPA);
}


double trejsy_NN_pion_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nl(w,q)+
		L33(E,w,q)*H33Nl(w,q)+
		2*L03(E,w,q)*H03Nl(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_pion_l(w,q,RPA);
}


double trejsy_NN_NN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nl(w,q)+
		L33(E,w,q)*H33Nl(w,q)+
		2*L03(E,w,q)*H03Nl(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_NN_l(w,q,RPA,nowy);
}


double trejsy_NN_NNN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00Nl(w,q)+
		L33(E,w,q)*H33Nl(w,q)+
		2*L03(E,w,q)*H03Nl(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_NNN_l(w,q,RPA);
}


double trejsy_DD_qel_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00D(w,q)+
		L33(E,w,q)*H33D(w,q)+
		2*L03(E,w,q)*H03D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_qel_l(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_D_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00D(w,q)+
		L33(E,w,q)*H33D(w,q)+
		2*L03(E,w,q)*H03D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_D_l(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_pion_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00D(w,q)+
		L33(E,w,q)*H33D(w,q)+
		2*L03(E,w,q)*H03D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_pion_l(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_NN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00D(w,q)+
		L33(E,w,q)*H33D(w,q)+
		2*L03(E,w,q)*H03D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_NN_l(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_NNN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00D(w,q)+
		L33(E,w,q)*H33D(w,q)+
		2*L03(E,w,q)*H03D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_NNN_l(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_ND_qel_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00ND(w,q)+
		L33(E,w,q)*H33ND(w,q)+
		2*L03(E,w,q)*H03ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/sqrt(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_qel_l(w,q,RPA);
}


double trejsy_ND_D_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00ND(w,q)+
		L33(E,w,q)*H33ND(w,q)+
		2*L03(E,w,q)*H03ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_D_l(w,q,RPA);
}


double trejsy_ND_pion_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00ND(w,q)+
		L33(E,w,q)*H33ND(w,q)+
		2*L03(E,w,q)*H03ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_pion_l(w,q,RPA);
}


double trejsy_ND_NN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00ND(w,q)+
		L33(E,w,q)*H33ND(w,q)+
		2*L03(E,w,q)*H03ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_NN_l(w,q,RPA,nowy);
}


double trejsy_ND_NNN_l (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L00(E,w,q)*H00ND(w,q)+
		L33(E,w,q)*H33ND(w,q)+
		2*L03(E,w,q)*H03ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_NNN_l(w,q,RPA);
}


double trejsy_NN_qel_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q,nu)*H12Nt(w,q)+
		L11(E,w,q)*H11Nt(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_qel_t(w,q,RPA);
}


double trejsy_NN_D_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q,nu)*H12Nt(w,q)+
		L11(E,w,q)*H11Nt(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_D_t(w,q,RPA);
}


double trejsy_NN_pion_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12Nt(w,q)+
		L11(E,w,q)*H11Nt(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_pion_t(w,q,RPA);
}


double trejsy_NN_NN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12Nt(w,q)+
		L11(E,w,q)*H11Nt(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_NN_t(w,q,RPA,nowy);
}


double trejsy_NN_NNN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12Nt(w,q)+
		L11(E,w,q)*H11Nt(w,q))*(mecM+sqrt(mecM2+q*q))/2/mecM*
		Pi_NN_NNN_t(w,q,RPA);
}


double trejsy_DD_qel_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12D(w,q)+
		L11(E,w,q)*H11D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_qel_t(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_D_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12D(w,q)+
		L11(E,w,q)*H11D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_D_t(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_pion_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12D(w,q)+
		L11(E,w,q)*H11D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_pion_t(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_NN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12D(w,q)+
		L11(E,w,q)*H11D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_NN_t(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_DD_NNN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12D(w,q)+
		L11(E,w,q)*H11D(w,q))*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*
		Pi_DD_NNN_t(w,q,RPA)*fperf*mecMD/sqrt(mecMD2+q*q);
}


double trejsy_ND_qel_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12ND(w,q)+
		L11(E,w,q)*H11ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_qel_t(w,q,RPA);
}


double trejsy_ND_D_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12ND(w,q)+
		L11(E,w,q)*H11ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_D_t(w,q,RPA);
}


double trejsy_ND_pion_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12ND(w,q)+
		L11(E,w,q)*H11ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_pion_t(w,q,RPA);
}


double trejsy_ND_NN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12ND(w,q)+
		L11(E,w,q)*H11ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_NN_t(w,q,RPA,nowy);
}


double trejsy_ND_NNN_t (double E, double w, double q, bool RPA, bool nu, bool nowy)
{
	return (L12(E,w,q, nu)*H12ND(w,q)+
		L11(E,w,q)*H11ND(w,q))*sqrt((mecM+sqrt(mecM2+q*q))/2/mecM*(mecMD+sqrt(mecMD2+q*q))/2/mecMD*fperf*mecMD/(mecMD2+q*q)*mecM/sqrt(mecM2+q*q))*
		Pi_ND_NNN_t(w,q,RPA);
}


double rozn_NN (double q, double w, double E, bool RPA, bool nu, bool nowy)
{
	return mecstala*(-mecVol/Pi)*q/Pi/32/E/E*
		(trejsy_NN_NN_l(E,w,q,RPA,nu,nowy)+trejsy_NN_NN_t(E,w,q,RPA,nu,nowy)+
		trejsy_DD_NN_l(E,w,q,RPA,nu,nowy)+trejsy_DD_NN_t(E,w,q,RPA,nu,nowy)+
		2*trejsy_ND_NN_l(E,w,q,RPA,nu,nowy)+2*trejsy_ND_NN_t(E,w,q,RPA,nu,nowy)
		);
}


double rozn_NNN (double q, double w, double E, bool RPA, bool nu, bool nowy)
{
	return mecstala*(-mecVol/Pi)*q/Pi/32/E/E*
		(trejsy_NN_NNN_l(E,w,q,RPA,nu,nowy)+trejsy_NN_NNN_t(E,w,q,RPA,nu,nowy)+
		trejsy_DD_NNN_l(E,w,q,RPA,nu,nowy)+trejsy_DD_NNN_t(E,w,q,RPA,nu,nowy)+
		2*trejsy_ND_NNN_l(E,w,q,RPA,nu,nowy)+2*trejsy_ND_NNN_t(E,w,q,RPA,nu,nowy)
		);
}


double qminus(double w, double E)//granice calkowania dla przekazu pedu
{
	return
		mecmax2(w, E-sqrt((E-w)*(E-w)-mecm2) );
	//fermilab hypothesis and corrections...
	;
}


double qplus(double w, double E)
{
	return
		mecmin2(
		sqrt(w*w-mecm2+2*E*(E-w)+2*E*sqrt((E-w)*(E-w)-mecm2))
		,meckf+sqrt(meckf2+w*w+2*w*mecEf) )
		;
}


///////////////////////////////////////////////////////////////////////////////////
void model_2body2 (nucleus t, double E, double w, double q, double Bin, particle &meclep, particle &nuc1, particle &nuc2, particle &nucini1, particle &nucini2, 
		   bool fsi, double poten)
{
	//it is assumed that neutrino direction is (0,0,1); perhaps should be relaxed...
	mecm=meclep.mass();
	mecm2=mecm*mecm;

	double muonmom = sqrt ((E-w)*(E-w) - mecm2);
	double cosmuon = (2.0*E*(E-w) - mecm2 - q*q + w*w)/2.0/E/muonmom;
	double phi = 2.0*Pi*los();

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

	double mass1=nuc1.mass();
	double mass2=nuc2.mass();

	do
	{
	  particle probe1 = t.get_nucleon ();
	  vec pos (probe1.x, probe1.y, probe1.z);
	  
	  double localfermi = t.localkf (probe1);
	  
	  nuc1.r.x=probe1.r.x;
	  nuc1.r.y=probe1.r.y;
	  nuc1.r.z=probe1.r.z;
	  
	  nuc2.r.x=probe1.r.x;
	  nuc2.r.y=probe1.r.y;
	  nuc2.r.z=probe1.r.z;
	  /*
	  nucleon[0] = t.get_nucleon (pos);
	  nucleon[1] = t.get_nucleon (pos);
	  
	  particle probe2 = t.get_nucleon ();
	  particle probe3 = t.get_nucleon ();
	  N1=probe1.p4();
	  N2=probe2.p4();
	  N2=probe3.p4();
	  */
		N1=rand_from_ball (localfermi);
		N2=rand_from_ball (localfermi);
		/*
		N1=spectral_choice (6, 6);
		N2=-N1;
		*/
	
		//N1=rand_from_ball (meckf);
		//N2=rand_from_ball (meckf);
		
		
		//cout<<E<<"  "<<w<<"  "<<q<<"  "<<N1<<"  "<<N2<<endl;
		length1= sqrt(mecM2 + N1.norm2());
		length2= sqrt(mecM2 + N2.norm2());//total energy

		N11=vect(N1, length1);
								 //Fermi energy is not subtracted
		N22=vect(N2, length2);

		vect sum=N11+N22;

		sumadd = sum+qqq;
	}

								 //to be able to make Lorentz boost and to "decay"
	while ( ( sumadd.t < sumadd.length() )   ||   ( sumadd*sumadd <(mass1+mass2)*(mass1+mass2) ) );
	//cout<<"po"<<endl;
	trans= sumadd.v();
	sumadd.boost2(trans);		 //boost to the CM frame

	Winv= sumadd.t;

	vec F1=rand_dir ()*sqrt( (Winv*Winv-mass1*mass1 -mass2*mass2)*(Winv*Winv-mass1*mass1 -mass2*mass2) -  4.0*mass1*mass1*mass2*mass2)/2.0/Winv;
	vec F2=-F1;

								 //hadron four momenta
	vect nuc11((Winv*Winv+mass1*mass1-mass2*mass2)/2.0/Winv,F1);
	nuc11.t= (Winv*Winv+mass1*mass1-mass2*mass2)/2.0/Winv;
	vect nuc22((Winv*Winv-mass1*mass1+mass2*mass2)/2.0/Winv ,F2);
	nuc22.t= (Winv*Winv-mass1*mass1+mass2*mass2)/2.0/Winv;

	nuc11.boost2(-trans);
	nuc22.boost2(-trans);		 //nucleons in the LAB frame
/*
	if (!fsi)
	{
		spowalniacz(Bin,nuc11);
		spowalniacz(Bin,nuc22);
	}

	if (fsi)
	{
		spowalniacz(-poten,nuc11);
		spowalniacz(-poten,nuc22);
	}
*/
	vec nucmom1 = vec (nuc11);
	vec nucmom2 = vec (nuc22);
	//cout<<"ped_x="<<cohleptonmom_x<<"   wynik="<<endl;
	meclep.set_momentum(kprim);
	nuc1.set_momentum(nucmom1);
	nuc2.set_momentum(nucmom2);
	//cout<<"dormo"<<endl;
	nucini1.set_momentum(N1);
	nucini2.set_momentum(N2);

}


///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void model_3body (nucleus t, double E, double w, double q, double Bin, particle &meclep, particle &nuc1, particle &nuc2, particle &nuc3, particle &nucini1, particle &nucini2, 
		  particle &nucini3, bool fsi, double poten)
{
	//it is assumed that neutrino direction is (0,0,1); perhaps should be relaxed...

	double muonmom = sqrt ((E-w)*(E-w) - mecm2);
	double cosmuon = (2.0*E*(E-w) - mecm2 - q*q + w*w)/2.0/E/muonmom;
	double phi = 2.0*Pi*los();

								 //momentum transfer
	vec qq(cos(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), sin(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon) , (E - muonmom*cosmuon));

								 //muon momentum
	vec kprim(-cos(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), -sin(phi)*muonmom*sqrt(1.0 - cosmuon*cosmuon), muonmom*cosmuon);
	vect qqq(qq,w);				 //muon 4-momentum

	vec N1, N2, N3;
	vect N11, N22, N33, sumadd;
	vec maksio_vec;

	double length1, length2, length3;
	double length;
	double Winv;

	vec trans;

	double mass1=nuc1.mass();
	double mass2=nuc2.mass();
	double mass3=nuc3.mass();

	do
	{
	  particle probe1 = t.get_nucleon ();
	  vec pos (probe1.x, probe1.y, probe1.z);
	  
	  double localfermi = t.localkf (probe1);
	  
	  nuc1.r.x=probe1.r.x;
	  nuc1.r.y=probe1.r.y;
	  nuc1.r.z=probe1.r.z;
	  
	  nuc2.r.x=probe1.r.x;
	  nuc2.r.y=probe1.r.y;
	  nuc2.r.z=probe1.r.z;
	  
	  nuc3.r.x=probe1.r.x;
	  nuc3.r.y=probe1.r.y;
	  nuc3.r.z=probe1.r.z;
	  
	  /*
	  nucleon[0] = t.get_nucleon (pos);
	  nucleon[1] = t.get_nucleon (pos);
	  
	  particle probe2 = t.get_nucleon ();
	  particle probe3 = t.get_nucleon ();
	  N1=probe1.p4();
	  N2=probe2.p4();
	  N2=probe3.p4();
	  */
		N1=rand_from_ball (localfermi);
		N2=rand_from_ball (localfermi);
		N3=rand_from_ball (localfermi);
	  
		//cout<<E<<"  "<<w<<"  "<<q<<"  "<<N1<<"  "<<N2<<endl;
		length1= sqrt(mecM2 + N1.norm2());
		length2= sqrt(mecM2 + N2.norm2());
		length3= sqrt(mecM2 + N3.norm2());

		N11=vect(N1, length1);
								 //Fermi energy is NOT subtracted
		N22=vect(N2, length2);
		N33=vect(N3, length3);

		vect sum=N11+N22+N33;

		sumadd = sum+qqq;
	}

								 //to be able to make Lorentz boost and to "decay"
	while ( ( sumadd.t < sumadd.length() )   ||   ( sumadd*sumadd < (mass1+mass2+mass3)*(mass1+mass2+mass3)) );
	//cout<<"po"<<endl;
	trans= sumadd.v();
	sumadd.boost2(trans);		 //boost to the CM frame

	Winv= sumadd.t;

	vect nuc11, nuc22, nuc33;

	threebody (Winv, mass1, mass2, mass3, nuc11, nuc22, nuc33);

	nuc11.boost2(-trans);
	nuc22.boost2(-trans);		 //nucleons in the LAB frame
	nuc33.boost2(-trans);
/*
	if (!fsi)
	{
		spowalniacz(Bin,nuc11);
		spowalniacz(Bin,nuc22);
		spowalniacz(Bin,nuc33);
	}

	if (fsi)
	{
		spowalniacz(-poten,nuc11);
		spowalniacz(-poten,nuc22);
		spowalniacz(-poten,nuc33);
	}
*/
	vec nucmom1 = vec (nuc11);
	vec nucmom2 = vec (nuc22);
	vec nucmom3 = vec (nuc33);
	//cout<<"ped_x="<<cohleptonmom_x<<"   wynik="<<endl;
	meclep.set_momentum(kprim);
	nuc1.set_momentum(nucmom1);
	nuc2.set_momentum(nucmom2);
	nuc3.set_momentum(nucmom3);
	//cout<<"dormo"<<endl;
	nucini1.set_momentum(N1);
	nucini2.set_momentum(N2);
	nucini3.set_momentum(N3);

}


///////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//      mecweight2 returns Xsec per nucleon for the Marteau model
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double mecweight2 (double E, bool nu, bool cc, nucleus& t, NSNWRO::params &p, particle &meclepton, particle &mecnucleon1,
		   particle &mecnucleon2, particle &mecnucleon3, particle &mecnucini1, 
		   particle &mecnucini2, particle &mecnucini3,
bool fsi, double poten, bool nowy)
{
	int pdg1=nu ? PDG::pdg_proton:PDG::pdg_neutron;
	int pdg2=nu ? PDG::pdg_neutron:PDG::pdg_proton;

	mecm=meclepton.mass();
	mecm2=mecm*mecm;
	double cut =0.0;
	double w = cut + (E-mecm-cut)*los();
	double q = qminus(w,E) + (qplus(w,E)-qminus(w,E))*los();
	vect nuc1, nuc2, nuc3;
	
// here is the isospin model; I assume that there is a fraction of p.mec_ratio_pp pn pairs in the initial state
	mecnucleon1.pdg = pdg1;		 //there must be at least one pdg1 in the final state
	
	double losso = los();
	 
	if (losso<p.mec_ratio_pp) // was 0.6
	{
		mecnucleon2.pdg = pdg1;
		mecnucini1.pdg = pdg1;
		mecnucini2.pdg = pdg2;
	}
	else
	{
		mecnucleon2.pdg = pdg2;
		mecnucini1.pdg = pdg2;
		mecnucini2.pdg = pdg2;
	}
	losso = los();
	if (losso<0.5)  // was losso<0.8
	{
		mecnucleon3.pdg = pdg1;
		mecnucini3.pdg = pdg1;
	}
	else
	{
		mecnucleon3.pdg = pdg2;
		mecnucini3.pdg = pdg2;
	}

	mecnucleon1.set_mass (PDG::mass (mecnucleon1.pdg));
	mecnucleon2.set_mass (PDG::mass (mecnucleon2.pdg));
	mecnucleon3.set_mass (PDG::mass (mecnucleon3.pdg));
	
	mecnucini1.set_mass (PDG::mass (mecnucini1.pdg));
	mecnucini2.set_mass (PDG::mass (mecnucini2.pdg));
	mecnucini3.set_mass (PDG::mass (mecnucini3.pdg));


	double weight2 = ( qplus(w,E)-qminus(w,E) )*(E-mecm-cut)*rozn_NN (q, w, E, true, nu, nowy)/12.0/1e38;
	if (weight2<0)
	{
		weight2 = 0;
	}
		
	double weight3 = ( qplus(w,E)-qminus(w,E) )*(E-mecm-cut)*rozn_NNN (q, w, E, true, nu, nowy)/12.0/1e38;
	if (weight3<0)
	{
		weight3 = 0;
	}
		//cout<<E<<"  "<<w<<"  "<<q<<endl;

	if ( weight2>(weight2+weight3)*los() )
	{							 //cout<<"dwa"<<endl;
		model_2body2 (t, E, w, q, 8, meclepton, mecnucleon1, mecnucleon2, mecnucini1, mecnucini2, fsi, poten);
		//cout<<"2body"<<endl;
		vec ped3 = vec (0,0,0);
		mecnucleon3.set_momentum(ped3);
	}
	else
	{							 //cout<<"trzy"<<endl;
		model_3body (t, E, w, q, 8, meclepton, mecnucleon1, mecnucleon2, mecnucleon3, mecnucini1, mecnucini2, mecnucini3, fsi, poten);
		//cout<<"3body"<<endl;
	}


	return weight2+weight3;

}
