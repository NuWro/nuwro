#include "qel_sigma.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "pdg.h"

using namespace std;

inline double min (double a, double b) {  return (a < b ? a : b);}
inline double max (double a, double b) {  return (a > b ? a : b);}
inline double pow2 (double x)          {  return x * x;}





/////////////////////////////////////////////////////////////////////////////////////
//                 semielastic neutrino nucleon scattering 
//                       I M P L E M E N T A T I O N
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
QEL::QEL (params &p)
/////////////////////////////////////////////////////////////////////////////////////
{  	
		
}


/////////////////////////////////////////////////////////////////////////////////////
QEL::QEL (double lepton_mass0, 
          double hadron_mass0, 
          bool anty0, 
          bool use_lambda0 
)
/////////////////////////////////////////////////////////////////////////////////////
       : m(lepton_mass0), 
         M (hadron_mass0), 
         anty(anty0),
         use_lambda (use_lambda0)
	 
  {
    mm=m*m;
    MM=M*M;
  }
  


/////////////////////////////////////////////////////////////////////////////////////
double QEL::sigma (double q2,int kind)
/////////////////////////////////////////////////////////////////////////////////////
//kind= 0 - cc 
//      1 - nc proton 
//      2 - nc neutron
{ 
    double F1,F2,Fa,Fp; 
    
    list(F1,F2)=f12(q2,kind);// Calculate Form Factors 
    list(Fa,Fp)=fap(q2,kind);
    
//cout<<F1<< ' '<<F2<<' '<<Fa<<' '<<Fp<<" new"<<endl;
    
    double x = q2 / MM ;
    double tau = - x/4 ;
    
    double faa=Fa*Fa;
    double fpp=Fp*Fp;
    double f11=F1*F1;
    double f22=F2*F2;
    
    double F1v=F1,F2v=F2;
      double A =
      (mm - q2) / (4 * MM) * ((4 - x) * faa
			      - (4 + x) * F1v * F1v
			      - x * pow2 (F2v) * (1 + x / 4)
			      - 4 * x * F1v * F2v
			      - mm / MM * (pow2 (F1v + F2v) +
					           pow2 (Fa + 2 * Fp) +
					           (x-4)*fpp
                               )
                );


    double B = -x * Fa * (F1v + F2v);

    double CC = (faa + f11 + tau * f22) / 4;


//  s - u
    double su = 4 * M * Enu + q2 - mm;
    double suM2 = su / MM;
    double ABC = A + suM2*( (anty? B : -B)   + suM2*CC );

    if(kind) // nc proton or nc neutron
	  return  (G*G*MM/8/Pi/Enu/Enu)*ABC;
    else
      return  (G*G*MM/8/Pi/Enu/Enu)*ABC *cos2thetac;
 }


// integration limits q^2
// factor=-1  -lower limit
// factor=+1  -upper limit
/////////////////////////////////////////////////////////////////////////////////////
double QEL::qq (double factor)
/////////////////////////////////////////////////////////////////////////////////////
  {
    double x;
    x = pow2 (2 * Enu * M - mm) - 4 * mm * MM;
    if (x < 0)
      x = 0;
    return -fabs (mm
		  + factor * Enu / (M + 2 * Enu) * sqrt (x)
		  - Enu * (2 * Enu * M + mm) / (M + 2 * Enu));
  }
/*  
/////////////////////////////////////////////////////////////////////////////////////
   double QEL::total_sigma(double E)
/////////////////////////////////////////////////////////////////////////////////////
   { Enu=E;
    return calg5(*this,&QEL::sigma,qq(-1),qq(+1),1e-41*cm2,10); 
   }
   */
/*
double QEL::sigmanc(double q2,int nukleon)
{ 

//double magneton =  4.71   ;
double sin2thetaw = 0.232  ;
//const static double M12_2=M12*M12;
int proton  = 1;
int neutron = -1;  // bylo -1 Czarek
int neutrino =1;
int antyneutrino = -1;
double lp = 1.793;
double ln= -1.913;

double Q2=-q2;
double M=M12;
double MM=M*M;
////////////////////////////////////////////////////////////////////////////
double Mv2=0.71*GeV2;
double g_a=1.267;
const double sin_2_theta_W = 0.23120;
double Ma=1.03*GeV;
double Ma2=Ma*Ma;

double Ge=1.0/pow(1+Q2/Mv2,2.0);
double Gm=4.71*Ge;

//nowe; BBBA05
double tau=Q2/4.0/MM;
double tau2=tau*tau;
double tau3=tau*tau2;
double tau4=tau2*tau2;

double Gep=(1-0.0577*tau)/(1+11.2*tau+13.6*tau2+33*tau3);
double Gmp=2.79*(1+0.15*tau)/(1+11.1*tau+19.6*tau2+7.54*tau3);
double Gen=(1.38*tau-0.214*tau2)/(1+8.51*tau+59.9*tau2+13.6*tau3+2.57*tau4);
double Gmn=-1.91*(1+1.82*tau)/(1+14.1*tau+20.7*tau2+69.7*tau3);


double Ge_B=Gep-Gen;
double Gm_B=Gmp-Gmn;


double Fa=-g_a/pow(1+Q2/Ma2,2.0);
double Fp=0;
double F2=(Gm-Ge)/(1+tau);
double F2p=(Gmp-Gep)/(1+tau);
double F2n=(Gmn-Gen)/(1+tau);

double F1=(Ge+Gm*tau)/(1+tau);
double F1p=(Gep+Gmp*tau)/(1+tau);
double F1n=(Gen+Gmn*tau)/(1+tau);

///////form faktory neutralne^
double MV2_NC = 0.843*0.843*GeV2;

//Trzy parametryzacje form faktorów dziwnych my korzystamy z pierwszej^
double F1s_0 = 0.53;
double F2s_0 = -0.40;
double delta_s = -0.21;
Ma =1.012*GeV;
Ma2=Ma*Ma;


double F1s= F1s_0/(1+tau)/(1+Q2/MV2_NC)/(1+Q2/MV2_NC);
double F2s= F2s_0/(1+tau)/(1+Q2/MV2_NC)/(1+Q2/MV2_NC);
double FAs= delta_s/(1+Q2/Ma2)/(1+Q2/Ma2);


//from faktroy dla protonu i neutronu ^

double F1p_nc = 0.5*(F1p*(1-4*sin_2_theta_W) - F1n - F1s);
double F2p_nc = 0.5*(F2p*(1-4*sin_2_theta_W) - F2n - F2s);
double F1n_nc = 0.5*(F1n*(1-4*sin_2_theta_W) - F1p - F1s);
double F2n_nc = 0.5*(F2n*(1-4*sin_2_theta_W) - F2p - F2s);

double FAn = 0.5*(-Fa + FAs );
double FAp = 0.5*( Fa + FAs );


double AA=0;
double BB=0;
double CC=0;

//czynniki kinematyczne w przekroju proton 
if(nukleon==proton)
{
AA = Q2/M2*( (1+tau)*FAp*FAp  -(1-tau)*F1p_nc*F1p_nc +tau*(1-tau)*F2p_nc*F2p_nc + 4*tau*F1p_nc*F2p_nc);
BB = Q2/M2*FAp*(F1p_nc+F2p_nc);
CC = 0.25*(FAp*FAp + F1p_nc*F1p_nc + tau*F2p_nc*F2p_nc );
}
//czynniki kinematyczne w przekroju neutron
if(nukleon==neutron)
{
AA = Q2/M2*( (1+tau)*FAn*FAn -(1-tau)*F1n_nc*F1n_nc +tau*(1-tau)*F2n_nc*F2n_nc + 4*tau*F1n_nc*F2n_nc);
BB = Q2/M2*FAn*(F1n_nc+F2n_nc);
CC = 0.25*(FAn*FAn + F1n_nc*F1n_nc + tau*F2n_nc*F2n_nc);
}

//w ulkladzie LAB z masa leptonu 0, bo neutrino koncowe

double su = 4*M12*En-Q2;
double suM2 = su / M2;
double ABC = AA + znak * BB * suM2 + CC * suM2 * suM2;
	   
//cout<<"nc="<<(G*G*M2/8/M_PI/En/En)*ABC/cm2<<endl;
//cin.get();

//return M2 * G * G *  / 8 / M_PI / En / En * ABC;

//cout<<"M2="<<M2<<endl;
//cout<<"ABC="<<ABC<<endl;
//cout<<"CR="<<(G*G*M2/8/M_PI/En/En)*(AA + znak* su/M2*BB + su*su/M2/M2*CC)<<endl;
//cin.get();
//return ABC*M*2*GeV2;	       

//cout<<"znak"<<znak<<endl;
return  (G*G*M2/8/M_PI/En/En)*(AA - znak* su/M2*BB + su*su/M2/M2*CC);

/////////////////////////////////////////////////////////////////////////////

}          

*/
