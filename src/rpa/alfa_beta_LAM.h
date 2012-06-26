#ifndef  _alfa_beta_LAM_h
#define _alfa_beta_LAM_h
#include <iostream>
#include <fstream>
#include <cmath>
#include "jednostki.h"
#include "stale_rpa.h"


namespace rpa
{
double logEf = log((kf+Ef)/Mef);

 
 
 double alfa(double q0, double qv)
 {                             
	double q2=q0*q0-qv*qv;
	double q4=q2*q2;
	double licznik        = fabs(q4 - 4*pow( q0*Ef - qv*kf, 2 ) );
	double mianownik = fabs(q4 - 4*pow(q0*Ef + qv*kf,  2));
	return log(licznik/mianownik );
 }
 double beta(double q0,double qv)
 {
	double q02 =q0*q0;
	double q2   =q02-qv*qv;
	double licznik =fabs((pow(q2 + 2*qv*kf, 2) - 4*q02*Ef2));
	double mianownik =fabs(pow(q2 - 2*qv*kf, 2) - 4*q02*Ef2 );
	return log(licznik/mianownik ) ;
 }
 double LAMBDA(double q0, double qv)
 {
	 double q2      =   q0*q0-qv*qv;
	 double qvkf    =   qv*kf;
	 double q0Ef   =  q0*Ef;
	 double q4Ef  = q2*q2*Ef;
	 double stal=sqrt(q2*(q2 -4*Mef2));
	 double q04Mef2 = 4*q0*Mef2;
	 double kfq2stal = kf*q2*stal;
	 double licznik1        = fabs( q4Ef - q04Mef2*(q0Ef - qvkf) - kfq2stal);
	 double Mianownik1 = fabs( q4Ef - q04Mef2*(q0Ef  - qvkf ) + kfq2stal);
	 double licznik2        = fabs( q4Ef - q04Mef2*(q0Ef + qvkf) - kfq2stal);
	 double Mianownik2 = fabs( q4Ef - q04Mef2*(q0Ef + qvkf) + kfq2stal);
	 
	if(q2==0)  
		return 0;
	else 
		return (log(licznik1/Mianownik1 ) 
				+log(licznik2/Mianownik2)  
				)/2/stal; 
  }
 
 
 int fey =0;   //Czesci Faynmanowskie sa rowne zero dla fey = 0;
 
 
 
 
 double ReHs(double q0, double qv)
 {
	 double q2 = q0*q0-qv*qv;
	 
	 double nawias = kf*Ef 
							  - ( 3*Mef2 - 0.5*q2 )*logEf 
				  - q0*(4*Mef2-q2)*alfa(q0,qv) /8/qv 
				  + Ef*(4*Mef2-q2)*beta(q0,qv)  /4/qv 
				  +pow(4*Mef2-q2,2)*LAMBDA(q0,qv)/4;
			 
	 double eta  = sqrt(1-4*Mef2/q2);		 
	 double rehs = M*M - 4*M*Mef + 13*Mef2/3- 4*q2/9 - (2*Mef2 -q2/3)*log(Mef/M)
				   -(4/3)*(Mef2-q2/4)*eta*atan(eta);
	 
	  
	 return (nawias + fey*3*rehs/2)/2/Pi/Pi;
 }
 
 double ReHl(double q0, double qv)
 {
	 double q2 = q0*q0-qv*qv;
	 double qv2= qv*qv;
	 double q02= q0*q0;
	 
	 double nawias= 2*kf*Ef/3 
							  - (qv2/6)*logEf 
				 +(q0/4/qv)*(Ef2 + (q02 - 3*qv2)/12 )*alfa(q0,qv) 
						  - Ef*(3*q2 + 4*Ef2)*beta(q0,qv)/24/qv
						 + qv2*(2*Mef2 + q2)*(4*Mef2 - q2)*LAMBDA(q0,qv)/12/q2 ;
	 
	 double eta  = sqrt(1 - 4*Mef2/q2);		 
	 
	 double rehl = -(q2/18 + q2*log(Mef/M)/3 +(q2+2*Mef2)*(eta*atan(eta)-1)/3);
	  
	 
	 return nawias*q2/Pi/Pi/qv2 + fey*rehl/(2*Pi*Pi);
 }
 
 double ReHt(double q0, double qv)
 {
	   double qv2=qv*qv;
	   double q02=q0*q0;
	double q2 = q02-qv2;
	   
	   double nawias= kf*Ef*(1+ 2*q02/qv2 )/3 
							 + (q2/3)*logEf
							 + (q0*q2/(4*qv2*qv))*((q02 + 3*qv2)/12 + qv2*Mef2/q2+ Ef2)*alfa(q0,qv)
							 + (Ef/qv)*((qv2*qv2-q02*q02)/8/qv2 - Mef2/2 - q2*Ef2/6/qv2)*beta(q0,qv)
							  -  (2*Mef2+q2)*(4*Mef2 - q2)*LAMBDA(q0,qv)/6;
	 
	   double eta  = sqrt(1-4*Mef2/q2);		 
	   double reht = q2/18 + q2*log(Mef/M)/3 +(q2+2*Mef2)*(eta*atan(eta)-1)/3;
	  
	  return nawias/2/Pi/Pi + fey*2*reht/(2*Pi*Pi);  // tutaj przemnozylem przez 2
 }

 double ReH0(double q0, double qv)
 {
	 double q2 = q0*q0-qv*qv;
	 double qv2=qv*qv;
	 
	 double nawias = qv*kf 
							 + 0.5*q0*Ef*alfa(q0,qv)
							  -(q2/8 + Ef2/2 + qv2*Mef2/2/q2)*beta(q0,qv);
	 
	 return nawias *Mef/2/Pi/Pi/qv;
 }

}
#endif
