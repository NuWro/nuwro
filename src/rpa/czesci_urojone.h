#ifndef _czesci_urojone_h_
#define _czesci_urojone_h_
#include<iostream>
#include<cmath>
#include "jednostki.h"
#include "stale_rpa.h"


namespace rpa
{
double Ming(double a,double b){ return (a<b ? a : b); }


double Maxg(double a, double b, double c)
{  
	double x = (a>b ? a : b); 
	return (x > c ? x : c) ;
}

double E_max(double q0,double qv)
{ 
	return Maxg( Mef , Ef - q0, 0.5*(- q0+qv * sqrt(1-4*Mef2/(q0*q0-qv*qv) )));
}

double E_(double q0,double qv)
       { return Ming(Ef, E_max(q0,qv));}

       
double E1(double q0,double qv)
{
	return Ef - E_(q0,qv); 
}

double E2(double q0, double qv)
{
	return (Ef*Ef -E_(q0,qv)*E_(q0,qv))/2;
}

double E3(double q0, double qv)
{
	return (Ef*Ef*Ef -E_(q0,qv)*E_(q0,qv)*E_(q0,qv))/3; 
}

double Im_H_vv_l(double q0,double qv)
{ 

	return ((q0*q0-qv*qv)/(2*Pi*qv*qv*qv))*
			( E3(q0,qv) + q0*E2(q0,qv) + (q0*q0-qv*qv)*E1(q0,qv)/4 ) ; 

}

double Im_H_vv_t(double q0,double qv)
{ 

	return  2*((q0*q0-qv*qv)/(4*Pi*qv*qv*qv))*   // dwojka
		 	(E3(q0,qv)+q0*E2(q0,qv)+(qv*qv*Mef2/(q0*q0-qv*qv)+0.25*(q0*q0+qv*qv))*E1(q0,qv)) ;
}

double Im_H_tt_l(double q0, double qv) 
{ 
	return   - ( ( q0*q0-qv*qv )/( 8*Pi*qv*qv*qv*M*M))*
   				( (q0*q0-qv*qv)*E3(q0,qv) + q0*(q0*q0-qv*qv)*E2(q0,qv)+
				( qv*qv*Mef2 +0.25*(q0*q0-qv*qv)*q0*q0 )*E1(q0,qv) );
}

double Im_H_tt_t(double q0, double qv)  
{ 

	return  ( (q0*q0-qv*qv)/(8*Pi*qv*qv*qv*M*M) )* // bylo 16 zamiast 8;
			( ( Mef2*qv*qv - 0.25*(q0*q0-qv*qv)*(q0*q0-qv*qv) )*E1(q0,qv)
			- q0*(q0*q0-qv*qv)*E2(q0,qv) -(q0*q0-qv*qv)*E3(q0,qv) ) ;
  
}

double Im_H_vt_l(double q0, double qv) 
{ 
	return  -( (q0*q0-qv*qv)*E1(q0,qv) *Mef )/(8*Pi*qv*M); 
}

double Im_H_vt_t(double q0, double qv)  
{ 
	return  -2*Im_H_vt_l(q0,qv);                             // dwojka
}

double Im_H_va(double q0, double qv)    
{ 
	return  ((q0*q0-qv*qv)/(8*Pi*qv*qv*qv))*(2*E2(q0,qv) + q0*E1(q0,qv));
}

double Im_H_a(double q0, double qv)  
{ 
	return (Mef2 * E1(q0,qv) )/(2*Pi*qv);
}
}
#endif

