#ifndef _FormFaktory_h_
#define _FormFaktory_h_
#include <iostream>
#include <cmath>
#include "jednostki.h"
#include "stale_rpa.h"

namespace rpa
{

const char FormFaktor = 'Z'; 

double MA_MiniBoonE = 1.20*GeV;                                                               



double G_A(double q0,double qv)                                         
{ 
	double q2 = q0*q0-qv*qv;

	double MA=1.20*GeV;                                                               
	return  -1.26/pow(1-q2/(MA*MA),2);                 
}                                                                



                                                                        
double G_E_1(double q0,double qv)                                      
{
	double MV2 = 0.71*GeV*GeV;                                                           
	return  1/( (1-(q0*q0-qv*qv)/MV2)*(1-(q0*q0-qv*qv)/MV2) ); 
}                                                            
                                                                        
double G_E(double q0, double qv)                  
{                                                           
	return  G_E_1(q0,qv) ;                                
}                                                            
                                                                        
double G_M(double q0,double qv){ return  G_E_1(q0,qv)*magneton ; }       
                                                                                                                                 
double F_1(double q0,double qv)                   
{                                                               
	return  ((q0*q0-qv*qv)*G_M(q0,qv) - 4*M12_2*G_E(q0,qv) )/        
		    ((q0*q0-qv*qv)-4*M12_2) ;
}                                                               
					                                
double F_2(double q0,double qv)
{                                                              
	return  4*M12_2*( G_M(q0,qv)-G_E(q0,qv) )/(4*M12_2-(q0*q0-qv*qv));    
}
                                                              
}
#endif
