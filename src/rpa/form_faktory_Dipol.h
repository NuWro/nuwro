// Postac Dipolowa Form Faktorow wedluf A. Bodek et. al... hep-ex/0308005 

#ifndef _form_faktory_Dipol_h_
#define _form_faktory_Dipol_h_
#include<iostream>
#include<cmath>
#include "jednostki.h"

namespace rpa
{
const double gA_Dipol   = -1.267;
const double MA_Dipol   = 1.00*GeV;
const double magneton   = 4.71;   
const double MV2_Dipol = 0.71*GeV*GeV; 

//const magneton_N = ;   // znajdz ????

//double magneton_proton   = 2.793*magneton_N;
//double magneton_neutron = -1.913*magneton_N;

double G_A_Dipol(double Q2)                                         
{                                                                
	return  gA_Dipol/(pow(1 + Q2/MA_Dipol/MA_Dipol,2));                 
}                                                                

double F_p_Dipol(double Q2)
{
	return 2*M12*G_A_Dipol(Q2)/(m_pi*m_pi + Q2); // mam ruznice z Bodek - M12*M12
}
                                                                        
double G_E_Dipol(double Q2)                                      
{                                                           
	return  1/pow(1+Q2/MV2_Dipol,2); 
}                                                            
                                                                        
double G_M_Dipol(double Q2)
{ 
	return  G_E_Dipol(Q2)*magneton; 
}       
                                                                                                                                 
double F_1_Dipol(double Q2)                   
{                                                               
	return  (Q2*G_M_Dipol(Q2) + 4*M12*M12*G_E_Dipol(Q2) )/(Q2+4*M12*M12) ;
}                                                               
					                                
double F_2_Dipol(double Q2)
{                                                              
	return 4*M12*M12*( G_M_Dipol(Q2)-G_E_Dipol(Q2) )/(4*M12*M12+Q2);    
}                                                              
}
#endif
