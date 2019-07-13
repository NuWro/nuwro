// Postac Dipolowa Form Faktorow wedluf A. Bodek et. al... hep-ex/0308005

#ifndef _form_faktory_Dipol_h_
#define _form_faktory_Dipol_h_
#include <iostream.h>
#include <math.h>
#include "jednostki.h"

namespace rpa
{
const long double gA_Dipol   = -1.267;
const long double MA_Dipol   = 1.00*GeV;
const long double magneton   = 4.71;   
const long double MV2_Dipol = 0.71*GeV*GeV; 

//const magneton_N = ;   // znajdz ????

//long double magneton_proton   = 2.793*magneton_N;
//long double magneton_neutron = -1.913*magneton_N;

long double G_A_Dipol(long double Q2)                                         
       {                                                                
       return  gA_Dipol/(pow(1 + Q2/MA_Dipol/MA_Dipol,2));                 
       }                                                                

long double F_p_Dipol(long double Q2)
{
return 2*M12*G_A_Dipol(Q2)/(m_pi*m_pi + Q2); // mam ruznice z Bodek - M12*M12
}
                                                                        
long double G_E_Dipol(long double Q2)                                      
           {                                                           
	    return  1/pow(1+Q2/MV2_Dipol,2); 
	   }                                                            
                                                                        
long double G_M_Dipol(long double Q2)
          { 
	  return  G_E_Dipol(Q2)*magneton; 
	  }       
                                                                                                                                 
long double F_1_Dipol(long double Q2)                   
        {                                                               
	return  (Q2*G_M_Dipol(Q2) + 4*M12*M12*G_E_Dipol(Q2) )/(Q2+4*M12*M12) ;
	}                                                               
					                                
long double F_2_Dipol(long double Q2)
        {                                                              
 return 4*M12*M12*( G_M_Dipol(Q2)-G_E_Dipol(Q2) )/(4*M12*M12+Q2);    
	}                                                              
}
#endif
