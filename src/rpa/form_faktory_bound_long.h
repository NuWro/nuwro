#ifndef _faktory_bound_h_
#define _faktory_bound_h_
#include <cmath>
#include "jednostki.h"

namespace rpa
{
const long double magneton_bound=4.71;

long double MA_bound =1.0*GeV;
long double MV2_bound=0.7*GeV*GeV;      

long double G_A_bound(long double Q2_)
       { 
       long double Q2 = Q2_/GeV/GeV;
       long double prop = 0.0054*pow(Q2,4) - 0.0361*pow(Q2,3) +0.1115*pow(Q2,2) -0.1239*Q2 + 0.9224;
       return  -1.26*prop/(pow(1+Q2_/MA_bound/MA_bound,2));
       }

long double G_E_bound(long double Q2)
           {
	    
	   return  1/pow(1+Q2/MV2_bound,2) ;
	   }


long double G_M_bound(long double Q2)
           { 
	   return  G_E_bound(Q2)*magneton_bound; 
	   }

long double F_1_bound(long double Q2_)
           { 
	long double Q2 = Q2_/GeV/GeV;
        long double prop = -0.0031*pow(Q2,4) + 0.0061*pow(Q2,3) +0.0386*pow(Q2,2) -0.0274*Q2 + 1.0842;
        return  prop*(Q2_*G_M_bound(Q2_) + 4*M12*M12*G_E_bound(Q2_) )/(Q2_+4*M12*M12) ;
	    }
					 
long double F_2_bound(long double Q2_)
        { 
	long double Q2 = Q2_/GeV/GeV;
        long double prop = -0.0056*pow(Q2,4) + 0.0229*pow(Q2,3) +0.0132*pow(Q2,2) -0.0366*Q2 + 0.9999;
        return  prop*4*M12*M12*( G_M_bound(Q2_)-G_E_bound(Q2_) )/(4*M12*M12 + Q2_);
	}

long double F_p_bound(long double Q2)
        { 
	return 2*G_A_bound(Q2) * M12 / ( m_pi*m_pi+Q2 );
	}
}
#endif
