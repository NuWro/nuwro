// Forma Faktory wziete z pracy: A. Bodek, H. Budd, and J. Arrington hep-ex/0308005 
// Zwane pozniej BBA - 2003 

#ifndef _from_faktory_BBA_2003_h_
#define _from_faktory_BBA_2003_h_
#include "jednostki.h"
#include <cmath>

namespace rpa
{
const double   MV2 = 0.71*GeV*GeV;
const double   gA_BBA_2003 = -1.267;
const double   MA_BBA_2003 = 1.0*GeV;

const double a2_p_E = 3.253/GeV/GeV; 
const double a2_p_M = 3.104/GeV/GeV;
const double a2_n_M = 3.043/GeV/GeV;

const double a4_p_E = 1.422/GeV/GeV/GeV/GeV; 
const double a4_p_M = 1.428/GeV/GeV/GeV/GeV;
const double a4_n_M = 0.8548/GeV/GeV/GeV/GeV;

const double a6_p_E = 0.08582/GeV/GeV/GeV/GeV/GeV/GeV; 
const double a6_p_M = 0.1112/GeV/GeV/GeV/GeV/GeV/GeV;
const double a6_n_M = 0.6806/GeV/GeV/GeV/GeV/GeV/GeV;

const double a8_p_E  = 0.3318/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV; 
const double a8_p_M = -0.006981/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV;
const double a8_n_M = -0.1287/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV;

const double a10_p_E  = - 0.09371/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV; 
const double a10_p_M = 0.0003705/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV;
const double a10_n_M = 0.008912/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV;

const double a12_p_E  = 0.01076/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV; 
const double a12_p_M = -0.7063e-05/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV/GeV;


const double G_p_E_BBA_2003_zero  =1;  ///  Znalezc
const double G_p_M_BBA_2003_zero  =magneton_proton; ///  Znalezc
const double G_n_M_BBA_2003_zero  =magneton_neutron; /// Znalezc


double G_A_BBA_2003(double Q2)
{
return gA_BBA_2003/pow(1 +Q2/MA_BBA_2003/MA_BBA_2003,2); 
}

double F_p_BBA_2003(double Q2)
{
return 2*M12*gA_BBA_2003/pow(1 + Q2/MA_BBA_2003/MA_BBA_2003,2)/(m_pi*m_pi+Q2); // mam roznice z Bodek M12*M12
}

double G_p_E_BBA_2003(double Q2)
{
return G_p_E_BBA_2003_zero/(1+a2_p_E*Q2 +a4_p_E*Q2*Q2 +a6_p_E*Q2*Q2*Q2 + a8_p_E*Q2*Q2*Q2*Q2
                                                 +a10_p_E*Q2*Q2*Q2*Q2*Q2 + a12_p_E*Q2*Q2*Q2*Q2*Q2*Q2);
}

double G_p_M_BBA_2003(double Q2)
{
return G_p_M_BBA_2003_zero/(1+a2_p_M*Q2 +a4_p_M*Q2*Q2 +a6_p_M*Q2*Q2*Q2 + a8_p_M*Q2*Q2*Q2*Q2
                                                 +a10_p_M*Q2*Q2*Q2*Q2*Q2 + a12_p_M*Q2*Q2*Q2*Q2*Q2*Q2);
}

double G_n_M_BBA_2003(double Q2)
{
return G_n_M_BBA_2003_zero/(1+a2_n_M*Q2 +a4_n_M*Q2*Q2 +a6_n_M*Q2*Q2*Q2 + a8_n_M*Q2*Q2*Q2*Q2
                                                 +a10_n_M*Q2*Q2*Q2*Q2*Q2 );
}

double G_n_E_BBA_2003(double Q2)
{
double a = 0.942;
double b = 4.61;
double tau = Q2/4/M12/M12;
return -magneton_neutron*a*tau/(1+b*tau)/pow(1+Q2/MV2,2);
}


double G_E_V_BBA_2003(double Q2)
{
return G_p_E_BBA_2003(Q2) -G_n_E_BBA_2003(Q2);
}

double G_M_V_BBA_2003(double Q2)
{
return G_p_M_BBA_2003(Q2)-G_n_M_BBA_2003(Q2);
}

double F_1_BBA_2003(double Q2)
{
return (G_E_V_BBA_2003(Q2) +Q2*G_M_V_BBA_2003(Q2)/4/M12/M12)/
           (1+Q2/4/M12/M12);
}

double F_2_BBA_2003(double Q2)
{
return (G_M_V_BBA_2003(Q2) -G_E_V_BBA_2003(Q2))/(1+Q2/4/M12/M12);
}

}
#endif                                                                                                                               

