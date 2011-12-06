// Forma Faktory wziete z pracy: Krzysztof M. Graczyk, Piotr Plonski, Robert Sulej JHEP 1009:053,2010

#ifndef _nnff_h_
#define _nnff_h_

#include "jednostki.h"
#include "pdg.h"
#include <cmath>



const double magneton_N_nnff       = 1;   // tego nie jestem pewien !!!!
const double magneton_proton_nnff  = 2.793*magneton_N_nnff;
const double magneton_neutron_nnff =-1.913*magneton_N_nnff;
static const double M12=(PDG::mass_proton+PDG::mass_neutron)/2;
//static const double sin2thetaW=0.23122;


// magneton_proton; ///  Znalezc
// magneton_neutron; /// Znalezc



 double tab_gen[7]  = { 10.19704, 2.36812, -1.144266,-4.274101,0.8149924,2.985524,-0.7864434};

 double tab_gmn[10] = { 3.19646, 2.565681, 6.441526, -2.004055, -0.2972361, 3.606737, -3.135199, 0.299523, 1.261638, 2.64747};

 double tab_gep[10] = { 3.930227, 0.1108384, -5.325479, -2.846154, -0.2071328, 0.8742101, 0.4283194, 2.568322, 2.577635, -1.185632};

 double tab_gmp[10] = {-2.862682, -1.560675, 2.321148, 0.1283189, -0.2803566, 2.794296, 1.726774, 0.861083, 0.4184286, -0.1526676};

 double tab_axial[10] = {-26.10885 ,1.823041, -8.391283,-7.737312, 15.27646, 0.3992788,-1.350184,-0.2021121,-2.870517, 3.879841};
 


double gA_nnff = -1.267;
double MA_nnff =  1.015*GeV;     // 1.03*GeV;	

double funkcja1( double Q2, double *w )
{

double q2 = Q2/GeV/GeV;

return w[4]/(1 + exp( - (q2*w[0] + w[1]) ) )+ w[5]/(1 + exp( - (q2*w[2] + w[3]) ))  + w[6];

}

double funkcja2( double Q2, double *w )
{
double q2 = Q2/GeV/GeV;

return w[6]/(1 + exp( - (q2*w[0] + w[1]) ) ) 
     + w[7]/(1 + exp( - (q2*w[2] + w[3]) )) 
     + w[8]/(1 + exp( - (q2*w[4] + w[5]) )) + w[9];

}

double G_p_E_nnff(double Q2)
{
return  funkcja2( Q2, tab_gep)/pow(1.0 + Q2/0.71/GeV/GeV,2);
}

double G_p_M_nnff(double Q2)
{
return funkcja2( Q2, tab_gmp)*magneton_proton_nnff/pow(1.0 + Q2/0.71/GeV/GeV,2);
}


double G_n_M_nnff(double Q2)
{
return funkcja2( Q2, tab_gmn)*magneton_neutron_nnff/pow(1.0 + Q2/0.71/GeV/GeV,2);
}

double G_n_E_nnff(double Q2)
{
return funkcja1( Q2, tab_gen);
}

/// Form Factory Charged Current


double G_E_nnff(double Q2)
{
return G_p_E_nnff(Q2) -G_n_E_nnff(Q2);
}

double G_M_nnff(double Q2)
{
return G_p_M_nnff(Q2) - G_n_M_nnff(Q2);
}

double F_1_p_nnff(double Q2)
{
return (G_p_E_nnff(Q2) +Q2*G_p_M_nnff(Q2)/4/M12/M12)/(1+Q2/4/M12/M12);
}

double F_1_n_nnff(double Q2)
{
return (G_n_E_nnff(Q2) +Q2*G_n_M_nnff(Q2)/4/M12/M12)/(1+Q2/4/M12/M12);
}

double F_2_n_nnff(double Q2)
{
return (G_n_M_nnff(Q2) -G_n_E_nnff(Q2))/(1+Q2/4/M12/M12);
}

double F_2_p_nnff(double Q2)
{
return (G_p_M_nnff(Q2) -G_p_E_nnff(Q2))/(1+Q2/4/M12/M12);
}

double F_1_nnff(double Q2)
{
return (G_E_nnff(Q2) +Q2*G_M_nnff(Q2)/4/M12/M12)/(1+Q2/4/M12/M12);
}

double F_2_nnff(double Q2)
{
return (G_M_nnff(Q2) -G_E_nnff(Q2))/(1+Q2/4/M12/M12);
}

/// Form Factory Dziwne

double G_E_s(double Q2)
{
return 0;
}

double G_M_s(double Q2)
{
return 0;
}


double F_1_s(double Q2)
{
return 0;
}

double F_2_s(double Q2)
{
return 0;
}

/// Form Factory Neutral Current


double G_E_p_NC_nnff(double Q2)
{
return  0.5*G_E_nnff(Q2) - 2*sin2thetaW*G_p_E_nnff(Q2) - 0.5*G_E_s(Q2);  /// Na razie bez dziwnego wkladu
}

double G_E_n_NC_nnff(double Q2)
{
return -0.5*G_E_nnff(Q2) - 2*sin2thetaW*G_n_E_nnff(Q2) - 0.5*G_E_s(Q2);  /// Na razie bez dziwnego wkladu
}


double G_M_p_NC_nnff(double Q2)
{
return 0.5*G_M_nnff(Q2) - 2*sin2thetaW*G_p_M_nnff(Q2) - 0.5*G_M_s(Q2);  /// Na razie bez dziwnego wkladu
}

double G_M_n_NC_nnff(double Q2)
{
return -0.5*G_M_nnff(Q2) - 2*sin2thetaW*G_n_M_nnff(Q2) - 0.5*G_M_s(Q2);  /// Na razie bez dziwnego wkladu
}

double F1_p_NC_nnff(double Q2){

return 0.5*F_1_nnff(Q2) - 2*sin2thetaW*F_1_p_nnff(Q2) - 0.5*F_1_s(Q2) ;
}

double F1_n_NC_nnff(double Q2){

return -0.5*F_1_nnff(Q2) - 2*sin2thetaW*F_1_n_nnff(Q2) - 0.5*F_1_s(Q2) ;
}

double F2_p_NC_nnff(double Q2){

return 0.5*F_2_nnff(Q2) - 2*sin2thetaW*F_2_p_nnff(Q2) - 0.5*F_2_s(Q2) ;
}

double F2_n_NC_nnff(double Q2){

return -0.5*F_2_nnff(Q2) - 2*sin2thetaW*F_2_n_nnff(Q2) - 0.5*F_2_s(Q2) ;
}


#endif                                                                                                                               
