//Oryginalnie porprawki bodka wprowadzone przez K.Graczyka
//Zmiany przez JN to: 
#include "jednostki.h"
#include "masses.h"
#include "parameters.h"

//double GeV=1000;
//double GeV2=GeV*GeV;


const double freez = 0.24 * GeV2;	// Granica zamrozenia funkcji strukturalnych

double ss (double Q2)			//to jest poczatek GRV94
{
  double a2 = 0.232 * 0.232 * GeV * GeV;
  double mi2 = 0.23 * GeV * GeV;
  if (Q2 <= freez)		//wszystko jest w MeV...
    return 0;
  else
    return log (log (Q2 / a2) / log (mi2 / a2));	//log naturalny
}



double uv (double xx, double Q2)	//u(val), dobre unormowanie!
{
  double s = ss (Q2);
  double s2 = s * s;

//do wzoru na u(val)
  double a = 0.59 - 0.024 * s;
  double b = 0.131 + 0.063 * s;
  double N = 2.284 + 0.802 * s + 0.055 * s2;
  double A = -0.449 - 0.138 * s - 0.076 * s2;
  double B = 0.213 + 2.669 * s - 0.728 * s2;
  double C = 8.854 - 9.135 * s + 1.979 * s2;
  double D = 2.997 + 0.753 * s - 0.076 * s2;
  return N * pow (xx, a) * (1 + A * pow (xx, b) + xx * (B + C * sqrt (xx))) * pow (1 - xx, D)/xx;
}

double dv (double xx, double Q2)	//dval
{
  double s = ss (Q2);
  double s2 = s * s;

  //do wzoru na d(val)


  double a = 0.376;
  double b = 0.486 + 0.062 * s;
  double N = 0.371 + 0.083 * s + 0.039 * s2;
  double A = -0.509 + 3.31 * s - 1.248 * s2;
  double B = 12.41 - 10.52 * s + 2.267 * s2;
  double C = 6.373 - 6.208 * s + 1.418 * s2;
  double D = 3.691 + 0.799 * s - 0.071 * s2;
  return N * pow (xx, a) * (1 + A * pow (xx, b) + xx * (B + C * sqrt (xx))) * pow (1 - xx, D)/xx;
}


double del (double xx, double Q2)	//dbar-ubar
{
  double s = ss (Q2);
  double s2 = s * s;

  //do wzoru na dbar-ubar


  double a = 0.409 - 0.005 * s;
  double b = 0.799 + 0.071 * s;
  double N = 0.082 + 0.014 * s + 0.008 * s2;
  double A = -38.07 + 36.13 * s - 0.656 * s2;
  double B = 90.31 - 74.15 * s + 7.645 * s2;
//  double C = 0;
  double D = 7.486 + 1.217 * s - 0.159 * s2;
  return N * pow (xx, a) * (1 + A * pow (xx, b) + xx * B) * pow (1 - xx, D);
}



double sb (double xx, double Q2)	//s=s(bar)
{
  double s = ss (Q2);
  double s2 = s * s;

//do wzoru na s=s(bar)

  double alfa = 0.914;
  double beta = 0.577;
  double a = 1.798 - 0.596 * s;
  double A = -5.548 + 3.669 * sqrt (s) - 0.616 * s;
  double B = 18.92 - 16.73 * sqrt (s) + 5.168 * s;
  double D = 6.379 - 0.35 * s + 0.142 * s2;
  double E = 3.981 + 1.638 * s;
  double Eprim = 6.402;

  return pow (s, alfa) / pow (log (1.0 / xx), a) * (1 + A * sqrt (xx) + B * xx) * pow (1 - xx, D) *
    exp (-E + sqrt (Eprim * pow (s, beta) * log (1.0 / xx)))/xx;
}


double udb (double xx, double Q2)	//ubar+dbar
{

  double s = ss (Q2);
  double s2 = s * s;

//do ubar+dbar

  double alfa = 1.451;
  double beta = 0.271;
  double a = 0.41 - 0.232 * s;
  double b = 0.534 - 0.457 * s;
  double A = 0.89 - 0.14 * s;
  double B = -0.981;
  double C = 0.32 + 0.683 * s;
  double D = 4.752 + 1.164 * s + 0.286 * s2;
  double E = 4.119 + 1.713 * s;
  double Eprim = 0.682 + 2.978 * s;

  return (pow (xx, a) * (A + B * xx + C * xx * xx) * pow (log (1.0 / xx), b) +
	  pow (s, alfa) * exp (-E + sqrt (Eprim * pow (s, beta) * log (1.0 / xx)))) * pow (1 - xx, D);
}


double dbar (double xx, double Q2)
{
  return (del (xx, Q2) + udb (xx, Q2)) / 2/xx;
}

double
ubar (double xx, double Q2)
{
  return (-del (xx, Q2) + udb (xx, Q2)) / 2/xx;
}

////////////////////////////////////////////////////////////////////
//////////////BY correntions

double R (double xx, double Q2)
{
  double theta = 1 + (12 * Q2 / (Q2 + GeV2)) * 0.015625 / (0.015625 + xx * xx);	//Calan-Gross----> #3

  if (Q2 < 0.35 * GeV2)

    return 3.207 * GeV2 * Q2 / (Q2 * Q2 + GeV2 * GeV2) * R (xx, 0.35 * GeV2);
  else
    return 0.0635 * theta / log (Q2 / 0.04 / GeV2) + 0.5747 * GeV2 / Q2
      - 0.3534 * GeV2 * GeV2 / (Q2 * Q2 + 0.09 * GeV2 * GeV2);

}

double
R_elektrony (double xx, double Q2)
{
  double theta = 1 + (12 * Q2 / (Q2 + GeV2)) * 0.015625 / (0.015625 + xx * xx);	//Calan-Gross----> #3

  if (Q2 < 0.3 * GeV2)
    return R_elektrony (xx, 0.35 * GeV2);
  else
    return 0.0635 * theta / log (Q2 / 0.04 / GeV2) + 0.5747 * GeV2 / Q2
      - 0.3534 * GeV2 * GeV2 / (Q2 * Q2 + 0.09 * GeV2 * GeV2);

}
///////////////////////////////////////////////////////////////////////////////////////////


/////////////////////NEUTRON///////////////////////
/////////Only GRV//////////////////////////////////
double F3_cc_nu_n_GRV94 (double xx, double Q2){return 2 * (uv (xx, Q2) + sb (xx, Q2) + ubar (xx, Q2) - dbar (xx, Q2)) ;}
double F2_cc_nu_n_GRV94 (double xx, double Q2){return 2 *xx* (uv (xx, Q2) + sb (xx, Q2) + ubar (xx, Q2) + dbar (xx, Q2)) ;}
double F1_cc_nu_n_GRV94 (double xx, double Q2){return F2_cc_nu_n_GRV94 (xx, Q2) / (2*xx);}


//kwarki w neutronie
//kwark d

double F3_cc_nu_n_d_grv (double xx, double Q2){return 2 * (uv (xx, Q2) + ubar (xx, Q2));}
double F2_cc_nu_n_d_grv (double xx, double Q2){return 2 * xx * (uv (xx, Q2) + ubar (xx, Q2));}
double F1_cc_nu_n_d_grv (double xx, double Q2){return F2_cc_nu_n_d_grv (xx, Q2)/2./xx;}

//kwark s 
double F3_cc_nu_n_sb_grv (double xx, double Q2){return 2 * (sb (xx, Q2));}
double F2_cc_nu_n_sb_grv (double xx, double Q2){return 2 * xx * (sb (xx, Q2));}
double F1_cc_nu_n_sb_grv (double xx, double Q2){return F2_cc_nu_n_sb_grv (xx, Q2) /2./xx;}

//kwark anty u
double F3_cc_nu_n_ubar_grv (double xx, double Q2){return 2 * ( -dbar (xx, Q2));}				
double F2_cc_nu_n_ubar_grv (double xx, double Q2){return 2 * xx * (  dbar (xx, Q2));}
double F1_cc_nu_n_ubar_grv (double xx, double Q2){return F2_cc_nu_n_ubar_grv (xx, Q2)/2./xx;}

//Funkcje struktury z poprawkami Bodka Yanga dla CC Neutron + neutrino 
double F2_cc_nu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez) {return prop * 2 * xw * (uv (xw, freez) + sb (xw, freez) + dbar (xw, freez) + ubar (xw, freez));}
          else{    return prop * 2 * xw * (uv (xw, Q2) + sb (xw, Q2) + dbar (xw, Q2) +  ubar (xw, Q2));}
}


double F1_cc_nu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_n (xx,Q2) /2./xw;
}


double F3_cc_nu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (uv (xw, freez) + sb (xw, freez) + ubar (xw, freez) - dbar (xw, freez));}
         else{    return prop * 2 * (uv (xw, Q2) + sb (xw, Q2) + ubar (xw, Q2) - dbar (xw, Q2));}
}

//////////////////////////////////////////////////////////////////////////////////
//kwark d
double F2_cc_nu_n_d_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * (uv (xw, freez) + ubar (xw, freez));}
         else{    return prop * 2 * xw * (uv (xw, Q2) +  ubar (xw, Q2));}
}


double F1_cc_nu_n_d_BY  (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_n_d_BY (xx,Q2) /(2*xw);
}


double F3_cc_nu_n_d_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (uv (xw, freez) + ubar (xw, freez) );}
             else{return prop * 2 * (uv (xw, Q2)  + ubar (xw, Q2) );}
}

// kwakr s
double F2_cc_nu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * (sb (xw, freez));}
  else{return prop * 2 * xw * (sb (xw, Q2));}
}


double F1_cc_nu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_n_sb_BY (xx,Q2) /2./xw;

}

double F3_cc_nu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (sb (xw, freez));}
         else{    return prop * 2 * (sb (xw, Q2));}
}
//antykwak u

double F2_cc_nu_n_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){ return prop * 2 * xw * (dbar (xw, freez));}
          else{    return prop * 2 * xw * (dbar (xw, Q2));}
}


double F1_cc_nu_n_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_n_ubar_BY (xx,Q2)/2./xw;
}


double F3_cc_nu_n_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * ( -dbar (xw, freez) );}
             else{return prop * 2 * ( -dbar (xw, Q2));}
}


///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////CC Neutron + Antyneutrino//////////////////////////////////////////
double F3_cc_anu_n_GRV94 (double xx, double Q2){return 2 * (dv (xx, Q2) - sb (xx, Q2) + dbar (xx, Q2) - ubar (xx, Q2));}				
double F2_cc_anu_n_GRV94 (double xx, double Q2){return 2 *xx* (dv (xx, Q2) + sb (xx, Q2) + ubar (xx, Q2) + dbar (xx, Q2));}
double F1_cc_anu_n_GRV94 (double xx, double Q2){return F2_cc_anu_n_GRV94 (xx, Q2)/2./xx;}
////////////kwarki
//kwark u

double F3_cc_anu_n_u_grv (double xx, double Q2){return 2 * (dv (xx, Q2) + dbar (xx, Q2));}				
double F2_cc_anu_n_u_grv (double xx, double Q2){return 2 * xx* (dv (xx, Q2) + dbar (xx, Q2));}
double F1_cc_anu_n_u_grv (double xx, double Q2){return F2_cc_anu_n_u_grv (xx, Q2)/2./xx;}

//antykwakt s
double F3_cc_anu_n_sb_grv (double xx, double Q2){return 2 * (-sb (xx, Q2));}
double F2_cc_anu_n_sb_grv (double xx, double Q2){return 2 * xx * ( sb (xx, Q2));}
double F1_cc_anu_n_sb_grv (double xx, double Q2){return F2_cc_anu_n_sb_grv (xx, Q2)/2./xx;}
//!!!!!!!!!!!
//antykwark d
double F3_cc_anu_n_dbar_grv (double xx, double Q2){return 2 * (-ubar (xx, Q2));}
double F2_cc_anu_n_dbar_grv (double xx, double Q2){return 2 * xx * ( ubar (xx, Q2));}
double F1_cc_anu_n_dbar_grv (double xx, double Q2){return F2_cc_anu_n_dbar_grv (xx, Q2)/2./xx;}


///////////BY correctiobns 

double F2_cc_anu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (dv (xw, freez) + sb (xw, freez) + ubar (xw, freez) + dbar (xw, freez));}
  else{return prop * 2 * xw * (dv (xw, Q2) + sb (xw, Q2) + ubar (xw, Q2) + dbar (xw, Q2));}
}

double F3_cc_anu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (dv (xw, freez) + sb (xw, freez) + dbar (xw, freez) - ubar (xw, freez));}
  else{return prop * 2 * (dv (xw, Q2) + sb (xw, Q2) + dbar (xw, Q2) - ubar (xw, Q2));}
}

double F1_cc_anu_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n (xx, Q2)/2./xw;
}





double F2_cc_anu_n_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (dv (xw, freez) + sb (xw, freez) + ubar (xw, freez) + dbar (xw, freez));}
             else{return prop * 2 * xw * (dv (xw, Q2)    + sb (xw, Q2)    + ubar (xw, Q2)    + dbar (xw, Q2));}
}

double F3_cc_anu_n_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (dv (xw, freez) - sb (xw, freez) + dbar (xw, freez) - ubar (xw, freez));}
             else{return prop * 2 * (dv (xw, Q2)    - sb (xw, Q2) + dbar (xw, Q2) - ubar (xw, Q2));}
}

double F1_cc_anu_n_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_BY (xx, Q2)/2./xw;
}


//kwark u
double F2_cc_anu_n_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  if (Q2 < freez){return prop * 2 * xw * (dv (xw, freez) + dbar (xw, freez));}
             else{return prop * 2 * xw * (dv (xw, Q2)    + dbar (xw, Q2));}
}


double F1_cc_anu_n_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
//  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_u_BY (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_u_BY (xx, Q2)/2/xw;
}

double F3_cc_anu_n_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (dv (xw, freez) + dbar (xw, freez) ) ;}
             else{return prop * 2 * (dv (xw, Q2)    + dbar (xw, Q2)) ;}
}

//antykwar s

double F2_cc_anu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (sb (xw, freez)) ;}
             else{return prop * 2 * xw * (sb (xw, Q2)) ;}
}


double F1_cc_anu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
//  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_sb_BY (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_sb_BY (xx, Q2) /2/xw;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double F3_cc_anu_n_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (-sb (xw, freez)) ;}
             else{return prop * 2 * (-sb (xw, Q2));}
}

//antykwar d
double F2_cc_anu_n_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  if (Q2 < freez){return prop * 2 * xw * (ubar (xw, freez));}
             else{return prop * 2 * xw * (ubar (xw, Q2) );}
}


double F1_cc_anu_n_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_n_dbar_BY (xx, Q2)/2./xw;
}

double F3_cc_anu_n_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * (- ubar (xw, freez));}
             else{return prop * 2 * (- ubar (xw, Q2));}
}






//////////////////////////////////////////////////////////////////////////////
//d gestosc kwarku dval  w protonie i uval w neutronie 
//dbar gestosc kwarku dsea i dsea_bar w protonie i usea i usea_bar w neutronie
//sb gestosc s lub sbar w protoni lub neutronie


/////////////////////////PROTON////////////////////////////////////////////////////////
/////////////////NEUTRINO///////////////////////////////////////////////////////////

double F3_cc_nu_p_GRV94 (double xx, double Q2){  return 2 * (dv (xx, Q2) + sb (xx, Q2) + dbar (xx, Q2) - ubar (xx, Q2));}
double F2_cc_nu_p_GRV94 (double xx, double Q2){  return 2 * xx* (dv (xx, Q2) + sb (xx, Q2) + dbar (xx, Q2) + ubar (xx, Q2));}
double F1_cc_nu_p_GRV94 (double xx, double Q2){  return F2_cc_nu_p_GRV94 (xx, Q2)/2./xx;}


/////////////////////////////


// Structure functions for CC neutrino interaction on proton with Bodek Yang corrections

double F2_cc_nu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (dv (xw, freez) + sb (xw, freez) + dbar (xw, freez) + ubar (xw, freez));}
             else{return prop * 2 * xw * (dv (xw, Q2)    + sb (xw, Q2)    + dbar (xw, Q2)    + ubar (xw, Q2)) ;}
}

double F1_cc_nu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_p (xx,Q2)/2./xw;
}

double F3_cc_nu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);
  
  if (Q2 < freez){return prop * 2 * (dv (xw, freez) + sb (xw, freez) + dbar (xw, freez) -  ubar (xw, freez)) ;}
             else{return prop * 2 * (dv (xw, Q2) + sb (xw, Q2) + dbar (xw, Q2) - ubar (xw, Q2));}
}


///////////////////////////////////////////////////////////////////
//Funkcje struktury dla kwarkow

//kwark d

double F3_cc_nu_p_d_grv (double xx, double Q2){  return 2 * (dv (xx, Q2) + dbar (xx, Q2) );}
double F2_cc_nu_p_d_grv (double xx, double Q2){  return 2 * xx * (dv (xx, Q2) + dbar (xx, Q2) );}
//double F1_cc_nu_p_d_grv (double xx, double Q2){  return F2_cc_nu_p_d_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_d_grv (double xx, double Q2){  return F2_cc_nu_p_d_grv (xx, Q2)/2./xx;}

//kwark sb

double F3_cc_nu_p_sb_grv (double xx, double Q2){  return 2 * (sb (xx, Q2)) ;}
double F2_cc_nu_p_sb_grv (double xx, double Q2){  return 2 * xx * (sb (xx, Q2));}
//double F1_cc_nu_p_sb_grv (double xx, double Q2){  return F2_cc_nu_p_sb_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_sb_grv (double xx, double Q2){  return F2_cc_nu_p_sb_grv (xx, Q2)/2./xx;}

//antykwark u 

double F3_cc_nu_p_ubar_grv (double xx, double Q2){  return 2 * (-ubar (xx, Q2)) ;}
double F2_cc_nu_p_ubar_grv (double xx, double Q2){  return 2 *xx* ( ubar (xx, Q2));}
//double F1_cc_nu_p_ubar_grv (double xx, double Q2){  return F2_cc_nu_p_ubar_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_ubar_grv (double xx, double Q2){  return F2_cc_nu_p_ubar_grv (xx, Q2)/2./xx;}

//poprawki bodka

//kwakr d
double F2_cc_nu_p_d_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (dv (xw, freez) + dbar (xw, freez) );}
             else{return prop * 2 * xw * (dv (xw, Q2) + dbar (xw, Q2));}
}

double F1_cc_nu_p_d_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_p_d_BY (xx,Q2)/2./xw;
}

double F3_cc_nu_p_d_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);
  
  if (Q2 < freez){return prop * 2 * (dv (xw, freez) + dbar (xw, freez)) ;}
             else{return prop * 2 * (dv (xw, Q2) + dbar (xw, Q2)) ;}
}


//kwark sb

double F2_cc_nu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * (sb (xw, freez) );}
             else{return prop * 2 * xw * (sb (xw, Q2) );}
}

double F1_cc_nu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_p_sb_BY (xx,Q2)/2./xw;
}

double F3_cc_nu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);
  
  if (Q2 < freez){return prop * 2 * ( sb (xw, freez) );}
             else{return prop * 2 * ( sb (xw, Q2) );}
}

//antykwark u

double F2_cc_nu_p_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez){return prop * 2 * xw * ( ubar (xw, freez));}
             else{return prop * 2 * xw * ( ubar (xw, Q2));}
}

double F1_cc_nu_p_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
//  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_p_ubar_BY (xx,Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_nu_p_ubar_BY (xx,Q2)/2./xw;
}

double F3_cc_nu_p_ubar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);
  
  if (Q2 < freez){return prop * 2 * ( - ubar (xw, freez));}
             else{return prop * 2 * ( - ubar (xw, Q2));}
}




///////////////////////ANTYNEUTRINO//////////////////////////////////////////////
double F3_cc_anu_p_GRV94 (double xx, double Q2){return 2 * (uv (xx, Q2) - sb (xx, Q2) + ubar (xx, Q2) - dbar (xx, Q2));}
double F2_cc_anu_p_GRV94 (double xx, double Q2){return 2 * xx * (uv (xx, Q2) + sb (xx, Q2) + ubar (xx, Q2) + dbar (xx, Q2));}
//double F1_cc_anu_p_GRV94 (double xx, double Q2){return F2_cc_anu_p_GRV94 (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_GRV94 (double xx, double Q2){return F2_cc_anu_p_GRV94 (xx, Q2)/2./xx;}


double F2_cc_anu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * (uv (xw, freez) + sb (xw, freez) + ubar (xw, freez) + dbar (xw, freez));}
             else{return prop * 2 * xw * (uv (xw, Q2) + sb (xw, Q2) + ubar (xw, Q2) + dbar (xw, Q2));}
}

double F1_cc_anu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_p (xx, Q2)/2./xw;
}

double F3_cc_anu_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * (uv (xw, freez) - sb (xw, freez) + ubar (xw, freez) - dbar (xw, freez));}
             else{return prop * 2 * (uv (xw, Q2)    - sb (xw, Q2)    + ubar (xw, Q2)    - dbar (xw, Q2));}
}

////////////kwarki
//!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!

//kwakr u
double F3_cc_anu_p_u_grv (double xx, double Q2){  return 2 * (uv (xx, Q2) + ubar (xx, Q2) );}
double F2_cc_anu_p_u_grv (double xx, double Q2){  return 2 * xx * (uv (xx, Q2) + ubar (xx, Q2) );}
//double F1_cc_anu_p_u_grv (double xx, double Q2){  return F2_cc_anu_p_u_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_u_grv (double xx, double Q2){  return F2_cc_anu_p_u_grv (xx, Q2)/2./xx;}

//kwark sb
double F3_cc_anu_p_sb_grv (double xx, double Q2){  return 2 * (-sb (xx, Q2) );}
double F2_cc_anu_p_sb_grv (double xx, double Q2){  return 2 * xx * ( sb (xx, Q2) );}
//double F1_cc_anu_p_sb_grv (double xx, double Q2){  return F2_cc_anu_p_sb_grv (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_sb_grv (double xx, double Q2){  return F2_cc_anu_p_sb_grv (xx, Q2)/2./xx;}

//antykwakr d
double F3_cc_anu_p_dbar_grv (double xx, double Q2){  return 2 * (-dbar (xx, Q2));}
double F2_cc_anu_p_dbar_grv (double xx, double Q2){  return 2 * xx * ( dbar (xx, Q2));}
//double F1_cc_anu_p_dbar_grv (double xx, double Q2){  return F2_cc_anu_p_dbar_grv (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_dbar_grv (double xx, double Q2){  return F2_cc_anu_p_dbar_grv (xx, Q2) /2./xx;}


////////kwarki z poprawkami bodka

//kwark u

double F2_cc_anu_p_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * (uv (xw, freez) + ubar (xw, freez));}
             else{return prop * 2 * xw * (uv (xw, Q2)    + ubar (xw, Q2) );}
}

double F1_cc_anu_p_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_p_u_BY (xx, Q2)/2./xw;
}

double F3_cc_anu_p_u_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * (uv (xw, freez) + ubar (xw, freez) ) ;}
             else{return prop * 2 * (uv (xw, Q2) + ubar (xw, Q2) ) ;}
}

//kwark sb

double F2_cc_anu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * ( sb (xw, freez)) ;}
         else{    return prop * 2 * xw * ( sb (xw, Q2)) ;}
}

double F1_cc_anu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_p_sb_BY (xx, Q2)/2./xw;
}

double F3_cc_anu_p_sb_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * ( -sb (xw, freez));}
             else{return prop * 2 * ( -sb (xw, Q2) ) ;}
}

//antykwark d
double F2_cc_anu_p_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * xw * ( dbar (xw, freez));}
             else{return prop * 2 * xw * ( dbar (xw, Q2));}
}

double F1_cc_anu_p_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
//  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_p_dbar_BY (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_cc_anu_p_dbar_BY (xx, Q2)/2./xw;
}

double F3_cc_anu_p_dbar_BY (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  
  if (Q2 < freez){return prop * 2 * ( - dbar (xw, freez));}
             else{return prop * 2 * ( - dbar (xw, Q2));}
}




//////////////////////////////////////////////////////////////////////////////////////////////////


const double g_V_u = 1 / 2. - 4 / 3. * sin_2_theta_W;
const double g_V_d = 2 / 3. * sin_2_theta_W - 1 / 2.;
const double g_A_u = 1 / 2.;
const double g_A_d = -1 / 2.;

//funkcje struktury dla NC by JN
//proton

double F3_nc_p_GRV94 (double xx, double Q2){return 2 * (g_V_u * g_A_u * uv (xx, Q2) + g_V_d * g_A_d * dv (xx, Q2));}

double F2_nc_p_GRV94 (double xx, double Q2)
{
  return xx*((g_V_u * g_V_u + g_A_u * g_A_u) * (uv (xx, Q2) + 2 * ubar (xx, Q2)) + (g_V_d * g_V_d + g_A_d * g_A_d) * (dv (xx, Q2) + 2 * dbar (xx, Q2) + 2 * sb (xx, Q2)));
}

double F1_nc_p_GRV94 (double xx, double Q2){return F2_nc_p_GRV94 (xx, Q2) /2./xx;}

///POPRAWKI BODKA WEDLUG POWYZSZEGO SCHEMATU/////////////////////////////////////////////////////////////////////

double F3_nc_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez) {return 2 * prop * (g_V_u * g_A_u * uv (xw, freez) + g_V_d * g_A_d * dv (xw, freez));}
           else   {return 2 * prop * (g_V_u * g_A_u * uv (xw, Q2)    + g_V_d * g_A_d * dv (xw, Q2));}
}

double F2_nc_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  if (Q2 < freez)
    {
      return prop * xw * (g_V_u * g_V_u + g_A_u * g_A_u) * (uv (xw, freez) + 2 * ubar (xw, freez))
           + (g_V_d * g_V_d + g_A_d * g_A_d) * (dv (xw, Q2) + 2 * dbar (xw, Q2) + 2 * sb (xw, freez));
    }

  else
    {
      return xw * prop * (g_V_u * g_V_u + g_A_u * g_A_u) * (uv (xw, Q2) + 2 * ubar (xw, Q2))
	+ (g_V_d * g_V_d + g_A_d * g_A_d) * (dv (xw, Q2) + 2 * dbar (xw, Q2) + 2 * sb (xw, Q2));
    }
}

double F1_nc_p (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_nc_p (xx, Q2) /2./xw;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//kwark d
double F3_nc_p_d (double xx, double Q2)
{
  return 2 * (g_V_d * g_A_d * dv (xx, Q2));
}

double F2_nc_p_d (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (dv (xx, Q2) + dbar (xx, Q2));
}

double
F1_nc_p_d (double xx, double Q2)
{
  return F2_nc_p_d (xx, Q2) /2./xx;
}

//kwark u
double F3_nc_p_u (double xx, double Q2)
{
  return 2 * (g_V_u * g_A_u * uv (xx, Q2)) ;
}

double F2_nc_p_u (double xx, double Q2)
{
  return xx*(g_V_u * g_V_u + g_A_u * g_A_u) * (uv (xx, Q2) + 2 * ubar (xx, Q2));
}

double F1_nc_p_u (double xx, double Q2)
{
  return F2_nc_p_u (xx, Q2) /2./xx;
}

//antykwark d
double F3_nc_p_dbar (double xx, double Q2)
{
  return 0;
}

double F2_nc_p_dbar (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (dbar (xx, Q2));
}

double F1_nc_p_dbar (double xx, double Q2)
{
  return F2_nc_p_dbar (xx, Q2) /2./xx;
}

//antykwark u 
double F3_nc_p_ubar (double xx, double Q2)
{
  return 0;
}

double F2_nc_p_ubar (double xx, double Q2)
{
  return xx*(g_V_u * g_V_u + g_A_u * g_A_u) * (ubar (xx, Q2));
}

double F1_nc_p_ubar (double xx, double Q2)
{
  return F2_nc_p_ubar (xx, Q2) /2./xx;
}

//kwark s 
double F3_nc_p_s (double xx, double Q2)
{
  return 0;
}

double F2_nc_p_s (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (sb (xx, Q2));
}

double F1_nc_p_s (double xx, double Q2)
{
  return F2_nc_p_s (xx, Q2)/2./xx;
}

//antykwark s 
double F3_nc_p_sbar (double xx, double Q2)
{
  return 0;
}

double F2_nc_p_sbar (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (sb (xx, Q2));
}

double F1_nc_p_sbar (double xx, double Q2)
{
  return F2_nc_p_sbar (xx, Q2) /2./xx;
}



//neutron
double F3_nc_n_GRV94 (double xx, double Q2)
{
  return 2 * (g_V_u * g_A_u * dv (xx, Q2) + g_V_d * g_A_d * uv (xx, Q2));
}

double F2_nc_n_GRV94 (double xx, double Q2)
{
  return xx*((g_V_u * g_V_u + g_A_u * g_A_u) * (dv (xx, Q2) + 2 * dbar (xx, Q2))+ (g_V_d * g_V_d + g_A_d * g_A_d) * (uv (xx, Q2) + 2 * ubar (xx, Q2) + 2 * sb (xx, Q2)));
}

double F1_nc_n_GRV94 (double xx, double Q2)
{
  return F2_nc_n_GRV94 (xx, Q2) /2./xx;
}

///POPRAWKI BODKA WEDLUG POWYZSZEGO SCHEMATU/////////////////////////////////////////////////////////////////////

double F3_nc_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  double prop = Q2 / (Q2 + 0.188 * GeV2);	

  if (Q2 < freez)
    {
      return 2 * prop * (g_V_u * g_A_u * dv (xw, freez) + g_V_d * g_A_d * uv (xw, freez));
    }
  else
    {
      return 2 * prop * (g_V_u * g_A_u * dv (xw, Q2) + g_V_d * g_A_d * uv (xw, Q2));
    }
}

double F2_nc_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	

  double prop = Q2 / (Q2 + 0.188 * GeV2);	
  if (Q2 < freez)
    {
      return prop * xw * (g_V_u * g_V_u + g_A_u * g_A_u) * (dv (xw, freez) + 2 * dbar (xw, freez))
	+ (g_V_d * g_V_d + g_A_d * g_A_d) * (uv (xw, freez) + 2 * ubar (xw, freez) + 2 * sb (xw, freez));
    }
  else
    {
      return prop * xw * (g_V_u * g_V_u + g_A_u * g_A_u) * (dv (xw, Q2) + 2 * dbar (xw, Q2))
	+ (g_V_d * g_V_d + g_A_d * g_A_d) * (uv (xw, Q2) + 2 * ubar (xw, Q2) + 2 * sb (xw, Q2));
    }
}


double F1_nc_n (double xx, double Q2)
{
  double xw = xx * (Q2 + 0.624 * GeV2) / (Q2 + 1.735 * GeV2 * xx);	
  return ((1 + 4 * M2 * xx * xx / Q2) / (1 + R (xx, Q2))) * F2_nc_n (xx, Q2) /2./xw;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//kwark d
double F3_nc_n_d (double xx, double Q2)
{
  return 2 * (g_V_d * g_A_d * uv (xx, Q2));
}

double F2_nc_n_d (double xx, double Q2)
{
  return xx*((g_V_d * g_V_d + g_A_d * g_A_d) * (uv (xx, Q2) + ubar (xx, Q2)));
}

double F1_nc_n_d (double xx, double Q2)
{
  return F2_nc_n_d (xx, Q2) /2./xx;
}

//kwark u
double F3_nc_n_u (double xx, double Q2)
{
  return 2 * (g_V_u * g_A_u * dv (xx, Q2)) ;
}

double F2_nc_n_u (double xx, double Q2)
{
  return xx*((g_V_u * g_V_u + g_A_u * g_A_u) * (dv (xx, Q2) + dbar (xx, Q2)));
}

double F1_nc_n_u (double xx, double Q2)
{
  return F2_nc_n_u (xx, Q2) /2./xx;
}

//antykwark d
double F3_nc_n_dbar (double xx, double Q2)
{
  return 0;
}

double F2_nc_n_dbar (double xx, double Q2)
{
  return xx*((g_V_d * g_V_d + g_A_d * g_A_d) * (ubar (xx, Q2)));
}

double F1_nc_n_dbar (double xx, double Q2)
{
  return F2_nc_n_dbar (xx, Q2) /2./xx;
}

//antykwark u
double F3_nc_n_ubar (double xx, double Q2)
{
  return 0;
}

double F2_nc_n_ubar (double xx, double Q2)
{
  return xx*(g_V_u * g_V_u + g_A_u * g_A_u) * (dbar (xx, Q2));
}

double F1_nc_n_ubar (double xx, double Q2)
{
  return F2_nc_n_ubar (xx, Q2) /2./xx;
}

//kwark s
double F3_nc_n_s (double xx, double Q2)
{
  return 0;
}

double F2_nc_n_s (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (sb (xx, Q2));
}

double F1_nc_n_s (double xx, double Q2)
{
  return F2_nc_n_s (xx, Q2) /2./xx;
}

//antykwark s
double F3_nc_n_sbar (double xx, double Q2)
{
  return 0;
}

double F2_nc_n_sbar (double xx, double Q2)
{
  return xx*(g_V_d * g_V_d + g_A_d * g_A_d) * (sb (xx, Q2));
}

double F1_nc_n_sbar (double xx, double Q2)
{
  return F2_nc_n_sbar (xx, Q2) /2./xx;
}


