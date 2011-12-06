#ifndef  _grv94_BODEK_h_
#define _grv94_BODEK_h_
//Oryginalnie porprawki bodka wprowadzone przez K.Graczyka
//Zmiany przez JN to: 
//#include "jednostki.h"
//#include "parameters.h"
//double GeV=1000;
//double GeV2=GeV*GeV;


//const double freez = 0.24 * GeV2;	// Granica zamrozenia funkcji strukturalnych

double ss (double Q2);

double uv (double xx, double Q2);

double dv (double xx, double Q2);

double del (double xx, double Q2);


double sb (double xx, double Q2);

double udb (double xx, double Q2);

double dbar (double xx, double Q2);

double ubar (double xx, double Q2);
////////////////////////////////////////////////////////////////////
//////////////BY correntions

double R (double xx, double Q2);
double R_elektrony (double xx, double Q2);
/////////////////////////////////////////////////////


/////////////////////NEUTRON///////////////////////
/////////Only GRV//////////////////////////////////
double F3_cc_nu_n_GRV94 (double xx, double Q2);
double F2_cc_nu_n_GRV94 (double xx, double Q2);
double F1_cc_nu_n_GRV94 (double xx, double Q2);

//kwarki w neutronie
//kwark d

double F3_cc_nu_n_d_grv (double xx, double Q2);

double F2_cc_nu_n_d_grv (double xx, double Q2);

double F1_cc_nu_n_d_grv (double xx, double Q2);

//kwark s 
double F3_cc_nu_n_sb_grv (double xx, double Q2);
double F2_cc_nu_n_sb_grv (double xx, double Q2);
double F1_cc_nu_n_sb_grv (double xx, double Q2);

//kwark anty u
double F3_cc_nu_n_ubar_grv (double xx, double Q2);
double F2_cc_nu_n_ubar_grv (double xx, double Q2);
double F1_cc_nu_n_ubar_grv (double xx, double Q2);

//Funkcje struktury z poprawkami Bodka Yanga dla CC Neutron + neutrino 
double F2_cc_nu_n (double xx, double Q2);

double F1_cc_nu_n (double xx, double Q2);

double F3_cc_nu_n (double xx, double Q2);
//////////////////////////////////////////////////////////////////////////////////
//kwark d
double F2_cc_nu_n_d_BY (double xx, double Q2);

double F1_cc_nu_n_d_BY  (double xx, double Q2);

double F3_cc_nu_n_d_BY (double xx, double Q2);
// kwakr s
double F2_cc_nu_n_sb_BY (double xx, double Q2);

double F1_cc_nu_n_sb_BY (double xx, double Q2);
double F3_cc_nu_n_sb_BY (double xx, double Q2);

double F2_cc_nu_n_ubar_BY (double xx, double Q2);

double F1_cc_nu_n_ubar_BY (double xx, double Q2);

double F3_cc_nu_n_ubar_BY (double xx, double Q2);

///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////CC Neutron + Antyneutrino//////////////////////////////////////////
double F3_cc_anu_n_GRV94 (double xx, double Q2);
double F2_cc_anu_n_GRV94 (double xx, double Q2);
double F1_cc_anu_n_GRV94 (double xx, double Q2);
////////////kwarki
//kwark u

double F3_cc_anu_n_u_grv (double xx, double Q2);
double F2_cc_anu_n_u_grv (double xx, double Q2);
double F1_cc_anu_n_u_grv (double xx, double Q2);

//antykwakt s
double F3_cc_anu_n_sb_grv (double xx, double Q2);
double F2_cc_anu_n_sb_grv (double xx, double Q2);
double F1_cc_anu_n_sb_grv (double xx, double Q2);
//!!!!!!!!!!!
//antykwark d
double F3_cc_anu_n_dbar_grv (double xx, double Q2);
double F2_cc_anu_n_dbar_grv (double xx, double Q2);
double F1_cc_anu_n_dbar_grv (double xx, double Q2);


///////////BY correctiobns 

double F2_cc_anu_n (double xx, double Q2);
double F3_cc_anu_n (double xx, double Q2);
double F1_cc_anu_n (double xx, double Q2);




double F2_cc_anu_n_BY (double xx, double Q2);
double F3_cc_anu_n_BY (double xx, double Q2);

double F1_cc_anu_n_BY (double xx, double Q2);


//kwark u
double F2_cc_anu_n_u_BY (double xx, double Q2);


double F1_cc_anu_n_u_BY (double xx, double Q2);
double F3_cc_anu_n_u_BY (double xx, double Q2);

//antykwar s

double F2_cc_anu_n_sb_BY (double xx, double Q2);


double F1_cc_anu_n_sb_BY (double xx, double Q2);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double F3_cc_anu_n_sb_BY (double xx, double Q2);
//antykwar d
double F2_cc_anu_n_dbar_BY (double xx, double Q2);
double F1_cc_anu_n_dbar_BY (double xx, double Q2);
double F3_cc_anu_n_dbar_BY (double xx, double Q2);





//////////////////////////////////////////////////////////////////////////////
//d gestosc kwarku dval  w protonie i uval w neutronie 
//dbar gestosc kwarku dsea i dsea_bar w protonie i usea i usea_bar w neutronie
//sb gestosc s lub sbar w protoni lub neutronie


/////////////////////////PROTON////////////////////////////////////////////////////////
/////////////////NEUTRINO///////////////////////////////////////////////////////////

double F3_cc_nu_p_GRV94 (double xx, double Q2);
double F2_cc_nu_p_GRV94 (double xx, double Q2);
double F1_cc_nu_p_GRV94 (double xx, double Q2);


/////////////////////////////


// Structure functions for CC neutrino interaction on proton with Bodek Yang corrections

double F2_cc_nu_p (double xx, double Q2);
double F1_cc_nu_p (double xx, double Q2);
double F3_cc_nu_p (double xx, double Q2);

///////////////////////////////////////////////////////////////////
//Funkcje struktury dla kwarkow

//kwark d

double F3_cc_nu_p_d_grv (double xx, double Q2);
double F2_cc_nu_p_d_grv (double xx, double Q2);
//double F1_cc_nu_p_d_grv (double xx, double Q2){  return F2_cc_nu_p_d_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_d_grv (double xx, double Q2);

//kwark sb

double F3_cc_nu_p_sb_grv (double xx, double Q2);
double F2_cc_nu_p_sb_grv (double xx, double Q2);
//double F1_cc_nu_p_sb_grv (double xx, double Q2){  return F2_cc_nu_p_sb_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_sb_grv (double xx, double Q2);
//antykwark u 

double F3_cc_nu_p_ubar_grv (double xx, double Q2);
double F2_cc_nu_p_ubar_grv (double xx, double Q2);
//double F1_cc_nu_p_ubar_grv (double xx, double Q2){  return F2_cc_nu_p_ubar_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_nu_p_ubar_grv (double xx, double Q2);

//poprawki bodka

//kwakr d
double F2_cc_nu_p_d_BY (double xx, double Q2);
double F1_cc_nu_p_d_BY (double xx, double Q2);
double F3_cc_nu_p_d_BY (double xx, double Q2);

//kwark sb

double F2_cc_nu_p_sb_BY (double xx, double Q2);
double F1_cc_nu_p_sb_BY (double xx, double Q2);
double F3_cc_nu_p_sb_BY (double xx, double Q2);
//antykwark u

double F2_cc_nu_p_ubar_BY (double xx, double Q2);
double F1_cc_nu_p_ubar_BY (double xx, double Q2);
double F3_cc_nu_p_ubar_BY (double xx, double Q2);


///////////////////////ANTYNEUTRINO//////////////////////////////////////////////
double F3_cc_anu_p_GRV94 (double xx, double Q2);
double F2_cc_anu_p_GRV94 (double xx, double Q2);
//double F1_cc_anu_p_GRV94 (double xx, double Q2){return F2_cc_anu_p_GRV94 (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_GRV94 (double xx, double Q2);


double F2_cc_anu_p (double xx, double Q2);
double F1_cc_anu_p (double xx, double Q2);
double F3_cc_anu_p (double xx, double Q2);
////////////kwarki
//!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!

//kwakr u
double F3_cc_anu_p_u_grv (double xx, double Q2);
double F2_cc_anu_p_u_grv (double xx, double Q2);
//double F1_cc_anu_p_u_grv (double xx, double Q2){  return F2_cc_anu_p_u_grv (xx, Q2)*(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_u_grv (double xx, double Q2);

//kwark sb
double F3_cc_anu_p_sb_grv (double xx, double Q2);
double F2_cc_anu_p_sb_grv (double xx, double Q2);
//double F1_cc_anu_p_sb_grv (double xx, double Q2){  return F2_cc_anu_p_sb_grv (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_sb_grv (double xx, double Q2);

//antykwakr d
double F3_cc_anu_p_dbar_grv (double xx, double Q2);
double F2_cc_anu_p_dbar_grv (double xx, double Q2);
//double F1_cc_anu_p_dbar_grv (double xx, double Q2){  return F2_cc_anu_p_dbar_grv (xx, Q2) *(xx*2*M12*M12*M12/Q2+1./2/xx);}
double F1_cc_anu_p_dbar_grv (double xx, double Q2);


////////kwarki z poprawkami bodka

//kwark u

double F2_cc_anu_p_u_BY (double xx, double Q2);
double F1_cc_anu_p_u_BY (double xx, double Q2);
double F3_cc_anu_p_u_BY (double xx, double Q2);
//kwark sb

double F2_cc_anu_p_sb_BY (double xx, double Q2);
double F1_cc_anu_p_sb_BY (double xx, double Q2);
double F3_cc_anu_p_sb_BY (double xx, double Q2);

//antykwark d
double F2_cc_anu_p_dbar_BY (double xx, double Q2);
double F1_cc_anu_p_dbar_BY (double xx, double Q2);
double F3_cc_anu_p_dbar_BY (double xx, double Q2);



//////////////////////////////////////////////////////////////////////////////////////////////////


//const double g_V_u = 1 / 2. - 4 / 3. * sin_2_theta_W;
//const double g_V_d = 2 / 3. * sin_2_theta_W - 1 / 2.;
//const double g_A_u = 1 / 2.;
//const double g_A_d = -1 / 2.;

//funkcje struktury dla NC by JN
//proton

double F3_nc_p_GRV94 (double xx, double Q2);
double F2_nc_p_GRV94 (double xx, double Q2);
double F1_nc_p_GRV94 (double xx, double Q2);

///POPRAWKI BODKA WEDLUG POWYZSZEGO SCHEMATU/////////////////////////////////////////////////////////////////////

double F3_nc_p (double xx, double Q2);
double F2_nc_p (double xx, double Q2);
double F1_nc_p (double xx, double Q2);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//kwark d
double F3_nc_p_d (double xx, double Q2);
double F2_nc_p_d (double xx, double Q2);
double F1_nc_p_d (double xx, double Q2);
//kwark u
double F3_nc_p_u (double xx, double Q2);
double F2_nc_p_u (double xx, double Q2);
double F1_nc_p_u (double xx, double Q2);

//antykwark d
double F3_nc_p_dbar (double xx, double Q2);
double F2_nc_p_dbar (double xx, double Q2);
double F1_nc_p_dbar (double xx, double Q2);
//antykwark u 
double F3_nc_p_ubar (double xx, double Q2);
double F2_nc_p_ubar (double xx, double Q2);
double F1_nc_p_ubar (double xx, double Q2);
//kwark s 
double F3_nc_p_s (double xx, double Q2);
double F2_nc_p_s (double xx, double Q2);
double F1_nc_p_s (double xx, double Q2);
//antykwark s 
double F3_nc_p_sbar (double xx, double Q2);
double F2_nc_p_sbar (double xx, double Q2);
double F1_nc_p_sbar (double xx, double Q2);


//neutron
double F3_nc_n_GRV94 (double xx, double Q2);
double F2_nc_n_GRV94 (double xx, double Q2);
double F1_nc_n_GRV94 (double xx, double Q2);
///POPRAWKI BODKA WEDLUG POWYZSZEGO SCHEMATU/////////////////////////////////////////////////////////////////////

double F3_nc_n (double xx, double Q2);
double F2_nc_n (double xx, double Q2);
double F1_nc_n (double xx, double Q2);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//kwark d
double F3_nc_n_d (double xx, double Q2);
double F2_nc_n_d (double xx, double Q2);
double F1_nc_n_d (double xx, double Q2);
//kwark u
double F3_nc_n_u (double xx, double Q2);
double F2_nc_n_u (double xx, double Q2);
double F1_nc_n_u (double xx, double Q2);
//antykwark d
double F3_nc_n_dbar (double xx, double Q2);
double F2_nc_n_dbar (double xx, double Q2);
double F1_nc_n_dbar (double xx, double Q2);
//antykwark u
double F3_nc_n_ubar (double xx, double Q2);
double F2_nc_n_ubar (double xx, double Q2);
double F1_nc_n_ubar (double xx, double Q2);
//kwark s
double F3_nc_n_s (double xx, double Q2);
double F2_nc_n_s (double xx, double Q2);
double F1_nc_n_s (double xx, double Q2);
//antykwark s
double F3_nc_n_sbar (double xx, double Q2);
double F2_nc_n_sbar (double xx, double Q2);
double F1_nc_n_sbar (double xx, double Q2);


#endif
