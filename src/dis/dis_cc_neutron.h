#ifndef _dis_cc_neutron_h_
#define _dis_cc_neutron_h_


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////NEUTRON//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double sigma_xy_cc_nu_n_Paschos(double E, double xx, double yy, double m );//troche inny wzor z pracy Paschosa

double cr_sec_cc_nu_n (double E, double W, double nu, double m);


double sigma_xy_cc_nu_n_grv(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_grv (double E, double W, double nu, double m);


double sigma_xy_cc_nu_n_grv(double E, double xx, double yy, double m,int ff);

double cr_sec_cc_nu_n_grv (double E, double W, double nu, double m, int ff);





/////////////////kwarki za pomoca funkcji struktury 

//kwark d
double sigma_xy_cc_nu_n_d_grv(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_d_grv (double E, double W, double nu, double m);
//kwark s

double sigma_xy_cc_nu_n_sb_grv(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_sb_grv (double E, double W, double nu, double m);

//antykwark u

double sigma_xy_cc_nu_n_ubar_grv(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_ubar_grv (double E, double W, double nu, double m);


///kwarki  poprawkami Bodka


double sigma_xy_cc_nu_n_BY(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_BY (double E, double W, double nu, double m);


//kwark d

double sigma_xy_cc_nu_n_d_BY(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_d_BY (double E, double W, double nu, double m);
//kwark s

double sigma_xy_cc_nu_n_sb_BY(double E, double xx, double yy, double m );

double cr_sec_cc_nu_n_sb_BY (double E, double W, double nu, double m);

//antykwark u

double sigma_xy_cc_nu_n_ubar_BY(double E, double xx, double yy, double m);

double cr_sec_cc_nu_n_ubar_BY (double E, double W, double nu, double m);



//krarki za pomoca pdf
//quark d valence + d sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
double sigma_xy_cc_nu_n_d(double E, double xx, double yy, double m);

double cr_sec_cc_nu_n_d (double E, double W, double nu, double m);

//quark s + sbar
double sigma_xy_cc_nu_n_sb(double E, double xx, double yy, double m);

double cr_sec_cc_nu_n_sb (double E, double W, double nu, double m);


//dbar
//in common convetion ubar_neutron = dbar(_proton)
double sigma_xy_cc_nu_n_ubar(double E, double xx, double yy, double m);

double cr_sec_cc_nu_n_ubar (double E, double W, double nu, double m);


// nu_bar + n --> mu_bar +X (CC, lepton)
// rozklad kwarkow sprawdzony 26.06.05

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double sigma_xy_cc_anu_n_Paschos(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n (double E, double W, double nu, double m);

double sigma_xy_cc_anu_n_grv(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_grv (double E, double W, double nu, double m);
/////////////////////////////////podzial na wklady od poszczegolnych czesci///////////////////////////
double sigma_xy_cc_anu_n_grv(double E, double xx, double yy, double m, int ff);
double cr_sec_cc_anu_n_grv (double E, double W, double nu, double m, int ff);
//////////////////////////////////////////////////////////////////////////////////////////////////////
//quark u valence + u sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
/////////////////kwarki za pomoca funkcji struktury 

//kwark d
double sigma_xy_cc_anu_n_u_grv(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_u_grv (double E, double W, double nu, double m);

double sigma_xy_cc_anu_n_sb_grv(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_sb_grv (double E, double W, double nu, double m);
//antykwark u

double sigma_xy_cc_anu_n_dbar_grv(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_dbar_grv (double E, double W, double nu, double m);
/// z poprwkami BY
//kwark u

double sigma_xy_cc_anu_n_BY(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_BY (double E, double W, double nu, double m);

double sigma_xy_cc_anu_n_u_BY(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_u_BY (double E, double W, double nu, double m);//antykwark s

double sigma_xy_cc_anu_n_sb_BY(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_sb_BY (double E, double W, double nu, double m);
//antykwark d

double sigma_xy_cc_anu_n_dbar_BY(double E, double xx, double yy, double m );
double cr_sec_cc_anu_n_dbar_BY (double E, double W, double nu, double m);

#endif
