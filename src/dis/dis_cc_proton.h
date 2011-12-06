#ifndef _dis_cc_proton_h_
#define _dis_cc_proton_h_

double sigma_xy_cc_nu_p_Paschos(double E, double xx, double yy, double m);
double cr_sec_cc_nu_p (double E, double W, double nu, double m);
double sigma_xy_cc_nu_p_grv(double E, double xx, double yy, double m);
double cr_sec_cc_nu_p_grv (double E, double W, double nu, double m);

double sigma_xy_cc_nu_p_BY(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_BY (double E, double W, double nu, double m);


//quark d valence + d sea(=dbar sea)
double sigma_xy_cc_nu_p_d(double E, double xx, double yy, double m);
double cr_sec_cc_nu_p_d (double E, double W, double nu, double m);
//quark sb
double sigma_xy_cc_nu_p_sb(double E, double xx, double yy, double m);
double cr_sec_cc_nu_p_sb (double E, double W, double nu, double m);

// antiquark u
double sigma_xy_cc_nu_p_ubar(double E, double xx, double yy, double m);
double cr_sec_cc_nu_p_ubar (double E, double W, double nu, double m);
//przekroje czynne z funkcjami struktury rozbitymina krawki

//kwrak d
double sigma_xy_cc_nu_p_d_grv(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_d_grv (double E, double W, double nu, double m);
double sigma_xy_cc_nu_p_d_BY(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_d_BY (double E, double W, double nu, double m);
//kwark sb

double sigma_xy_cc_nu_p_sb_grv(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_sb_grv (double E, double W, double nu, double m);
double sigma_xy_cc_nu_p_sb_BY(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_sb_BY (double E, double W, double nu, double m);
//antykwark u

double sigma_xy_cc_nu_p_ubar_grv(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_ubar_grv (double E, double W, double nu, double m);
double sigma_xy_cc_nu_p_ubar_BY(double E, double xx, double yy, double m);

double cr_sec_cc_nu_p_ubar_BY (double E, double W, double nu, double m);


// nu_bar + p --> mu_bar +X (CC, lepton)
//////DO SPRAWDZENIA ROZK?AD KWARKOW!!!!!!!!!!!
//kwarki sprawdzone 26.06.05

double sigma_xy_cc_anu_p_Paschos(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p (double E, double W, double nu, double m);
double sigma_xy_cc_anu_p_grv(double E, double xx, double yy, double m);
double cr_sec_cc_anu_p_grv (double E, double W, double nu, double m);
double sigma_xy_cc_anu_p_BY(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p_BY (double E, double W, double nu, double m);
//quark u valence + u sea(=ubar sea)
double sigma_xy_cc_anu_p_u(double E, double xx, double yy, double m);

///przekroje czynne z funkcjami struktury rozbitymina krawki

//kwrak u
double sigma_xy_cc_anu_p_u_grv(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p_u_grv (double E, double W, double nu, double m);
double sigma_xy_cc_anu_p_u_BY(double E, double xx, double yy, double m);
double cr_sec_cc_anu_p_u_BY (double E, double W, double nu, double m);
//kwark sb

double sigma_xy_cc_anu_p_sb_grv(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p_sb_grv (double E, double W, double nu, double m);
double sigma_xy_cc_anu_p_sb_BY(double E, double xx, double yy, double m);
double cr_sec_cc_anu_p_sb_BY (double E, double W, double nu, double m);
//antykwark d

double sigma_xy_cc_anu_p_dbar_grv(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p_dbar_grv (double E, double W, double nu, double m);
double sigma_xy_cc_anu_p_dbar_BY(double E, double xx, double yy, double m);

double cr_sec_cc_anu_p_dbar_BY (double E, double W, double nu, double m);



#endif
