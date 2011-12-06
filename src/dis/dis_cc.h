#ifndef _dis_cc_h_
#define _dis_cc_h_


//Cross section for inelastic scattering on proton
// nu + p --> mu +X (CC, lepton)

double sigma_xy_cc_nu_p_Paschos(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);
  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(xx*F3_cc_nu_p(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(xx*F1_cc_nu_p(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(xx*F2_cc_nu_p(xx,Q2))
		       )/xx
;}

double cr_sec_cc_nu_p (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_Paschos(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

double sigma_xy_cc_nu_p_grv(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(F3_cc_nu_p_GRV94(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(F1_cc_nu_p_GRV94(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(F2_cc_nu_p_GRV94(xx,Q2))
		       )/xx
;}

double cr_sec_cc_nu_p_grv (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_grv(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}



//quark d valence + d sea(=dbar sea)
double sigma_xy_cc_nu_p_d(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*(dv(xx,Q2)+dbar(xx,Q2)) +
		       (xy*yy+m*m*yy/2/M12/E)*(dv(xx,Q2)+dbar(xx,Q2))+
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*(dv(xx,Q2)+dbar(xx,Q2))
		       )/xx

;}

double cr_sec_cc_nu_p_d (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_d(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//sum of quark s and sbar
double sigma_xy_cc_nu_p_sb(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;

  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*sb(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*sb(xx,Q2)+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*sb(xx,Q2)//ok
		       )/xx
;}

double cr_sec_cc_nu_p_sb (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_sb(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


// antiquark u
double sigma_xy_cc_nu_p_ubar(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;


  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(-1)*2*ubar(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*ubar(xx,Q2)+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*ubar(xx,Q2)//ok
		       )/xx


;}

double cr_sec_cc_nu_p_ubar (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_ubar(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


// nu_bar + p --> mu_bar +X (CC, lepton)
//////DO SPRAWDZENIA ROZK?AD KWARKOW!!!!!!!!!!!
//kwarki sprawdzone 26.06.05

double sigma_xy_cc_anu_p_Paschos(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*((-1)*xx*F3_cc_anu_p(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(xx*F1_cc_anu_p(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(xx*F2_cc_anu_p(xx,Q2))
		       )/xx
;}


double cr_sec_cc_anu_p (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_Paschos(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


double sigma_xy_cc_anu_p_grv(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*((-1)*F3_cc_anu_p_GRV94(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(F1_cc_anu_p_GRV94(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(F2_cc_anu_p_GRV94(xx,Q2))
		       )/xx
;}


double cr_sec_cc_anu_p_grv (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_grv(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//quark u valence + u sea(=ubar sea)
double sigma_xy_cc_anu_p_u(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(-2)*(uv(xx,Q2)+ubar(xx,Q2)) +
		       (xy*yy+m*m*yy/2/M12/E)*(uv(xx,Q2)+ubar(xx,Q2))+
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*(uv(xx,Q2)+ubar(xx,Q2))
		       )/xx

;}

double cr_sec_cc_anu_p_u (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_u(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//sum of quark s and sbar
double sigma_xy_cc_anu_p_sb(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;

  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(-2)*sb(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*sb(xx,Q2)+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*sb(xx,Q2)//ok
		       )/xx
;}

double cr_sec_cc_anu_p_sb (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_sb(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


// antiquark d
double sigma_xy_cc_anu_p_dbar(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;


  return stala*M12*E/pi*((xy-xx*yy*yy/2-yy*m*m/4/M12/E)*(-2)*dbar(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*(-dbar(xx,Q2))+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*(-ubar(xx,Q2))//ok
		       )/xx


;}

double cr_sec_cc_anu_p_dbar (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_dbar(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////NEUTRON//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double sigma_xy_cc_nu_n_Paschos(double E, double xx, double yy, double m )//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(xx*F3_cc_nu_n(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(xx*F1_cc_nu_n(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(xx*F2_cc_nu_n(xx,Q2))
		       )/xx
;}

double cr_sec_cc_nu_n (double E, double W, double nu, double m)//przekroj w zmiennych W, nu
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_n_Paschos(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


double sigma_xy_cc_nu_n_grv(double E, double xx, double yy, double m )//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(F3_cc_nu_n_GRV94(xx,Q2))+
		                         (xy*yy+m*m*yy/2/M12/E)*(F1_cc_nu_n_GRV94(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(F2_cc_nu_n_GRV94(xx,Q2))
		       )/xx
;}

double cr_sec_cc_nu_n_grv (double E, double W, double nu, double m)//przekroj w zmiennych W, nu
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_n_grv(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}



//quark d valence + d sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
double sigma_xy_cc_nu_n_d(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*(uv(xx,Q2)-dbar(xx,Q2)) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*(uv(xx,Q2)+dbar(xx,Q2))+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*(uv(xx,Q2)+dbar(xx,Q2))//ok
		       )/xx
;}

double cr_sec_cc_nu_n_d (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_n_d(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//quark s + sbar

double sigma_xy_cc_nu_n_sb(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*sb(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*sb(xx,Q2)+//ok
		       (1-yy-M12*xx*yy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*2*xx*sb(xx,Q2)//ok
		       )/xx
;}

double cr_sec_cc_nu_n_sb (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_n_sb(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


//dbar
//in common convetion ubar_neutron = dbar(_proton)
double sigma_xy_cc_nu_n_ubar(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*ubar(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*ubar(xx,Q2)+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*ubar(xx,Q2)//ok
		       )/xx

;}

double cr_sec_cc_nu_n_ubar (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_n_ubar(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


// nu_bar + n --> mu_bar +X (CC, lepton)
//////DO SPRAWDZENIA ROZK?AD KWARKOW!!!!!!!!!!!
// rozklad kwarkow sprawdzony 26.06.05

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double sigma_xy_cc_anu_n_Paschos(double E, double xx, double yy, double m )//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*((-1)*xx*F3_cc_anu_n(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(xx*F1_cc_anu_n(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(xx*F2_cc_anu_n(xx,Q2))
		       )/xx
;}

double cr_sec_cc_anu_n (double E, double W, double nu, double m)//przekroj w zmiennych W, nu
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_n_Paschos(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


double sigma_xy_cc_anu_n_grv(double E, double xx, double yy, double m )//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*((-1)*F3_cc_anu_n_GRV94(xx,Q2))+
		       (xy*yy+m*m*yy/2/M12/E)*(F1_cc_anu_n_GRV94(xx,Q2))+
		       (1-yy-M12*xy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*(F2_cc_anu_n_GRV94(xx,Q2))
		       )/xx
;}

double cr_sec_cc_anu_n_grv (double E, double W, double nu, double m)//przekroj w zmiennych W, nu
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_n_grv(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}





//quark u valence + u sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
double sigma_xy_cc_anu_n_u(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*2*(-dv(xx,Q2)+ubar(xx,Q2)) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*(dv(xx,Q2)+ubar(xx,Q2))+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*(dv(xx,Q2)+ubar(xx,Q2))//ok
		       )/xx
;}

double cr_sec_cc_anu_n_u (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_n_u(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//quark s + sbar

double sigma_xy_cc_anu_n_sb(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(-2)*sb(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*sb(xx,Q2)+//ok
		       (1-yy-M12*xx*yy/2/E-m*m/E/E/4-m*m/2/M12/E/xx)*2*xx*sb(xx,Q2)//ok
		       )/xx
;}

double cr_sec_cc_anu_n_sb (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_n_sb(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}


//dbar
//in common convetion ubar_neutron = dbar(_proton)
double sigma_xy_cc_anu_n_dbar(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2-yy*m*m/4/M12/E)*(-2)*dbar(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2/M12/E)*dbar(xx,Q2)+//ok
		       (1-yy-M12*xy/2/E-pow(m/E,2.0)/4-m*m/2/M12/E/xx)*2*xx*dbar(xx,Q2)//ok
		       )/xx

;}

double cr_sec_cc_anu_n_dbar (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_n_dbar(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}







#endif
