#include "jednostki.h"
#include "masses.h"
#include "grv94_bodek.h"
#include "parameters.h"

//Cross section for inelastic scattering on proton
// nu + p --> mu +X (CC, lepton)

double sigma_xy_cc_nu_p_Paschos(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;

  double propagator = 1./kwad(1.+ Q2/MW/MW);
  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p(xx,Q2))+
		                          (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p(xx,Q2))+
    	            (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p(xx,Q2)));
}

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

  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_GRV94(xx,Q2))+
                 		          (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_GRV94(xx,Q2))+
 	            (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_GRV94(xx,Q2)));
}


double cr_sec_cc_nu_p_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_grv(E,xx,yy,m)*W/(M12*E*nu);
}


double sigma_xy_cc_nu_p_BY(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p (xx,Q2))+
		       (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p (xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p(xx,Q2))   );
}


double cr_sec_cc_nu_p_BY (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_BY (E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}



//quark d valence + d sea(=dbar sea)
double sigma_xy_cc_nu_p_d(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*2*(dv(xx,Q2)+dbar(xx,Q2)) +
		       (xy*yy+m*m*yy/2./M12/E)*(dv(xx,Q2)+dbar(xx,Q2))+
		       (1-yy-M12*xy/2./E-pow(m/E,2.0)/4.-m*m/2./M12/E/xx)*2*xx*(dv(xx,Q2)+dbar(xx,Q2))
		       )/xx;
}

double cr_sec_cc_nu_p_d (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_d(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//quark sb
double sigma_xy_cc_nu_p_sb(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;

  return stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*2*sb(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2./M12/E)*sb(xx,Q2)+//ok
		       (1-yy-M12*xy/2./E-pow(m/E,2.0)/4.-m*m/2./M12/E/xx)*2*xx*sb(xx,Q2)//ok
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


  return stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-1)*2*ubar(xx,Q2) +//ok
		       (xy*yy+m*m*yy/2./M12/E)*ubar(xx,Q2)+//ok
		       (1-yy-M12*xy/2./E-pow(m/E,2.0)/4.-m*m/2./M12/E/xx)*2*xx*ubar(xx,Q2)//ok
		       )/xx


;}

double cr_sec_cc_nu_p_ubar (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_ubar(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

///przekroje czynne z funkcjami struktury rozbitymina krawki

//kwrak d
double sigma_xy_cc_nu_p_d_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_d_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_d_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_d_grv(xx,Q2))   );
}


double cr_sec_cc_nu_p_d_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_d_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_nu_p_d_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_d_BY(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_d_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_d_BY(xx,Q2))   );
}


double cr_sec_cc_nu_p_d_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_d_BY(E,xx,yy,m)*W/(M12*E*nu);
}

//kwark sb

double sigma_xy_cc_nu_p_sb_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_sb_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_sb_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_sb_grv(xx,Q2))   );
}


double cr_sec_cc_nu_p_sb_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_sb_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_nu_p_sb_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_sb_BY(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_sb_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_sb_BY(xx,Q2))   );
}


double cr_sec_cc_nu_p_sb_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_sb_BY(E,xx,yy,m)*W/(M12*E*nu);
}

//antykwark u

double sigma_xy_cc_nu_p_ubar_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_ubar_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_ubar_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_ubar_grv(xx,Q2))   );
}


double cr_sec_cc_nu_p_ubar_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_ubar_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_nu_p_ubar_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(F3_cc_nu_p_ubar_BY(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_nu_p_ubar_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_nu_p_ubar_BY(xx,Q2))   );
}


double cr_sec_cc_nu_p_ubar_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_nu_p_ubar_BY(E,xx,yy,m)*W/(M12*E*nu);
}




// nu_bar + p --> mu_bar +X (CC, lepton)
//////DO SPRAWDZENIA ROZK?AD KWARKOW!!!!!!!!!!!
//kwarki sprawdzone 26.06.05

double sigma_xy_cc_anu_p_Paschos(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p(xx,Q2))+
		                                (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p(xx,Q2))+
		          (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p(xx,Q2)));
}


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

  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_GRV94(xx,Q2))+
		       (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_GRV94(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_GRV94(xx,Q2))
		       );
		       
}


double cr_sec_cc_anu_p_grv (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_grv(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

double sigma_xy_cc_anu_p_BY(double E, double xx, double yy, double m)//troche inny wzor z pracy Paschosa
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p(xx,Q2))+
		       (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p(xx,Q2))
		       )/xx;
}


double cr_sec_cc_anu_p_BY (double E, double W, double nu, double m)
{ double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_BY(E,xx,yy,m)*W/(M12*E*nu)//jakobian
;}

//quark u valence + u sea(=ubar sea)
double sigma_xy_cc_anu_p_u(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  
  return stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-2)*(uv(xx,Q2)+ubar(xx,Q2)) +
		       (xy*yy+m*m*yy/2./M12/E)*(uv(xx,Q2)+ubar(xx,Q2))+
		       (1-yy-M12*xy/2./E-pow(m/E,2.0)/4.-m*m/2./M12/E/xx)*2*xx*(uv(xx,Q2)+ubar(xx,Q2))
		       );
		       
}



///przekroje czynne z funkcjami struktury rozbitymina krawki

//kwrak u
double sigma_xy_cc_anu_p_u_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_u_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_u_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_u_grv(xx,Q2))   );
}


double cr_sec_cc_anu_p_u_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_u_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_anu_p_u_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_u_BY(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_u_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_u_BY(xx,Q2))   );
}


double cr_sec_cc_anu_p_u_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_u_BY(E,xx,yy,m)*W/(M12*E*nu);
}

//kwark sb

double sigma_xy_cc_anu_p_sb_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_sb_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_sb_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_sb_grv(xx,Q2))   );
}


double cr_sec_cc_anu_p_sb_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_sb_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_anu_p_sb_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_sb_BY(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_sb_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_sb_BY(xx,Q2))   );
}


double cr_sec_cc_anu_p_sb_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_sb_BY(E,xx,yy,m)*W/(M12*E*nu);
}

//antykwark d

double sigma_xy_cc_anu_p_dbar_grv(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_dbar_grv(xx,Q2))+
                        	           (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_dbar_grv(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_dbar_grv(xx,Q2))   );
}


double cr_sec_cc_anu_p_dbar_grv (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_dbar_grv(E,xx,yy,m)*W/(M12*E*nu);
}

double sigma_xy_cc_anu_p_dbar_BY(double E, double xx, double yy, double m)
{
  double Q2=2*M12*E*xx*yy;
  double xy=xx*yy;
  double ME4=4*M12*E;
  double propagator = 1./kwad(1.+Q2/MW/MW);

  return   propagator*stala*M12*E/pi*((xy-xy*yy/2.-yy*m*m/4./M12/E)*(-F3_cc_anu_p_dbar_BY(xx,Q2))+
                         	            (xy*yy+m*m*yy/2./M12/E)*(F1_cc_anu_p_dbar_BY(xx,Q2))+
		       (1-yy-M12*xy/2./E-m*m/E/E/4.-m*m/2./M12/E/xx)*(F2_cc_anu_p_dbar_BY(xx,Q2)) );
}


double cr_sec_cc_anu_p_dbar_BY (double E, double W, double nu, double m)
{
  double xx=(M2-W*W+2*M12*nu)/(2*M12*nu);
  double yy=nu/E;
  return sigma_xy_cc_anu_p_dbar_BY(E,xx,yy,m)*W/(M12*E*nu);
}

