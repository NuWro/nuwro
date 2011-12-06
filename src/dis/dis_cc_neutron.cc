#include "jednostki.h"
#include "masses.h"
#include "grv94_bodek.h"
#include "parameters.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////NEUTRON//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double
sigma_xy_cc_nu_n_Paschos (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2. - yy * m * m / 4 / M12 / E) * (F3_cc_nu_n (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_cc_nu_n (xx, Q2)) 
     + (1 - yy -  M12 *  xy / (2 *  E) -  m * m / (E * E *  4) 
        - m * m / (2 * M12 * E * xx)) * (F2_cc_nu_n (xx, Q2)));
}

double
cr_sec_cc_nu_n (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_Paschos (E, xx, yy, m) * W / (M12 * E * nu);
}


double
sigma_xy_cc_nu_n_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_GRV94 (xx, Q2)) 
     + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_GRV94 (xx, Q2)) 
     + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
        * (F2_cc_nu_n_GRV94 (xx, Q2)));
}

double
cr_sec_cc_nu_n_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_grv (E, xx, yy, m) * W / (M12 * E * nu);
}


double
sigma_xy_cc_nu_n_grv (double E, double xx, double yy, double m, int ff)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);
  double czesc = 0;

  if (ff == 1)
    {
      czesc =
	(xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_GRV94 (xx, Q2));
    }
  if (ff == 2)
    {
      czesc =
	(1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. -
	 m * m / 2. / M12 / E / xx) * (F2_cc_nu_n_GRV94 (xx, Q2));
    }
  if (ff == 3)
    {
      czesc =
	(xy - xy * yy / 2. -
	 yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_GRV94 (xx, Q2));
    }

  return propagator * stala * M12 * E / pi * czesc;

}

double
cr_sec_cc_nu_n_grv (double E, double W, double nu, double m, int ff)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_grv (E, xx, yy, m, ff) * W / (M12 * E * nu);
}





/////////////////kwarki za pomoca funkcji struktury 

//kwark d
double
sigma_xy_cc_nu_n_d_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_d_grv (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_d_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_nu_n_d_grv (xx, Q2)));
}

double
cr_sec_cc_nu_n_d_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_d_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

//kwark s

double
sigma_xy_cc_nu_n_sb_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_sb_grv (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_sb_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_nu_n_sb_grv (xx, Q2)));
}

double
cr_sec_cc_nu_n_sb_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_sb_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

//antykwark u

double
sigma_xy_cc_nu_n_ubar_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_ubar_grv (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_ubar_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_nu_n_ubar_grv (xx, Q2)));
}

double
cr_sec_cc_nu_n_ubar_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_ubar_grv (E, xx, yy, m) * W / (M12 * E * nu);
}


///kwarki  poprawkami Bodka


double
sigma_xy_cc_nu_n_BY (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n (xx, Q2)) 
	    + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n (xx, Q2)) 
	    + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	    * (F2_cc_nu_n (xx, Q2)));
}

double
cr_sec_cc_nu_n_BY (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_BY (E, xx, yy, m) * W / (M12 * E * nu);
}


//kwark d

double
sigma_xy_cc_nu_n_d_BY (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_d_BY (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_d_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_nu_n_d_BY (xx, Q2)));
}

double
cr_sec_cc_nu_n_d_BY (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_d_BY (E, xx, yy, m) * W / (M12 * E * nu);
}

//kwark s

double
sigma_xy_cc_nu_n_sb_BY (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_sb_BY (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_sb_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_nu_n_sb_BY (xx, Q2)));
}

double
cr_sec_cc_nu_n_sb_BY (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_sb_BY (E, xx, yy, m) * W / (M12 * E * nu);
}

//antykwark u

double
sigma_xy_cc_nu_n_ubar_BY (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi 
         * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * (F3_cc_nu_n_ubar_BY (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_nu_n_ubar_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_nu_n_ubar_BY (xx, Q2)));
}

double
cr_sec_cc_nu_n_ubar_BY (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_ubar_BY (E, xx, yy, m) * W / (M12 * E * nu);
}



//krarki za pomoca pdf
//quark d valence + d sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
double
sigma_xy_cc_nu_n_d (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;

  return stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * 2 * (uv (xx, Q2) - dbar (xx, Q2)) 
	 + (xy * yy + m * m * yy / 2. / M12 / E) * (uv (xx, Q2) + dbar (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - pow (m / E, 2.0) / 4. - m * m / 2. / M12 / E / xx) 
	   * 2 * xx * (uv (xx, Q2) + dbar (xx, Q2))) / xx;
}

double
cr_sec_cc_nu_n_d (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_d (E, xx, yy, m) * W / (M12 * E * nu);
}

//quark s + sbar
double
sigma_xy_cc_nu_n_sb (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;

  return stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) * 2 * sb (xx, Q2) 
         + (xy * yy + m * m * yy / 2. / M12 / E) * sb (xx, Q2) 
         + (1 - yy - M12 * xx * yy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * 2 * xx * sb (xx, Q2)) / xx;
}

double
cr_sec_cc_nu_n_sb (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_sb (E, xx, yy, m) * W / (M12 * E * nu);
}


//dbar
//in common convetion ubar_neutron = dbar(_proton)
double
sigma_xy_cc_nu_n_ubar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;

  return stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * 2 * ubar (xx, Q2) + (xy * yy + m * m * yy / 2. / M12 / E) * ubar (xx, Q2) 
	 + (1 - yy - M12 * xy / 2. / E - pow (m / E, 2.0) / 4. - m * m / 2. / M12 / E / xx) 
	   * 2 * xx * ubar (xx, Q2)) / xx;
}

double
cr_sec_cc_nu_n_ubar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_nu_n_ubar (E, xx, yy, m) * W / (M12 * E * nu);
}


// nu_bar + n --> mu_bar +X (CC, lepton)
// rozklad kwarkow sprawdzony 26.06.05

//Cross section for inelastic scattering on neutron
// nu + n --> mu +X (CC)

double
sigma_xy_cc_anu_n_Paschos (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) * (F2_cc_anu_n (xx, Q2)));
}

double
cr_sec_cc_anu_n (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_Paschos (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


double
sigma_xy_cc_anu_n_grv (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_GRV94 (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_GRV94 (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_anu_n_GRV94 (xx, Q2)));
}

double
cr_sec_cc_anu_n_grv (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

/////////////////////////////////podzial na wklady od poszczegolnych czesci///////////////////////////
double
sigma_xy_cc_anu_n_grv (double E, double xx, double yy, double m, int ff)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  double czesc = 0;


  if (ff == 1)
    {
      czesc =
	(xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_GRV94 (xx, Q2));
    }
  if (ff == 2)
    {
      czesc =
	(1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. -
	 m * m / 2. / M12 / E / xx) * (F2_cc_anu_n_GRV94 (xx, Q2));
    }
  if (ff == 3)
    {
      czesc =
	(xy - xy * yy / 2. -
	 yy * m * m / 4. / M12 / E) * (-F3_cc_anu_n_GRV94 (xx, Q2));
    }

  return propagator * stala * M12 * E / pi * czesc;
}

double
cr_sec_cc_anu_n_grv (double E, double W, double nu, double m, int ff)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_grv (E, xx, yy, m, ff) * W / (M12 * E * nu);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//quark u valence + u sea(=ubar sea) - due to convetion u_neutron=d (d_proton=d)
/////////////////kwarki za pomoca funkcji struktury 

//kwark d
double
sigma_xy_cc_anu_n_u_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_u_grv (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_u_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_anu_n_u_grv (xx, Q2)));
}

double
cr_sec_cc_anu_n_u_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_u_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

//kwark s

double
sigma_xy_cc_anu_n_sb_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_sb_grv (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_sb_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_anu_n_sb_grv (xx, Q2)));
}

double
cr_sec_cc_anu_n_sb_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_sb_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

//antykwark u

double
sigma_xy_cc_anu_n_dbar_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_dbar_grv (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_dbar_grv (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_anu_n_dbar_grv (xx, Q2)));
}

double
cr_sec_cc_anu_n_dbar_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_dbar_grv (E, xx, yy, m) * W / (M12 * E * nu);
}

/// z poprwkami BY
//kwark u

double
sigma_xy_cc_anu_n_BY (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_BY (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	   * (F2_cc_anu_n_BY (xx, Q2)));
}

double
cr_sec_cc_anu_n_BY (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_BY (E, xx, yy, m) * W / (M12 * E * nu);
}




double
sigma_xy_cc_anu_n_u_BY (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_u_BY (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_u_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_anu_n_u_BY (xx, Q2)));
}

double
cr_sec_cc_anu_n_u_BY (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_u_BY (E, xx, yy, m) * W / (M12 * E * nu);
}


//antykwark s

double
sigma_xy_cc_anu_n_sb_BY (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_sb_BY (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_sb_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_anu_n_sb_BY (xx, Q2)));
}

double
cr_sec_cc_anu_n_sb_BY (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_sb_BY (E, xx, yy, m) * W / (M12 * E * nu);
}

//antykwark d

double
sigma_xy_cc_anu_n_dbar_BY (double E, double xx, double yy, double m)	//troche inny wzor z pracy Paschosa
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MW / MW);

  return propagator * stala * M12 * E / pi * ((xy - xy * yy / 2. - yy * m * m / 4. / M12 / E) 
         * (-F3_cc_anu_n_dbar_BY (xx, Q2)) + (xy * yy + m * m * yy / 2. / M12 / E) * (F1_cc_anu_n_dbar_BY (xx, Q2)) 
	 + (1 - yy - M12 * xy / 2. / E - m * m / E / E / 4. - m * m / 2. / M12 / E / xx) 
	 * (F2_cc_anu_n_dbar_BY (xx, Q2)));
}

double
cr_sec_cc_anu_n_dbar_BY (double E, double W, double nu, double m)	//przekroj w zmiennych W, nu
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_cc_anu_n_dbar_BY (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}
