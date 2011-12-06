#include "grv94_bodek.h"
//#include "dis_cr_sec.h"
//#include "fragmentation_nc.h"
#include "jednostki.h"
#include "masses.h"
#include "parameters.h"


//07/07/2005
//Cross section for inelastic scattering on proton
// nu + p --> nu +X (NC, lepton)


double
sigma_xy_nc_nu_p_Paschos (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);


  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p (xx, Q2)) + (1 - yy -
								  M12 * xy /
								  2 / E -
								  m * m / E /
								  E / 4 -
								  m * m / 2 /
								  M12 / E /
								  xx) *
     (F2_nc_p (xx, Q2)));
}

double
cr_sec_nc_nu_p (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_Paschos (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

/////////////////tylk grv/////////////////////////////////////////////////
double
sigma_xy_nc_nu_p_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * (F3_nc_p_GRV94 (xx,
						  Q2)) + (xy * yy +
							  m * m * yy / 2 /
							  M12 / E) *
     (F1_nc_p_GRV94 (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				 m * m / E / E / 4 -
				 m * m / 2 / M12 / E / xx) *
     (F2_nc_p_GRV94 (xx, Q2)));
}

double
cr_sec_nc_nu_p_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_grv (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//kwark uval + usea
double
sigma_xy_nc_nu_p_u (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_u (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_u (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_p_u (xx, Q2)));
}

double
cr_sec_nc_nu_p_u (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_u (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark dval + dsea
double
sigma_xy_nc_nu_p_d (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_d (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_d (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_p_d (xx, Q2)));
}

double
cr_sec_nc_nu_p_d (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_d (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//antykwark u 
double
sigma_xy_nc_nu_p_ubar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_ubar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_ubar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_p_ubar (xx, Q2)));
}

double
cr_sec_nc_nu_p_ubar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_ubar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwar d
double
sigma_xy_nc_nu_p_dbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_dbar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_dbar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_p_dbar (xx, Q2)));
}

double
cr_sec_nc_nu_p_dbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_dbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark s
double
sigma_xy_nc_nu_p_s (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_s (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_s (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_p_s (xx, Q2)));
}

double
cr_sec_nc_nu_p_s (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_s (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwark s
double
sigma_xy_nc_nu_p_sbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_p_sbar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_p_sbar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_p_sbar (xx, Q2)));
}

double
cr_sec_nc_nu_p_sbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_p_sbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}



// nu_bar + p --> nu_bar +X (NC, lepton)


double
sigma_xy_nc_anu_p_Paschos (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p (xx,
						   Q2)) + (xy * yy +
							   m * m * yy / 2 /
							   M12 / E) *
     (F1_nc_p (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			   m * m / 2 / M12 / E / xx) * (F2_nc_p (xx, Q2)));
}

double
cr_sec_nc_anu_p (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_Paschos (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}



double
sigma_xy_nc_anu_p_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_GRV94 (xx,
							 Q2)) + (xy * yy +
								 m * m * yy /
								 2 / M12 /
								 E) *
     (F1_nc_p_GRV94 (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				 m * m / E / E / 4 -
				 m * m / 2 / M12 / E / xx) *
     (F2_nc_p_GRV94 (xx, Q2)));
}

double
cr_sec_nc_anu_p_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_grv (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//kwark uval + usea
double
sigma_xy_nc_anu_p_u (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_u (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_p_u (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_p_u (xx,
								     Q2)));
}

double
cr_sec_nc_anu_p_u (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_u (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark dval + dsea
double
sigma_xy_nc_anu_p_d (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_d (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_p_d (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_p_d (xx,
								     Q2)));
}

double
cr_sec_nc_anu_p_d (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_d (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//antykwark u 
double
sigma_xy_nc_anu_p_ubar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_ubar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_p_ubar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_p_ubar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_p_ubar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_ubar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwar d
double
sigma_xy_nc_anu_p_dbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_dbar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_p_dbar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_p_dbar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_p_dbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_dbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark s
double
sigma_xy_nc_anu_p_s (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_s (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_p_s (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_p_s (xx,
								     Q2)));
}

double
cr_sec_nc_anu_p_s (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_s (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwark s
double
sigma_xy_nc_anu_p_sbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_p_sbar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_p_sbar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_p_sbar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_p_sbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_p_sbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////NEUTRON//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// nu + n --> nu +X (NC, lepton)

double
sigma_xy_nc_nu_n_Paschos (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n (xx, Q2)) + (1 - yy -
								  M12 * xy /
								  2 / E -
								  m * m / E /
								  E / 4 -
								  m * m / 2 /
								  M12 / E /
								  xx) *
     (F2_nc_n (xx, Q2)));
}

double
cr_sec_nc_nu_n (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_Paschos (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//kwark uval + usea
double
sigma_xy_nc_nu_n_u (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_u (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_u (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_n_u (xx, Q2)));
}

/////////////////tylko grv
double
sigma_xy_nc_nu_n_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * (F3_nc_n_GRV94 (xx,
						  Q2)) + (xy * yy +
							  m * m * yy / 2 /
							  M12 / E) *
     (F1_nc_n_GRV94 (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				 m * m / E / E / 4 -
				 m * m / 2 / M12 / E / xx) *
     (F2_nc_n_GRV94 (xx, Q2)));
}

double
cr_sec_nc_nu_n_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_grv (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}







///////////////////////





double
cr_sec_nc_nu_n_u (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_u (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark dval + dsea
double
sigma_xy_nc_nu_n_d (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_d (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_d (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_n_d (xx, Q2)));
}

double
cr_sec_nc_nu_n_d (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_d (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//antykwark u 
double
sigma_xy_nc_nu_n_ubar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_ubar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_ubar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_n_ubar (xx, Q2)));
}

double
cr_sec_nc_nu_n_ubar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_ubar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwar d
double
sigma_xy_nc_nu_n_dbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_dbar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_dbar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_n_dbar (xx, Q2)));
}

double
cr_sec_nc_nu_n_dbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_dbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark s
double
sigma_xy_nc_nu_n_s (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_s (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_s (xx, Q2)) + (1 - yy -
								    M12 * xy /
								    2 / E -
								    m * m /
								    E / E /
								    4 -
								    m * m /
								    2 / M12 /
								    E / xx) *
     (F2_nc_n_s (xx, Q2)));
}

double
cr_sec_nc_nu_n_s (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_s (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwark s
double
sigma_xy_nc_nu_n_sbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 - yy * m * m / 4 / M12 / E) * (F3_nc_n_sbar (xx, Q2)) +
     (xy * yy + m * m * yy / 2 / M12 / E) * (F1_nc_n_sbar (xx, Q2)) + (1 -
								       yy -
								       M12 *
								       xy /
								       2 / E -
								       m * m /
								       E / E /
								       4 -
								       m * m /
								       2 /
								       M12 /
								       E /
								       xx) *
     (F2_nc_n_sbar (xx, Q2)));
}

double
cr_sec_nc_nu_n_sbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_nu_n_sbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}



// nu_bar + n --> nu_bar +X (NC, lepton)


double
sigma_xy_nc_anu_n_Paschos (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n (xx,
						   Q2)) + (xy * yy +
							   m * m * yy / 2 /
							   M12 / E) *
     (F1_nc_n (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			   m * m / 2 / M12 / E / xx) * (F2_nc_n (xx, Q2)));
}

double
cr_sec_nc_anu_n (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_Paschos (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}



double
sigma_xy_nc_anu_n_grv (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  double propagator = 1. / kwad (1. + Q2 / MZ / MZ);

  return propagator * stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_GRV94 (xx,
							 Q2)) + (xy * yy +
								 m * m * yy /
								 2 / M12 /
								 E) *
     (F1_nc_n_GRV94 (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				 m * m / E / E / 4 -
				 m * m / 2 / M12 / E / xx) *
     (F2_nc_n_GRV94 (xx, Q2)));
}

double
cr_sec_nc_anu_n_grv (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_grv (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//kwark uval + usea
double
sigma_xy_nc_anu_n_u (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_u (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_n_u (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_n_u (xx,
								     Q2)));
}

double
cr_sec_nc_anu_n_u (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_u (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark dval + dsea
double
sigma_xy_nc_anu_n_d (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_d (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_n_d (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_n_d (xx,
								     Q2)));
}

double
cr_sec_nc_anu_n_d (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_d (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}

//antykwark u 
double
sigma_xy_nc_anu_n_ubar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_ubar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_n_ubar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_n_ubar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_n_ubar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_ubar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwar d
double
sigma_xy_nc_anu_n_dbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_dbar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_n_dbar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_n_dbar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_n_dbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_dbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//kwark s
double
sigma_xy_nc_anu_n_s (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_s (xx,
						     Q2)) + (xy * yy +
							     m * m * yy / 2 /
							     M12 / E) *
     (F1_nc_n_s (xx, Q2)) + (1 - yy - M12 * xy / 2 / E - m * m / E / E / 4 -
			     m * m / 2 / M12 / E / xx) * (F2_nc_n_s (xx,
								     Q2)));
}

double
cr_sec_nc_anu_n_s (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_s (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}


//antykwark s
double
sigma_xy_nc_anu_n_sbar (double E, double xx, double yy, double m)
{
  double Q2 = 2 * M12 * E * xx * yy;
  double xy = xx * yy;
  double ME4 = 4 * M12 * E;
  return stala * M12 * E / pi *
    ((xy - xy * yy / 2 -
      yy * m * m / 4 / M12 / E) * ((-1) * F3_nc_n_sbar (xx,
							Q2)) + (xy * yy +
								m * m * yy /
								2 / M12 / E) *
     (F1_nc_n_sbar (xx, Q2)) + (1 - yy - M12 * xy / 2 / E -
				m * m / E / E / 4 -
				m * m / 2 / M12 / E / xx) * (F2_nc_n_sbar (xx,
									   Q2)));
}

double
cr_sec_nc_anu_n_sbar (double E, double W, double nu, double m)
{
  double xx = (M2 - W * W + 2 * M12 * nu) / (2 * M12 * nu);
  double yy = nu / E;
  return sigma_xy_nc_anu_n_sbar (E, xx, yy, m) * W / (M12 * E * nu)	//jakobian
    ;
}
