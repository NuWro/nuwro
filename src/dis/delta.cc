#include <iostream>
#include <cmath>
#include "parameters.h"
#include "pdg_name.h"
#include "jednostki.h"
#include "masses.h"
#include "LeptonMass.h"
#include <cstdlib>
using namespace std;

double Masa = M12;

double MM = Masa * Masa;

int znak;

const int neutrino = 1;

const int antyneutrino = -1;

const int CC = 1;
const int NC = -1;

double masa_leptonu;

int W_Z;			// Decyduje czy mamy prady naladowane czy neutralne
 //int Form_Faktory;

//double sin2thetaW = sin_2_theta_W;
double masa_rezonansu;

double W_limit (1000 * GeV);	// --- jesli nie jest ustalna to az 1000*GeV !!!

int form_faktory_wektorowe (int FFset, double Q2, double W, int NC_czy_CC,
                            double *c3v, double *c4v, double *c5v, double *c6v)
{
  double funkcja;

  double MW = 0.84 * GeV;
  double MW2 = MW * MW;


  switch (FFset)
    {
    case 1:
      funkcja = 1.0 / (1 + Q2 / MW2) / (1 + Q2 / MW2);
      break;			//Graczyk Sobczyk, from Olga

    case 2:
      funkcja = 1.95 / (1 + Q2 / MW2) / (1 + Q2 / MW2) / (1 + Q2 / 4 / MW2);
      break;			//Olga

    case 3:
      funkcja = 1.95 / (1 + Q2 / MW2) / (1 + Q2 / MW2) / (1 + Q2 / 4 / MW2);
      break;			//Olga

    case 4:
      funkcja = 1.95 / (1 + Q2 / MW2) / (1 + Q2 / MW2) / (1 + Q2 / 4 / MW2);
      break;			//Olga

    case 5:
      funkcja = 1.95 / (1 + Q2 / MW2) / (1 + Q2 / MW2) / (1 + Q2 / 4 / MW2);
      break;			//Olga

    case 6:
      funkcja = 2.05 / pow (1 + Q2 / 0.54 / GeV2, 2);
      break;			//Alvarez Ruso et al; Jarek's choice

    default:
      funkcja = 1.0 / (1 + Q2 / MW2) / (1 + Q2 / MW2);
      break;
    }


  if (NC_czy_CC == CC)
    {

switch (FFset)
    {
    case 1:
      *c3v = 2.13*funkcja/(1+Q2/4.0/MW2);
      *c4v = - 1.51*funkcja/(1+Q2/4.0/MW2);
      *c5v = 0.48*funkcja/(1+Q2/0.776/MW2);
      *c6v = 0;
      break;			//Graczyk Sobczyk

    case 2:
      *c3v = funkcja;
      *c4v = - funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 3:
     *c3v = funkcja;
      *c4v = - funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 4:
     *c3v = funkcja;
      *c4v = - funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 5:
      *c3v = funkcja;
      *c4v = - funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 6:
      *c3v = funkcja;
      *c4v = - funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Alvarez Ruso et al; Jarek's choice

//case 7 is missing

    default:
      *c3v = 2.13*funkcja/(1+Q2/4.0/MW2);
      *c4v = - 1.51*funkcja/(1+Q2/4.0/MW2);
      *c5v = 0.48*funkcja/(1+Q2/0.776/MW2);
      *c6v = 0;
      break;
    }
      
    }
  if (NC_czy_CC == NC)
    {
      //cerr<<"neutralne"<<endl;
switch (FFset)
    {
    case 1:
      *c3v = (1 - 2 * sin2thetaW) *2.13*funkcja/(1+Q2/4.0/MW2);
      *c4v = - (1 - 2 * sin2thetaW) *1.51*funkcja/(1+Q2/4.0/MW2);
      *c5v = (1 - 2 * sin2thetaW) *0.48*funkcja/(1+Q2/0.776/MW2);
      *c6v = 0;
      break;			//Graczyk Sobczyk

    case 2:
      *c3v = (1 - 2 * sin2thetaW) *funkcja;
      *c4v = - (1 - 2 * sin2thetaW) *funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 3:
     *c3v = (1 - 2 * sin2thetaW) *funkcja;
      *c4v = - (1 - 2 * sin2thetaW) *funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 4:
     *c3v = (1 - 2 * sin2thetaW) *funkcja;
      *c4v = - (1 - 2 * sin2thetaW) *funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 5:
      *c3v = (1 - 2 * sin2thetaW) *funkcja;
      *c4v = - (1 - 2 * sin2thetaW) *funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Olga

    case 6:
      *c3v = (1 - 2 * sin2thetaW) *funkcja;
      *c4v = - (1 - 2 * sin2thetaW) *funkcja*Masa/W;
      *c5v = 0;
      *c6v = 0;
      break;			//Alvarez Ruso et al; Jarek's choice

//case 7 is missing

    default:
      *c3v = (1 - 2 * sin2thetaW) *2.13*funkcja/(1+Q2/4.0/MW2);
      *c4v = - (1 - 2 * sin2thetaW) *1.51*funkcja/(1+Q2/4.0/MW2);
      *c5v = (1 - 2 * sin2thetaW) *0.48*funkcja/(1+Q2/0.776/MW2);
      *c6v = 0;
      break;
    }

    }

  return 0;
}



int form_faktory_aksjalne (int FFset, double delta_axial_mass, double delta_C5A, double Q2, double W,
                           double *c3a, double *c4a, double *c5a, double *c6a)
{

  double MA;			// 1.00*GeV;
  double MA2;
  double C5A0;
  double funkcja3, funkcja4, funkcja5, funkcja6;

double a3A, b3A, c3A, d3A, e3A, a4A, b4A, c4A, d4A, e4A, a5A, b5A, c5A, d5A, e5A, a6A, b6A, c6A, d6A, e6A;

  //double  c4(-0.3), c5(1.2);
  double a = -1.21;
  double b = 2 * GeV2;

  switch (FFset)
    {

    case 1:
	MA = delta_axial_mass * GeV;
	MA2=MA*MA;
	C5A0=delta_C5A;

      funkcja5 = C5A0 / (1 + Q2 / MA2) / (1 + Q2 / MA2);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;

    
    case 2:
      MA = 1.05 * GeV;
      MA2 = MA * MA;

      funkcja5 = 1.2 / (1 + Q2 / MA2) / (1 + Q2 / MA2) / (1 + Q2 / MA2 / 3);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;			//Olga 2.12

    case 3:
      MA = 0.84 * GeV;
      MA2 = MA * MA;

      funkcja5 = 1.2 / (1 + Q2 / MA2) / (1 + Q2 / MA2) / (1 + Q2 / MA2 / 3);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;			//Olga 2.12

    case 4:
      MA = 1.05 * GeV;
      MA2 = MA * MA;

      funkcja5 = 1.2 / (1 + Q2 / MA2) / (1 + Q2 / MA2) / (1 + 2 * Q2 / MA2);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;

    case 5:
      MA = 0.95 * GeV;
      MA2 = MA * MA;

      funkcja5 =
	1.2 / (1 + Q2 / MA2) / (1 + Q2 / MA2) / (1 + Q2 / MA2 / 3) / (1 +
								      Q2 /
								      MA2 /
								      3);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;

    case 6:
      MA = 1.0 * GeV;
      MA2 = MA * MA;

      funkcja3 = 0;
      funkcja4 =
	-0.3 * (1 + a * Q2 / (b + Q2)) / (1 + Q2 / MA2) / (1 + Q2 / MA2);
      funkcja5 =
	1.2 * (1 + a * Q2 / (b + Q2)) / (1 + Q2 / MA2) / (1 + Q2 / MA2);
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;

	case 7:

//c3A

a3A               = 0.0350791 * GeV * GeV * GeV * GeV * GeV * GeV;
b3A               = -0.0460531 * GeV * GeV * GeV * GeV;
c3A               = 0.697158 * GeV * GeV * GeV * GeV;
d3A               = -0.0159286 * GeV * GeV;
e3A               = 0.430003;

funkcja3 = (a3A + b3A*Q2)/(1 + c3A*Q2 + d3A*Q2*Q2 + e3A*Q2*Q2*Q2);

//c4a

a4A               = -0.262166 * GeV * GeV * GeV * GeV * GeV * GeV;
b4A             = 0.546486 * GeV * GeV * GeV * GeV;
c4A               = 2.22715 * GeV * GeV * GeV * GeV;
d4A               = -2.48587 * GeV * GeV;
e4A               = 5.71295;

funkcja4 = (a4A + b4A*Q2)/(1 + c4A*Q2+ d4A*Q2*Q2 + e4A*Q2*Q2*Q2);

//c5a

a5A               = 0.933472 * GeV * GeV * GeV * GeV * GeV * GeV;
b5A               = -1.05174  * GeV * GeV * GeV * GeV;
c5A               = 1.50965 * GeV * GeV * GeV * GeV;
d5A               = 0.539646 * GeV * GeV;
e5A               = 1.04352;

funkcja5 = (a5A + b5A*Q2)/(1 + c5A*Q2+ d5A*Q2*Q2 + e5A*Q2*Q2*Q2);

//c6a

a6A = 0.0466997 * GeV * GeV * GeV * GeV * GeV * GeV;
b6A               = 0.286002 * GeV * GeV * GeV * GeV;
c6A               = 0.604646 * GeV * GeV * GeV * GeV;
d6A               = 0.423329 * GeV * GeV;
e6A               = -0.0392785;

funkcja6 = (a6A + b6A*Q2)/(1 + c6A*Q2+ d6A*Q2*Q2 + e6A*Q2*Q2*Q2);

break;

    default:
      MA = 1.05 * GeV;
      MA2 = MA * MA;

      funkcja5 = 1.2 / (1 + Q2 / MA2) / (1 + Q2 / MA2) / (1 + Q2 / MA2 / 3);
      funkcja3 = 0;
      funkcja4 = -0.25 * funkcja5;
      funkcja6 = funkcja5 * MM / (Q2 + m_pi * m_pi);
      break;			//Olga 2.12

    }

  *c3a = funkcja3;

  *c5a = funkcja5;

  *c4a = funkcja4;

  *c6a = funkcja6;

  return 0;
}


int tensor_hadronowy (int FFset, double delta_axial_mass, double delta_C5A, double E, double Q2, double W,
                      double *W1, double *W2, double *W3, double *W4, double *W5, double *W6)
/* Ta funkcja zwraca funkcje struktury (patrz notattki k.g.), ktore sa lorentzowsko niezmiennicze */
{

  double Md = masa_rezonansu;

  double Md2 = Md * Md;

  double qp = (Q2 + W * W - MM) / 2;

  double c3v, c4v, c5v, c6v;
  double c3a, c4a, c5a, c6a;

  form_faktory_wektorowe (FFset, Q2, W, W_Z, &c3v, &c4v, &c5v, &c6v);

  form_faktory_aksjalne (FFset, delta_axial_mass, delta_C5A, Q2, W, &c3a, &c4a, &c5a, &c6a);

// Wklady tensorowe z praca: O. Lalakulich, E. Paschos, Phys. ReV. D 71, 074003 (2005) //

  double c3v2 = c3v * c3v;
  double c4v2 = c4v * c4v;
  double c5a2 = c5a * c5a;
  double c4a2 = c4a * c4a;
  double c6a2 = c6a * c6a;


  double c3vc4v = c3v * c4v;
  double c4ac5a = c4a * c5a;
  double c3vc4a = c3v * c4a;
  double c3vc5a = c3v * c5a;
  double c4vc4a = c4v * c4a;
  double c4vc5a = c4v * c5a;
  double c4ac6a = c4a * c6a;
  double c5ac6a = c5a * c6a;


  double qp_M_2MdM = qp + MM - 2 * Masa * Md;

  double V1 =
    2.0 * (c3v2 / MM / Md2) * (pow (qp - Q2, 2) * (qp + MM) +
			       Md2 * (qp * qp + Q2 * MM + Q2 * Masa * Md)) +
    2.0 * (c4v2 / MM / MM) * pow (qp - Q2,
				  2) * (qp + MM - Masa * Md) +
    2.0 * (c3vc4v / MM / Masa / Md) * (qp - Q2) * ((qp - Q2) * qp_M_2MdM +
						   Md2 * qp) +
    2.0 * ((c4a2 / MM / MM) * pow (qp - Q2, 2) + c5a2 +
	   2 * (c4ac5a / MM) * (qp - Q2)) * (qp + MM + Masa * Md);

  *W1 = -V1 / 2 / Masa / Md;

  double V2 =
    2 * (c3v2 / Md2) * Q2 * (qp + MM + Md2) + 2 * (c4v2 / MM) * Q2 * (qp +
								      MM -
								      Masa *
								      Md) +
    2 * (c3vc4v / Masa / Md) * Q2 * (qp + pow (Md - Masa, 2)) +
    2 * (c5a2 * MM / Md2 + c4a2 * Q2 / MM) * (qp + MM + Masa * Md);

  *W2 = V2 / 2 / Masa / Md;

  double V3 =
    (4 / Md) * (-(c3vc4a / Masa) * (qp - Q2) - c3vc5a * Masa) * (2 * Md2 +
								 2 * Masa *
								 Md + Q2 -
								 qp) +
    4 * (qp - Q2) * ((-c4vc4a / MM) * (qp - Q2) - c4vc5a);

  *W3 = V3 / 2 / Masa / Md;


  double V4 = 2 * (c3v2 / Md2) * ((2 * qp - Q2) * (qp + MM) - Md2 * (MM + Masa * Md)) + 2 * (c4v2 / MM) * (2 * qp - Q2) * (qp + MM - Masa * Md) + 2 * (c3vc4v / Masa / Md) * ((2 * qp - Q2) * qp_M_2MdM + qp * Md2)	// <-- w ost czlonie nie ma razy dwa
    + 2 * (c5a2 * MM / Md2 + (c4a2 / MM) * (2 * qp - Q2) +
	   (c6a2 / MM / Md2) * (pow (Q2 - qp, 2) + Q2 * Md2) + 2 * c4ac5a -
	   2 * c4ac6a * qp / MM - 2 * (c5ac6a / Md2) * (Md2 + Q2 -
							qp)) * (qp + MM +
								Masa * Md);

  *W4 = V4 / Masa / Md / 2;

  double V5 =
    2 * (c3v2 / Md2) * qp * (qp + MM + Md2) + 2 * (c4v2 / MM) * qp * (qp +
								      MM -
								      Masa *
								      Md) +
    2 * (c3vc4v / Masa / Md) * qp * (qp + pow (Md - Masa, 2)) +
    2 * ((c4a2 / MM) * qp + c5a2 * MM / Md2 + c4ac5a - c4ac6a * Q2 / MM -
	 (c5ac6a / Md2) * (Q2 - qp)) * (qp + MM + Masa * Md);

  *W5 = V5 / Masa / Md;

  *W6 = 0;

  return 100;

}

double qw (double W)
{
  double l = pow (W * W + MM - m_pi * m_pi, 2) - 4 * W * W * MM;

  return sqrt (l) / 2 / W;
}


double dif_cross_q0_W (int FFset, double delta_axial_mass, double delta_C5A, double E, double q0, double W, double Meff)
{
  double Q2 = 2 * Meff * q0 + Meff * Meff - W * W;

  double mm2 = masa_leptonu * masa_leptonu;

  double W1, W2, W3, W4, W5, W6;

  double qp = (Q2 + W * W - MM) / 2;

  double kp = E * Masa;		// w ukladzie LAB

  masa_rezonansu = 1232 * MeV;

  tensor_hadronowy (FFset, delta_axial_mass, delta_C5A, E, Q2, W, &W1, &W2, &W3, &W4, &W5, &W6);

  double czynnik =
    (2 * W1 * (-Q2 - mm2) +
     2 * (W2 / MM) * (2 * kp * (kp - qp) - 0.5 * MM * (Q2 + mm2)) -
     2 * znak * (W3 / MM) * (Q2 * kp - 0.5 * qp * (Q2 + mm2)) +
     2 * (W4 / MM) * mm2 * (Q2 + mm2) / 2 -
     (W5 / MM) * mm2 * kp) * masa_rezonansu;
  // mnoze przez 2 Md aby odzyskac normalizacje Paschosa 


  double Gamma_0 = 120 * MeV;

  double Gamma_D =
    Gamma_0 * pow (qw (W) / qw (masa_rezonansu), 3) * masa_rezonansu / W;

  double GAMMA =
    Gamma_D * masa_rezonansu /
    (pow (W * W - masa_rezonansu * masa_rezonansu, 2) +
     masa_rezonansu * masa_rezonansu * Gamma_D * Gamma_D) / Pi;

  double przekroj_Q2_W =
    G * G * cos_2_theta_C * GAMMA * czynnik * W / Masa / E / E / 4 / Pi;

  if (W_limit * W_limit > W * W && W * W > pow (Masa + m_pi, 2))
    return 2 * Meff * przekroj_Q2_W; // d(Q2) = 2*Meff*d(omega)
  else
    return 0;

}


// Ponizsza ma za argumenty Energia=Energia neutrina, q0=przekaz energii k_0 - {k'}_O m_l =masa leptonu
// Zwraca ona 8 liczb
// Przekrojow czynnych NC, 
// 4 pierwsze: rozpraszanie neutrin: 
//             produkcja proton  pi_0 :  *NC_n_p_p0   
//             produkcja neutron pi_+ :  *NC_n_n_pi_plus  
//             produkcja neutron pi_0 :  *NC_n_n_p0
//             produkcja proton  po_- :  *NC_n_p_pi_minus
// 4 nastepne: rozpraszanie antyneutrin: 
//             produkcja proton  pi_0 :  *NC_an_p_p0   
//             produkcja neutron pi_+ :  *NC_an_n_pi_plus  
//             produkcja neutron pi_0 :  *NC_an_n_p0
//             produkcja proton  po_- :  *NC_an_p_pi_minus
int Przekroje_Czynne_q0_W_NC (int FFset, double delta_axial_mass, double delta_C5A,
                              double Energia, double q0, double W, double Meff, int FF,
                              double *NC_n_p_p0, double *NC_n_n_pi_plus, double *NC_n_n_p0, double *NC_n_p_pi_minus,
                              double *NC_an_p_p0, double *NC_an_n_pi_plus,double *NC_an_n_p0, double *NC_an_p_pi_minus)
{

  masa_rezonansu = 1232 * MeV;

  W_Z = NC;

  double czynnik_n, czynnik_an;

  znak = neutrino;

  masa_leptonu = 0;

  czynnik_n = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, Energia, q0, W, Meff);

  *NC_n_p_p0 = (2 / 9.) * czynnik_n;
  *NC_n_n_pi_plus = (1 / 9.) * czynnik_n;
  *NC_n_n_p0 = (2 / 9.) * czynnik_n;
  *NC_n_p_pi_minus = (1 / 9.) * czynnik_n;

  znak = antyneutrino;

  czynnik_an = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, Energia, q0, W, Meff);

  *NC_an_p_p0 = (2 / 9.) * czynnik_an;
  *NC_an_n_pi_plus = (1 / 9.) * czynnik_an;
  *NC_an_n_p0 = (2 / 9.) * czynnik_an;
  *NC_an_p_pi_minus = (1 / 9.) * czynnik_an;

  return 0;
}

// Ponizsza ma za argumenty Energia=Energia neutrina, q0=przsekaz energii k_0 - {k'}_O m_l =masa leptonu
// Zwraca ona 6 liczb
// Przekrojow czynnych CC, 
// 3 pierwsze: rozpraszanie neutrin: 
//             produkcja proton  pi_+ :  *CC_n_p_plus   
//             produkcja neutron pi_+ :  *CC_n_n_p_plus  
//             produkcja proton  pi_0 :  *CC_n_p_p0
// 3 nastepne: rozpraszanie antyneutrin: 
//             produkcja proton  pi_+ :  *CC_an_p_plus   
//             produkcja neutron pi_+ :  *CC_an_n_p_plus  
//             produkcja proton  pi_0 :  *CC_an_p_p0
int Przekroje_Czynne_q0_W_CC (int FFset, double delta_axial_mass, double delta_C5A,
                              double Energia, double q0, double W, double Meff, double m_l, int FF,
                              double *CC_n_p_plus, double *CC_n_n_p_plus, double *CC_n_p_p0,
                              double *CC_an_p_plus, double *CC_an_n_p_plus, double *CC_an_p_p0)
{

  masa_rezonansu = 1232 * MeV;

  W_Z = CC;
  masa_leptonu = m_l;

  znak = neutrino;

  double liczba_a = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, Energia, q0, W, Meff);

  *CC_n_p_plus = liczba_a;
  *CC_n_n_p_plus = (1 / 9.) * liczba_a;
  *CC_n_p_p0 = (2 / 9.) * liczba_a;


  znak = antyneutrino;

  double liczba_an = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, Energia, q0, W, Meff);

  *CC_an_p_plus = liczba_an;
  *CC_an_n_p_plus = (1 / 9.) * liczba_an;
  *CC_an_p_p0 = (2 / 9.) * liczba_an;

  return 0;
}


//function by JN
double cr_sec_delta (int FFset, double delta_axial_mass, double delta_C5A,
                     double E, double W, double nu, double Meff,
                     int lepton_in, int nukleon_in, int nukleon_out, int meson_out, bool current)
{

  double przekroj = 0, wynik = 0;

  masa_rezonansu = 1232 * MeV;

  switch (current)
    {
    case true:
      W_Z = CC;
      masa_leptonu = lepton_mass (abs (lepton_in), current);
      if (lepton_in > 0)
	{
	  znak = neutrino;
	}
      if (lepton_in < 0)
	{
	  znak = antyneutrino;
	}
      przekroj = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, E, nu, W, Meff);
////////////?????????????????????????????????????????????sprawdzic czy w porzadku sa reguly izospuinowe//////////////////////
//neutrino
      if (lepton_in > 0 && nukleon_in == proton && nukleon_out == proton
	  && meson_out == piplus)
	{
	  wynik = przekroj;
	}			//nu+proton --> mu + proton + piplus
      if (lepton_in > 0 && nukleon_in == neutron && nukleon_out == neutron
	  && meson_out == piplus)
	{
	  wynik = (1 / 9.) * przekroj;
	}			//nu+neutron --> mu + neutron + piplus
      if (lepton_in > 0 && nukleon_in == neutron && nukleon_out == proton
	  && meson_out == pizero)
	{
	  wynik = (2 / 9.) * przekroj;
	}			//nu+neutron --> mu + proton + pizero
//antyneutrino
      if (lepton_in < 0 && nukleon_in == neutron && nukleon_out == neutron
	  && meson_out == piminus)
	{
	  wynik = przekroj;
	}			//nubar+proton --> mubar + neutron + piminus
      if (lepton_in < 0 && nukleon_in == proton && nukleon_out == neutron
	  && meson_out == pizero)
	{
	  wynik = (2 / 9.) * przekroj;
	}			//nubar+neutron --> mubar + neutron + pizero
      if (lepton_in < 0 && nukleon_in == proton && nukleon_out == proton
	  && meson_out == piminus)
	{
	  wynik = (1 / 9.) * przekroj;
	}			//nu+neutron --> mu + proton + pizero
      break;

    case false:
      W_Z = NC;
      masa_leptonu = 0;
      if (lepton_in > 0)
	{
	  znak = neutrino;
	}
      if (lepton_in < 0)
	{
	  znak = antyneutrino;
	}
      przekroj = dif_cross_q0_W (FFset, delta_axial_mass, delta_C5A, E, nu, W, Meff);
      //neutrino & antyneutrino
      if (nukleon_in == proton && nukleon_out == proton
	  && meson_out == pizero)
	{
	  wynik = (2 / 9.) * przekroj;
	}			//nu+proton  --> nu + proton  + pizero
      if (nukleon_in == proton && nukleon_out == neutron
	  && meson_out == piplus)
	{
	  wynik = (1 / 9.) * przekroj;
	}			//nu+proton  --> nu + neutron + piplus
      if (nukleon_in == neutron && nukleon_out == neutron
	  && meson_out == pizero)
	{
	  wynik = (2 / 9.) * przekroj;
	}			//nu+neutron --> nu + neutron + pizero
      if (nukleon_in == neutron && nukleon_out == proton
	  && meson_out == piminus)
	{
	  wynik = (1 / 9.) * przekroj;
	}			//nu+neutron --> nu + proton  + piminus
      break;
    }

// cout<<"wynik="<<wynik*1e38/cm2;
  return wynik / cm2 * 1e38;
}


/*main()
{

cr_sec_delta(10000, 1400, 300, nu_mu, proton,proton,piplus,true);
cr_sec_delta(1000, 1400, 300, -14, 0,0,0,true);


}*/
