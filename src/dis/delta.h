#ifndef _delta_h_
#define _delta_h_


int form_faktory_wektorowe (int Form_Faktory, double Q2, double W, int NC_czy_CC,
                            double *c3v, double *c4v, double *c5v, double *c6v);

int form_faktory_aksjalne (int Form_Faktory, double delta_axial_mass, double delta_C5A, double Q2, double W,
                           double *c3a, double *c4a, double *c5a, double *c6a);

int tensor_hadronowy (int Form_Faktory, double delta_axial_mass, double delta_C5A, double E, double Q2, double W,
                      double *W1, double *W2, double *W3, double *W4, double *W5, double *W6);

double qw (double W);

double dif_cross_q0_W (int Form_Faktory, double delta_axial_mass, double delta_C5A, double E, double q0, double W, double Meff);


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
int Przekroje_Czynne_q0_W_NC (int Form_Faktory, double delta_axial_mass, double delta_C5A,
                              double Energia, double q0, double W, double Meff, int FF,
                              double *NC_n_p_p0, double *NC_n_n_pi_plus, double *NC_n_n_p0, double *NC_n_p_pi_minus,
                              double *NC_an_p_p0, double *NC_an_n_pi_plus, double *NC_an_n_p0, double *NC_an_p_pi_minus);

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
int Przekroje_Czynne_q0_W_CC (int Form_Faktory, double delta_axial_mass, double delta_C5A,
                              double Energia, double q0, double W, double Meff, double m_l, int FF,
                              double *CC_n_p_plus, double *CC_n_n_p_plus, double *CC_n_p_p0,
                              double *CC_an_p_plus, double *CC_an_n_p_plus, double *CC_an_p_p0);

//function by JN
double cr_sec_delta (int FFset, double delta_axial_mass, double delta_C5A,
                     double E, double W, double nu, double Meff,
                     int lepton_in, int nukleon_in, int nukleon_out, int meson_out, bool current);

#endif
