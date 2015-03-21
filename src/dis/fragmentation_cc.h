#ifndef _fragmentation_cc_h_
#define _fragmentation_cc_h_

//choosing the hit parton and spectator diquark
//written by J.Nowak


/////////////////////Hadronization for neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////

int hit_parton_cc_nu_p(double E, double W, double nu, double m);
int transfer_cc_nu_p(int hit_parton, double W, double nu);
int frag_parton_cc_nu_p(int hit_parton, int transfer);

int spectator_diquark1_cc_nu_p(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_cc_nu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);
void hadronization_cc_nu_p(double E, double W, double nu, double m);

/////////////////////Hadronization for anti-neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////
///////////////do sprawdzenia kat cabibo

int hit_parton_cc_anu_p(double E, double W, double nu, double m);
int transfer_cc_anu_p(int hit_parton, double W, double nu);
int frag_parton_cc_anu_p(int hit_parton, int transfer);

int spectator_diquark1_cc_anu_p(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_cc_anu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);


void hadronization_cc_anu_p(double E, double W, double nu, double m);





/////////////////////////NEUTRON/////////////////////////////////////////

int hit_parton_cc_nu_n(double E, double W, double nu, double m);
int transfer_cc_nu_n(int hit_parton, double W, double nu);
int frag_parton_cc_nu_n(int hit_parton, int transfer);

int spectator_diquark1_cc_nu_n(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_cc_nu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);
void hadronization_cc_nu_n(double E, double W, double nu, double m);

/////////////////////Hadronization for anti-neutrinos////////////////////////
/////////////////////////Neutron/////////////////////////////////////////
///////////////do sprawdzenia kat cabibo

int hit_parton_cc_anu_n(double E, double W, double nu, double m);
int transfer_cc_anu_n(int hit_parton, double W, double nu);
int frag_parton_cc_anu_n(int hit_parton, int transfer);

int spectator_diquark1_cc_anu_n(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_cc_anu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);


void hadronization_cc_anu_n(double E, double W, double nu, double m);


#endif
