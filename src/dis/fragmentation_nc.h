#ifndef _fragmentation_nc_h_
#define _fragmentation_nc_h_



//choosing the hit parton and spectator diquark
//written by J.Nowak


/////////////////////Hadronization for neutrinos////////////////////////
/////////////////////////PROTON/////////////////////////////////////////

int hit_parton_nc_nu_p(double E, double W, double nu, double m);
int transfer_nc_nu_p(int hit_parton);
int frag_parton_nc_nu_p(int hit_parton, int transfer);

int spectator_diquark1_nc_nu_p(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_nc_nu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);

void hadronization_nc_nu_p(double E, double W, double nu, double m);
///antineurtinos proton

int hit_parton_nc_anu_p(double E, double W, double nu, double m);
int transfer_nc_anu_p(int hit_parton);
int frag_parton_nc_anu_p(int hit_parton, int transfer);

int spectator_diquark1_nc_anu_p(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_nc_anu_p(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);
void hadronization_nc_anu_p(double E, double W, double nu, double m);


///////////NEUTRON///////////////////////

int hit_parton_nc_nu_n(double E, double W, double nu, double m);
int transfer_nc_nu_n(int hit_parton);
int frag_parton_nc_nu_n(int hit_parton, int transfer);
int spectator_diquark1_nc_nu_n(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_nc_nu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);
void hadronization_nc_nu_n(double E, double W, double nu, double m);
///antineurtinos proton

int hit_parton_nc_anu_n(double E, double W, double nu, double m);
int transfer_nc_anu_n(int hit_parton);
int frag_parton_nc_anu_n(int hit_parton, int transfer);

int spectator_diquark1_nc_anu_n(int hit_parton, int transfer, int frag_parton, double W);

int spectator_meson_nc_anu_n(int hit_parton, int transfer, int frag_parton, int spectator_diquark1, double W);
void hadronization_nc_anu_n(double E, double W, double nu, double m);


#endif
