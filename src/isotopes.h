#ifndef _ISOTOPE_H_
#define _ISOTOPE_H_

/*
 N    Z   A  EL    O     MASS EXCESS           BINDING ENERGY/A    BETA-DECAY ENERGY         ATOMIC MASS
                           (keV)                    (keV)                (keV)                (micro-u)
*/


struct isotope
{
	const int N;
	const int Z;
	const int A;
	const char* EL;
	const char*  O;
	const double masss_excess;           
	const double masss_excess_error;
	const double binding_energy;
	const double binding_energy_error;
	const double beta_decay_energy;
	const double beta_decay_energy_error;
	const int    atomic_mass;
	const double atomic_mass_error;
};

extern isotope isotopes[];
extern const int N_isotopes;
isotope * isotope_find(int p, int n);


#endif
