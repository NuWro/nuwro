#ifndef _sfevent_h_
#define _sfevent_h_
#include "event1.h"
#include "nucleus.h"
#include "params.h"

// TODO: can't include sf/CSpectralFunc.h because of double M2 global
class CSpectralFunc;

//! generate kinematics and calculate cross section using SF
double sfevent(params &par, event &e, nucleus &t);

//! check if SF exists for a given nucleus (FG will be used if False)
bool has_sf(nucleus &t, int method);

//! get removal energy (for given nucleon momentum p)
double get_E(CSpectralFunc *sf, double p);

//! check if the interaction occurs on correlated pair
bool is_src(double p, double E, int Z, int N, bool is_on_p);

//! return transparency (C) for given Q2 [GeV^2]
double transparency(double Q2);

//! return real part of the potential for given kinetic energy [MeV]
double potential_real(double Tk);

//! return random energy shift according to distribution given by Fq
double random_omega();

//! return Couloumb correction to final charged lepton energy (opposite sign for nu and nubar)
double coulomb_correction(bool is_anti, int p, int n);

//! return Couloumb correction to neutron energy levels
double coulomb_correction_neutron(int p, int n)


#endif
