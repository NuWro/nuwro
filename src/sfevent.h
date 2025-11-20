#ifndef _sfevent_h_
#define _sfevent_h_
#include "event1.h"
#include "nucleus.h"
#include "params.h"

class CSpectralFunc;

// generate kinematics and calculate cross section using SF
double sfevent(params &par, event &e, nucleus &t);
// check if SF exists for a given nucleus (FG will be used if False)
bool has_sf(nucleus &t, int method);
// check if the interaction occurs on correlated pair
bool is_src(double p, double E, int Z, int N, bool is_on_p);
// return random energy shift according to distribution given by Fq
double random_omega();
// return Couloumb correction to neutron energy levels
double coulomb_correction_neutron(int p, int n);
// Approximate kinetic energy
double tPPrime_approx(const double eK, const double cosOfScattAngle, const bool EM, const bool NC, const double m_leptSq, const double m_nucl); 


#endif
