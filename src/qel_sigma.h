/////////////////////////////////////////////////////////////////////////////////////
// cala informacja o kwazielastycznym rozpraszaniu neutrin na nukleonach 
// na podstawie Llewelyn-Smith 
// (C) C. Juszczak 2001
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _QEL_h_
#define _QEL_h_

/// semielastic neutrino nucleon scattering cross section  (Llewelyn-Smith)
double qel_sigma (double Enu,double q2, int kind, bool anty, double m, double M);
 
#endif
