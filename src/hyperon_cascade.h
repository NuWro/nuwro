#include "particle.h"

#ifndef _hyperon_cascade_h_
#define _hyperon_cascade_h_

void hyperon_exp_xsec(double E, double Plab, double sigma[], int hyp_state);

void hyperon_state(int hyp_state, double sigma[], int &ij, particle p[]);

#endif
