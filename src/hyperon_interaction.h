#include "vect.h"

#ifndef _hyperon_interaction_h_
#define _hyperon_interaction_h_

double Hyperon_Interaction(double Q2, double E_nu, int h, vect k1, vect p1, vect k2, vect p2,bool anti);

double Singh_Model(double Q2, double E_nu, int h , vect k1 , vect p1 , vect k2, vect p2,bool neu);

#endif
