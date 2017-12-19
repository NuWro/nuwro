#ifndef _resevent2_h_
#define _resevent2_h_

#include "event1.h"
#include "params.h"

//! cross section reduction to remove contribution from pion-less delta decay
double pdd_red(double energy);

//! get the value of binding energy according to setting from params
double get_binding_energy(const params& p, const vec& momentum);

void resevent2(params& p, event& e, bool cc);

#endif
