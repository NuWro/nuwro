#ifndef _resevent2_h_
#define _resevent2_h_

#include "event1.h"
#include "params.h"

//! in SPP language: 0 -> pi+, 1 -> pi0, 2 -> pi-
enum { pip, pi0, pim } spp_code;

//! cross section reduction to remove contribution from pion-less delta decay
double pdd_red(double energy);

void resevent2(params& p, event& e, bool cc);

#endif
