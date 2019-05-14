#ifndef _resevent_hybrid_h_
#define _resevent_hybrid_h_

#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"

// Generate RES event
void resevent_hybrid(params& p, event& e, bool cc);

// Double-differential cross section in CMS without LAB factors
double hybrid_dsdQ2dW (res_kinematics* kin);

#endif
