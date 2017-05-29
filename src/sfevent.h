#ifndef _sfevent_h_
#define _sfevent_h_
#include "event1.h"
#include "nucleus.h"
#include "params.h"

//! check if SF exists for a given nucleus (FG will be used if False)
bool has_sf(nucleus &t, int method);

//! generate kinematics and calculate cross section for CC QEL event using SF
double sfevent2cc(params &p, event &e, nucleus &t);
//! generate kinematics and calculate cross section for NC EL event using SF
double sfevent2nc(params &p, event &e, nucleus &t);

//! generate kinematics and calculate cross section using SF
double sfevent(params &par, event &e, nucleus &t);

#endif
