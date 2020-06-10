#ifndef E_EL_EVENT_H
#define E_EL_EVENT_H

#include "params.h"
#include "nucleus.h"
#include "event1.h"

using namespace NUWRO;

double e_el_event(params&p, event & e, nucleus &t, bool nc);
double e_el_event2(params&p, event & e, nucleus &t, bool nc);
double e_el_event2orig(params&p, event & e, nucleus &t, bool nc);

#endif
