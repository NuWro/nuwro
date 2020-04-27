#ifndef _lepevent_h_
#define _lepevent_h_

#include "params.h"
#include "event1.h"
#include "pdg.h"
#include "nu_e_el_sigma.h"
#include "jednostki.h"

// Generate LEP event
double lepevent(params& p, event& e); //, bool cc);

void kinemat_out(int kind_, int switch_sigma_, double &m_prime_, particle &lept_out_, particle &neut_out_);

#endif
