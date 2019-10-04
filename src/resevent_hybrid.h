#ifndef _resevent_hybrid_h_
#define _resevent_hybrid_h_

#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"

// Cuts in Q2 and W for the hybrid model
const double Q2max_hybrid=1.91*GeV2;
const double Q2min_hybrid=0.01*GeV2;
const double Q2spc_hybrid=0.10*GeV2;
const int Q2bin_hybrid=(Q2max_hybrid-Q2min_hybrid)/Q2spc_hybrid+1;
const double  Wmax_hybrid=1495*MeV;
const double  Wmin_hybrid=1000*MeV;
const double  Wspc_hybrid=   5*MeV;
const int Wbin_hybrid=(Wmax_hybrid-Wmin_hybrid)/Wspc_hybrid+1;

// Generate RES event
void resevent_hybrid(params& p, event& e, bool cc);

// Double-differential cross section in CMS without LAB factors
// In terms of tabularized hadronic responses
double hybrid_dsdQ2dW(res_kinematics* kin, int channel);

// Four-fold differential cross section in CMS without LAB factors
// In terms of ABCDE decomposition
double hybrid_dsdQ2dWdOm(res_kinematics* kin, int channel, vect final_pion);

#endif