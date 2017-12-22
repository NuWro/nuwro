#ifndef _resevent2_h_
#define _resevent2_h_

#include <TPythia6.h>
#include "event1.h"
#include "params.h"

//! in SPP language: 0 -> pi+, 1 -> pi0, 2 -> pi-
enum { pip, pi0, pim } spp_code;

//! map PDG code to SPP code
inline int pdg2spp(const int pdg) {
  switch (pdg) {
    case PDG::pdg_piP:
      return pip;
    case PDG::pdg_pi:
      return pi0;
    case -PDG::pdg_piP:
      return pim;
    default:
      cerr << "[ERROR] not valid pion pdg";
      return -1;
  };
}

//! cross section reduction to remove contribution from pion-less delta decay
double pdd_red(double energy);

//! set up PYTHIA6 parameters
TPythia6* get_pythia();

void resevent2(params& p, event& e, bool cc);

#endif
