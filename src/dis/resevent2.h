#ifndef _resevent2_h_
#define _resevent2_h_

#include <TPythia6.h>
#include "event1.h"
#include "params.h"
#include "res_kinematics.h"
#include "nucleus.h"

//! in SPP language: 0 -> pi+, 1 -> pi0, 2 -> pi-
typedef enum { pip, pi0, pim } spp_code;

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

//! set up PYTHIA6 parameters
TPythia6* get_pythia();

//! get id-th particle from pythia's particle list
particle get_pythia_particle(Pyjets_t* pythia_particles, const int particle_id, res_kinematics kin);

//! save Pythia particles
void save_pythia_particles(event& e, Pyjets_t* pythia_particles, const int nof_particles, const res_kinematics& kin);

//! generate RES event
void resevent2(params& p, event& e, nucleus& t, bool cc);

#endif
