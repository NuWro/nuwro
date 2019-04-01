#ifndef _mecevent_common_h_
#define _mecevent_common_h_

#include "pdg.h"
#include "particle.h"

// Size of the randomization region in energy transfer
extern double width_q0;

// Lepton mass
extern double ml;
extern double ml2;

// Reaction threshold or binding energy
extern double Bmec;

// Antinetrino
extern int ap;

// Nucleus heavier than nitrogen
extern int nucl;

// Pauli Blocking?
extern bool PB;

// How many repetitions in kinematic sampling
static int calls_max=50;

// Nucleon mass
const double MN = (PDG::mass_proton);// + PDG::mass_neutron) / 2.0;
const double MN2 = MN * MN;

////////////////////////////////////////

// This sets all final nucleon isospin behaviour (T.G.)
void mec_do_cc (particle *p, double ratio);

#endif
