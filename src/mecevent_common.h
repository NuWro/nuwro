#ifndef _mecevent_common_h_
#define _mecevent_common_h_

#include "pdg.h"
#include "particle.h"

// Size of the randomization region in energy transfer
extern double width_q0;

// Cut in momentum transfer for Valencia model
const double qmax=1.2*GeV;

// Lepton mass
extern double ml;
extern double ml2;

// Efffecitive CC binding energies in MeV: Carbon-neutrino, Carbon-antineutrino
// Oxygen-neutrino, Oxygen-antineutrino
static double qvalues[6]={16.827*MeV,13.880*MeV,14.906*MeV,10.931*MeV,13.809*MeV,1.822*MeV};
//static double qvalues[6]={16.827*MeV,13.880*MeV,14.906*MeV,10.931*MeV,14.906*MeV,10.931*MeV};

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
