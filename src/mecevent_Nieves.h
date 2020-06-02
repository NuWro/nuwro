#ifndef _mecevent_Nieves_h_
#define _mecevent_Nieves_h_

#include <cmath>

#include "particle.h"
#include "jednostki.h"
#include "params.h"
#include "beam.h"
#include "pdg.h"
#include "nucleus.h"
#include "event1.h"

// Global variables and common functions in Jakub's implementation
#include "mecevent_common.h"

// Data file with hadronic tensor elements
#include "Nieves_MEC.h"
using namespace NSNWRO;

// Cut in momentum transfer for Valencia model
const double qmax_Nieves=1.2*GeV;

// Efffecitive CC binding energies in MeV: Carbon-neutrino, Carbon-antineutrino
// Oxygen-neutrino, Oxygen-antineutrino
static double qvalues_Nieves[6]={16.827*MeV,13.880*MeV,14.906*MeV,10.931*MeV,13.809*MeV,1.822*MeV};
//static double qvalues_Nieves[6]={16.827*MeV,13.880*MeV,14.906*MeV,10.931*MeV,14.906*MeV,10.931*MeV};

// Generate event
void mecevent_Nieves (params & p, event & e, nucleus & t, bool cc);

// Generate lepton kinematics and calculate the cross section
double Nieves_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t,
                              double mec_central, double mec_smearing, double binding,
                              int ile_pb, double sampling);

// Double-differential cross section
double Nieves_dsdEdc (double E, double q0, double ct);

#endif
