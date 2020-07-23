#ifndef _mecevent_SuSA_h_
#define _mecevent_SuSA_h_

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
#include "SuSA_MEC.h"
using namespace NUWRO;

// Cut in momentum transfer for SuSA model
const double qmax_SuSA=2*GeV;

// Efffecitive CC binding energies in MeV: Carbon-neutrino, Carbon-antineutrino
// Oxygen-neutrino, Oxygen-antineutrino
static double qvalues_SuSA[4]={40.*MeV,40.*MeV,32.*MeV,32.*MeV};

// Generate event
void mecevent_SuSA (params & p, event & e, nucleus & t, bool cc);

// Generate lepton kinematics and calculate the cross section
double SuSA_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t,
                            double mec_central, double mec_smearing, double binding,
                            int ile_pb, double sampling);

// Double-differential cross section
double SuSA_dsdEdc (double E, double q0, double ct);

#endif
