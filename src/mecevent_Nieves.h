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

// Generate event
void mecevent_Nieves (params & p, event & e, nucleus & t, bool cc);

// Generate lepton kinematics and calculate the cross section
double Nieves_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t,
                              double mec_central, double mec_smearing, double binding,
                              int ile_pb, double sampling);

// Double-differential cross section
double dsdEdc (double E, double q0, double ct);

#endif
