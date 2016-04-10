#ifndef _mecevent_Nieves_h_
#define _mecevent_Nieves_h_

//#include <cassert>
#include "particle.h"
#include "jednostki.h"
//#include <fstream>
//#include <sstream>
#include <cmath>
#include "params.h"
#include "beam.h"
#include "pdg.h"
#include "nucleus.h"
#include "event1.h"

//data file with response function
#include "Nieves_MEC.h"

//size of the randomization region in energy transfer
extern double width_q0;

//cut in momentum transfer for Valencia model
const double qmax=1.2*GeV;
//lepton mass
extern double ml;
extern double ml2;
//efffecitive CC binding energies in MeV: Carbon-neutrino, Carbon-antineutrino
//Oxygen-neutrino, Oxygen-antineutrino
static double qvalues[6]={16.827*MeV,13.880*MeV,14.906*MeV,10.931*MeV,13.809*MeV,1.822*MeV};
//reaction threshold or binding energy
extern double Bmec;
//antinetrino
extern int ap;
//nucleus hesvier, than nitrogen
extern int nucl;
//Pauli Blocking?
extern bool PB;

//double-differential cross section
double dsdEdc(double E, double q0, double ct);
//how many repetitions in kinematic sampling
static int calls_max=50;
const double MN = (PDG::mass_proton);// + PDG::mass_neutron) / 2.0;
const double MN2 = MN * MN;

//because the cross section is defined by muon kinematics.
double Nieves_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t, double mec_central, double mec_smearing, 
			      double binding, int ile_pb, double sampling);
//This sets all final nucleon behaviour. Thanks to T.G.
void Nieves_do_cc (particle *p, double ratio);

#endif
