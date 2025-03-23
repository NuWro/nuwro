#ifndef _mecevent_2020Valencia_h_
#define _mecevent_2020Valencia_h_


#include <cmath>

#include "particle.h"
#include "jednostki.h"
#include "params.h"
#include "beam.h"
#include "pdg.h"
#include "nucleus.h"
#include "event1.h"
#include "mecevent_Nieves.h"

// Global variables and common functions in Jakub's implementation
#include "mecevent_common.h"


// Efffecitive CC binding energies in MeV: Carbon-neutrino, Carbon-antineutrino
// Oxygen-neutrino, Oxygen-antineutrino
static double E_corr[6]={16.827*MeV,13.878*MeV,14.906*MeV,10.931*MeV,13.809*MeV,1.822*MeV};


// Flags for event topolgy in Valencia 2020 
extern bool flag_2p2h;
extern bool flag_2p2h_pn;
extern bool flag_3p3h;


int ConvertCoordinate_2D_to_1D(int n, int m);
double Bilinear_Interpolation(double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22);
void mecevent_2020Valencia (params & p, event & e, nucleus & t, bool cc);
double Valencia2020_kin_and_weight_2p2h (double E, double *individual_dsdqdw, particle &meclep, particle *inc_nucleon_2p2h, particle *out_nucleon_2p2h, bool &flag_2p2h_pn, nucleus &t, double mec_central, double mec_smearing, double binding, int ile_PB, double* sampling, int* strength);
double Valencia2020_kin_and_weight_3p3h ( double E, double* individual_dsdqdw, particle &meclep, particle *inc_nucleon_3p3h, particle *out_nucleon_3p3h, nucleus &t, double binding, int ile_PB);
double Valencia2020_dsdEdc(double E, int outgoing_pair, double q0, double Ep, double ct, nucleus &T);
void Isospin_model_2p2h_2020Valencia (particle *in_p, particle *out_p, double ratio_pp, double ratio_np, bool &flag_2p2h_pn);
void Isospin_model_3p3h_2020Valencia (particle *in_p, particle *out_p, nucleus &T_nucleus);
void Generate_nucleon_kinematics_2p2h (particle *inc_nucleon_mec, particle *out_nucleon_mec, bool &flag_2p2h_pn, nucleus &T, double &xsec_result, vect qqq, double mec_central, double mec_smearing, double binding, int ile_PB, double* sampling, int* strength);
void Generate_nucleon_kinematics_3p3h( particle *inc_nucleon_mec, particle *out_nucleon_mec, nucleus &T, double &xsec_result, vect qqq, int ile_PB);
void threebodydecay (double W, double m1, double m2, double m3, vect &p1, vect &p2, vect &p3);


#endif
