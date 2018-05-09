#ifndef _crossection_h_
#define _crossection_h_

#include "DM.h"
//#include "units.h"
#include "../jednostki.h"
#include "constants.h"
#include "util2.h"
#include "calga.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                 CHANNELS:                          ///////////////////////////////
///////     1: e + p -> e + p + \pi^\0                                     ///////////////////////////////
///////     2: e + n -> e + p + \pi^-                                      ///////////////////////////////
///////     3: e + p -> e + n + \pi^+                                      ///////////////////////////////
///////     4: e + n -> e + n + \pi^0                                      ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

extern int chan; // set to choose the channel you want
extern int PV;  //
extern int FP; //
//how many integration subregions
extern int no1; //precision of integration
extern int ffset;//which form factors to use
extern bool selfenergy;

double dsigma_domega_dEprime_depi2_(double Epi);
//double dsigma_domega_dEprime_depi_(double Epi);
double dsigma_domega_dEprime_(double E, double q0, double cosine, double *CV);
double dsigma_dQ2_dW_(double E, double Q2, double W, double *CV);
//laboratory electron energy, outgoing electron 4-momentum, hadronic CMS pion solid angle w.r.t hadronic cms delta direction

//Done in laboratory frame: electron coming along the z-axis, any outgoing electron angle and energy any direction of nucleon
double dsigma_dq0_dOmegal_dOmegaCMS(double E, double mass, D4V<double> lprimemu, double costhetacms, double phicms, D4V<double> pn);

//Done in laboratory frame: electron coming along the z-axis, any outgoing electron angle and energy any direction of nucleon+ Oset Delta selfenergy
//with additional fermi momentum and density
double dsigma_dq0_dOmegal_dOmegaCMS_Oset(double E, double mass, D4V<double> lprimemu, double costhetacms, double phicms, D4V<double> pn, double KF, double rhorho0);
#endif
