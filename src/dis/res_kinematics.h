#ifndef _RES_KINEMATICS_H
#define _RES_KINEMATICS_H

#include "event1.h"

//! store kinematics used by RES
struct res_kinematics {
  static const double Wmin;  //!< invariant mass threshold

  res_kinematics(const event& e);  //!< initialize basic kinmatics

  particle neutrino;  //!< initial neutrino
  particle target;    //!< target nucleon

  double lepton_mass;  //!< outgoing lepton mass
  double lepton_mass2; //!< outgoing lepton mass squared
};

//! get the value of binding energy according to setting from params
double get_binding_energy(const params& p, const vec& momentum);

#endif