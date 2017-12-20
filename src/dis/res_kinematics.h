#ifndef _RES_KINEMATICS_H
#define _RES_KINEMATICS_H

#include "event1.h"

//! store kinematics used by RES
struct res_kinematics {
  static const double Wmin;              //!< invariant mass threshold
  static const double avg_nucleon_mass;  //!< average nucleon mass

  res_kinematics(const event& e);  //!< initialize basic kinmatics

  bool generate_kinematics(const double& res_dis_cut);  //!< set the rest of kinematics

  bool is_above_threshold();  //!< check it neutrino energy is above thresholds

  particle neutrino;  //!< initial neutrino
  particle target;    //!< target nucleon
  particle lepton;    //!< outgoing lepton (note: not in LAB frame)

  double lepton_mass;   //!< outgoing lepton mass
  double lepton_mass2;  //!< outgoing lepton mass squared

  double effective_mass;   //!< effective mass of a nucleon
  double effective_mass2;  //!< effective mass squared

  double W;   //!< invariant mass
  double W2;  //!< invariant mass squared

  vect q;  //!< four-momentum transfer

  vec hadron_speed;  //!< the speed of hadronic system

  double jacobian;  //!< integration over z
};

//! get the value of binding energy according to setting from params
double get_binding_energy(const params& p, const vec& momentum);

inline double pow2(double x) { return x * x; }

#endif