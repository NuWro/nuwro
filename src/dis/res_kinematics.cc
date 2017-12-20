#include "res_kinematics.h"

const double res_kinematics::Wmin = 1080;  // TODO: it is not exactly pion mass + nucleon mass

res_kinematics::res_kinematics(const event& e) : neutrino(e.in[0]), target(e.in[1]) {
  // final lepton mass = 0 for NC or corresponding lepton mass (nu PDG - 1)
  lepton_mass = e.flag.cc * PDG::mass(abs(neutrino.pdg) - 1);
  lepton_mass2 = lepton_mass * lepton_mass;
}