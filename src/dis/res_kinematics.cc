#include "res_kinematics.h"
#include "jednostki.h"

const double res_kinematics::Wmin = 1080;  // TODO: it is not exactly pion mass + nucleon mass
const double res_kinematics::avg_nucleon_mass = (PDG::mass_proton + PDG::mass_neutron) / 2.0;

res_kinematics::res_kinematics(const event &e) : neutrino(e.in[0]), target(e.in[1]) {
  // final lepton mass = 0 for NC or corresponding lepton mass (nu PDG - 1)
  lepton_mass = e.flag.cc * PDG::mass(abs(neutrino.pdg) - 1);
  lepton_mass2 = lepton_mass * lepton_mass;

  // subtract binding energy from nucleon energy inside nucleus
  target.t -= get_binding_energy(e.par, target.p());

  // boost to the bound nucleon rest frame
  neutrino.boost(-target.v());

  // effective nucleon mass depends on binding energy
  effective_mass = min(sqrt(target.p4() * target.p4()), avg_nucleon_mass);
  effective_mass2 = effective_mass * effective_mass;
}

double get_binding_energy(const params &p, const vec &momentum) {
  switch (p.nucleus_target) {
    case 0:  // free nucleon
      return 0;
    case 1:  // (global) Fermi gas
      return p.nucleus_E_b;
    case 2:  // local Fermi gas TODO: why 0?
      return 0;
    case 3:  // Bodek-Ritchie
      return 0;
    case 4:  // effective spectral function
      return binen(momentum, p.nucleus_p, p.nucleus_n);
    case 5:  // deuterium
      return deuter_binen(momentum);
    case 6:  // deuterium with constant binding energy
      return p.nucleus_E_b;
    default:
      return 0;
  }
}
