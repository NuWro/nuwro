#include "res_xsec.h"
// #include "dis2res.h"
// #include "dis_cr_sec.h"

res_xsec::res_xsec(res_kinematics &kin, const bool cc) {
  // determine indices for SPP table (see singlepion.cc)
  j = kin.neutrino.pdg < 0;
  k = not cc;
  l = kin.target.pdg != PDG::pdg_proton;

  // total electric charge of the pion-nucleon system
  final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  // the contribution to the cross section coming from DIS
  // from_dis = max(0.0, cr_sec_dis(kin.neutrino.E(), kin.W, kin.q.t, kin.neutrino.pdg, kin.target.pdg, cc));
}