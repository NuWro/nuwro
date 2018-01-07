#include "res_xsec.h"
#include "dis2res.h"
#include "dis_cr_sec.h"

res_xsec::res_xsec(res_kinematics &kin, const bool cc)
    : dis_pip(0), dis_pi0(0), dis_pim(0), delta_pip(0), delta_pi0(0), delta_pim(0), dis_spp(0), delta_spp(0) {
  // determine indices for SPP table (see singlepion.cc)
  j = kin.neutrino.pdg < 0;
  k = not cc;
  l = kin.target.pdg != PDG::pdg_proton;

  // total electric charge of the pion-nucleon system
  final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  // the contribution to the cross section coming from DIS
  from_dis = max(0.0, cr_sec_dis(kin.neutrino.E(), kin.W, kin.q.t, kin.neutrino.pdg, kin.target.pdg, cc));
}