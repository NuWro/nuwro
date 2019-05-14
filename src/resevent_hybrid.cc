/*
This function calculates RES events from the Hybrid model
*/

#include "resevent_hybrid.h"
#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"
#include "particle.h"
#include "vect.h"
#include "dis/LeptonMass.h"


void resevent_hybrid(params &p, event &e, bool cc) {      // free nucleon only!
  e.weight = 0;              // if kinematically forbidden

  res_kinematics kin(e);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(0.)) return;

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  // final state particles
  particle final_pion, final_nucleon;

  final_pion.set_pdg_and_mass( PDG::pdg_piP );
  final_nucleon.set_pdg_and_mass( PDG::pdg_proton );

  // produces 4-momenta of final pair: nucleon + pion
  kin2part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion);

  // Omega_pi^* was chosen in hadronic CMS
  
  // calculate the cross section with only CMS, additional LAB factors come later
  e.weight = hybrid_dsdQ2dW(&kin);

  // boost back to LAB frame
  final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
  final_nucleon.p4() = final_nucleon.boost(kin.target.v());

  final_pion.p4() = final_pion.boost(kin.hadron_speed);
  final_pion.p4() = final_pion.boost(kin.target.v());

  // the cross section needs a factor (lepton momentum)^-2 in LAB
  // the cross section needs a jacobian: dw = dQ2/2M
  e.weight /= final_lepton.momentum2() * 2 * res_kinematics::avg_nucleon_mass;

  // save final state hadrons
  e.out.push_back(final_pion);
  e.out.push_back(final_nucleon);

  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
}

double hybrid_dsdQ2dW (res_kinematics *kin)
{
  double result = 0.;

  // double CS = 0;
  // for (int i = 0 ; i < 5 ;  i++)
  // {
  //   q11 = CS_table->table[index_QQ][index_W][i];
  //   q12 = CS_table->table[index_QQ+1][index_W][i];
  //   q21 = CS_table->table[index_QQ][index_W +1][i];
  //   q22 = CS_table->table[index_QQ+1][index_W+1][i];
  //
  //   CS += Lepton_factors[i]*(
  //   q11 * x2x_y2y +
  //   q21 * xx1_y2y +
  //   q12 * x2x_yy1 +
  //   q22 * xx1_yy1 );
  // }
  // CS/pow(El_inc,2) - E_lepton in LAB

  return result;
}
