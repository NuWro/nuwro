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
#include "hybrid_RES.h"
#include "hybrid/hybrid_gateway.h"


void resevent_hybrid(params &p, event &e, bool cc) {      // free nucleon only!
  e.weight = 0;           // if kinematically forbidden

  res_kinematics kin(e);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(1500)) return;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  // final state particles
  particle final_pion, final_nucleon;

  // selection of the final state

  // specify the total electric charge of the pion-nucleon system
  int j = kin.neutrino.pdg < 0;
  int k = not cc;
  int l = kin.target.pdg != PDG::pdg_proton;
  int final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  double xsec_pip=0, xsec_pi0=0, xsec_pim=0, xsec_inclusive=0; // cross sections

  switch (final_charge) { // calculate the cross section with only the CMS variables
    case  2:  // pi+ + proton (nu_11)
      xsec_pip = hybrid_dsdQ2dW(&kin, 11);
      xsec_inclusive = xsec_pip;
      if ( not (xsec_inclusive > 0) ) return;
      final_pion.set_pdg_and_mass( PDG::pdg_piP );
      final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
      break;
    case  1:  // pi+ + neutron (nu_22) or pi0 + proton (nu_21)
      xsec_pip = hybrid_dsdQ2dW(&kin, 22);
      xsec_pi0 = hybrid_dsdQ2dW(&kin, 21);
      xsec_inclusive = xsec_pip + xsec_pi0;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pip / xsec_inclusive > frandom() ) // random selection
      {
        final_pion.set_pdg_and_mass( PDG::pdg_piP );
        final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
      }
      else
      {
        final_pion.set_pdg_and_mass( PDG::pdg_pi );
        final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
      }
      break;
    case  0:  // pi0 + neutron (anu_11) or pi- + proton (anu_12)
      xsec_pi0 = hybrid_dsdQ2dW(&kin, 21); // anu_11 has the tables of nu_21
      xsec_pim = hybrid_dsdQ2dW(&kin, 22); // anu_12 has the tables of nu_22
      xsec_inclusive = xsec_pi0 + xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pi0 / xsec_inclusive > frandom() ) // random selection
      {
        final_pion.set_pdg_and_mass( PDG::pdg_pi );
        final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
      }
      else
      {
        final_pion.set_pdg_and_mass( -PDG::pdg_piP );
        final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
      }
      break;
    case -1:  // pi- + neutron (anu_22)
      xsec_pim = hybrid_dsdQ2dW(&kin, 11); // anu_22 has the tables of nu_11
      xsec_inclusive = xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      final_pion.set_pdg_and_mass( -PDG::pdg_piP );
      final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
      break;
    default:
      cerr << "[WARNING]: Reaction charge out of range\n";
  };

  // produces 4-momenta of final pair: nucleon + pion
  kin2part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion);

  // Omega_pi^* was chosen in hadronic CMS

  // 2d cross section only
  e.weight = xsec_inclusive;

  // use 3d cross section
  e.weight = hybrid_dsdQ2dWdcth(&kin, 11, final_pion);

  // use 4d cross section
  e.weight = hybrid_dsdQ2dWdOm(&kin, 11, final_pion);

  // the cross section needs a factor (initial lepton momentum)^-2 in LAB
  // the cross section needs a jacobian: dw = dQ2/2M
  e.weight /= e.in[0].E() * e.in[0].E() / 2 / res_kinematics::avg_nucleon_mass / kin.jacobian;
  // coupling
  e.weight *= G*G*cos2thetac/2;
  // units
  e.weight /= cm2;

  // boost back to LAB frame
  final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
  final_nucleon.p4() = final_nucleon.boost(kin.target.v());

  final_pion.p4() = final_pion.boost(kin.hadron_speed);
  final_pion.p4() = final_pion.boost(kin.target.v());

  // save final state hadrons
  e.out.push_back(final_pion);
  e.out.push_back(final_nucleon);

  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
}

double hybrid_dsdQ2dW(res_kinematics *kin, int channel)
{
  double *hybrid_grid;
  switch (channel)
  {
    case 11:
      hybrid_grid = hybrid_grid_11;
      break;

    case 21:
      hybrid_grid = hybrid_grid_21;
      break;

    case 22:
      hybrid_grid = hybrid_grid_22;
      break;

    default:
      cerr << "[WARNING]: no hybrid inclusive grid!\n";
  }

  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // we build leptonic tensor in CMS with q along z, see JES paper App. A
  vect kl_inc_lab = kin->neutrino.p4();   //
  vect kl_lab     = kin->lepton.p4();     //  They are not in LAB, but in target rest frame!
  vect q_lab      = kin->q;               //

  double v = q_lab.length() / (q_lab[0] + res_kinematics::avg_nucleon_mass);
  double g = 1 / sqrt(1 - v*v);
  double c = (kl_inc_lab.length() - kl_lab[3]) / q_lab.length();
  double s = sqrt(pow(kl_lab.length(),2) - pow(kl_lab[3],2)) / q_lab.length();

  vect kl_inc (g*(kl_inc_lab[0] - v*kl_inc_lab.length()*c), kl_inc_lab.length()*s,
               0, g*(-v*kl_inc_lab[0] + kl_inc_lab.length()*c));
  vect kl     (g*(kl_lab[0] - v*(kl_inc_lab.length()*c - q_lab.length())), kl_inc_lab.length()*s,
               0, g*(-v*kl_lab[0] + kl_inc_lab.length()*c - q_lab.length()));

  double kl_inc_dot_kl = kl_inc * kl;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // interpolating the nuclear tensor elements
  double w00=0;
  double w03=0;
  double w33=0;
  double w11=0;
  double w12=0;

  //if( Q2 < Q2min_hybrid ) Q2 = Q2min_hybrid;
  //if( Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  if( Q2 >= Q2min_hybrid && Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  {
    int m=int((Q2-Q2min_hybrid)/Q2spc_hybrid);
    int n=int((W-Wmin_hybrid)/Wspc_hybrid);
    int pos=(n*Q2bin_hybrid+m)*5;

    int a=pos;
    int b=a+5;
    int c=pos+Q2bin_hybrid*5;
    int d=c+5;

    double H1=W-Wmin_hybrid-n*Wspc_hybrid;
    double H2=Q2-Q2min_hybrid-m*Q2spc_hybrid;

    w00=(H2*(H1*hybrid_grid[d]+(Wspc_hybrid-H1)*hybrid_grid[b])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c]+(Wspc_hybrid-H1)*hybrid_grid[a])/Wspc_hybrid)/Q2spc_hybrid;
    w03=(H2*(H1*hybrid_grid[d+1]+(Wspc_hybrid-H1)*hybrid_grid[b+1])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+1]+(Wspc_hybrid-H1)*hybrid_grid[a+1])/Wspc_hybrid)/Q2spc_hybrid;
    w33=(H2*(H1*hybrid_grid[d+2]+(Wspc_hybrid-H1)*hybrid_grid[b+2])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+2]+(Wspc_hybrid-H1)*hybrid_grid[a+2])/Wspc_hybrid)/Q2spc_hybrid;
    w11=(H2*(H1*hybrid_grid[d+3]+(Wspc_hybrid-H1)*hybrid_grid[b+3])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+3]+(Wspc_hybrid-H1)*hybrid_grid[a+3])/Wspc_hybrid)/Q2spc_hybrid;
    w12=(H2*(H1*hybrid_grid[d+4]+(Wspc_hybrid-H1)*hybrid_grid[b+4])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+4]+(Wspc_hybrid-H1)*hybrid_grid[a+4])/Wspc_hybrid)/Q2spc_hybrid;
  }

  double l00 = (2.*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  double l03 = (-kl_inc[0]*kl[3] + -kl[0]*kl_inc[3]);
  double l33 = (2*kl_inc[3]*kl[3] + kl_inc_dot_kl);
  double l11 = (2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  double l22 = (2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  double l12 = (kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);

  result = l00*w00 + 2*l03*w03 + l33*w33 + 0.5*(l11+l22)*w11 - 2*l12*w12;

  return result;
}

double hybrid_dsdQ2dWdcth(res_kinematics* kin, int channel, vect final_pion)
{
  // placeholders
  double costh[1];
  int    params[4];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // specify the params (note strange order!)
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (kin->neutrino.pdg > 0)); // helicity
  params[3] = int(channel/10);                     // nucleon
  params[1] = channel % 10;                        // decay

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result = ABCDE[0][0];
  result *= pion_momentum / pow(2*Pi,4);
  result *= 4 * Pi; // Phase space!

  return result;
}

double hybrid_dsdQ2dWdOm(res_kinematics* kin, int channel, vect final_pion)
{
  // placeholders
  double costh[1];
  int    params[4];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // specify the params (note strange order!)
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (kin->neutrino.pdg > 0)); // helicity
  params[3] = int(channel/10);                     // nucleon
  params[1] = channel % 10;                        // decay

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result = ABCDE[0][0] + ABCDE[0][1]*cos(pion_phi) + ABCDE[0][2]*cos(2*pion_phi)
                       + ABCDE[0][3]*sin(pion_phi) + ABCDE[0][4]*sin(2*pion_phi);
  result *= pion_momentum / pow(2*Pi,4);
  result *= 4 * Pi; // Phase space!

  return result;
}
