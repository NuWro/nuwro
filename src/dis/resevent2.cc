// clang-format off
////////////////This function calculates RES event
////////////////RES region is defined in the hadronic mass from 1080 to the parameter res_dis_cut introduced in params.txt
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////The following model of the #Sect is postulated
///////// XSect (nu + N -> mu + N' + pi)/dW d omega  = [XSect (nu + N -> mu + N' + pi)/dW d omega ]_Delta * alfadelta (W) +
//////////                                             [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_SPP* SPP(W) * alfados(W) +
//////////                                             [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_nonSPP
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// The main idea behind is that DIS contribution simulates non-resonant part
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// The functions alfa and beta are a priori arbitrary but they have to satisfy the conditions
/////////// alfadis(res_dis_cut)=1, ....alfadelta(res_dis_cut=0,............alfadis (1080)=0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////alfadis and alfadelta are defined in alfa.cc
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////In principle for each channel the functions should be different. They should be fitted to experimental data
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////Channels are labelled by 4 integers corresponding to: nu/anu   cc/nc   proton/neutron    pi+/pi0/pi-
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// clang-format on

#include <TMCParticle.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "LeptonMass.h"
#include "alfa.h"
#include "charge.h"
#include "delta.h"
#include "dis2res.h"
#include "dis_cr_sec.h"
#include "event1.h"
#include "fragmentation.h"
#include "grv94_bodek.h"
#include "jednostki.h"
#include "masses.h"
#include "parameters.h"
#include "params.h"
#include "pauli.h"
#include "pdg_name.h"
#include "resevent2.h"
#include "res_kinematics.h"
#include "singlepion.h"
#include "vect.h"

extern "C" int pycomp_(const int *);
extern double SPP[2][2][2][3][40];

double pdd_red(double energy) {
  if (energy >= 1000)
    return 0.85;
  else if (energy > 750)
    return 0.65 + energy * 0.05 / 250.0;
  else  // if (en<=750)
    return 0.2 + energy * 0.2 / 250.0;
}

TPythia6 *get_pythia() {
  TPythia6 *pythia71 = new TPythia6();

  //////////////////////////////////////////////
  //      Setting Pythia parameters
  //      Done by Jaroslaw Nowak
  //////////////////////////////////////////////

  // stable pi0
  pythia71->SetMDCY(pycomp_(&pizero), 1, 0);

  pythia71->SetMSTU(20, 1);  // advirsory warning for unphysical flavour switch off
  pythia71->SetMSTU(23, 1);  // It sets counter of errors at 0
  pythia71->SetMSTU(26, 0);  // no warnings printed

  // PARJ(32)(D=1GeV) is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet
  // parton system
  pythia71->SetPARJ(33, 0.1);

  // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below
  // which the fragmentation of a parton system is stopped and two final hadrons formed.
  pythia71->SetPARJ(34, 0.5);
  pythia71->SetPARJ(35, 1.0);

  // PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point
  // of the fragmentation. Strongly corlated with PARJ(33-35)
  pythia71->SetPARJ(37, 1.);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  // MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and
  // thus allow a cluster to decay rather than collaps
  pythia71->SetMSTJ(18, 3);  // do not change

  //////////////////////////////////////////////
  //      End of setting Pythia parameters
  //////////////////////////////////////////////

  return pythia71;
}

void resevent2(params &p, event &e, bool cc) {
  e.weight = 0;  // if kinmetically forbidden

  res_kinematics kin(e);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(p.res_dis_cut)) return;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  // determine indices for SPP table (see singlepion.cc)
  const int j = kin.neutrino.pdg < 0;
  const int k = not cc;
  const int l = kin.target.pdg != PDG::pdg_proton;

  // total electric charge of the pion-nucleon system
  const int finalcharge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  // the contribution to the cross section coming from DIS
  const double fromdis = max(0.0, cr_sec_dis(kin.neutrino.E(), kin.W, kin.q.t, kin.neutrino.pdg, kin.target.pdg, cc));

  if (not kin.is_above_pythia_threshold() || fromdis == 0) {
    // contributions from DIS
    const double dis_pip = fromdis * SPP[j][k][l][pip][0] * betadis(j, k, l, pip, kin.W, p.bkgrscaling);
    const double dis_pi0 = fromdis * SPP[j][k][l][pi0][0] * betadis(j, k, l, pi0, kin.W, p.bkgrscaling);
    const double dis_pim = fromdis * SPP[j][k][l][pim][0] * betadis(j, k, l, pim, kin.W, p.bkgrscaling);

    // contributions from Delta
    double delta_pip = 0, delta_pi0 = 0, delta_pim = 0;

    switch (finalcharge) {
      case 2:  // pi+ + proton
        delta_pip = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_proton, PDG::pdg_piP, cc) *
                    alfadelta(j, k, l, pip, kin.W);
        break;
      case 1:  // pi+ + neutron or pi0 + proton
        delta_pip = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_neutron, PDG::pdg_piP, cc) *
                    alfadelta(j, k, l, pip, kin.W);
        delta_pi0 = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_proton, PDG::pdg_pi, cc) *
                    alfadelta(j, k, l, pi0, kin.W);
        break;
      case 0:  // pi0 + neutron or pi- + proton
        delta_pi0 = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_neutron, PDG::pdg_pi, cc) *
                    alfadelta(j, k, l, pi0, kin.W);
        delta_pim = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_proton, -PDG::pdg_piP, cc) *
                    alfadelta(j, k, l, pim, kin.W);
        break;
      case -1:  // pi- + neutron
        delta_pim = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                 kin.neutrino.pdg, kin.target.pdg, PDG::pdg_neutron, -PDG::pdg_piP, cc) *
                    alfadelta(j, k, l, pim, kin.W);
        break;
      default:
        cerr << "[WARNING]: charge out of rangen\n";
    };

    // we arrived at the overall strength !!!
    double total = dis_pip + dis_pi0 + dis_pim + delta_pip + delta_pi0 + delta_pim;

    // save cross section in appropriate units
    e.weight = total * 1e-38 * kin.jacobian;

    // final state particles
    particle final_pion, final_nucleon;

    // contributions from different pions to xsec
    const double pip_fraction = (dis_pip + delta_pip) / total;
    const double pi0_fraction = (dis_pi0 + delta_pi0) / total;
    const double pim_fraction = (dis_pim + delta_pim) / total;

    // randomly select final state pion
    double rand01 = frandom();

    if (pip_fraction > rand01)
      final_pion.set_piP();
    else if (pip_fraction + pi0_fraction > rand01)
      final_pion.set_pi();
    else
      final_pion.set_piM();

    // determine isospin of a final nucleon
    if (nukleon_out_(kin.W, kin.neutrino.pdg, kin.target.pdg, final_pion.pdg, cc) == PDG::pdg_proton)
      final_nucleon.set_proton();
    else
      final_nucleon.set_neutron();

    // produces 4-momenta of final pair: nucleon + pion
    kin2part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion);

    // boost back to LAB frame
    final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
    final_nucleon.p4() = final_nucleon.boost(kin.target.v());

    final_pion.p4() = final_pion.boost(kin.hadron_speed);
    final_pion.p4() = final_pion.boost(kin.target.v());

    // save final state hadrons
    e.out.push_back(final_pion);
    e.out.push_back(final_nucleon);
  } else  // the algorithm starts from the production of PYTHIA event
  {
    TPythia6 *pythia71 = get_pythia();

    int nof_particles = 0;     // number of particles in the final state
    Pyjets_t *pythiaParticle;  // pythia particles placeholder

    // force less than 5 particles in the final state (lepton, target, nucleon, pion?)
    while (nof_particles < 5) {
      hadronization(kin.neutrino.E(), kin.W, kin.q.t, kin.lepton_mass, kin.neutrino.pdg, kin.target.pdg, cc);
      pythiaParticle = pythia71->GetPyjets();
      nof_particles = pythia71->GetN();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////There are three different possible outcomes:
    ///////////     a) spp event nof_particles=5
    ///////////     b) more inelastic event nof_particles>5, typically 7
    ///////////     c) single kaon production; this causes technical complications because nof_particles=5 also in this
    /// case
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (nof_particles == 5 && (pythiaParticle->K[1][3] == 211 || pythiaParticle->K[1][4] == 211 ||
                               pythiaParticle->K[1][3] == 111 || pythiaParticle->K[1][4] == 111 ||
                               pythiaParticle->K[1][3] == -211 || pythiaParticle->K[1][4] == -211))  // spp condition
    {
      int t;
      int pion_pdg;

      if (pythiaParticle->K[1][3] == 211 || pythiaParticle->K[1][4] == 211)  // the second part
      {
        t = 0;
        pion_pdg = 211;
      }

      if (pythiaParticle->K[1][3] == 111 || pythiaParticle->K[1][4] == 111) {
        t = 1;
        pion_pdg = 111;
      }

      if (pythiaParticle->K[1][3] == -211 || pythiaParticle->K[1][4] == -211) {
        t = 2;
        pion_pdg = -211;
      }
      // cout<<pythiaParticle->K[1][3]<<" "<< pythiaParticle->K[1][4]<<endl;
      double dis_spp = fromdis * betadis(j, k, l, t, kin.W, p.bkgrscaling);  // dis contribution

      int nukleoncharge = finalcharge + t - 1;  // the charge of final nucleon

      int nukleon2;

      if (nukleoncharge == 1) {
        nukleon2 = 2212;
      } else {
        nukleon2 = 2112;
      }

      // delta contribution
      // if (SPPF (j,k,l,t,W)<0.01)
      // cout<<SPPF (j,k,l,t,W)<<" "<<j<<" "<<k<<" "<<l<<" "<<t<<" "<<W<<" "<<endl;

      double delta_spp = cr_sec_delta(p.delta_FF_set, p.pion_axial_mass, p.pion_C5A, kin.neutrino.E(), kin.W, kin.q.t,
                                      kin.neutrino.pdg, kin.target.pdg, nukleon2, pion_pdg, cc) /
                         SPPF(j, k, l, t, kin.W) * alfadelta(j, k, l, t, kin.W);
      // cout<<delta_spp<<endl;

      // approximate implementation of pionless delta decays
      if ((p.nucleus_p + p.nucleus_n) > 7) {  // cout<<delta_spp<<"  ";
        // double ennergy = e.in[0].t;
        // double rescale = pdd_red (ennergy);
        delta_spp *= pdd_red(e.in[0].t);
        // cout<<ennergy<<"  "<<rescale<<"  "<<delta_spp<<endl;
      }
      // approximate implementation of pionless delta decays

      double spp_strength = dis_spp + delta_spp;
      // cout<<"spp "<<W<<" "<<nu<<" "<<spp_strength<<endl;
      e.weight = spp_strength * 1e-38 * kin.jacobian;

      double reldis = dis_spp / spp_strength;
      double reldelta = delta_spp / spp_strength;

      double los = frandom();

      if (reldis > los)  // disevent
      {
        for (int i = 0; i < nof_particles; i++) {
          particle part;
          part.t = pythiaParticle->P[3][i] * GeV;
          part.x = pythiaParticle->P[0][i] * GeV;
          part.y = pythiaParticle->P[1][i] * GeV;
          part.z = pythiaParticle->P[2][i] * GeV;
          rotation(part, kin.q);

          part = part.boost(kin.hadron_speed);  // correct direction ???
          part = part.boost(kin.target.v());

          part.ks = pythiaParticle->K[0][i];
          part.pdg = pythiaParticle->K[1][i];
          part.orgin = pythiaParticle->K[2][i];

          e.temp.push_back(part);
          if (part.ks == 1)  // condition for a real particle in the final state
          {
            e.out.push_back(part);
          }
        }
      } else  // deltaevent
      {
        vect finnuk, finpion;

        kin.neutrino.boost(-kin.hadron_speed);  // a boost from nu-N CMS to the hadronic CMS
        kin.lepton.boost(-kin.hadron_speed);    // a boost from nu-N CMS to the hadronic CMS
        kin4part(kin.neutrino, kin.lepton, kin.W, nukleon2, pion_pdg, finnuk, finpion,
                 p.delta_angular);  // produces 4-momenta of final pair: nucleon + pion with density matrix information
        e.weight *= angrew;         // reweight according to angular correlation

        // kin2part (W, nukleon2, pion, finnuk, finpion);	//produces 4-momenta of the final pair: nucleon + pion

        kin.neutrino.boost(kin.hadron_speed);  // a boost back to the nu-N CMS frame
        kin.neutrino.boost(kin.target.v());    // a boost back to tha LAB frame

        kin.lepton.boost(kin.hadron_speed);  // a boost back to the nu-N CMS frame
        kin.lepton.boost(kin.target.v());    // a boost back to tha LAB frame

        finnuk = finnuk.boost(kin.hadron_speed);
        finnuk = finnuk.boost(kin.target.v());

        finpion = finpion.boost(kin.hadron_speed);
        finpion = finpion.boost(kin.target.v());

        particle ppion(finpion);
        particle nnukleon(finnuk);

        ppion.pdg = pion_pdg;
        nnukleon.pdg = nukleon2;

        e.out.push_back(ppion);
        e.out.push_back(nnukleon);
      }

    }
    // end of spp case

    else  // more inelastic final state or single kaon production

    {  // cout<<"inel "<<W<<" "<<nu<<endl;
      e.weight = fromdis * 1e-38 * kin.jacobian;

      for (int i = 0; i < nof_particles; i++) {
        particle part;
        part.t = pythiaParticle->P[3][i] * GeV;
        part.x = pythiaParticle->P[0][i] * GeV;
        part.y = pythiaParticle->P[1][i] * GeV;
        part.z = pythiaParticle->P[2][i] * GeV;
        rotation(part, kin.q);

        part = part.boost(kin.hadron_speed);  // correct direction ???
        part = part.boost(kin.target.v());

        part.ks = pythiaParticle->K[0][i];
        part.pdg = pythiaParticle->K[1][i];
        part.orgin = pythiaParticle->K[2][i];

        e.temp.push_back(part);
        if (part.ks == 1)  // condition for a real particle in the final state
        {
          e.out.push_back(part);
        }
      }
    }
    // end of more inelastic part or single kaon production

    delete pythia71;
  }  // end of W>1210 &&  !fromdis==0

  // E above threshold
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
  for (int j = 0; j < e.out.size(); j++) cout << e.out[j] << "\n";
}