#ifndef _RES_XSEC_H
#define _RES_XSEC_H

#include "res_kinematics.h"

//! store xsec / spp parameters and xsec contributions from different channels
struct res_xsec {
  int j, k, l;       //!< indices for SPP table (see singlepion.cc)
  int final_charge;  //!< total electric charge of the pion-nucleon system

  double dis_pip, dis_pi0, dis_pim;        //!< contribution to RES from DIS
  double delta_pip, delta_pi0, delta_pim;  //!< contribution to RES from Delta

  double dis_spp, delta_spp;  //!< total contribution from RES and DIS

  double from_dis;  //!< the strength of DIS background

  res_xsec(res_kinematics &kin, const bool cc);
};

#endif