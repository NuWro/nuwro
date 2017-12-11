#include "../event1.h"
#include "../ff.h"
#include "../nucleus.h"
#include "../params.h"
#include "../qel_sigma.h"
#include "../rpa_2013.h"
#include "rewparams.h"

static double E_bind(params &p, nucleus &t, particle &N0) {
  /*
      Binding energy of nucleon N0 inside target nucleus t
  */
  switch (p.nucleus_target) {
  case 0:
    return 0;
  case 1:
    return p.nucleus_E_b;
  case 2:
    return t.Ef(N0) + p.kaskada_w;
  case 3:
    return 0;
  case 4:
    return binen(N0.p(), p.nucleus_p, p.nucleus_n);
  case 5:
    return deuter_binen(N0.p()); // deuterium
  case 6:
    return p.nucleus_E_b; // deuterium like Fermi gas
  default:
    return 0;
  }
}

double calcQEL(event &e, params &p, nucleus &t) {
  /*
      Returns number proportional to quasi elastic cross section.

      Depends on global rew parameters mainly through form factor
      functions f12() and fap() used in qel_sigma() and rpa_ratio()
  */

  if (!e.flag.qel)
    return 1;

  double weight = 1;

  particle &nu = e.in[0];
  particle &N0 = e.in[1];

  particle &lepton = e.out[0];
  particle &N1 = e.out[1];

  int kind =
      e.flag.nc *
      (N0.proton() ? 1
                   : 2); /// process type: 0 - cc, 1 - nc proton, 2 - nc neutron

  double _E_bind = E_bind(e.par, t, N0);

  weight *= qel_sigma(nu.t, e.q2(), kind, nu.pdg < 0, lepton.mass(), N0.mass());

  bool new_ver = true;

  switch (p.qel_rpa) { //     qv   ,  q0            ,  E         , nu_pdg,
                       //     lepton_mass  ,  Meff     , kF   , version
  case 1:
    weight *= ratio_rpa(e.qv(), e.q0() - _E_bind, nu.t - _E_bind, nu.pdg,
                        lepton.mass(), N1.mass(), t.kF(), new_ver);
    break;
  case 3:
    weight *= ratio_rpa(e.qv(), e.q0() - _E_bind, nu.t - _E_bind, nu.pdg,
                        lepton.mass(), t.Mf(), t.kF(), new_ver);
    break;
  }

  return weight;
}
