#include <fstream>

#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "params.h"
#include "sf/CSFOptions.h"
#include "sf/CSpectralFunc.h"
#include "sf/GConstants.h"

static double mevtofm = 8e6 / Pi2 / 4;

static inline double pow2(double x) { return x * x; }

// method=1 - grid (from Benhar)
// method=2 - sum of gaussians
bool has_sf(nucleus &t, int method) {
  switch (1000 * t.Z() + t.N()) {
    case 6006:
      return method == 1;
    case 8008:
      return method == 1 || method == 2;
    case 20020:
      return method == 2;
    case 18022:
      return method == 1 || method == 2;
    case 26030:
      return method == 1;
    default:
      return false;
  }
}

double sfevent2cc(params &par, event &e, nucleus &t) {
  particle l0 = e.in[0];  // neutrino
  particle N0 = e.in[1];  // initial nucleon
  particle l1;
  particle N1;
  particle N2;

  N1.r = N0.r;
  N2.r = N0.r;
  if ((l0.pdg > 0 and N0.pdg == pdg_proton) or
      (l0.pdg < 0 and N0.pdg == pdg_neutron))
    return 0;  // no CC interaction possible on this nucleon

  if (l0.pdg < 0) {
    N1.set_neutron();
    N2.set_neutron();
    l1.pdg = l0.pdg + 1;
  } else {
    N1.set_proton();
    N2.set_proton();
    l1.pdg = l0.pdg - 1;
  }

  double m = 0;
  switch (abs(l1.pdg)) {
    case pdg_e:
      m = mass_e;
      break;
    case pdg_mu:
      m = mass_mu;
      break;
    case pdg_tau:
      m = mass_tau;
      break;
  }

  l1.set_mass(m);

  double mm = m * m;
  double M = N1.m();
  double MM = M * M;

  e.in[1] = N0;

  //	cout<<"Options to create"<<endl;
  CSFOptions options(par, 1, N0.pdg == pdg_proton,
                     l0.pdg < 0);  // czy proton // czy antyneutrino
                                   //	cout<<"sf about to create"<<endl;
  CSpectralFunc *sf = options.get_SF();
  //	cout<<"sf created"<<endl;
  const double pBlock = sf->get_pBlock();
  //	cout<<"sf used"<<endl;

  const double p = sf->MomDist()->generate();

  N0.set_momentum(rand_dir() * p);
  //    cout<<"ola"<<endl;
  double EBEN;
  do {
    EBEN = sf->generateE(p);
  } while (!(EBEN == EBEN));

  vect s = l0 + N0;

  s.t = l0.E() + N0.mass() - EBEN;

  // below approximate seperation of the correlated part in order to introduce a
  // correlated (spectator) nucleon
  // MF and corr momentum distributions are compared and a fraction of corr is
  // evaluated
  // corr contribution is assumed to correspond to larger values of the E
  // argument

  bool corr;
  corr = false;

  if (t.p == 8 && t.n == 8 && par.sf_method == 1)  // Benhar SF; basically for
                                                   // protons but taken the same
                                                   // for neutrons
  {
    if (p < 85 && EBEN > 63) corr = true;

    if (p > 85 && p < 320 && EBEN > (73.4 - 0.167 * p)) corr = true;

    if (p > 320 && p < 390 && EBEN > 19.1) corr = true;

    if (p > 395) corr = true;
  }

  if (t.p == 6 && t.n == 6 && par.sf_method == 1)  // Benhar SF; basically for
                                                   // protons but taken the same
                                                   // for neutrons
  {
    if (p < 330 && EBEN > (52.27 + 0.00428 * p - 0.0004618 * p * p))
      corr = true;

    if (p > 330) corr = true;
  }

  if (t.p == 26 && t.n == 30 && par.sf_method == 1)  // Benhar SF; basically for
                                                     // protons but taken the
                                                     // same for neutrons
  {
    if (p < 335 && EBEN > (60.41 + 0.004134 * p - 0.0004343 * p * p))
      corr = true;

    if (p > 335) corr = true;
  }

  if (t.p == 18 && t.n == 22 && par.sf_method == 1 &&
      N0.pdg == pdg_proton)  // Argon protons
  {
    if (p < 230 &&
        EBEN > (49.11 - 0.08305 * p + 0.0008781 * p * p + 1.045e-7 * p * p * p -
                8.312e-9 * p * p * p * p))
      corr = true;

    if (p > 230 && p < 395 && EBEN > (-52.97 + 0.8571 * p - 0.001696 * p * p))
      corr = true;

    if (p > 395) corr = true;
  }

  if (t.p == 18 && t.n == 22 && par.sf_method == 1 &&
      N0.pdg == pdg_neutron)  // Argon neutrons
  {
    if (p < 225 &&
        EBEN > (50.03 - 0.0806 * p + 0.0006774 * p * p + 1.717e-6 * p * p * p -
                1.236e-8 * p * p * p * p))
      corr = true;

    if (p > 225 && p < 395 && EBEN > (-17.23 + 0.6314 * p - 0.001373 * p * p))
      corr = true;

    if (p > 395) corr = true;
  }

  vec mom(N0.x, N0.y, N0.z);
  N2.set_momentum(-mom);

  double ss = s * s;

  if (ss < pow2(M + m)) return 0;

  vec v = s.v();

  /// here we do ::decay(s,N1,l1) by hand

  double pcms = sqrt(0.25 * pow2(ss + mm - MM) / ss - mm);
  vec dircms = rand_dir();
  vec p1 = pcms * dircms;
  l1.set_momentum(p1);
  N1.set_momentum(-p1);
  l1.boost(v);
  N1.boost(v);

  const double omega = l0.E() - l1.E();

  const double omegaTil = N1.E() - N0.E();

  if (false)
    if (omegaTil < 0) return 0;
  if (false)
    if (omega < 0) return 0;

  if (par.pauli_blocking) {
    if (par.sf_pb == 0 and N1.momentum() < pBlock)
      return 0;
    else if (par.sf_pb == 1 and N1.momentum() < t.localkf(N1))
      return 0;
    else if (par.sf_pb == 2 and
             frandom() <
                 sf->MomDist()->Tot(N1.momentum()) / sf->MomDist()->Tot(0))
      return 0;
  }

  vect q4til = N1 - N0;
  const double q4til2 = q4til * q4til;
  const double p4k4 = l0 * N0;
  if (false)
    if (p4k4 < 0 or q4til2 > 0) return 0;

  double pp = pcms * pcms;

  double vol = 4 * pi * pp;

  double gamma = 1 / sqrt(1 - v * v);
  double graddelta = (l1.v() - N1.v()).length();
  double surfscale =
      sqrt(1 + (1 - pow2(v * dircms) / (v * v)) * (gamma * gamma - 1));
  double val = G * G * cos2thetac / 8 / pi / pi * vol *
               (surfscale / graddelta) / (l1.E() * l0.E() * N0.E() * N1.E()) *
               options.evalLH(q4til * q4til, l0 * N0, l1 * N0, q4til * N0,
                              l0 * q4til, l1 * q4til, l0 * l1);

  e.weight = val / cm2;
  N0.t = N0.mass() - EBEN;
  e.in[1] = N0;
  e.out.push_back(l1);
  e.out.push_back(N1);

  if (corr && (l0.t - l1.t - N1.Ek() - N2.Ek()) > 14) e.out.push_back(N2);

  // cout<<val/cm2<<endl<<endl;
  return val;
}

double sfevent2nc(params &par, event &e, nucleus &t) {
  //	cout<<"sf created"<<endl;
  //  cout<<"sf used"<<endl;

  particle l0 = e.in[0];  // neutrino
  particle N0 = e.in[1];  // initial nucleon
  particle l1 = l0;
  particle N1 = N0;
  particle N2;

  if (N0.pdg == pdg_proton)
    N2.set_neutron();
  else
    N2.set_proton();

  double m = 0;
  double mm = m * m;
  double M = N1.m();
  double MM = M * M;

  e.in[1] = N0;

  CSFOptions options(par, 0, N0.pdg == pdg_proton,
                     l0.pdg < 0);  // czy proton // czy antyneutrino
  CSpectralFunc *sf = options.get_SF();
  double pBlock = sf->get_pBlock();

  const double p = sf->MomDist()->generate();

  N0.set_momentum(rand_dir() * p);

  double EBEN;
  do {
    EBEN = sf->generateE(p);
  } while (!(EBEN == EBEN));

  vect s = l0 + N0;

  s.t = l0.E() + N0.mass() - EBEN;

  bool corr;
  corr = false;

  if (t.p == 8 && t.n == 8 && par.sf_method == 1) {
    if (p < 85 && EBEN > 63) corr = true;

    if (p > 85 && p < 320 && EBEN > (73.4 - 0.167 * p)) corr = true;

    if (p > 320 && p < 390 && EBEN > 19.1) corr = true;

    if (p > 395) corr = true;
  }

  if (t.p == 6 && t.n == 6 && par.sf_method == 1)  // Benhar SF; basically for
                                                   // protons but taken the same
                                                   // for neutrons
  {
    if (p < 330 && EBEN > (52.27 + 0.00428 * p - 0.0004618 * p * p))
      corr = true;

    if (p > 330) corr = true;
  }

  if (t.p == 26 && t.n == 30 && par.sf_method == 1)  // Benhar SF; basically for
                                                     // protons but taken the
                                                     // same for neutrons
  {
    if (p < 335 && EBEN > (60.41 + 0.004134 * p - 0.0004343 * p * p))
      corr = true;

    if (p > 335) corr = true;
  }

  if (t.p == 18 && t.n == 22 && par.sf_method == 1 &&
      N0.pdg == pdg_proton)  // Argon protons
  {
    if (p < 230 &&
        EBEN > (49.11 - 0.08305 * p + 0.0008781 * p * p + 1.045e-7 * p * p * p -
                8.312e-9 * p * p * p * p))
      corr = true;

    if (p > 230 && p < 395 && EBEN > (-52.97 + 0.8571 * p - 0.001696 * p * p))
      corr = true;

    if (p > 395) corr = true;
  }

  if (t.p == 18 && t.n == 22 && par.sf_method == 1 &&
      N0.pdg == pdg_neutron)  // Argon neutrons
  {
    if (p < 225 &&
        EBEN > (50.03 - 0.0806 * p + 0.0006774 * p * p + 1.717e-6 * p * p * p -
                1.236e-8 * p * p * p * p))
      corr = true;

    if (p > 225 && p < 395 && EBEN > (-17.23 + 0.6314 * p - 0.001373 * p * p))
      corr = true;

    if (p > 395) corr = true;
  }

  vec mom(N0.x, N0.y, N0.z);
  N2.set_momentum(-mom);

  double ss = s * s;

  if (ss < pow2(M + m)) return 0;

  vec v = s.v();

  /// here we do ::decay(s,N1,l1) by hand

  double pcms = sqrt(0.25 * pow2(ss + mm - MM) / ss - mm);
  vec dircms = rand_dir();
  vec p1 = pcms * dircms;
  l1.set_momentum(p1);
  N1.set_momentum(-p1);
  l1.boost(v);
  N1.boost(v);

  const double omega = l0.E() - l1.E();

  const double omegaTil = N1.E() - N0.E();

  if (false)
    if (omegaTil < 0) return 0;
  if (false)
    if (omega < 0) return 0;

  if (par.pauli_blocking) {
    if (par.sf_pb == 0 and N1.momentum() < pBlock)
      return 0;
    else if (par.sf_pb == 1 and N1.momentum() < t.localkf(N1))
      return 0;
    else if (par.sf_pb == 2 and
             frandom() < sf->MomDist()->Tot(N1.momentum() * mevtofm))
      return 0;
  }

  vect q4til = N1 - N0;
  const double q4til2 = q4til * q4til;
  const double p4k4 = l0 * N0;
  if (false)
    if (p4k4 < 0 or q4til2 > 0) return 0;

  double pp = pcms * pcms;

  double vol = 4 * pi * pp;

  double gamma = 1 / sqrt(1 - v * v);
  double graddelta = (l1.v() - N1.v()).length();
  double surfscale =
      sqrt(1 + (1 - pow2(v * dircms) / (v * v)) * (gamma * gamma - 1));
  double val = G * G / 8 / pi / pi * vol * (surfscale / graddelta) /
               (l1.E() * l0.E() * N0.E() * N1.E()) *
               options.evalLHnc(q4til * q4til, l0 * N0, l1 * N0, N0 * q4til,
                                l0 * q4til, l1 * q4til, l0 * l1);

  e.weight = val / cm2;
  N0.t = N0.mass() - EBEN;
  e.in[1] = N0;
  e.out.push_back(l1);
  e.out.push_back(N1);

  if (corr && (l0.t - l1.t - N1.Ek() - N2.Ek()) > 14) e.out.push_back(N2);
  return val;
}
