#include <fstream>
#include <random>

#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "params.h"
#include "sf/CSFOptions.h"
#include "sf/CSpectralFunc.h"
#include "sf/GConstants.h"
#include "sfevent.h"

using namespace NUWRO;
static inline double pow2(double x) { return x * x; }

// main function to generate event using spectral function
double sfevent(params &par, event &e, nucleus &t) {
  // references to initial particles (for convenience)
  particle &l0 = e.in[0];  // incoming neutrino
  particle &N0 = e.in[1];  // target nucleon

  // flags used to set up final state configuration
  const bool is_anti = l0.pdg < 0;                  // true for anti-neutrino
  const bool is_on_n = N0.pdg == pdg_neutron;       // true for target neutron
  const bool is_cc_possible = is_anti xor is_on_n;  // true for nu+n and nubar+p

  if (e.flag.cc and not is_cc_possible) return 0;  // CC is not possible

  particle l1;  // outgoing lepton
  particle N1;  // outgoing nucleon
  particle N2;  // nucleon spectator

  N2.r = N1.r = N0.r;  // final nucleons position = target nucleon position

  // outoing nucleon isospin
  // CC on proton  -> neutron (xor = 1)
  // CC on neutron -> proton  (xor = 0)
  // NC on proton  -> proton  (xor = 0)
  // NC on neutron -> neutron (xor = 1)
  is_on_n xor e.flag.cc ? N1.set_neutron() : N1.set_proton();

  // outgoing lepton pdg
  // NC = the same as incoming neutrino
  // CC nu = neutrino pdg - 1
  // CC nubar = neutrino pdg + 1
  l1.pdg = l0.pdg;
  if (e.flag.cc) l1.pdg += is_anti ? 1 : -1;

  // spectator isospin (assuming pn pairs for SRC)
  is_on_n ? N2.set_proton() : N2.set_neutron();

  const double m = mass(l1.pdg);  // outgoing lepton mass
  const double M = N1.mass();     // outgoing nucleon mass
  const double m2 = m * m;        // lepton mass squared
  const double M2 = M * M;        // nucleon mass squared

  l1.set_mass(m);  // set outgoing lepton mass

  CSFOptions options(par, e.flag.cc, !is_on_n, is_anti);  // SF configuration
  CSpectralFunc *sf = options.get_SF();                   // create spectral function

  // target nucleon momentum (p) and removal energy (E) generated according to probability distribution given by SF
  const double p = sf->MomDist()->generate();  // target nucleon momentum
  double E = get_E(sf, p);                     // removal energy

  // if the interaction occurs on neutron apply Coulomb correction to energy levels
  if (par.sf_coulomb and is_on_n) E += coulomb_correction_neutron(par.nucleus_p, par.nucleus_n);

  // set target nucleon momentum randomly from Fermi sphere
  N0.set_momentum(rand_dir() * p);
  // set opposite momentum for nucleon spectator
  N2.set_momentum(-N0.p());

  // s mandelstam-like
  vect s = l0 + N0;
  s.t = l0.E() + N0.mass() - E;
  const double s2 = s * s;

  if (s2 < pow2(M + m)) return 0;  // check if kinematics possible

  const vec v = s.v();  // the velocity of cms frame

  // kinematics in cms
  const double mom_cms = sqrt(0.25 * pow2(s2 + m2 - M2) / s2 - m2);
  const vec dir_cms = rand_dir();
  // set lepton and nucleon momenta (in cns)
  l1.set_momentum(mom_cms * dir_cms);
  N1.set_momentum(-l1.p());
  // boost to lab frame
  l1.boost(v);
  N1.boost(v);

  // check Pauli blocking
  if (par.pauli_blocking) {
    if (par.sf_pb == 0 and N1.momentum() < sf->get_pBlock())
      return 0;
    else if (par.sf_pb == 1 and N1.momentum() < t.localkf(N1))
      return 0;
    else if (par.sf_pb == 2 and frandom() < sf->MomDist()->Tot(N1.momentum()) / sf->MomDist()->Tot(0))
      return 0;
  }

  // four-momentum transfer
  vect q = N1 - N0;
  // sphere volume in cms
  const double vol = 4 * pi * mom_cms * mom_cms;
  // gradient for Dirac delta when integrating over k'
  const double graddelta = (l1.v() - N1.v()).length();
  // surface scaling when going from lab (elipsoide) to cms (sphere)
  const double surfscale = sqrt(1 - pow2(v * dir_cms)) / sqrt(1 - v * v);
  // cross section
  const double common = G * G / 8 / pi / pi * vol * (surfscale / graddelta) / (l1.E() * l0.E() * N0.E() * N1.E());
  const double val =
      e.flag.cc ? common * cos2thetac * options.evalLH(q * q, l0 * N0, l1 * N0, q * N0, l0 * q, l1 * q, l0 * l1)
                : common * options.evalLHnc(q * q, l0 * N0, l1 * N0, N0 * q, l0 * q, l1 * q, l0 * l1);

  double q0_shift = 0.0;  // energy transfer shift due to FSI and/or Coulomb correction

  // apply Couloumb corrections for charged leptons
  if (par.sf_coulomb and e.flag.cc) q0_shift += coulomb_correction(is_anti, par.nucleus_p, par.nucleus_n);

  if (par.sf_fsi and par.nucleus_p == 6 and par.nucleus_n == 6) {
    // apply FSI as described in: A. Ankowski et al, PRD91 (2015) 033005
    // express knock-out nucleon kinetic energy in terms of beam energy and scattering angle (eq. 7)
    const double Ek = e.in[0].E();
    const double x = 1 - l1.p().z / l1.momentum();
    const double Tk = Ek * Ek * x / (M + Ek * x);

    // energy transfer shift (as defined in eq. 3)
    if( Tk < 299.088 )
      q0_shift += potential_real(Tk);  // real part of optical potential
    // apply folding function smearing (eq. 2)
    if (frandom11() > sqrt(transparency(2 * M * Tk))) {
      // repeat until energy transfer > 0
      // loop stopped after 100 tries (although it should never happen)
      int n_tries = 0;
      while (n_tries++ < 100) {
        // calculate total shift
        const double shift = q0_shift + random_omega();
        if (l1.E() - shift < l0.E()) {
          // accept random omega
          q0_shift = shift;
          break;
        }
      }
    }
  }

  // modify lepton kinetic energy or xsec = 0 if not possible
  if (l1.Ek() > q0_shift)
    l1.set_energy(l1.E() - q0_shift);
  else
    return 0;

  // modify nucleon kinetic energy or xsec = 0 if not possible
  if (N0.mass() > E + q0_shift)
    N0.t = N0.mass() - E - q0_shift;
  else
    return 0;

  e.weight = val / cm2;

  // push final state particles
  // N0.t = N0.mass() - E;
  e.in[1] = N0;
  e.out.push_back(l1);
  e.out.push_back(N1);

  // add a spectator if on correlated pair
  if (par.sf_method == 1 and is_src(p, E, t.p, t.n, !is_on_n) and (l0.t - l1.t - N1.Ek() - N2.Ek()) > 14)
    e.out.push_back(N2);

  return val;
}

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

// get removal energy for given momentum p
double get_E(CSpectralFunc *sf, double p) {
  // TODO: do we need this loop? or generateE(p) should be modified?
  double E;

  do {
    E = sf->generateE(p);
  } while (!(E == E));

  return E;
}

// determine if scattering occured on correlated pair of nucleons
bool is_src(double p, double E, int Z, int N, bool is_on_p) {
  // oxygen
  if (Z == 8 and N == 8) {
    if (p < 85 and E > 63) return true;
    if (p > 85 and p < 320 and E > (73.4 - 0.167 * p)) return true;
    if (p > 320 and p < 390 and E > 19.1) return true;
    if (p > 395) return true;
  }
  // carbon
  // Benhar SF; basically for protons but taken the same for neutrons
  if (Z == 6 and N == 6) {
    if (p < 330 and E > (52.27 + 0.00428 * p - 0.0004618 * p * p)) return true;
    if (p > 330) return true;
  }
  // iron
  // Benhar SF; basically for protons but taken the same for neutrons
  if (Z == 26 and N == 30) {
    if (p < 335 and E > (60.41 + 0.004134 * p - 0.0004343 * p * p)) return true;
    if (p > 335) return true;
  }
  // argon
  if (Z == 18 and N == 22) {
    if (is_on_p) {  // protons
      if (p < 230 and E > (49.11 - 0.08305 * p + 0.0008781 * p * p + 1.045e-7 * p * p * p - 8.312e-9 * p * p * p * p))
        return true;

      if (p > 230 and p < 395 and E > (-52.97 + 0.8571 * p - 0.001696 * p * p)) return true;

      if (p > 395) return true;
    } else {  // neutrons
      if (p < 225 and E > (50.03 - 0.0806 * p + 0.0006774 * p * p + 1.717e-6 * p * p * p - 1.236e-8 * p * p * p * p))
        return true;

      if (p > 225 and p < 395 and E > (-17.23 + 0.6314 * p - 0.001373 * p * p)) return true;

      if (p > 395) return true;
    }
  }

  return false;
}

// return the value of transparency for given Q2 (used to determine FSI)
double transparency(double Q2) {
  // parametrization for Carbon: O. Benhar et al. Phys.Rev. D72 (2005) 053005

  Q2 /= 1000.0;  // please note units [GeV^2 *1000]

  // constant for Q2 > 1000
  if (Q2 > 1000) return 0.5792642140468227;

  // fit (polynomial)
  static const double coeff[] = {7.71692837e-01, -2.77751361e-04, 2.24980171e-06, -1.11358859e-08,
                                 1.98862243e-11, -1.50900788e-14, 4.17699547e-18};

  double T = coeff[0];
  double x = Q2;

  for (int i = 1; i < 7; i++) {
    T += coeff[i] * x;
    x *= Q2;
  }

  return T;
}

// real part of the potential which modifies energy transfer
double potential_real(double Tk) {
  // parametrization for Carbon: A. Ankowski et al, PRD91 (2015) 033005
  // fit (polynomial)
  static const double coeff[] = {-3.76929648e+01, 4.35269313e-01, -2.59678634e-03, 9.55434214e-06,
                                 -2.15373898e-08, 3.07501687e-11, -2.83810998e-14, 1.69043802e-17,
                                 -6.27515290e-21, 1.32038136e-24, -1.20270294e-28};

  double V = coeff[0];
  double x = Tk;

  for (int i = 1; i < 11; i++) {
    V += coeff[i] * x;
    x *= Tk;
  }

  return V;
}

// gaussian fit to energy transfer shift
double random_omega() {
  // for |q| = 1 GeV: O. Benhar PRC 87 (2013) 024606
  // fit to gauss
  static std::default_random_engine generator;
  static std::normal_distribution<double> distribution(5.43264624e-04, 8.88774322e+01);
  return distribution(generator);
}

// Coulomb correction for outgoing charged lepton
double coulomb_correction(bool is_anti, int p, int n) {
  double shift = 0.0;  // correction in MeV

  if (p == 6 and n == 6) shift = 3.5;  // carbon

  return is_anti ? -shift : shift;
}

// Coulomb correction to the neutron energy levels
double coulomb_correction_neutron(int p, int n) {
  switch (1000 * p + n) {
    case 6006:
      return 2.8;  // carbon
    default:
      return 0;
  }
}
