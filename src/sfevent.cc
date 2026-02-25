#include <fstream>
#include <vector>
#include <iostream>
#include <random>
#include <array>
#include <map>
#include <tuple>
#include <limits>
#include <algorithm>

#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "params.h"

#include "sf/CSFOptions.h"
#include "sf/CSpectralFunc.h"
#include "sf/CSpectralFunctions.h"
#include "sf/GConstants.h"

#include "sfevent.h"

static inline double pow2(double x) { return x * x; }
typedef std::tuple<TargetNucleus, bool, bool> SFKey;

// Generate kinematics and calculate cross section using SF
double sfevent(params &par, event &e, nucleus &t)
{

  static std::map<SFKey, CSpectralFunctions*> sfCache;

    particle &l0 = e.in[0];  // incoming lepton
    particle &N0 = e.in[1];  // target nucleon

  // Get target nucleus composition
  const int Z = par.nucleus_p;
  const int N = par.nucleus_n;

  // Interaction flags
  const bool isAnti = l0.pdg < 0;              // true for anti-neutrino
  const bool onNeutron = N0.pdg == pdg_neutron;// true for target neutron
  const bool ccAllowed = isAnti ^ onNeutron;   // true for nu+n and nubar+p
  const bool isElectron = (l0.pdg == pdg_e);   // true for electron interactions
  const bool isNC = e.flag.nc;                 // true for neutral current intercations
  const bool doPB = par.pauli_blocking == 1;   // true for enabled pauli blocking
  const bool doRecoil = par.sf_recoil == 1;    // true for nuclear recoil
  const bool doCoulomb = par.sf_Coulomb == 1;  // true for coulomb distortion
  const bool doFSI = par.FSI_on == 1;          // true for FSI
  const bool doSRC = par.sf_src == 1;          // true for SRC

  int method = par.sf_method;
  bool isCorrelated = false;
  double pot( 0.0 ), pot_cascade( 0.0 ), prefactor( 1.0 ), xs( 1.0 );
  double factor = isAnti ? -1 : 1 ;

  if ((e.flag.cc && !ccAllowed) || method == 0) return 0; // If charged current (CC) interaction isn't possible, exit

    particle l1;  // outgoing lepton
    particle N1;  // outgoing struck nucleon
    particle N2;  // spectator second nucleon

  N2.r = N1.r = N0.r; // final nucleons position = target nucleon position

  // outoing nucleon isospin
  // CC on proton  -> neutron (xor = 1)
  // CC on neutron -> proton  (xor = 0)
  // NC on proton  -> proton  (xor = 0)
  // NC on neutron -> neutron (xor = 1)
  onNeutron ^ e.flag.cc ? N1.set_neutron() : N1.set_proton();

  // outgoing lepton pdg
  // NC = the same as incoming neutrino or electron
  // CC nu = neutrino pdg - 1
  // CC nubar = neutrino pdg + 1
  l1.pdg = l0.pdg;
  if (e.flag.cc) l1.pdg += isAnti ? 1 : -1;

  // spectator isospin (assuming only pn pairs for NN SRC)
  onNeutron ? N2.set_proton() : N2.set_neutron();

  const double mi = mass( l0.pdg );       // incoming lepton mass
  const double m = mass( l1.pdg );        // outgoing lepton mass
  const double M = N1.mass();             // outgoing nucleon mass // check?
  const double mi2 = mi * mi;             // incoming lepton mass squared
  const double m2 = m * m;                // outgoing lepton mass squared
  const double MSq = M * M;               // outgoing nucleon mass squared
  l1.set_mass( m );                       // set outgoing lepton mass

   // Set the target nucleus
   TargetNucleus target = ( Z == carbonZ && N == carbonN ) ? C12_newSF :
                          ( Z == oxygenZ && N == oxygenN ) ? O16_SF :
                          ( Z == argonZ  && N == argonN )  ? Ar40_SF :
                          ( Z == ironZ   && N == ironN )   ? Fe56_SF :
                          TargetNucleus::Unsupported;
   if (target == TargetNucleus::Unsupported) return 0;

   const bool isAr40 = (target == Ar40_SF);

    // Create the spectral function object for sf_method = 1
    SFKey key = std::make_tuple(target, doSRC, doFSI);
    CSpectralFunctions *sf; // CSpectralFunctions instance

    auto it = sfCache.find( key );
    if (it != sfCache.end()) sf = it->second;
    else
    {
          sf = new CSpectralFunctions(target, doSRC, doFSI);
          sfCache[ key ] = sf;
    }

  // Apply Coulomb distortion effects to the incoming charged lepton
  double averCE = ( doCoulomb ) ? factor * sf->get_CoulombAvEnergy() : 0;
  double maxCE = ( doCoulomb ) ? factor * sf->get_CoulombMaxEnergy() : 0;

  particle l0eff = l0;

  if (isElectron) l0eff.set_energy( l0.E() + averCE );

  CSFOptions SF( par, e.flag.cc, !onNeutron, isAnti ); // Create cross-section object

  // Draw the missing momentum and removal energy according to probability distribution given by SF based on MF or correlation
  double p(0.0), E(0.0);

  // new implementation of SF and nuclear corrections by RWIK DHARMAPAL BANERJEE, 2025
  if (method == 1)
  {
    if (doSRC && isAr40)
    {
        double corrFraction = isAnti ? sf->get_corrProtonFraction() : sf->get_corrNeutronFraction();
        if (frandom11() < corrFraction)
        {
            // Sample from correlated spectral function
            isCorrelated = true;
            if (onNeutron)
            {
                p = sf->generateNeutronMomentum_corr();
                E = sf->generateNeutronRemovalEnergy_corr(p);
            }
            else
            {
                p = sf->generateProtonMomentum_corr();
                E = sf->generateProtonRemovalEnergy_corr(p);
            }
        }
        else
        {
            // Sample from mean-field spectral function
            isCorrelated = false;
            if (onNeutron)
            {
                p = sf->generateNeutronMomentum_MF();
                E = sf->generateNeutronRemovalEnergy_MF(p);
            }
            else
            {
                p = sf->generateProtonMomentum_MF();
                E = sf->generateProtonRemovalEnergy_MF(p);
            }
        }
     }
     else
     {
        // Use total spectral function
        p = onNeutron ? sf->generateNeutronMomentum() : sf->generateProtonMomentum();
        E = onNeutron ? sf->generateNeutronRemovalEnergy(p) : sf->generateProtonRemovalEnergy(p);
     }
   }
   else if (method == 2 || method == 3) // older SF implementations
   {
    CSpectralFunc *sf_ = SF.get_SF();
    p = sf_->MomDist()->generate();
    E = sf_->generateE(p);
    if (onNeutron) E += coulomb_correction_neutron(Z, N);
  }

  N0.set_momentum( rand_dir() * p ); // set target nucleon momentum randomly from Fermi sphere

  if (isAr40 && method == 1)  // Generate the momentum of the spectator nucleon for argon
  {
   const vec spectatorNucleonMomentum = onNeutron
        ? sf->generateAdditionalProtonsMomentum( p, E, N0.p().dir() )
        : sf->generateAdditionalNeutronsMomentum( p, E, N0.p().dir() );

    N2.set_momentum( spectatorNucleonMomentum );
  }
  else
    N2.set_momentum( -N0.p() ); // Back-to-back approximation for other targets

  // Apply Nuclear recoil
  const double targNucleusMass = sf->get_targetMass();          // target nucleus mass
  const double residualNucleusMass =  targNucleusMass - M + E ; // Residual (& excited) nucleus mass after one nucleon knock-out
  const double residualNucleusEnergy = std::sqrt( pow2( residualNucleusMass ) + pow2( p ) ) ; // residual nucleus energy

  vect s = l0eff + N0; // s Mandelstam-like
  s.t = doRecoil
        ? l0eff.E() + targNucleusMass - residualNucleusEnergy
        : l0eff.E() + M - E;         // adjust the initial energy
  const double s2 = s * s;
  if (s2 < pow2( M + m )) return 0;  // Unphysical invariant mass

  const vec v = s.v();               // the velocity of cms frame
  const double mom_cms = sqrt( 0.25 * pow2( s2 + m2 - MSq ) / s2 - m2 );
  const vec dir_cms = rand_dir();

  l1.set_momentum( mom_cms * dir_cms );  // set lepton momenta (in cms)
  N1.set_momentum( -l1.p() );            // set nucleon momenta (in cms)

  // Boost back to lab frame
  l1.boost( v );
  N1.boost( v );

  // Calculate the focusing factor
  const double E_l0 = l0.E();
  const double kSq = pow2( l0.momentum() );
  const double kSq_eff_max = pow2( E_l0 + maxCE ) - mi2;
  const double focusingFactorIncomingSq = ( isElectron ) ? ( kSq_eff_max / kSq ) : 1;
  
  const double E_l1 = ( isAnti && l1.Ek() + maxCE < 0) ? l1.E() - maxCE : l1.E();
  const double kPrimeSq_eff_max = pow2( E_l1 + maxCE ) - m2;
  const double kPrimeSq_eff = pow2( E_l1 + averCE ) - m2;
  const double focusingFactorOutgoingSq = kPrimeSq_eff_max / kPrimeSq_eff;
  
  if ((kSq <= 0) || (kPrimeSq_eff <= 0)) return 0;

  const double focusingFactorSq = ( focusingFactorIncomingSq * focusingFactorOutgoingSq );
  
   // FSI
   bool isTransparent = true;

   if (doFSI)
   {
       const double Tk = tPPrime_approx( e.in[0].E(), l1.p().z / l1.momentum(), isElectron, isNC, m2, M ) ;
       const double rOP = sf->eval_realOP( Tk );                      // Real part of the optical potential
       pot = ( N1.pdg == pdg_neutron ) ? rOP - factor * averCE: rOP ; // correct for proton and neutron
       int tidx = par.sf_transparency_table_idx - 1;  // 1→0, 2→1

      isTransparent = ( frandom11() <= sf->eval_sqrtOfTransparency(Tk, par.sf_transparency_scale, tidx) );
      if (!isTransparent)
      {
          const double shift = pot + random_omega();
          if( l1.E() - l0eff.E() > shift ) return 0;
          else pot = shift;
       }
   }
   
   double U = std::min(0.0, pot);
   
  // Pauli blocking   
  if (doPB && par.sf_pb != 0 &&
     ((par.sf_pb == 1 && sf->eval_particleSF(N1.pdg, N1.momentum()) == 0.0) ||
      (par.sf_pb == 2 && N1.momentum() < t.localkf(N1))))
  return 0;

  vect  q = N1 - N0; // four-momentum transfer (hadronic side)
  vect qq = l0 - l1; // four-momentum transfer (leptonic side)

  const double sphereArea = 4 * Pi * mom_cms * mom_cms;                         
  const double graddelta = ( l1.v() - N1.v() ).length();                        // Gradient for Dirac delta when integrating over k'
  const double surfscale = sqrt( 1 - pow2( v * dir_cms ) ) / sqrt( 1 - v * v ); // Surface scaling when going from lab (elipsoide) to cms (sphere)
  const double jacobian =  sphereArea * ( surfscale / graddelta );
  const double energyDenominator = ( l1.E() * l0eff.E() * N0.E() * N1.E() );
  const double phasespaceJac = jacobian / energyDenominator;

  // Evaluate cross-section
  if (isElectron)
  {
    prefactor = focusingFactorSq / pow2( reciprocalAlpha ) / ( qq*qq ) / ( qq*qq ) * phasespaceJac;
    xs = prefactor * SF.evalLHel( q*q, l0eff*N0, l1*N0, q*N0, l0eff*q, l1*q, l0eff*l1 );
  }
  else
  {
    prefactor =  focusingFactorSq * pow2( GF ) / 8 / Pi2 * phasespaceJac;
    xs = e.flag.cc ? prefactor * cos2ThetaC * SF.evalLH( q*q, l0*N0, l1*N0, q*N0, l0*q, l1*q, l0*l1 )
                   : prefactor * SF.evalLHnc( q*q, l0*N0, l1*N0, N0*q, l0*q, l1*q, l0*l1 );
  }

   // Apply nuclear effects
   double en_thr = isAnti ? pot : averCE + pot;
   if (l1.Ek() > en_thr) l1.set_energy( l1.E() - averCE - pot );
   else return 0;

   if(!onNeutron) {
     if( N1.Ek() > -averCE ) N1.set_energy( N1.E() + averCE );
     else return 0;
   }

   //fixing a bug in initial nucleon energy
   if (N0.mass() > E) N0.t = N0.mass() - E;
   else return 0;

  // Cut on the excitation energy
  static const double minimalEX( 0.0 * MeV );
  const double coefASq( targNucleusMass * minimalEX + 0.5 * ( pow2(minimalEX) - m2 ) + pow2(l0eff.E()) );
  const double coefB ( l0eff.E() + targNucleusMass );
  const double coefC ( l0eff.E() * ( l1.p().z / l1.momentum() ) );
  const double radicand( pow2( coefASq - coefB * l0eff.E() ) - ( pow2(coefB) - pow2(coefC) ) * m2 );
  const double root( radicand < 0.0 ? 0.0 : sqrt( radicand ) );
  const double minimalOmegaExact( ( coefASq * coefB - pow2(coefC) * l0eff.E() - coefC * root ) / ( pow2(coefB) - pow2(coefC) ) );
  if (l0eff.E() - l1.E() < minimalOmegaExact) return 0;
  
    // Approximated SRC (except argon SF)
    if ((method == 1 && !isAr40) || (method > 1))
    {
      if (doSRC && is_src(p, E, Z, N, !onNeutron) && (l0eff.t - l1.t - N1.Ek() - N2.Ek()) > 14)
        isCorrelated = true;
      else
        isCorrelated = false;
    }

  e.weight = xs / cm2;                   // Add cross-section as the weight
  
  e.in[1] = N0;                          // Update target nucleon state
  e.out.push_back(l1);                   // Add outgoing lepton to event.
  e.out.push_back(N1);                   // Struck nucleon added to out
  if (isCorrelated && N2.momentum() > 1e-6) {e.out.push_back(N2);} // Handle spectator nucleon

  e.averageCE = averCE;
  e.flag.isTransparent = isTransparent;
  e.optical_potential  = U; 
  e.flag.isCorrelated  = isCorrelated;  

  // Apply acceptance cut for electron scattering
  if (isElectron)
  {
    double cosTheta = l1.p().z / l1.momentum();
    if (cosTheta < (par.el_costh_lab - par.el_costh_del) || cosTheta > (par.el_costh_lab + par.el_costh_del))
    {
      e.weight = 0;
      return 0;
     }
    else
    {
      e.weight /= 2 * par.el_costh_del;
      xs /= 2 * par.el_costh_del;
     }
   }

  return xs;
}

// Approximate kinetic energy
double tPPrime_approx( const double eK, const double cosOfScattAngle, const bool EM, const bool NC, const double m_leptonSq, const double m_nucleon ) {
    const double x = 1.0 - cosOfScattAngle;
    if (EM || NC) return pow2(eK) * x / ( m_nucleon + eK * x );
    const double coefASq ( pow2(eK) - 0.5 * m_leptonSq );
    const double coefB ( eK + m_nucleon );
    const double coefC ( eK * cosOfScattAngle );
    const double radicand = pow2(coefASq - coefB * eK) - (pow2(coefB) - pow2(coefC)) * m_leptonSq;
    const double root( radicand < 0.0 ? 0.0 : sqrt( radicand ) );
    const double tPPrime_approx( ( coefASq * coefB - pow2(coefC) * eK - coefC * root ) / ( pow2(coefB) - pow2(coefC) ) );

    return tPPrime_approx;
}

bool has_sf(nucleus &t, int method)
{
  const int key = 1000 * t.Z() + t.N();
  switch (key)
  {
    case CARBON:
      return method == 2;
    case OXYGEN:
      return method == 3;
    case CALCIUM:
      return method == 3;
    case ARGON:
      return method == 2 || method ==3;
    default:
      return false;
  }
}

// Approximated SRC
bool is_src(double p, double E, int Z, int N, bool is_on_p)
{
  // carbon
  // Benhar SF; basically for protons but taken the same for neutrons
  if (Z == 6 && N == 6)
  {
    if (p < 330 && E > (52.27 + 0.00428 * p - 0.0004618 * p * p)) return true;
    if (p > 330) return true;
  }
  // oxygen
  // Benhar SF; basically for protons but taken the same for neutrons  
  if (Z == 8 && N == 8)
  {
    if (p < 85 && E > 63) return true;
    if (p > 85 && p < 320 && E > (73.4 - 0.167 * p)) return true;
    if (p > 320 && p < 390 && E > 19.1) return true;
    if (p > 395) return true;
  }
  // argon
  if (Z == 18 && N == 22)
  {
  if (is_on_p) {  // protons
    if (p < 230 && E > (49.11 - 0.08305 * p + 0.0008781 * p * p + 1.045e-7 * p * p * p - 8.312e-9 * p * p * p * p))
      return true;

    if (p > 230 && p < 395 && E > (-52.97 + 0.8571 * p - 0.001696 * p * p)) return true;

    if (p > 395) return true;
  }
  else
  {  // neutrons
    if (p < 225 && E > (50.03 - 0.0806 * p + 0.0006774 * p * p + 1.717e-6 * p * p * p - 1.236e-8 * p * p * p * p))
      return true;

    if (p > 225 && p < 395 && E > (-17.23 + 0.6314 * p - 0.001373 * p * p)) return true;

    if (p > 395) return true;
  }
  }  
  // iron
  // Benhar SF; basically for protons but taken the same for neutrons
  if (Z == 26 && N == 30)
  {
    if (p < 335 && E > (60.41 + 0.004134 * p - 0.0004343 * p * p)) return true;
    if (p > 335) return true;
  }

  return false;
}

// Gaussian fit to folding function
double random_omega()
{
  static std::default_random_engine generator;
  static std::normal_distribution<double> distribution( -7.00637e-15, 88.3146 );
  return distribution( generator );
}

// Coulomb correction to the neutron energy levels
double coulomb_correction_neutron(int p, int n)
{
  const int key = 1000 * p + n;
  switch (key)
  {
    case CARBON:
      return carbon11Mass - boron11Mass + nMass - pMass - eMass;  // carbon
    case OXYGEN:
      return oxygen15Mass - nitrogen15Mass + nMass - pMass - eMass; // oxygen
    case ARGON:
      return (argon39Mass + nMass - argon40Mass) - (scandium47Mass + pMass + eMass - titanium48Mass); // argon
    case CALCIUM:
      return (argon39Mass + nMass - argon40Mass) - (scandium47Mass + pMass + eMass - titanium48Mass); // calcium, taken same as argon
    default:
      return 0;
  }
}
