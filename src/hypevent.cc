#include "jednostki.h"
#include <cassert>          
#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "event1.h"
#include "hypevent.h"
#include "kinematics.h"
#include "pdg.h"
#include "nucleus.h"
#include <cstdlib>
#include "hiperon_sigma.h"
#include "hyperon_interaction.h"
#include "scatter.h"
#include "dis/LeptonMass.h"
#define LOCALKF localkf_O


double hypevent(params&p, event &e, nucleus &t)
{
  /*
   * 
   * This is the hyperon production code. Only for antineutrino scattering.
   * 
   * p - global nuwro parameters 
   * e - event object to modify 
   * t - nucleus on which scattering occurs
   * returns the scattering cross section, produces the kinematics
   * and stores the outgoing particles in the e.out vector
   *
   * Based on the code by Chris Thorpe received on 18.06.19,
   * the previous code has also been exploited
   * 
   * by KN, Aug 2019
   *
   */

  e.flag.hyp=true;
  e.flag.cc=true;
  e.weight=0;

  if(e.in[0].pdg>0) 
    return 0; // this dynamics is for antineutrinos only

  particle nu = e.in[0];  // initial antineutrino
  particle N0 = e.in[1];  // initial nucleon
  particle lepton;        // outgoing antilepton
  particle hyperon;       // created hyperon
  // same coordinates for initial/final particles
  lepton.r=N0.r;
  hyperon.r=N0.r;

  // adjust the lepton flavour
  lepton.pdg=nu.pdg-(nu.pdg>0 ? 1 :-1);

  int h=0; // hyperon produced: 1,2,3
  switch(N0.pdg)
  {
    case PDG::pdg_proton: // Lambda or Sigma_0 (cross sections should add up)
      if(p.hyp_lambda)
      {
        hyperon.pdg=PDG::pdg_Lambda;
        h = 1;
      }
      else if(p.hyp_sigma_zero)
      {
        hyperon.pdg=PDG::pdg_Sigma;
        h = 2;
      }
      else
        return 0;
      break;

    case PDG::pdg_neutron: 
      if(p.hyp_sigma_minus)
      {
        hyperon.pdg=PDG::pdg_SigmaM;
        h = 3;
      }
      else
        return 0;
      break;
  }

  // set final particle masses
  hyperon.set_mass(PDG::mass(hyperon.pdg));
  lepton.set_mass(PDG::mass(lepton.pdg));

  double _E_bind=0; //binding energy

  /*
  # 0 is free target; 
  # 1 is Fermi gas; 
  # 2 is local Fermi gas; 
  # 3 is Bodek-Ritchie; 
  # 4 is "effective" spectral function (carbon or oxygen); 
  # 5 is deuterium; 
  # 6 is deuterium with constant binding energy nucleus_E_b (for tests only!)
  */

  // if nucleus conains more than one nucleon, calculate binding energy
  if(t.p + t.n > 1)
  { 
    switch(p.nucleus_target)
    {
      case 0: _E_bind = 0;        break;
      case 1: _E_bind = p.nucleus_E_b; break;
      case 2: _E_bind = t.Ef(N0) + p.kaskada_w;break; //temporary
      case 3: _E_bind = bodek_binding_energy(N0, t.p , t.n); break;
      case 4: _E_bind = binen (N0.p(), p.nucleus_p, p.nucleus_n); break;
      case 5: _E_bind = deuter_binen (N0.p());break; //deuterium 
      case 6: _E_bind = p.nucleus_E_b;      break; //deuterium like Fermi gas
      default: _E_bind=0;
    }
  }

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  N0_Eb.t -= _E_bind;

  double xsec = 0;
  double jakobian=0;  // the value will be set by kinematics generator

  // DEPRECATED
  // double q2 = qel_kinematics(_E_bind, nu, N0, lepton, hyperon, jakobian)
  // double q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, hyperon, jakobian); // simplest choice for hiperon
   double q2 = scatter_2 (nu, N0_Eb, lepton, hyperon);

  if(q2==0) return 0; //indicates interaction is forbidden by kinematics

  vect nu4 = nu;
  
  if(p.hyp_effmass) nu4.boost (-N0_Eb.v());  // go to target frame
  else nu4.boost (-N0.v());  

  double Enu0=nu4.t;    // neutrino energy in target frame   

  // DEPRECATED
  // OLD CZAREK's APPROACH
  // xsec=jakobian*hiperon_sigma(Enu0,-q2,lepton.mass(),hyperon.pdg);
  
  // Four momenta of particles involved
  vect v1(nu);
  vect v2(N0_Eb);
  vect v3(lepton);
  vect v4(hyperon);

  //boost these to CMS for calculations
  vec vcms = (vect(nu) + vect(N0_Eb)).v ();

  v1.boost(-vcms);
  v2.boost(-vcms);
  v3.boost(-vcms);
  v4.boost(-vcms);

  // Generate vector with direction of Lambda in CMS

  vec cms_dir = vec(v3)/vec(v3).length();

  /////////////////////////////////////////////////////////
  // Generate Cross Sections
  /////////////////////////////////////////////////////////

  double kin = v1.length(); //incoming neutrino momentum
  double kout = v3.length(); //outgoing lepton momentum  

  double Nuc_mass;
  double Hyp_mass;

  // Switch to use effective masses of particles in cross section calculation
  if(p.hyp_effmass){
     Nuc_mass = sqrt(N0_Eb*N0_Eb);
     Hyp_mass = sqrt(hyperon*hyperon);
  }
  else {
     Nuc_mass = N0.mass();
     Hyp_mass = hyperon.mass();
  }

  //DEPRECATED
  //double dif = Hyperon_Interaction(-q2,Enu0,h,v1,v2,v3,v4,true);
  //double dif = Singh_Model(-q2,Enu0,h,v1,v2,v3,v4,true);
  double dif = Singh_Model2(-q2,Enu0,h,Nuc_mass,Hyp_mass,lepton.mass(),true);

  // Various coefficients in cross section
  double M2 = Nuc_mass*Nuc_mass;
  double pf = G*G*(1-cos2thetac)/(8*Pi*Enu0*Enu0*M2);

  // Transform from Q2 to costheta
  jakobian = 4*kin*kout; 
  xsec = dif*pf*jakobian;

  // Sigma0 cross section calculation
  if(h==1 && p.hyp_sigma_zero)
  {
    double xsec2=0;

    // take new particles
    particle lepton2 = lepton;
    particle hyperon2 = hyperon;
    hyperon2.pdg = PDG::pdg_Sigma;
    hyperon2.set_mass(PDG::mass(hyperon2.pdg));

    //generate cms blob
    vect cms_Eb = vect(N0_Eb) + vect(nu);

    // Try to solve kinematics for new hyperon, returns false if forbidden
    if( rescale_momenta(cms_Eb,cms_dir,lepton2,hyperon2) )
    {
       v3 = vect(lepton2);
       v4 = vect(hyperon2);

       v3.boost(-vcms);
       v4.boost(-vcms);

       kout = v3.length();

       // calculate q2
       vect p13 = nu - lepton2;
       double q22 = p13 * p13;

       if(p.hyp_effmass)
       {
          Hyp_mass = sqrt(hyperon2*hyperon2);
       }
       else
       {
          Hyp_mass = hyperon2.mass();
       }

       // Calculate the cross section
       //DEPRECATED
       //dif = Singh_Model(-q22,Enu0,2,v1,v2,v3,v4,true);
       dif = Singh_Model2(-q22,Enu0,2,Nuc_mass,Hyp_mass,lepton.mass(),true);

       jakobian = 4*kin*kout; 
       // Prefactor pf only depends on initial state kinematics
       xsec2 = dif*pf*jakobian;
    }

    // choose a proper channel
    if (frandom() < xsec2/(xsec + xsec2))
    {
      lepton = lepton2;
      hyperon = hyperon2;
    }

    // cross sections sum up
    xsec += xsec2;
  }

  if(p.flux_correction)
  {
  // take into account the neutrino flux and nucleon proper time 
  // corrections to the cross section on the whole nucleus
  // int qel_relat=0;

    double vv,vz;
    vv = N0_Eb.v().length ();
    vz = N0_Eb.v() * nu.v().dir(); // this is general
    xsec *= sqrt ( (1 - vz) * (1 - vz) );
  }

  // Set potential 
  if(p.nucleus_target && t.p + t.n > 1) hyperon.set_fermi(t.hyp_BE(hyperon.r.length(),hyperon.pdg));
  else hyperon.set_fermi(0);

  // Check if hyperon can be generated in potential (performed at start of cascade)
  // Interaction didn't happen if not
  if( hyperon.Ek() + hyperon.his_fermi < 0 ){
     return 0;
  }

  e.temp.push_back(lepton);
  e.temp.push_back(hyperon);
  e.out.push_back(lepton);
  e.out.push_back(hyperon);

  e.weight=xsec/cm2;

  return e.weight*cm2;

}
