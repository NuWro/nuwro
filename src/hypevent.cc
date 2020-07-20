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

  int h; // hyperon produced: 1,2,3
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
      case 0: _E_bind=0;        break;
      case 1: _E_bind= p.nucleus_E_b; break;
      case 2: _E_bind= t.Ef(N0) + p.kaskada_w;break; //temporary
      case 3: _E_bind=0;              break;         //temporary
      case 4: _E_bind = binen (N0.p(), p.nucleus_p, p.nucleus_n);
           //in the future it is possible to add SF for other nuclei as well  
           //cout<<ped<<"  "<<_E_bind<<endl;//SF
           break;
      case 5: _E_bind= deuter_binen (N0.p());break; //deuterium 
      case 6: _E_bind= p.nucleus_E_b;      break; //deuterium like Fermi gas
      default: _E_bind=0;
    }
  }

  // to force zero binding energy
  // _E_bind=0; 

  particle N0_Eb = N0; // nucleon with 4 momentum adjusted for binding energy
  N0_Eb.t -= _E_bind;

  //generate hyperon binding energy as a function of 
  //local nuclear density

  double Y_Eb=0;

  if(p.hyp_Eb && p.nucleus_target && t.p + t.n > 1) Y_Eb = t.hyp_BE(hyperon.r.length());

  double xsec = 0;
  double jakobian=0;  // the value will be set by kinematics generator

  //need to ensure the sigma is also gnerated going in the same direction in the cms for 
  //comparison of differentials, to do this generate a unit vector in the cms and boost into the
  //modified cms (defined in scatter.cc). Modified cms is the same frame for both sigma
  //and lambda production, record the direction of this unit vector in the modified cms

  vec cms_dir; //scattering angle in bound CMS

  // double q2_BE = scatter2_with_BE_SC(nu,N0_Eb,lepton,hyperon_Eb,Y_Eb);
  double q2 = scatter2_with_BE(nu,N0_Eb,lepton,hyperon,Y_Eb,cms_dir);

  // DEPRECATED
  // double q2 = qel_kinematics(_E_bind, nu, N0, lepton, hyperon, jakobian)
  // double q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, hyperon, jakobian); // simplest choice for hiperon
  // double q2 = scatter_2 (nu, N0_Eb, lepton, hyperon);

  if(q2==0) return 0; //indicates interaction is forbidden by kinematics

  vect nu4 = nu;
  nu4.boost (-N0_Eb.v());  // go to target frame
  double Enu0=nu4.t;    // neutrino energy in target frame   

  // DEPRECATED
  // OLD CZAREK's APPROACH
  // xsec=jakobian*hiperon_sigma(Enu0,-q2,lepton.mass(),hyperon.pdg);

  // NEW CHRIS's APPROACH
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

  /////////////////////////////////////////////////////////
  // Generate Cross Sections
  /////////////////////////////////////////////////////////

  double kin = v1.length(); //incoming neutrino momentum
  double kout = v3.length(); //outgoing lepton momentum  
  double rs = sqrt((v1+v2)*(v1+v2));

  //old version
  //double dif = Hyperon_Interaction(-q2,Enu0,h,v1,v2,v3,v4,true);
  
  double dif = Singh_Model(-q2,Enu0,h,v1,v2,v3,v4,true);

  double M2 = v2*v2;

  double pf = G*G*(1-cos2thetac)/(8*Pi*Enu0*Enu0*M2);
  jakobian = 4*kin*kout; 

  xsec = dif*pf*jakobian;
    
  /////////////////////////////////////////////////////////////////
  // for testing of the sigma zero cross section calculation using the 
  // code below, set xsec to zero
  // xsec=0;
  // if the interaction was on proton, we calculated Lambda and we still was Sigma_0
  // we need to calcualte the cross section and weight the choice
  // check if sigma production is allowed by kinematics

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

    // bool rescale = rescale_momenta(cms_Eb,cms_dir,lepton2,hyperon2,Y_Eb);

    //reuturns false if sigma0 prod forbidden by kinematics

    //try to solve kinematics for new hyperon, returns false if forbidden
    if( rescale_momenta(cms_Eb,cms_dir,lepton2,hyperon2,Y_Eb) )
    {
      v3 = vect(lepton2);
      v4 = vect(hyperon2);

      v3.boost(-vcms);
      v4.boost(-vcms);

      kout = v3.length();

      // calculate q2
      vect p13 = nu - lepton2;
      double q22 = p13 * p13;

      // calculate the cross section
      dif = Singh_Model(-q22,Enu0,2,v1,v2,v3,v4,true);
 
      jakobian = 4*kin*kout; 
      //prefactor pf is the same
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

  //return final particles to their mass shells
  hyperon.set_mass(PDG::mass(hyperon.pdg));

  e.temp.push_back(lepton);
  e.temp.push_back(hyperon);
  e.out.push_back(lepton);
  e.out.push_back(hyperon);
  e.weight=xsec/cm2;

  return e.weight*cm2;
}
