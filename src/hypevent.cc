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

  switch(N0.pdg)
  {
    case PDG::pdg_proton: // Lambda or Sigma_0 (cross sections should add up)
      if(p.hyp_lambda)
        hyperon.pdg=PDG::pdg_Lambda;
      else
        hyperon.pdg=PDG::pdg_Sigma;
      break;

    case PDG::pdg_neutron: 
      if(p.hyp_lambda)
        return 0;
      else
        hyperon.pdg=PDG::pdg_SigmaM;
      break;
  }

  // set final particle masses
  hyperon.set_mass(PDG::mass(hyperon.pdg));
  lepton.set_mass(PDG::mass(lepton.pdg));

  double E_bind=0;  // binding energy set to 0 // SUBTRACT PROPER ONE!
  double xsec = 0;
  //double jakobian=0;  // the value will be set by kinematics generator

  // DEPRECATED
  // double q2 = qel_kinematics(E_bind, nu, N0, lepton, hyperon, jakobian)
  // double q2 = czarek_kinematics2(E_bind, nu, N0, lepton, hyperon, jakobian); // simplest choice for hiperon
  double q2 = scatter_2 (nu, N0, lepton, hyperon);

  if(q2==0) return 0;

  vect nu4 = nu;
  nu4.boost (-N0.v());  // go to target frame
  double Enu0=nu4.t;    // neutrino energy in target frame   

  // DEPRECATED
  // OLD CZAREK's APPROACH
  // xsec=jakobian*hiperon_sigma(Enu0,-q2,lepton.mass(),hyperon.pdg);

  // NEW CHRIS's APPROACH
  vect v1(nu);
  vect v2(N0);
  vect v3(lepton);
  vect v4(hyperon);
  //boost these to CMS for calculations
  vec vcms = (vect(nu) + vect(N0)).v ();
  v1.boost(-vcms);
  v2.boost(-vcms);
  v3.boost(-vcms);
  v4.boost(-vcms);

  double kin = v1.length(); //incoming neutrino momentum
  double kout = v3.length(); //outgoing lepton momentum

  double rs = sqrt((v1+v2)*(v1+v2)); //CMS energy sqrt(s)
  double pf = (1./(4*v1*v2))*(kout/(4*rs*4*Pi*Pi))*G*G*(1-cos2thetac)/2;
  double dif = Hyperon_Interaction(-q2,Enu0,1,v1,v2,v3,v4);
  xsec = 4*dif*pf*Pi;

  if(p.flux_correction)
  {
  // take into account the neutrino flux and nucleon proper time 
  // corrections to the cross section on the whole nucleus
  // int qel_relat=0;
    double vv,vz;
    vv = N0.v().length ();
    vz = N0.v() * nu.v().dir(); // this is general
    xsec *= sqrt ( (1 - vz) * (1 - vz) );
  }

  e.temp.push_back(lepton);
  e.temp.push_back(hyperon);
  e.out.push_back(lepton);
  e.out.push_back(hyperon);
  e.weight=xsec/cm2;

  return e.weight*cm2;
}


