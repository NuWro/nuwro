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


 
  particle N0_Eb = N0;
  N0_Eb.t -= _E_bind;

  double xsec = 0; double xsec2 = 0;
  //double jakobian=0;  // the value will be set by kinematics generator

  // DEPRECATED
  // double q2 = qel_kinematics(_E_bind, nu, N0, lepton, hyperon, jakobian)
  // double q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, hyperon, jakobian); // simplest choice for hiperon
  double q2 = scatter_2 (nu, N0_Eb, lepton, hyperon);


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
double dif = Hyperon_Interaction(-q2,Enu0,h,v1,v2,v3,v4);
  xsec = 4*dif*pf*Pi;

 
  //for testing of the sigma zero cross section calculation using the 
  //code below, set xsec to zero
  //xsec=0;

  // if the interaction was on proton, we calculated Lambda and we still was Sigma_0
  // we need to calcualte the cross section and weight the choice
  
  //check if sigma production is allowed by kinematics
  if( h == 1 && p.hyp_sigma_zero && rs > PDG::mass(PDG::pdg_Sigma)+lepton.mass())
  {
    // find the lepton unit vector
    vec unit_vec = v3 / v3.length();
    
    //3 momentum of outgoing sigma minus   
  double pp = cms_momentum2(rs*rs,lepton.mass2(),PDG::mass(PDG::pdg_Sigma)*PDG::mass(PDG::pdg_Sigma));   

 // take new particles
    particle lepton2 = lepton;
    particle hyperon2 = hyperon;
    hyperon2.pdg = PDG::pdg_Sigma;
    hyperon2.set_mass(PDG::mass(hyperon2.pdg));
    // adjust momentum
    lepton2.set_momentum (unit_vec * sqrt(pp));
    hyperon2.set_momentum (-unit_vec * sqrt(pp));
    // boost to LAB
    lepton2.boost (vcms);
    hyperon2.boost (vcms);

    // calculate q2
    vect p13 = nu - lepton2;
    double q22 = p13 * p13;

    // calculate the cross section

    //4 momenta of outgoing lepton and hyperon in cms
    v3 = vect(sqrt(lepton2.mass()*lepton2.mass()+pp),unit_vec*sqrt(pp));
    v4 = vect(sqrt(hyperon2.mass()*hyperon2.mass()+pp),-unit_vec*sqrt(pp));

    kout = sqrt(pp);

    pf = (1./(4*v1*v2))*(kout/(4*rs*4*Pi*Pi))*G*G*(1-cos2thetac)/2;
    dif  = Hyperon_Interaction(-q22,Enu0,2,v1,v2,v3,v4);

    xsec2 = 4*dif*pf*Pi;


    // range in q2
    double range = 4 * nu.E() * lepton.momentum();
    double range2 = 4 * nu.E() * lepton2.momentum();
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

  e.temp.push_back(lepton);
  e.temp.push_back(hyperon);
  e.out.push_back(lepton);
  e.out.push_back(hyperon);
  e.weight=xsec/cm2;

    return e.weight*cm2;
}


