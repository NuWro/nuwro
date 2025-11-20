#include "jednostki.h"
#include "kinsolver.h"
#include <cassert>
#include "particle.h"
#include "qel_sigma.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "event1.h"
#include "qelevent.h"
#include "kinematics.h"
#include "pdg.h"
#include "nucleus.h"
#include "nucleus_data.h"
#include "dis/LeptonMass.h"
#include <cstdlib>

#define LOCALKF localkf_O

#include "rpa_2013.h"
#include "hyperon_interaction.h"

////////////////////////////////////////////////////////////////////////
double qelevent1(params&p, event & e, nucleus &t,bool nc)
////////////////////////////////////////////////////////////////////////
{
  e.weight=0;

  particle nu=e.in[0];               // neutrino (or electron)
  particle N0=e.in[1];               // initial nucleon
  particle lepton;
  particle N1;

  N1.r=N0.r;
  lepton.r=N0.r;

  int kind=0; // 0 - cc //  1 - nc proton // 2 - nc neutron
  if(nc)
  {
    lepton=nu;
    N1=N0;
    kind=(N0.pdg==pdg_proton?1:2);
  }
  else if((nu.pdg>0 && N0.pdg==PDG::pdg_proton) ||( nu.pdg<0 && N0.pdg==PDG::pdg_neutron))
  {
    return 0;
  }
  else
  {
    lepton.pdg=nu.pdg-(nu.pdg>0 ? 1 :-1);
    lepton.set_mass(PDG::mass(lepton.pdg));
    N1.pdg=(nu.pdg>0 ? PDG::pdg_proton : PDG::pdg_neutron);//zmiana JN
    switch(p.qel_rpa)
    {
      case 2:
      case 3:
        N0.set_mass(t.Mf());
        N1.set_mass(t.Mf());
        break;
      default:
        N1.set_mass(PDG::mass(N1.pdg));
        break;
    }
  }

  if (nu.pdg==11 && N0.pdg==pdg_proton)
    kind=10;
  if (nu.pdg==11 && N0.pdg==pdg_neutron)
    kind=11;//will be used when FF are selected in the file ff.cc

  double _E_bind=0; //binding energy

  if(t.A() > 1)
  {
    switch(p.nucleus_target)
    {
      case 0: _E_bind = 0; break; // free nucleon
      case 1: _E_bind = p.nucleus_E_b; break; // GFG
      case 2: _E_bind = t.Ef(N0) + p.kaskada_w; break; // LFG
      case 3: _E_bind = bodek_binding_energy(N0, t.p, t.n); break; // Bodek-Ritchie
      case 4: _E_bind = binen (N0.p(), p.nucleus_p, p.nucleus_n); break; // effective SF
      case 5: _E_bind = deuter_binen (N0.p()); break; // deuterium
      case 6: _E_bind = 0; break; // effective potential
      default: _E_bind= 0;
    }
  }

  vect aa;
  aa = vect (N0);
  aa.t-=_E_bind;
  double dlu=sqrt( (aa + vect(nu))*(aa+vect(nu)) );
  if (dlu<= (lepton.mass() + N1.mass() ) )
  {
    e.weight=0;
    return 0;
  }

  // cross section (is 0 until the reaction occurs)
  double xsec = 0;
  double q2,jakobian;

  // parameter qel_kinematics is no longer used :)

  if (!(p.nucleus_target==6)) // generic case
    q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, N1,jakobian);
  else // effective momentum dependent potential
    q2=momentum_dependent_potential_kinematics(nu,N0,lepton,N1,jakobian);

  vect nu4 = nu;
  nu4.boost (-N0.v ());  // go to nucleon rest frame
  double Enu0 = nu4.t;   // neutrino energy in target frame
  xsec = jakobian * qel_sigma(Enu0, q2, kind, nu.pdg<0, lepton.mass(), N0.mass());

  /*
  /////////////////////////////////////////////////////////
  // Aligarh Model for QEL (includes second class current)
  /////////////////////////////////////////////////////////

  vect v1(nu); //neutrino
  vect v2(N0); //initial nucleon adjusted for binding energy
  vect v3(lepton); //final lepton
  vect v4(N1); //final nucleon

  vec vcms =(vect(nu) + vect(N0)).v ();

  v1.boost(-vcms);
  v2.boost(-vcms);
  v3.boost(-vcms);
  v4.boost(-vcms);

  double kin = v1.length();
  double kout = v3.length();
  double M2 = v2*v2;

  double pf = G*G*cos2thetac/(8*Pi*Enu0*Enu0*M2);

  bool anti;
  if(nu.pdg > 0){anti = false;}
  else if(nu.pdg < 0){anti = true;}
  else {return 0;}

  xsec = pf*Singh_Model(-q2,Enu0,kind-11,v1,v2,v3,v4,anti)*jakobian;
  */

  // now take into account the neutrino flux and nucleon proper time
  // corrections to the cross section on the whole nucleus
  // int qel_relat=0;

  if(p.flux_correction)
  {
    double vv,vz;
    vv = N0.v().length ();
    vz = N0.v() * nu.v().dir(); // this is general
    xsec *= sqrt ( (1 - vz) * (1 - vz) );
  }

  // (KN) TODO: can N1 have wrong mass here?
  switch(p.qel_rpa)
  {
    case 2:
    case 3:
      N1.set_mass(PDG::mass(N1.pdg)); // go back to real mass
    //N1.set_energy(e.in[1].t+max(nu.t-lepton.t,0.)); // recover the energy conservation
  }

  e.temp.push_back(lepton);
  e.temp.push_back(N1);
  e.out.push_back(lepton);
  e.out.push_back(N1);
  e.weight=xsec/cm2;

  if(!nc && (p.nucleus_target == 1 || p.nucleus_target == 2) && t.A() > 1)
  {
    bool new_ver=true;
    switch(p.qel_rpa)
    {                    //          qv  ,   q0          ,    E        , nu_pdg, lepton_mass  ,    Meff  ,   kF  , version
      case 1:e.weight *= ratio_rpa(e.qv(), e.q0()-_E_bind, nu.t-_E_bind, nu.pdg, lepton.mass(), N1.mass(), t.kF(), new_ver);break;
      case 3:e.weight *= ratio_rpa(e.qv(), e.q0()-_E_bind, nu.t-_E_bind, nu.pdg, lepton.mass(),    t.Mf(), t.kF(), new_ver);break;
    }
  }

  if( p.pauli_blocking) // inlined mypauli_qel
    if(t.pauli_blocking_old(N1, N0.length()))
      e.weight = 0;

  // selection of events for electron scattering using acceptance information
  if (nu.pdg==11)
  {
    double kosine=lepton.z/lepton.momentum();
    if ( kosine < (p.el_costh_lab-p.el_costh_del) || kosine > (p.el_costh_lab+p.el_costh_del) )
      e.weight=0;
    else
      e.weight*=Pi2*8.0/137.03599908/137.03599908/G/G/q2/q2 / (2*p.el_costh_del);
      // change of propagators in weak and em processes
      // KN: the output should be the differential in costh!
      //     (so the user doesn't have to remember what the width was)
  }

  return e.weight*cm2;
}
