#include "kinsolver.h"
#include <cassert>          
#include "particle.h"
//#include "qelcc.h"
//#include "qelnc.h"
#include "jednostki.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "scatter.h"

#define LOCALKF localkf_O

//////////////////////////////////////////////////////////////////////////////////////////
double jakob(vect nu,vect N0,vect lepton)
//////////////////////////////////////////////////////////////////////////////////////////
{
  vec vcms = (vect(nu) + vect(N0)).v ();
  if(vcms*vcms>=1)
     return 0;
  else
  {
    nu.boost (-vcms);
    lepton.boost (-vcms);	
    return vec (nu).length () * vec (lepton).length () * 4;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
double bodek_kinematics(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double & jakobian)
//////////////////////////////////////////////////////////////////////////////////////////
{
  static double M12=0.5*(PDG::mass_proton+PDG::mass_neutron);
  double Ma = 56 * M12;	// Argon nucleus mass 
  double Ma1 = Ma - M12 + Eb;  // mass of the spectator nucleus
  N0.t = 56 * M12 - sqrt ( Ma1*Ma1 + N0.momentum2 () );
  //if(N0*N0<0) return 0;
  double q2 = scatter_2 (nu, N0, lepton, N1);
  jakobian=jakob(nu,N0,lepton);

  return q2;
}

//////////////////////////////////////////////////////////////////////////////////////////
double bodek_binding_energy(particle N0, int A)
//////////////////////////////////////////////////////////////////////////////////////////
{
  static double M12=0.5*(PDG::mass_proton+PDG::mass_neutron);
  double Ma  = A * M12;
  double Ma1 = (A-1) * M12;

  double Eb = sqrt(Ma1*Ma1 + N0.momentum2()) - Ma + sqrt(M12*M12 + N0.momentum2());

  return Eb;
}


//////////////////////////////////////////////////////////////////////////////////////////
double czarek_kinematics(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double &jakobian)
//////////////////////////////////////////////////////////////////////////////////////////
{
  N0.t -= Eb;
  double q2 = scatter_2 (nu, N0, lepton, N1);
  jakobian = jakob(nu,N0,lepton);
  return q2;
}

// This version os faster due to enhancement of the forward direction
// in the decay 2 function
//////////////////////////////////////////////////////////////////////////////////////////
double czarek_kinematics2(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double &jakobian)
//////////////////////////////////////////////////////////////////////////////////////////
{
  N0.t -= Eb;
  if(not decay2(nu, nu+N0, lepton, N1,jakobian))
    return 0;
  jakobian *=jakob(nu,N0,lepton)/2;
  vect q=nu-lepton;
  return q*q;
}


double kF_potencjal=200*MeV;
double kF_P=kF_potencjal;
double kF_N=kF_potencjal;

//////////////////////////////////////////////////////////////////////////////////////////
double V(double p,double kf_t)
//////////////////////////////////////////////////////////////////////////////////////////
{
  double kF =  kf_t; // zeby nie mieszac w dalszych czesciach programu!!!
  double kF2= kF*kF;

  double p2 =   p*p;
  double p4 = p2*p2;

  double const a = 206.;
  double const a2=  a*a;

  double const b =   582.;
  double const c =  -322.;
  double const c4=c*c*c*c;

  double const d = 289.;
  double const d3=d*d*d;

  double const e =  442.; 
  double const e3= e*e*e;

  return (-1)*a2*kF*kF*(kF+b)/(c4+e3*kF+d3*p2/kF+p4);
}


// finds the cms momentum and sets the outgoing particles 
//////////////////////////////////////////////////////////////////////////////////////////
double 
find_momentum(double E_lab,
              vec cms_p_dir,
	      vec cmsspeed,
	      particle nu,
	      particle N0,
	      particle & lepton, 
	      particle &N1)
//////////////////////////////////////////////////////////////////////////////////////////
{ 
//  cout<< "Finding momentum ..."<<endl;
  kinsolver k(E_lab,cms_p_dir,cmsspeed, nu, N0, lepton, N1);
  double p=k.findmomentum();
//  cout<< "Found momentum p="<<p<<endl;
  return p;
}


//////////////////////////////////////////////////////////////////////////////////////////
double momentum_dependent_potential_kinematics(
          particle nu,
	  particle N0, 
	  particle &lepton, 
	  particle &N1, 
          double & jakobian)
//////////////////////////////////////////////////////////////////////////////////////////
{    particle l1=lepton;
     particle N2=N1;
   
    vec vcms=(vect(nu)+vect(N0)).v();
    if(vcms*vcms>1) throw("faster than light!!!!");
    // cout<<nu<<N0<<"CMSSPPED"<<vcms<<endl;
    
    // Initial energy in the lab
    double E_lab=nu.t+N0.t+V(N0.momentum(),kF_N);

    // decay direction in cms
    vec dir1=rand_dir();
#ifdef DEBUG_dE
    cout<<"dir1="<<dir1<<endl;
#endif    
    double p1=find_momentum(E_lab, dir1, vcms, nu, N0, lepton, N1);
#ifdef DEBUG_dE
    cout<<"DIF="<<nu+N0-lepton-N1<<endl;			    
#endif
    if(p1<0) return 0;			    
    double q2a=q2(nu,lepton);
    
    vec dir2=(dir1+vcms.dir()/100).dir();
    double p2=find_momentum(E_lab,dir2, vcms, nu, N0, l1, N2);
    if(p2<0) return 0;			    
    vect cc=nu;
    cc.boost(-vcms); 
    double q2b=q2(nu,l1);
    
    double jakobian1=(q2a-q2b)/(cos(dir1,vec(cc))-cos(dir2,vec(cc)))*2;
    
#ifdef DEBUG_jakobian
    double jakobian2=jakob(nu,N0,lepton);
    cout<<"jakobian1/jakobian2="<<jakobian1/jakobian2<<endl;
#endif

    jakobian=jakobian1;
    return q2a;			    
} 




