#include "jednostki.h"
#include <cassert>          
#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "event1.h"
#include "hipevent.h"
#include "kinematics.h"
#include "pdg.h"
#include "nucleus.h"
#include <cstdlib>
#include "hiperon_sigma.h"
#define LOCALKF localkf_O

using namespace NUWRO;


double hipevent(params&p, event & e, nucleus &t,bool lambda)
{
   /*
    * 
    * This is hiperon production code. Only for antineutrino scattering.
    * 
    * p - global nuwro parameters 
    * e - event object to modify 
    * t - nucleus on which scattering occurs
    * lambda = true  // Lambda production
    * lambda = false // Sigma production
    * returns the scattering cross section, produces the kinematics
    * and stores the outgoing particles in the e.out vector
    * 
    */ 
    
    e.flag.hip=true;
    e.flag.cc=true;
    e.weight=0;

    if(e.in[0].pdg>0) 
	return 0; // This Dynamics is for antineytrinos only

    particle nu=e.in[0];	// initial antineutrino	
    particle N0=e.in[1];	// initial nucleon
    particle lepton;   		// outgoing antilepton
    particle N1;       		// created hiperon
    N1.r=N0.r;
    lepton.r=N0.r;
    if(nu.pdg>0)      // only anty neutrino possible
    	return 0;
	
    lepton.pdg=nu.pdg-(nu.pdg>0 ? 1 :-1);
    
    switch(N0.pdg)
    {
	case PDG::pdg_proton: // Lambda or Sigma_0 (cross sections should add up)
	    if(lambda)
		    N1.pdg=PDG::pdg_Lambda; // (1) in 0606235
	    else
		    N1.pdg=PDG::pdg_Sigma;  // (2) in 0606235
	    break;

	case PDG::pdg_neutron: 
	    if(lambda)
		return 0;
	    else
		N1.pdg=PDG::pdg_SigmaM;     // (3) in 0606235
	    break;
    } 

    N1.set_mass(PDG::mass(N1.pdg));
    lepton.set_mass(PDG::mass(lepton.pdg));

    double E_bind=0;	// binding energy set to 0 
    double xsec = 0;		
    double jakobian=0;  // the value will be set by kinematics generator
    //double q2 = qel_kinematics(E_bind, nu, N0, lepton,N1, jakobian)
    double q2 = czarek_kinematics2(E_bind, nu, N0, lepton, N1,jakobian); // simplest choice for hiperon

    if(q2==0)
    	return 0;
    vect nu4 = nu;      
    nu4.boost (-N0.v());  // go to target frame
    double Enu0=nu4.t;     // neutrino energy in target frame   
    
    xsec=jakobian*hiperon_sigma(Enu0,-q2,lepton.mass(),N1.pdg); 

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
    e.temp.push_back(N1);
    e.out.push_back(lepton);
    e.out.push_back(N1);
    e.weight=xsec/cm2;
    
    return e.weight*cm2;
}


