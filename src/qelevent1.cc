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

//double qelm;

//static double E_b(int choice,ped,

////////////////////////////////////////////////////////////////////////
double qelevent1(params&p, event & e, nucleus &t,bool nc)
////////////////////////////////////////////////////////////////////////
{
    e.flag.qel=true;
    e.flag.nc=nc;
    e.flag.cc=!nc;
    e.flag.dis=false;
    e.weight=0;

    particle nu=e.in[0];		           // neutrino
    particle N0=e.in[1];		           // initial nucleon
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

	double _E_bind=0;	//binding energy

	/*
	# 0 is free target; 
	# 1 is Fermi gas; 
	# 2 is local Fermi gas; 
	# 3 is Bodek-Ritchie; 
	# 4 is "effective" spectral function (carbon or oxygen); 
	# 5 is deuterium; 
	# 6 is deuterium with constant binding energy nucleus_E_b (for tests only!)
	*/
	
	vec ped=N0.p();
	switch(p.nucleus_target)
	{
		case 0: _E_bind=0; 				break;
		case 1: _E_bind= p.nucleus_E_b;	break;
		case 2: _E_bind= t.Ef(N0) + p.kaskada_w;              break; //temporary
		case 3: _E_bind=0;              break; //temporary
		case 4: _E_bind = binen (ped, p.nucleus_p, p.nucleus_n);
				 //in the future it is possible to add SF for other nuclei as well	
				 //cout<<ped<<"  "<<_E_bind<<endl;//SF
				 break;
		case 5: _E_bind= deuter_binen (ped);break; //deuterium 
		case 6: _E_bind= p.nucleus_E_b; 	break; //deuterium like Fermi gas
		default: _E_bind=0;
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
    
//    int qel_kinematics=0;  // force standard behaviour
    
//    q2=qel_kin_gen(qel_kinematics,_E_bind, nu, N0,lepton, N1, jakobian);
    switch(p.qel_kinematics)
    { 
    case  0:    // Subtract the binding energy from the nucleon energy
		        // cout<<"jestem w switch qel-kinematics"<<endl;
		        
//                q2 = czarek_kinematics(_E_bind, nu, N0, lepton, N1,jakobian);
                q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, N1,jakobian);

//              czarek_kinematics2 is a few times faster due to qadratic 
//              enhancement of the forward direction in decay2
//              it also gives smaller errors
		        // cout<<"Q2=  "<<-q2<<endl;
                break;
		
    case  1:    // preservation of energy momentum for the whole system 
                // the nucleon is off shell 
		        // binding energy used for spectator nucleus mass calculation
                q2 = bodek_kinematics (_E_bind, nu, N0, lepton, N1,jakobian);
 	            break;
    case  2:    // effective mass trick kinematics the most doubtful    
                // first lower the nucleon mass down the binging energy
	            N0.set_mass((mass_proton+mass_neutron)/2-27*MeV);  
		        // next do the usual kinematics 
		        // preserving energy and momentum
		        // (not taking into account the spactator nucleus)
                q2 = czarek_kinematics(_E_bind  // should be 0 here??? 
		                      , nu, N0, lepton, N1, jakobian);		
	         
		        // use this mass also for dynamics????
		        // should it appear in the form factors
		        break;
    case  3:    // use momentum dependent potential V(p)
                // V(p) is on both sides of the equation but with different p
                // N1.set_energy(N1.energy()+V(N1.momentum(),kF_P));
                q2=momentum_dependent_potential_kinematics(nu,N0,lepton,N1,jakobian);
                break;	       	  
    case  4:    // use momentum dependent potential V(p)
                // V(p) is on both sides of the equation but with different p
	            // then adjust the kinetic energy 
                 q2=momentum_dependent_potential_kinematics(nu,N0,lepton,N1,jakobian);
	             N1.set_energy(N1.energy()+V(N1.momentum(),t.localkf(N1)));
                 break;	       	  
    default: 
                 cerr<<"Kinematics code: '"<<p.qel_kinematics<<"' is invalid."<<endl;	  
	             exit(9);
    }  
    if (q2 != 0)  		           // there was scattering
		if(p.qel_kinematics==3)
			N1.set_energy(N1.energy()+V(N1.momentum(), t.localkf(N1)));
      
      
	vect nu4 = nu;
	nu4.boost (-N0.v ());  // go to nucleon rest frame 
	double Enu0=nu4.t;     // neutrino energy in target frame   
	xsec=jakobian* qel_sigma (Enu0, q2, kind, nu.pdg<0, lepton.mass(), N0.mass()); 
	
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


	switch(p.qel_rpa) 
	{
		 case 2:
		 case 3:
			N1.set_mass(PDG::mass(N1.pdg));// go back to real mass 
			//N1.set_energy(e.in[1].t+max(nu.t-lepton.t,0.)); // recover the energy conservation 
	}
		
	e.temp.push_back(lepton);
	e.temp.push_back(N1);
	e.out.push_back(lepton);
	e.out.push_back(N1);
	e.weight=xsec/cm2;

	if(!nc)
	{
		bool new_ver=true;
		switch(p.qel_rpa)
 		{                         //     qv   ,  q0            ,  E         , nu_pdg, lepton_mass  ,  Meff     , kF   , version                               
			case 1:e.weight*= ratio(e.qv(), e.q0()-_E_bind, nu.t-_E_bind, nu.pdg, lepton.mass(), N1.mass(), t.kF(), new_ver);break;
			case 3:e.weight*= ratio(e.qv(), e.q0()-_E_bind, nu.t-_E_bind, nu.pdg, lepton.mass(),    t.Mf(), t.kF(), new_ver);break;
		} 
	}
	if( p.pauli_blocking) 	/// inlined mypauli_qel
		if(t.pauli_blocking_old (N1, N0.length() ) ) 
			e.weight = 0;
		
	return e.weight*cm2;
}
