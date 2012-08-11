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

#include "rpa_lib.h"
//double qelm;

//static double E_b(int choice,ped,

////////////////////////////////////////////////////////////////////////
double qelevent1(params&p, event & e, nucleus &t,bool nc)
////////////////////////////////////////////////////////////////////////
  { //cout<<"jestem w qelevent"<<endl;
    e.flag.qel=true;
    e.flag.nc=nc;
    e.flag.cc=!nc;
    e.flag.dis=false;
    e.weight=0;

    particle nu=e.in[0];		           // neutrino
    particle N0=e.in[1];		           // initial nucleon
    particle lepton;  
    particle N1;
    int kind=0; // 0 - cc //  1 - nc proton // 2 - nc neutron
    if(nc)
    {
    	lepton=nu;
     	N1=N0;
     	kind=(N0.pdg==pdg_proton?1:2);
    }
    else 
     if((nu.pdg>0 && N0.pdg==PDG::pdg_proton) ||( nu.pdg<0 && N0.pdg==PDG::pdg_neutron))
       return 0;
     else  
     { 
       lepton.pdg=nu.pdg-(nu.pdg>0 ? 1 :-1);
       lepton.set_mass(PDG::mass(lepton.pdg));
       N1.pdg=(nu.pdg>0 ? PDG::pdg_proton : PDG::pdg_neutron);//zmiana JN
       N1.set_mass(PDG::mass(N1.pdg));
     } 
	double qelmlep=lepton.mass();
    //cout<<"masa leptonu= "<<lepton.mass()<<endl;

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
	case 2: _E_bind=0;              break; //temporary
	case 3: _E_bind=0;              break; //temporary
	case 4: _E_bind = binen (ped, p.nucleus_p, p.nucleus_n);	
 		     //cout<<ped<<"  "<<_E_bind<<endl;//SF
	         break;
	case 5: _E_bind= deuter_binen (ped);break; //deuterium 
	case 6: _E_bind= p.nucleus_E_b; 	break; //deuterium like Fermi gas
	default: _E_bind=0;//in the future it is possible to add SF for other nuclei as well
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


    int qel_cc_japan=0;
    QEL qel(lepton.mass(),
            //N0.mass(),
            (PDG::mass_proton+PDG::mass_neutron)/2, 
            lepton.pdg<0,   //antyneutrino
            qel_cc_japan);
    
    // cross section (is 0 until the reaction occurs)   
    double x = 0;		
    double q2,jakobian;  
    
    int qel_kinematics=0;  // force standard behaviour
    
    switch(qel_kinematics)
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
	            N0.set_mass(M12-27*MeV);  
		        // next do the usual kinematics 
		        // preserving energy and momentum
		        // (not taking into account the spactator nucleus)
                q2 = czarek_kinematics(_E_bind  //chyba powinno by�0 bo ju odj�e
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
	             N1.set_energy(N1.energy()+V(N1.momentum(),t.localkf(N1.pdg,N1.r.length())));
                 break;	       	  
    default: 
                 cerr<<"Kinematics code: '"<<qel_kinematics<<"' is invalid."<<endl;	  
	             exit(1);
    }  
//    cout<<"q2="<<q2<<endl;
    /* TA CZESC BYLA OK , ZROBIAONA PRZEZ CZARKA
    if (q2 != 0  		   // there was scattering
        &&                         // and
        N1.momentum () >= kF   // Pauli-blocking did not occur	        )
       )
    */
    //!!!!!!!!ZMIANA JARKA PED FERMIEGO
    //e.flag.q2=q2;
    //e.flag.W=N1.mass();
    if (q2 != 0  		           // there was scattering
//        &&                       // and
//        !t.pauli_blocking(N1)    // Pauli-blocking did not occur	        
       )

    if(qel_kinematics==3)
      N1.set_energy(N1.energy()+V(N1.momentum(),t.localkf(N1.pdg,N1.r.length())));
      
      
    vect nu4 = nu;
    nu4.boost (-N0.v ());     // calculate in  nucleon rest frame
    qel.set_energy (nu4.t);       // neutrino energy in nucleon rest frame
    x=qel.sigma (q2,kind) *jakobian;
	
     // now take into account the neutrino flux and nucleon proper time 
     // corrections to the cross section on the whole nucleus
     // int qel_relat=0;
     
     if(p.flux_correction)
     {
		double vv,vz;
		vv = N0.v().length ();
		vz = N0.v() * nu.v().dir(); // this is general
		x *= sqrt ( (1 - vz) * (1 - vz) );
      }  

     e.temp.push_back(lepton);
     e.temp.push_back(N1);
     e.out.push_back(lepton);
     e.out.push_back(N1);
     e.weight=x/cm2;

     if(!nc)
		 if(p.qel_rpa)
		 {  //nucleus_data *data=best_data(t.p,t.n);
			rpa::configure(e.in[0].t, 0, e.in[0].pdg, p.qel_rpa==2, 0,t.kF(),t.Mf());
			double rpa_frac=rpa::ratio_rpa_fg(e.in[0].t,e.q0(),e.qv());
			
			if(rpa_frac==rpa_frac) //assert no nan
				e.weight*=rpa_frac;
		 } 
	//cout<<"kFP="<<kF_P<<endl;

	//cout<<"CR="<<e.weight<<endl;
	//cout<<(e.flag.nc ? "nc " : "cc " )<<x/cm2<<endl;
     return x;
  }

