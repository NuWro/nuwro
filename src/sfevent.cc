#include <fstream>

#include "event1.h"
#include "params.h"
#include "nucleus.h"
#include "sf/GConstants.h"
#include "sf/CSFOptions.h"
#include "sf/CSpectralFunc.h"
#include "kinematics.h"

static inline double pow2(double x)
{return x*x;
}

bool has_sf(nucleus &t)
{  switch(1000*t.Z()+t.N())
   {
	   case  6006:
	   case  8008:
	   case 18022:
	   case 20020:
	   case 26030:
	     return true;
	   default: 
	     return false;
   }
}

double sfevent2(params&par, event & e, nucleus &t)
{   


    particle l0=e.in[0];		           // neutrino
    particle N0=e.in[1];		           // initial nucleon
    particle l1;  
    particle N1;

	if( (l0.pdg>0 and N0.pdg==pdg_proton)
	    or (l0.pdg<0 and N0.pdg==pdg_neutron)
	  )
	   return 0; // no CC interaction possible on this nucleon

	if(l0.pdg<0) 
	  {N1.set_neutron();
	   l1.pdg=l0.pdg+1;
      }
	else 
	  {N1.set_proton();
	   l1.pdg=l0.pdg-1;
      }
		
	double m = 0;
	switch(abs(l1.pdg))
	{case pdg_e:   m = mass_e; break;
	 case pdg_mu:  m = mass_mu; break;
	 case pdg_tau: m = mass_tau; break;
    }
	
	l1.set_mass(m);
	
	double mm=m*m;
	double M=N1.m();
	double MM=M*M;
	
	e.in[1]=N0;
	
//	cout<<"Options to create"<<endl;
	CSFOptions options(par,1,N0.pdg==pdg_proton,l0.pdg<0); // czy proton // czy antyneutrino
//	cout<<"sf about to create"<<endl;
	CSpectralFunc* sf = options.get_SF();
//	cout<<"sf created"<<endl;
	const double pBlock = sf->get_pBlock() ;
//	cout<<"sf used"<<endl;
	
        
	const double p = sf->MomDist()->generate() ;

	N0.set_momentum(rand_dir()*p);		
//    cout<<"ola"<<endl;
    double EBEN;
    do{
          EBEN=sf->generateE(p);
    } while(!(EBEN==EBEN));
     
    vect s=l0+N0;
        
    s.t=l0.E()+N0.mass()-EBEN;
        
    double ss=s*s;
 
    if(ss<pow2(M+m))
        return 0;

    vec v=s.v();    
 

/// here we do decay(s,N1,l1) by hand    

    double pcms=sqrt(0.25*pow2(ss+mm-MM)/ss-mm);
    vec dircms=rand_dir();
    vec p1=pcms*dircms;
    l1.set_momentum(p1);
    N1.set_momentum(-p1);
    l1.boost(v);
    N1.boost(v); 
      
       
	const double omega= l0.E() - l1.E();
		
	const double omegaTil = N1.E() - N0.E() ; 
		

    if(false)
	if ( omegaTil < 0)
		return 0;
    if(false)
    if( omega < 0 )   
		return 0;
	if ( par.pauli_blocking and  N1.momentum()<pBlock  )
		return 0;
 	    
    vect q4til=N1-N0;
	const double q4til2= q4til*q4til;
	const double p4k4= l0*N0 ;
    if(false)
	if ( p4k4 < 0 or q4til2 > 0 ) 
			return 0; 

    double pp=pcms*pcms;

    double vol=4*pi*pp;

    double gamma=1/sqrt(1-v*v);
    double graddelta= (l1.v()-N1.v()).length();
    double surfscale=sqrt(1+(1-pow2(v*dircms)/(v*v))*(gamma*gamma-1));
	double	val = G*G*cos2thetac/8/pi/pi           
	          *vol *(surfscale/graddelta)
	          /(l1.E()*	
	            l0.E()*
	            N0.E()*
	            N1.E()
	           )
			 *options.evalLH(q4til*q4til,l0*N0,l1*N0,q4til*N0,l0*q4til,l1*q4til,l0*l1);
		
    e.weight=val/cm2;
    N0.t=N0.mass()-EBEN;
    e.in[1]=N0;   
	e.out.push_back(l1);
	e.out.push_back(N1);
  // cout<<val/cm2<<endl<<endl;
   return val;
}


double sfevent2nc(params&par, event & e, nucleus &t)
{   
//	cout<<"sf created"<<endl;
  //  cout<<"sf used"<<endl;
	

    particle l0=e.in[0];		           // neutrino
    particle N0=e.in[1];		           // initial nucleon
    particle l1=l0;  
    particle N1=N0;

	
	double m = 0;
	double mm=m*m;
	double M=N1.m();
	double MM=M*M;
	
	e.in[1]=N0;
	
	CSFOptions options(par,0,N0.pdg==pdg_proton,l0.pdg<0); // czy proton // czy antyneutrino
	CSpectralFunc* sf = options.get_SF();
	double pBlock = sf->get_pBlock() ;
        
	const double p = sf->MomDist()->generate() ;

	N0.set_momentum(rand_dir()*p);		

    double EBEN;
    do{
          EBEN=sf->generateE(p);
    } while(!(EBEN==EBEN));
     
    vect s=l0+N0;
        
    s.t=l0.E()+N0.mass()-EBEN;
        
    double ss=s*s;
 
    if(ss<pow2(M+m))
        return 0;

    vec v=s.v();    
 

/// here we do decay(s,N1,l1) by hand    

    double pcms=sqrt(0.25*pow2(ss+mm-MM)/ss-mm);
    vec dircms=rand_dir();
    vec p1=pcms*dircms;
    l1.set_momentum(p1);
    N1.set_momentum(-p1);
    l1.boost(v);
    N1.boost(v); 
      
       
	const double omega= l0.E() - l1.E();
		
	const double omegaTil = N1.E() - N0.E() ; 
		

    if(false)
	if ( omegaTil < 0)
		return 0;
    if(false)
    if( omega < 0 )   
		return 0;
	if ( par.pauli_blocking and  N1.momentum()<pBlock  )
		return 0;
 	    
    vect q4til=N1-N0;
	const double q4til2= q4til*q4til;
	const double p4k4= l0*N0 ;
    if(false)
	if ( p4k4 < 0 or q4til2 > 0 ) 
			return 0; 

    double pp=pcms*pcms;

    double vol=4*pi*pp;

    double gamma=1/sqrt(1-v*v);
    double graddelta= (l1.v()-N1.v()).length();
    double surfscale=sqrt(1+(1-pow2(v*dircms)/(v*v))*(gamma*gamma-1));
	double	val = G*G/8/pi/pi           
	          *vol *(surfscale/graddelta)
	          /(l1.E()*	
	            l0.E()*
	            N0.E()*
	            N1.E()
	           )
			 *options.evalLHnc(q4til*q4til,l0*N0,l1*N0,N0*q4til,l0*q4til,l1*q4til,l0*l1);
		
    e.weight=val/cm2;
    N0.t=N0.mass()-EBEN;
    e.in[1]=N0;   
	e.out.push_back(l1);
	e.out.push_back(N1);
   return val;
}
