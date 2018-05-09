#include "e_el_event.h"
#include "params.h"
#include "nucleus.h"
#include "event1.h"
#include "espp/crossection.h"
#include "kinematics.h"


double e_spp_event(params&p, event & e, nucleus &t, bool nc)
{
	//
	//  electron  elastic scattering 
	//
	
	e.weight=0; // event impossible before setting non zero cross section

	particle l0=e.in[0];  // in electron
	particle N0=e.in[1];  // in nucleon
	particle l1=l0;       // out lepton (good for nc)
	particle Delta=N0;    // intermadiate delta 
	particle N1=N0;       // nucleon from delta decay
	particle pion(pdg_pi);// pion from delta decay
	

	double _E_bind=0;	// binding energy				
	
	_E_bind=t.Ef(); // calculate binding energy 
  
	vect p4= l0+N0;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////                                 CHANNELS:                          ///////////////////////////////
	///////     1: e + p -> e + p + \pi^\0                                     ///////////////////////////////
	///////     2: e + n -> e + p + \pi^-                                      ///////////////////////////////
	///////     3: e + p -> e + n + \pi^+                                      ///////////////////////////////
	///////     4: e + n -> e + n + \pi^0                                      ///////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(frandom()<0.5)
	{
		N1.set_proton();
		if(N0.pdg==pdg_proton)
		{
			chan=1;
			pion.set_pi();
		}
		else
		{
			chan=2;
			pion.set_piM();
		}
	}
	else
	{
		N1.set_neutron();
		if(N0.pdg==pdg_proton)
		{
			chan=3;
			pion.set_piP();
		}
		else
		{
			chan=4;
			pion.set_pi();
		}
	}



	double Wmax=sqrt(p4*p4)-l0.mass();
	double Wmin=N1.mass()+pion.mass();
	
	double W=Wmin+frandom()*(Wmax-Wmin);
	Delta.set_mass(W);

	double jakobian;  
	double q2 = czarek_kinematics2(_E_bind, l0, N0, l1, Delta,jakobian);
	if(q2==0)
	{
		return 0;
	}
	
	
	vect nu4 = l0;         // electron fourmonentum 
	nu4.boost( -N0.v() );  // go to nucleon rest frame 
	double Ee0=nu4.t;      // electron energy in target nucleon frame   

    ffset=p.e_spp_ff_set;
	e.weight= jakobian * dsigma_dQ2_dW_( Ee0, -q2,  W, NULL)*(Wmax-Wmin);
	
	double dz=0.01;
	double ctl=cos(p.eel_theta_lab*Pi/180);
	if(fabs( costhetalab(l0,l1) - cos(p.eel_theta_lab*Pi/180)) > dz)
		e.weight*=1e-20;
	else	
		e.weight/=(min(1.0,ctl+dz)-max(-1.0,ctl-dz))/2; // this is the allowed/all ratio in costhetalab
	

	e.weight*=2;    // because channels  1 and 3 for proton
	                // and  2 and 4 for neutron 
	                // contribute to the cross section in 
	                // additive way
	
	Delta.decay(N1,pion);
	
	e.temp.push_back(l1); 
	e.temp.push_back(Delta);
	
	e.out.push_back(l1); // now e.q0(), e.qv() can be used 
	e.out.push_back(N1);
	e.out.push_back(pion);
        
	if( p.pauli_blocking && t.pauli_blocking_old (N1, N0.momentum() ) )  
		e.weight = 0;
		
	return e.weight*cm2;
	
}


vec turn_by_theta(vec v, double theta)
{
	/// Generate random unit vector w 
	/// such that angle betweeen v and w is theta

	vec a=v.dir();
	vec d=rand_dir();
	vec c=d-(d*a)*a; 
	vec p=c.dir();
	return a*cos(theta)+p*sin(theta);
}


double e_spp_event2(params&p, event & e, nucleus &t, bool nc)
{
	//
	//  electron  elastic scattering 
	//
	
	e.weight=0; // event impossible before setting non zero cross section

	particle l0=e.in[0];  // in electron
	particle N0=e.in[1];  // in nucleon
	particle l1=l0;       // out lepton (good for nc)
	particle Delta=N0;    // intermadiate delta 
	particle N1=N0;       // nucleon from delta decay
	particle pion(pdg_pi);// pion from delta decay
	

	double _E_bind=0;	// binding energy				
	
	_E_bind=t.Ef(); // calculate binding energy 
  
	vect p4= l0+N0;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////                                 CHANNELS:                          ///////////////////////////////
	///////     1: e + p -> e + p + \pi^\0                                     ///////////////////////////////
	///////     2: e + n -> e + p + \pi^-                                      ///////////////////////////////
	///////     3: e + p -> e + n + \pi^+                                      ///////////////////////////////
	///////     4: e + n -> e + n + \pi^0                                      ///////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(frandom()<0.5)
	{
		N1.set_proton();
		if(N0.pdg==pdg_proton)
		{
			chan=1;
			pion.set_pi();
		}
		else
		{
			chan=2;
			pion.set_piM();
		}
	}
	else
	{
		N1.set_neutron();
		if(N0.pdg==pdg_proton)
		{
			chan=3;
			pion.set_piP();
		}
		else
		{
			chan=4;
			pion.set_pi();
		}
	}
	//p4.t-=_E_bind;


	double Wmax=sqrt(p4*p4)-l0.mass();
	double Wmin=N1.mass()+pion.mass();
	if(Wmax<=Wmin) 
		return 0;
	double W=Wmin+frandom()*(Wmax-Wmin);
	Delta.set_mass(W);
	//~ cout<<"DElta="<<Delta<<endl;
	double M2=W*W;
	double m2=l0.mass2();
	vec k=turn_by_theta(l0,p.eel_theta_lab*Pi/180);
        k=k.dir();
        double E=p4.t;
        double E2=E*E;
        vec P=vec(p4);	
        double p2=P*P;
        double kp=k*P;
        double kp2=kp*kp;
        
        /// Energy preservation 
        /// E==sqrt(m2+x^2)+sqrt(m2+p*p+x^2-2*x*(p*k))
        
        double delta=E2*E2-2*E2*m2+4*kp2*m2+m2*m2 -2*E2*M2+M2*M2
                     -2*E2*p2-2*m2*p2+2*M2*p2+p2*p2;
        
        
//        cout<<delta<<endl;
        
        if(delta<0) 
		return 0;
	
	double A=0.5/(-E2+kp2);
	double B=(-E2-m2+M2+p2)*kp;
	
	double x[2]={ A*(B-E*sqrt(delta)),
	              A*(B+E*sqrt(delta))};
	
	double X=0;
	for(int i=0;i<2;i++)
		if(x[i]>0)
		{
			X=x[i];
			//~ cout<<x[i]<<endl;
			l1.set_momentum(k*x[i]);
			//~ cout<<"k="<<k<<endl;
			//~ cout<<"W="<<W<<endl;
			//~ cout<<"X="<<X<<endl;
			//~ cout<<"l1="<<l1<<endl;
			Delta.set_momentum(P-l1.p());
			
			vect p5=p4-l1-Delta;
	//		if(fabs(p5.t)<100)
				//~ cout<<"OO "<<p4<<endl
				    //~ <<l1<<endl<<Delta<<endl;
				//~ cout<<"p5="<<p5<<endl;
			//cout<<sqrt(l1*l1)<<endl;
			//~ cout<<thetalab(l0,l1)*180/Pi<<endl;
		}
	
	int ile=(x[0]>0)+(x[1]>0);
	
	if(ile==2 and frandom()<0.5)
		X=x[0];
	double jakobian=1;  

	vect q=l0-l1;   // energy momentum transfer
	
	
	vect nu4 = l0;         // electron fourmonentum 
	nu4.boost( -N0.v() );  // go to nucleon rest frame 
	double Ee0=nu4.t;      // electron energy in target nucleon frame   


	Delta.decay(N1,pion);

	//~ e.weight= jakobian * dsigma_dQ2_dW_( Ee0, -q*q,  W, NULL)*(Wmax-Wmin);
	//~ e.weight*=4*Pi*X*X/Pi*ile;

	vect N1p4=N1;
	N1p4.boost(-nu4.v());
	ffset=p.e_spp_ff_set;
	e.weight= jakobian * dsigma_dQ2_dW_( Ee0, -q*q,  W, NULL)*(Wmax-Wmin);
	e.weight*=4*Pi*X*X/Pi*ile;
	
	//~ double costhetacms=costhetalab(N1p4.v(),Delta.v());
	//~ double phicms=frandom()*2*Pi;
	//~ e.weight=dsigma_dq0_dOmegal_dOmegaCMS(
		//~ l0.E(), 
	        //~ D4V<double>(l1.t,l1.x,l1.y,l1.z), 
		//~ costhetacms, 
		//~ phicms, 
		//~ D4V<double>(N0.t,N0.x,N0.y,N0.z), 
		//~ NULL);
//~ 
	

	e.weight*=2;    // because channels  1 and 3 for proton
	                // and  2 and 4 for neutron 
	                // contribute to the cross section in 
	                // additive way


	
	e.temp.push_back(l1); 
	e.temp.push_back(Delta);
	
	e.out.push_back(l1); // now e.q0(), e.qv() can be used 
	e.out.push_back(N1);
	e.out.push_back(pion);
        
	if( p.pauli_blocking && t.pauli_blocking_old (N1, N0.momentum() ) )  
		e.weight = 0;
		
	return e.weight*cm2;
	
}

double e_spp_event3(params&p, event & e, nucleus &t, bool nc)
{
	//
	//  electron  elastic scattering 
	//
	
	e.weight=0; // event impossible before setting non zero cross section

	particle l0=e.in[0];  // in electron
	particle N0=e.in[1];  // in nucleon
	particle l1=l0;       // out lepton (good for nc)
	particle N1=N0;       // nucleon from delta decay
	particle pion(pdg_pi);// pion from delta decay

	double theta= fabs(p.eel_theta_lab*Pi/180);
	vec k1=turn_by_theta(l0.p(),theta);         // out mom dir
	l1.set_momentum(k1);
	
	vect sum=l0+N0;


	double El1max=sum.t-N1.mass()-pion.mass();
	double El1min=l1.mass();

	l1.set_energy(El1min+frandom()*(El1max-El1min));
	
	vect Delta=sum-l1;    // intermediate Delta 
	double W2=Delta*Delta;
	if(Delta.t<=0 or sqrt(W2)<=N1.mass()+pion.mass())
		return 0;
	if(not ::decay(Delta,N1,pion))
		return 0;

	if(frandom()<0.5)
	{
		N1.set_proton();
		if(N0.pdg==pdg_proton)
		{
			chan=1;
			pion.set_pi();
		}
		else
		{
			chan=2;
			pion.set_piM();
		}
	}
	else
	{
		N1.set_neutron();
		if(N0.pdg==pdg_proton)
		{
			chan=3;
			pion.set_piP();
		}
		else
		{
			chan=4;
			pion.set_pi();
		}
	}
	vec vD=Delta.v();
	vect piD=pion;
	vect l0D=l0;
	piD.boost(-vD);
	l0D.boost(-vD);
	double costhetacms=costhetalab(vD,piD.v());	
	
	vec p1n=vecprod(vec(piD),vD).dir();
	vec l0n=vecprod(vec(l0D),vD).dir();
	double cosphi=p1n*l0n;
	double sinphi=vecprod(p1n,l0n)*vD.dir();

	double phicms=atan2(sinphi,cosphi);
	selfenergy=p.delta_selfenergy;
	ffset=p.e_spp_ff_set;
	
	/*FOR THE DELTA SELF-ENERGY
	# 0 is free target; 
	# 1 is Fermi gas; 
	# 2 is local Fermi gas; 
	# 3 is Bodek-Ritchie;  INACTIVE (free Delta)
	# 4 is "effective" spectral function (carbon or oxygen); INACTIVE (free Delta)
	# 5 is deuterium; INACTIVE (free Delta)
	# 6 is deuterium with constant binding energy nucleus_E_b (for tests only!) INACTIVE (free Delta)
	*/
	
	double KF=0;          ///< local Fermi momentum
	double rho_rel=0;///< local nuclear matter density relative to rho_0
	
	if(selfenergy)
	switch(p.nucleus_target)
	{
		case 1: 
			KF=t.kF();          
			rho_rel=0.75;
			break;
		case 2:
			KF=t.localkf (N0);
			rho_rel=t.density (N0.r.length())/rho_0*fermi3;
			break;
		default: break;
	}
	
	//std::cerr<<KF<<" "<<rho_rel<<"\n\n";
	e.weight=dsigma_dq0_dOmegal_dOmegaCMS_Oset(
			l0.E(),
			e.in[0].mass(),
			D4V<double>(l1.t,l1.x,l1.y,l1.z), 
			costhetacms, 
			phicms, 
			D4V<double>(N0.t,N0.x,N0.y,N0.z),
			KF,
			rho_rel
			); //in cm2 per
	
	e.weight*=4*Pi; 
	e.weight*=4*Pi; 
	e.weight*=El1max-El1min;

	e.weight*=2;    // because channels  1 and 3 for proton
	                // and  2 and 4 for neutron 
	                // contribute to the cross section in 
	                // additive way


	if(e.weight>0)
	{
		e.out.push_back(l1);
		e.out.push_back(N1);
		e.out.push_back(pion);
	}	
	return e.weight*cm2;
}
