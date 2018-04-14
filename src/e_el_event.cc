#include "e_el_event.h"
#include "params.h"
#include "nucleus.h"
#include "event1.h"
#include "e_el_sigma.h"
#include "kinematics.h"
#include "ff.h"

double LH_(particle k, particle p, particle kprim, vect pprim);


double e_el_event(params&p, event & e, nucleus &t, bool nc)
{
    //
    //  electron  elastic scattering 
    //
    
    e.weight=0; // event impossible before setting non zero cross section

    particle l0=e.in[0];  // in electron
    particle N0=e.in[1];  // in nucleon
    particle l1=l0;       // out lepton (good for nc)
    particle N1=N0;       // out nucleon (good for nc)

    double _E_bind=0;  
			    
    _E_bind=t.Ef(); //calculate binding energy 
		    // some other options are also possible here

    double q2,jakobian;  // kinematics jacobian in q2 for crossection calcualtion

    switch(p.qel_kinematics) // use the same parameter as for qel_kinematics
    { 
	case 0: 
	    q2 = czarek_kinematics2(_E_bind, l0, N0, l1, N1,jakobian); 
	    break;
        case 1: 
	    // preservation of energy momentum for the whole system
            // the nucleon is off shell
            // binding energy used for spectator nucleus mass calculation
            q2 = bodek_kinematics(_E_bind, l0, N0, l1, N1, jakobian);
            break;
        case 2: 
	    // effective mass trick kinematics the most doubtful
            // first lower the nucleon mass down the binging energy
            N0.set_mass(N0.mass() - 27 * MeV);
            // next do the usual kinematics
            // preserving energy and momentum
            // (not taking into account the spactator nucleus)
            q2 = czarek_kinematics(_E_bind // should be 0 here???
                                   , l0, N0, l1, N1, jakobian);

            // use this mass also for dynamics????
            // should it appear in the form factors
            break;
        case 3: 
	    // use momentum dependent potential V(p) 
            // V(p) is on both sides of the equation but with different p
            // N1.set_energy(N1.energy()+V(N1.momentum(),kF_P));
            q2 = momentum_dependent_potential_kinematics(l0, N0, l1, N1, jakobian);
            break;
        case 4: 
	    // use momentum dependent potential V(p)
            // V(p) is on both sides of the equation but with different p
            // then adjust the kinetic energy 
            q2 = momentum_dependent_potential_kinematics(l0, N0, l1, N1, jakobian);
            N1.set_energy(N1.energy() + V(N1.momentum(), t.localkf(N1)));
            break;
        default:
            cerr << "Kinematics code: '" << p.qel_kinematics << "' is invalid." << endl;
            exit(9);

    }
    if (q2 == 0)  // no interaction 
    {
        return 0;
    }

    vect q=l0-l1;   // energy momentum transfer
	

	vect nu4 = l0;         // electron fourmonentum 
	nu4.boost( -N0.v() );  // go to nucleon rest frame 
	double Ee0=nu4.t;     // electron energy in target nucleon frame   
	
//	cout<<l0.t<<' '<<Ee0<<endl;
	int kind= (N0.pdg==pdg_proton ? 10 :11);
	// calculate the cross section	
	double dz=p.eel_dz;
	e.weight= jakobian *e_el_sigma (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; 
	double ctl=cos(p.eel_theta_lab*Pi/180);
	if(fabs( costhetalab(l0,l1) - cos(p.eel_theta_lab*Pi/180)) > dz)
		e.weight*=1e-20;
	else	
		e.weight/=(min(1.0,ctl+dz)-max(-1.0,ctl-dz))/2; // this is the allowed/all ratio in costhetalab
//	cout<<e.weight<<endl;

	if(p.flux_correction)
	{
		// take into account the neutrino flux and nucleon proper time 
		// corrections to the cross section on the whole nucleus
        // According to CJ this is proportional to projectile momentum 
        // divided by target energy
        vec v=N0.v();
        double cor=vec(nu4).length()/l0.momentum()*sqrt(1-v*v);
        e.weight*=cor;
	} 

		
	e.temp.push_back(l1); 
	e.temp.push_back(N1);
	e.out.push_back(l1); // now e.q0(), e.qv() can be used 
	e.out.push_back(N1);
        
	if( p.pauli_blocking && t.pauli_blocking_old (N1, N0.momentum() ) )  
		e.weight = 0;
		
	return e.weight*cm2;
	
}

inline static vec turn_by_theta(vec v, double theta)
{
	/// Generate random unit vector w 
	/// such that angle betweeen v and w is theta

	vec a=v.dir();
	vec d=rand_dir();
	vec c=d-(d*a)*a; 
	vec p=c.dir();
	return a*cos(theta)+p*sin(theta);
}



inline static double pow2 (double x)          
{  
    return x * x;
}

static double alpha = 1.0/137;

/// elastic electron - nucleon scattering cross section 

double e_el_M2_4 ( double Ee, // electron energy in the target frame
                    double q2, // 
                    int kind,  // process type: 10 - ep, 11 - en elastic scattering
                    bool anty, // true for positrons
                    double m,  // lepton mass
                    double M   // nucleon (effective) mass
                  )
{ 
   /*
    * elastic electron - nucleon scattering |Mif|^2 averaged over spin states
    */

    double F1,F2;
    list(F1,F2)=f12(q2,kind);
	    
//	double s = pow2(Ee+M)-pow2(Ee);
	double s = m*m + 2*Ee*M + M*M;   // 
   
    double t = q2;	
    double m2=m*m;
    double m4=m2*m2;
    double M2=M*M;  
    double M4=M2*M2;
    double ABC2  =    
    (4*(4*F1*F2*M2*t*(2*m2 + t) 
    + 2*pow(F1,2)*M2*(2*m4 + 2*M4 + 4*m2*(M2 - s) - 4*M2*s + 2*pow(s,2) + 2*s*t + pow(t,2)) + 
     pow(F2,2)*t*(-m4 - M4 + 2*M2*(s + t) - s*(s + t) + m2*(2*M2 + 2*s + t))))/M2; 
     // nie podzielone przez 4, ale usrednione po spinach

    // this is e^4/q^4*LH 
    double e2=alpha*4*Pi;
    return e2*e2*ABC2/t/t/4; // this is the |Mfi|^2 
}



double e_el_event2(params&p, event & e, nucleus &t, bool nc)
{
	/*
	 * electron  elastic scattering  
	 * this is meant to be a faster implementation
	 * but should give IDENTICAL results
	 */
	
	e.weight=0; // event impossible before setting non zero cross section

	particle l0=e.in[0];  // in electron
	particle N0=e.in[1];  // in nucleon
    
    //l0.set_mass(0);
    //~ N0.set_momentum(vec(0,0,0));
    //N0.set_proton();
    
	particle l1=l0;       // out lepton (good for nc)
	particle N1=N0;       // out nucleon (good for nc)

    
    
	double m2=l0.mass2(); /// == l1.mass();
	double M2=N0.mass2();
  	
    vect p4=l0+N0;
    double s=p4*p4;
    
    double F=4*sqrt(pow2(l0*N0)-m2*M2); // Lorentz invariant flux
    double pstar2=(s-pow2(l0.m()-N0.m()))*(s-pow2(l0.m()+N0.m()))/4/s;

	double _E_bind=0;	//binding energy
				//calculate binding energy 
                
	_E_bind=t.Ef(); 
    
    if(t.A()>4)  //subtract the binding energy from total energy
    {
        p4.t-=_E_bind;
    }

    vec k=turn_by_theta(l0,p.eel_theta_lab*Pi/180);
    k=k.dir();  // direction of the outgoing electron
    double E=p4.t;
    double E2=E*E;
    vec P=vec(p4);	//total momentum
    double p2=P*P;
    double kp=k*P;
    double kp2=kp*kp;
    
    /// Energy preservation 
    /// E==sqrt(m2+x^2)+sqrt(M2+p*p+x^2-2*x*(p*k))
    /// x is l1 momemenum 
    /// x can be calculated from the energy momentum preservation 
    /// one obtains quadratic equation with 
    
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
    int ile=0;
	for(int i=0;i<2;i++)
	if(x[i]>0)
	{   
        ile++;
		X=x[i];
		//~ cout<<x[i]<<endl;
		l1.set_momentum(k*x[i]);
		
		N1.set_momentum(P-l1.p());//forece momentum conservation
		
		vect p5=p4-l1-N1;
		if(fabs(p5.t)>0.001) // test energy-momentum conservation
			cout<<p5<<endl;

//CC	cout<<p5<< "-----"<<vect(N0)+l0-l1-N1<<endl;
		//~ //cout<<sqrt(l1*l1)<<endl;
		//cout<<thetalab(l0,l1)*180/Pi<<endl;
	}
	
	
    double df=X/sqrt(m2 + X*X) 
             -(P*k-X)/sqrt(M2 + P*P - 2*P*k*X + X*X); 
    //~ cout<<df<<endl;    
    //~ 
    //~ df=X/l1.E() - (N1.p()*k)/N1.E();//equivalent expression
    //~ df=X/l1.E() - N1.v()*k;//equivalent expression
    //~ df=(l1.v() - N1.v())*k;//equivalent expression
    //~ cout<<df<<endl<<endl;
                	   
  
	double jakobian=1;// kinematics jacobian in q2 for crossection calcualtion

	vect q=l0-l1;     // energy momentum transfer
    double q2=q*q;
    double tau=-q*q/4/M2;
	

	vect nu4 = l0;         // electron fourmomentum 
	nu4.boost( -N0.v() );  // go to nucleon rest frame 
	double Ee0=nu4.t;     // electron energy in target nucleon frame   
	
    int kind= (N0.pdg==pdg_proton ? 10 :11);
	
	
    if(p.eel_alg=="new") // equivelent to faster
    {
        // as shown in ewro-doc (5)
        double e2=4*Pi/137;
        double tens=e_el_M2_4 (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; 
        double lh=LH_(l0,N0,l1,N1)*e2*e2/q2/q2*4/cm2;
        e.weight=tens;
//        cout<<tens<<endl<<lh<<endl<<tens/lh<<endl<<endl;
        e.weight*=X*X;  // the numerator
        e.weight/=16*F*Pi*Pi*l1.t*df*N1.t; //denominator
        e.weight*=4*Pi; // integration area in dOmega        
        //~ cout<<"new "<< e.weight<<endl;
    }
    if(p.eel_alg=="newdf") // new + de forest
    {
        // as shown in ewro-doc (5)
        double e2=4*Pi/137;
//        double tens=e_el_M2_4 (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; 
        double lh=LH_(l0,N0,l1,N1)*e2*e2/q2/q2*4/cm2;
        e.weight=lh;
//        cout<<tens<<endl<<lh<<endl<<tens/lh<<endl<<endl;
        e.weight*=X*X;  // the numerator
        e.weight/=16*F*Pi*Pi*l1.t*df*N1.t; //denominator
        e.weight*=4*Pi; // integration area in dOmega        
        //~ cout<<"new "<< e.weight<<endl;
    }

    if(p.eel_alg=="new0")
    {
        //M.A.Thomson slides page 37 (bottom of the page) m1=0
        e.weight= e_el_M2_4 (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; //LH
        e.weight/=64*Pi*Pi;
//        e.weight/= pow2(N0.mass() +l0.E() * (1-cos(l0.p(),l1.p())));
        e.weight*= pow2(l1.E()/l0.E()/N0.mass());// equivalent to previous line
        e.weight*=4*Pi;        
        //~ cout<<"new0 "<< e.weight<<' '<<tensN<<endl;
    }

    
    if(p.eel_alg=="new1")
    {
        //M.A.Thomson slides page 38 (bottom of the page) m1>0
        e.weight= e_el_M2_4 (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; //LH
        e.weight/=64*Pi*Pi;
        e.weight/=l0.momentum()*N0.mass();
        e.weight*=l1.momentum2();
        e.weight/=l1.momentum()*(l0.E()+N0.mass()) - l1.E()*l0.momentum()*cos(l0.p(),l1.p());
        e.weight*=4*Pi;        
        //~ cout<<"new1 "<< e.weight<<' '<<tensN<<endl;
    }
    double r1,r2,r3;
    r1=e.weight;
    
    if(p.eel_alg=="faster")  // = 4*Pi*Pi * poprzednie podejścia
    {   // should be equivalent to orig 
        e.weight= e_el_sigma (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; 
        e.weight*=4*Pi*X*X/Pi*ile;
        //~ cout<<"faster "<< e.weight<<endl;
    }
    r2=e.weight; // = 4*Pi*Pi * poprzednie podejścia
	
    if(p.eel_alg=="Rosenbluth")	
    /// Rosenbluth formula 
    {
        double F1,F2;
        list(F1,F2)=f12(q*q,kind);
        
        double ct2=pow2(cos(p.eel_theta_lab/2*Pi/180));
        double st2=pow2(sin(p.eel_theta_lab/2*Pi/180));
        double st4=st2*st2;
        e.weight=alpha*alpha*l1.t
                /(4*l0.t*l0.t*l0.t*st4)
                *(ct2*(F1*F1+tau*F2*F2)+2*tau*st2*pow2(F1+F2))
                /cm2; //nuwro-doc wzorek (28)
        /// exactly from the script 
        e.weight*=4*Pi;// obszar całkowania Monte Carlo   
        //~ cout<<"E="<<l0.t<<"  E'="<<l1.t<<" sin^2(theta/2)="<<st2<<endl;          
//        cout<<"rosenbluth 1 "<< e.weight<<endl;
        
        //~ double GE=F1-tau*F2;
        //~ double GM=F1+F2;
        //~ double tensR=(ct2*(GE*GE+tau*GM*GM)/(1+tau)+2*tau*st2*GM*GM)*st4;
        //~ e.weight=alpha*alpha*l1.t
                //~ /(4*l0.t*l0.t*l0.t*st4)
                //~ *(ct2*(GE*GE+tau*GM*GM)/(1+tau)+2*tau*st2*GM*GM)
                //~ /cm2; //wzorek (M.A. Thompson) //identyczny wynik jak (28)
        //~ e.weight*=4*Pi;// obszar całkowania Monte Carlo             
        //~ ///end Rosenbluth          
        //~ cout<<"rosenbluth 2 "<< e.weight<<' '<<tensR<<endl;
        //~ 
        //~ e.weight=dsigma_dOmega_ksiazka (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2;
        //~ e.weight*=4*Pi;
//~ 
        //~ cout<<"rosenbluth K "<< e.weight<<endl;
//~ 
        //~ e.weight=e_el_sigma(Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass(),4)/cm2;
        //~ e.weight*=4*Pi;
//~ 
        //~ cout<<"KG as is "<< e.weight<<endl;
    } 
    r3=e.weight;   
    //~ cout<<"r12="<<r1/r2<<"   r2/1="<<r2/r1<<endl;
    //~ cout<<"r13="<<r1/r3<<"   r3/1="<<r3/r1<<endl;
    //~ cout<<endl;



	if(p.flux_correction)
	{
		// take into account the neutrino flux and nucleon proper time 
		// corrections to the cross section on the whole nucleus

		double cor1= 1 - N0.v() * l0.v().dir(); // this is good for neutrino beam

        // According to CJ this is proportional to projectile momentum 
        // divided by target energy
        
//        double cor=vec(nu4).length()/vec(l0).length()*sqrt(1-N0.v().norm2()); // this good is for any projectile mass
        double cor=vec(nu4).length()/vec(l0).length()*N0.mass()/N0.E(); // this good is for any projectile mass
        //~ cout<<setprecision(10);
        //~ cout<<costhetalab(l0,N0)<<endl;
        //~ cout<< cor<<endl;
        //~ cout<<cor1<<endl<<endl;   // co1,cor have common 5-6 digits for electrons
        e.weight*=cor; // we use the more accurate result
	} 

        
	e.temp.push_back(l1); 
	e.temp.push_back(N1);
	e.out.push_back(l1); // now e.q0(), e.qv() can be used 
	e.out.push_back(N1);
        
	if( p.pauli_blocking && t.pauli_blocking_old (N1, N0.momentum() ) )  
		e.weight = 0;
//	exit(0);	
	return e.weight*cm2;
	
}

double e_el_event2orig(params&p, event & e, nucleus &t, bool nc)
{
	/*
	 * electron  elastic scattering  
	 * this is meant to be a faster implementation
	 *  but should give IDENTICAL results
	 */

	e.weight=0; // event impossible before setting non zero cross section

	particle l0=e.in[0];  // in electron
	particle N0=e.in[1];  // in nucleon
	particle l1=l0;       // out lepton (good for nc)
	particle N1=N0;       // out nucleon (good for nc)

	double m2=l0.mass2(); /// == l1.mass();
	double M2=N0.mass2();
	double _E_bind=0;	//binding energy
				//calculate binding energy 

	if(t.A()>4)
	    _E_bind=t.Ef(); 

	vect p4=l0+N0;
	p4.t-=_E_bind;

	vec k=turn_by_theta(l0,p.eel_theta_lab*Pi/180);
        k=k.dir();
        double E=p4.t;
        double E2=E*E;
        vec P=vec(p4);	
        double p2=P*P;
        double kp=k*P;
        double kp2=kp*kp;

        /// Energy preservation 
        /// E==sqrt(m2+x^2)+sqrt(M2+p*p+x^2-2*x*(p*k))
        /// x is l1 momemenum

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
		
		N1.set_momentum(P-l1.p());
		
		vect p5=p4-l1-N1;
//		if(fabs(p5.t)<100)
			//~ cout<<p4<<endl
			    //~ <<l1+N1<<endl;
			//~ cout<<p5<<endl;
		//cout<<sqrt(l1*l1)<<endl;
		//~ cout<<thetalab(l0,l1)*180/Pi<<endl;
	}
	
	int ile=(x[0]>0)+(x[1]>0);
	
	if(ile==2 and frandom()<0.5)
		X=x[0];

	vect q=l0-l1;   // energy momentum transfer
	

	vect nu4 = l0;         // electron fourmonentum 
	nu4.boost( -N0.v() );  // go to nucleon rest frame 
	double Ee0=nu4.t;     // electron energy in target nucleon frame   
	
	int kind= (N0.pdg==pdg_proton ? 10 :11);
	
	
	// calculate the cross section	
	double dz=0.01;
	e.weight= e_el_sigma (Ee0, q*q, kind, l0.pdg<0, l1.mass(), N0.mass())/cm2; 
	e.weight*=4*Pi*X*X/Pi*ile;

	if(p.flux_correction)
	{
		// take into account the neutrino flux and nucleon proper time 
		// corrections to the cross section on the whole nucleus
		// int qel_relat=0;
		double vv,vz;
		vv = N0.v().length ();  // why not used ??
		vz = N0.v() * l0.v().dir(); // this is general
		e.weight *= sqrt ( (1 - vz) * (1 - vz) );
	} 
	
		
	e.temp.push_back(l1); 
	e.temp.push_back(N1);
	e.out.push_back(l1); // now e.q0(), e.qv() can be used 
	e.out.push_back(N1);
        
	if( p.pauli_blocking && t.pauli_blocking_old (N1, N0.momentum() ) )  
		e.weight = 0;
		
	return e.weight*cm2;
	
}


double LH_(particle k, particle p, particle kprim, vect pprim)
{
    vect q=k-kprim;
    vect q_=pprim-p;
    const double q_2=q_*q_;
    const double pk=p*k;
    const double pK=p*kprim;
    const double pq_=p*q_;
    const double kq_=k*q_;
    const double Kq_=kprim*q_;
    const double kK=k*kprim;

    double M2=pprim*pprim;
    
	const double qM=q_2/M2;
	const double tau=-qM/4.0;
	
    double a1= 2.0*M2*kK;
    double a2= 2.0*pk*pK-M2*kK+pk*Kq_+pK*kq_-kK*pq_;
    double a3= pK*kq_-pk*Kq_;
    double a4= kK*q_2-2.0*kq_*Kq_;

    bool m_qel_new=true;
	
    if (k.pdg%2==0) // neutrino scattering 
	{
		double f1,f2,fa,fp;
        int kind=0;//neutrino cc
        kind=1; //neutrino proton nc
        kind=2; //neutrino neutron nc
        kind=(k.pdg!=kprim.pdg ? 0   // cc
                               : (p.pdg==pdg_proton? 1 :2)
             );
        list(f1,f2)=f12(q_2,kind);
        list(fa,fp)=fap(q_2,kind);

		const double f11= f1*f1;
		const double f22= f2*f2;
		const double faa= fa*fa;
		const double fpp= fp*fp;

		const double ff= f1 + f2;

		const double h1= faa*(1.0+tau) + ff*ff*tau ;
		const double h2= faa + f11 + f22*tau ;
		const double h3= (k.pdg>0 ? 2.0*fa*ff :-2.0*fa*ff);
		const double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 + fa*fp - fpp*tau;

		const double lh= 2.0*(a1*h1 + a2*h2 + a3*h3 + a4 * h4) ;

		return  lh;
	
	}
	else // electron scattering
	{
        
        double f1,f2;
        list(f1,f2)=f12(q_2,(p.pdg==pdg_proton?10:11));

        double f11= f1*f1 ;
        double f22= f2*f2 ;

        double ff= f1 + f2 ;

        double h1= ff*ff*tau ;
        double h2= f11 + f22*tau ;
        double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 ;

        double lh= 2.0*(a1*h1 +  a2*h2 + a4*h4 );
        
        return lh;
	}
}


