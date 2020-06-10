#include <fstream>
#include <sstream>
#include <cassert>
#include <limits>
#include <cmath>
#include "jednostki.h"
#include "nucleus.h"
#include "event1.h"
using namespace NUWRO;

//      Coherent single pion production
//      Implementation of the Rein-Sehgal model
//      Nucl. Phys. B223 (1983) 29-44 
//      implemented by C. Juszczak 2011

static const double maxdouble=numeric_limits<double>::max();

static double const incl[]={0,0,
	0.2,25  ,   0.3,56.5,  	0.4, 17.0, 	0.6, 8.0,   	1.0,25.5,  
	1.3,20.5, 	1.5,26.0,   1.8,23.0,  	2.0, 24.0,   	5.0,22.0,
	maxdouble,22.0};

static double const tot[]={0,0,
	0.25,109,  0.6,26.0,  0.7,29.0, 0.8,27.0, 1.0,38.0,
	1.2,31.5,  1.5, 36.0, 2.0,30.5, 5.0,26.0, maxdouble,26.0};

double sigma(const double * t, double E)
{   E/=GeV;
	while(E>t[2]) t+=2;
	return (t[1]+(E-t[0])/(t[2]-t[0])*(t[3]-t[1]))*millibarn;
}
//Berger-Sehgal tables
static double const A_BS[]={0,0,
	0.076,11600.0,  0.080, 14700.0,  0.100, 18300.0, 0.148, 21300.0, 0.162, 22400.0,
	0.226,16400.0,  0.486, 5730.0, 0.584, 4610.0, 0.662,4570.0, 0.776, 4930.0, 0.870, 5140.0, maxdouble,5140.0};
	
static double const B_BS[]={0,116.0,
	0.076,116.0,  0.080, 109.0,  0.100, 89.8, 0.148, 91.0, 0.162, 89.2,
	0.226,80.8,  0.486, 54.6, 0.584, 55.2, 0.662,58.4, 0.776, 60.5, 0.870, 62.2, maxdouble,62.2};
//like sigma, no milbarn
double sigma1(const double * t, double E)
{   E/=GeV;
	while(E>t[2]) t+=2;
	return (t[1]+(E-t[0])/(t[2]-t[0])*(t[3]-t[1]));
}


inline double pow2(double x){return x*x;}

////////////////////////////////////////////////////////////////////////
void
cohevent_cj (params & par, event & e, nucleus & jadro, bool cc, bool fast)
////////////////////////////////////////////////////////////////////////
{
  e.weight=0;
  particle nu = e.in[0];
  particle lepton(nu);
  particle pion;
 
	if(jadro.A()<4)
		return; // no coherent scattering if A<4
  if (not cc)  pion.pdg = 111;//pi0
  else if (nu.pdg > 0)
    { //CC pi+ production by neutrino
      lepton.pdg --;
      pion.pdg = 211;//pi+
    }
  else
    { //CC pi- production by antyneutrino
      lepton.pdg ++;
      pion.pdg = -211;//pi-
    }
    double mlep=PDG::mass (lepton.pdg);
	lepton.set_mass (mlep);
    double mpi=PDG::mass (pion.pdg);
	mpi=138*MeV;
	pion.set_mass (mpi);
    
    int A = jadro.A();
	double E=nu.t;
	if(E<mpi) return;
	static const double M=0.5*(PDG::mass_proton+PDG::mass_neutron);
	double mpi2=mpi*mpi;
	const double fpi=0.93*mpi;
//	static const double fpi=130*MeV;;
	static const double r2=0.2;
	static const double ma2=1*GeV2;

    static const double r0=1.*fermi;
    static const double C0=r0*r0/3.;
	double b=C0*pow(A,2./3.);

    double z1=frandom();
	double t2=-log (1-z1)/ b;
	double t=sqrt(t2);
	double w1=1/b;
	
    double z2=frandom();
    double q0=mpi+(E-mpi)*z2*z2;
    double y=q0/E;
    double w2=(E-mpi)/E*2*z2;

    double p1=sqrt(pow2(E-q0)-lepton.mass2());
    double p2=sqrt(q0*q0-mpi2);
    double k=nu.momentum();
    //k=E;
    
    double z3=frandom();
    double qmin=max(abs(p1-k),abs(p2-t));
    double qmax=min(p1+k,p2+t);
    if(qmax<=qmin)return;
    
    double q=qmin+(qmax-qmin)*z3;
    double Q2=-q0*q0+q*q;
    double x=Q2/2/M/q0;
    double w3=1./2/M/q0*2*q*(qmax-qmin);
    
	double Fabs=exp(-9*pow(A,1./3.)/16/Pi/r0/r0*sigma(incl,q0));

	static const double C1=G*G*M/2/Pi/Pi/16/Pi*(1+r2);             //constants from formula (10)
	double value=C1*A*A*fpi*fpi*(E-q0)*pow2(sigma(tot,q0))*pow2(ma2/(ma2+Q2))*Fabs; //(10)

	double jacob=w1*w2*w3;
	e.weight=value/A*jacob/cm2;
	double Q2min=mlep*mlep*y*(1-y);
    if(cc)
    {//   double Q2min=mlep*y*(1-y);
		double Corr = pow2(1.0 - 0.5 * Q2min / (Q2 + mpi2))
			      + y/4 * Q2min * (Q2 - Q2min) / pow2(Q2 + mpi2);
		e.weight*=2*Corr;
    }
//    cout<<"new "<<mlep<<' '<<E<<' '<<Q2<<' '<<Q2min<<endl;
    if ((par.coh_mass_correction) && 
       ((y > 1.0 - mlep / E) || (Q2 < Q2min)))
        e.weight = 0.0;
    vec lep=rand3(vec(nu),p1,q);
    vec pio=rand3(vec(nu)-lep,p2,t);
    
    lepton.set_momentum(lep);
    pion.set_momentum(pio);
    
	if(e.weight>0)
	{  //e.weight/=1.5;//the last multiplicative factor is taken from Nuance MC
	   e.out.push_back (lepton);
	   e.out.push_back (pion);
	 }
}

////////////////////////////////////////////////////////////////////////
void
cohevent_cj_kin (params & par, event & e, nucleus & jadro, bool cc, bool fast)
////////////////////////////////////////////////////////////////////////
{
  e.weight=0;
  particle nu = e.in[0];
  particle lepton(nu);
  particle pion;
 
	if(jadro.A()<4)
		return; // no coherent scattering if A<4
  if (not cc)  pion.pdg = 111;//pi0
  else if (nu.pdg > 0)
    { //CC pi+ production by neutrino
      lepton.pdg --;
      pion.pdg = 211;//pi+
    }
  else
    { //CC pi- production by antyneutrino
      lepton.pdg ++;
      pion.pdg = -211;//pi-
    }
    double mlep=PDG::mass (lepton.pdg);
	lepton.set_mass (mlep);
    double mpi=PDG::mass (pion.pdg);
	mpi=138*MeV;
	pion.set_mass (mpi);
    
    int A = jadro.A();
	double E=nu.t;
	if(E<mpi) return;
	static const double M=0.5*(PDG::mass_proton+PDG::mass_neutron);
	double mpi2=mpi*mpi;
	const double fpi=0.93*mpi;
//	static const double fpi=130*MeV;;
	static const double r2=0.2;
	static const double ma2=1.00*GeV2;

    static const double r0=1.00*fermi;
    static const double C0=r0*r0/3.;
	double b=C0*pow(A,2./3.);

    double z1=frandom();
	double t2=-log (1-z1)/ b;
	double t=sqrt(t2);
	double w1=1/b;
	
    double z2=frandom();
    double q0=mpi+(E-mpi)*z2*z2;
    double y=q0/E;
    double w2=(E-mpi)/E*2*z2;

    double p1=sqrt(pow2(E-q0)-lepton.mass2());
    double p2=sqrt(q0*q0-mpi2);
    double k=nu.momentum();
    //k=E;
    
    double z3=frandom();
    double qmin=max(abs(p1-k),abs(p2-t));
    double qmax=min(p1+k,p2+t);
    if(qmax<=qmin)return;
    
    double q=qmin+(qmax-qmin)*z3;
    double Q2=-q0*q0+q*q;
    double x=Q2/2/M/q0;
    double w3=1./2/M/q0*2*q*(qmax-qmin);
    
	double Fabs=exp(-9*pow(A,1./3.)/16/Pi/r0/r0*sigma(incl,q0));

	static const double C1=G*G*M/2/Pi/Pi/16/Pi*(1+r2);             //constants from formula (10)
	double value=C1*A*A*fpi*fpi*q0*((2*E-q0)*(2*E-q0)-q*q)/4/E/q*pow2(sigma(tot,q0))*pow2(ma2/(ma2+Q2))*Fabs; //(10)

	double jacob=w1*w2*w3;
	e.weight=value/A*jacob/cm2;
	double Q2min=mlep*mlep*y*(1-y);
    if(cc)
    {//   double Q2min=mlep*y*(1-y);
		double Corr = pow2(1.0 - 0.5 * Q2min / (Q2 + mpi2))
			      + y/4 * Q2min * (Q2 - Q2min) / pow2(Q2 + mpi2);
		e.weight*=2*Corr;
    }
//    cout<<"new "<<mlep<<' '<<E<<' '<<Q2<<' '<<Q2min<<endl;
    if ((par.coh_mass_correction) && 
       ((y > 1.0 - mlep / E) || (Q2 < Q2min)))
        e.weight = 0.0;
    vec lep=rand3(vec(nu),p1,q);
    vec pio=rand3(vec(nu)-lep,p2,t);
    
    lepton.set_momentum(lep);
    pion.set_momentum(pio);
    
	if(e.weight>0)
	{  //e.weight/=1.5;//the last multiplicative factor is taken from Nuance MC
	   e.out.push_back (lepton);
	   e.out.push_back (pion);
	 }
}
//Berger-Sehgal
////////////////////////////////////////////////////////////////////////
void
cohevent_bs (params & par, event & e, nucleus & jadro, bool cc, bool fast)
////////////////////////////////////////////////////////////////////////
{
  e.weight=0;
  particle nu = e.in[0];
  particle lepton(nu);
  particle pion;
 
	if(jadro.A()<4)
		return; // no coherent scattering if A<4
  if (not cc)  pion.pdg = 111;//pi0
  else if (nu.pdg > 0)
    { //CC pi+ production by neutrino
      lepton.pdg --;
      pion.pdg = 211;//pi+
    }
  else
    { //CC pi- production by antyneutrino
      lepton.pdg ++;
      pion.pdg = -211;//pi-
    }
    double mlep=PDG::mass (lepton.pdg);
	lepton.set_mass (mlep);
    double mpi=PDG::mass (pion.pdg);
	mpi=138*MeV;
	pion.set_mass (mpi);
    
    int A = jadro.A();
	double E=nu.t;
	if(E<mpi) return;
	static const double M=0.5*(PDG::mass_proton+PDG::mass_neutron);
	double mpi2=mpi*mpi;
	const double fpi=0.93*mpi;

	static const double r2=0.2;
	static const double ma2=1.00*GeV2;

    static const double r0=1.00*fermi;
    static const double C0=r0*r0/3.;
	double b=C0*pow(12.0,2./3.);

    double z1=frandom();
	double t2=-log (1-z1)/ b;
	double t=sqrt(t2);
	double w1=1/b;
	
    double z2=frandom();
    double q0=mpi+(E-mpi)*z2*z2;
    double y=q0/E;
    double w2=(E-mpi)/E*2*z2;

    double p1=sqrt(pow2(E-q0)-lepton.mass2());
    double p2=sqrt(q0*q0-mpi2);
    double k=nu.momentum();
    
    double z3=frandom();
    double qmin=max(abs(p1-k),abs(p2-t));
    double qmax=min(p1+k,p2+t);
    if(qmax<=qmin)return;
    
    double q=qmin+(qmax-qmin)*z3;
    double Q2=-q0*q0+q*q;
    double x=Q2/2/M/q0;
    double w3=1./2/M/q0*2*q*(qmax-qmin);
    
	

	
    vec lep=rand3(vec(nu),p1,q);
    vec pio=rand3(vec(nu)-lep,p2,t);
    
    lepton.set_momentum(lep);
    pion.set_momentum(pio);
    
    double tpi=sqrt(pio*pio+mpi2)-mpi;
    
    double a_bs=sigma1(A_BS,tpi);
    double b_bs=sigma1(B_BS,tpi);

	double value=M*q0*G*G/(8*Pi*Pi*E)*fpi*fpi*((2*E-q0)*(2*E-q0)-q*q)/q*pow2(ma2/(ma2+Q2))*a_bs*pow(1-z1,b_bs/b/GeV2-1.0);

	
	//For some reason the following is less efficient in spite of using b_bs to sample t insted of global b....
    /*double z2=frandom();
    double q0=mpi+(E-mpi)*z2*z2;
    double y=q0/E;
    double w2=(E-mpi)/E*2*z2;
    
    double tpi=q0-mpi;
    
    double a_bs=sigma1(A_BS,tpi);
    double b_bs=sigma1(B_BS,tpi)/GeV2;
    
    double z1=frandom();
	double t2=-log (1-z1)/ b_bs;
	double t=sqrt(t2);
	double w1=1/b_bs;

    double p1=sqrt(pow2(E-q0)-lepton.mass2());
    double p2=sqrt(q0*q0-mpi2);
    double k=nu.momentum();
    
    double z3=frandom();
    double qmin=max(abs(p1-k),abs(p2-t));
    double qmax=min(p1+k,p2+t);
    if(qmax<=qmin)return;
    
    double q=qmin+(qmax-qmin)*z3;
    double Q2=-q0*q0+q*q;
    double x=Q2/2/M/q0;
    double w3=1./2/M/q0*2*q*(qmax-qmin);
    
	

	
    vec lep=rand3(vec(nu),p1,q);
    vec pio=rand3(vec(nu)-lep,p2,t);
    
    lepton.set_momentum(lep);
    pion.set_momentum(pio);
    
    

	double value=M*q0*G*G/(8*Pi*Pi*E)*fpi*fpi*((2*E-q0)*(2*E-q0)-q*q)/q*pow2(ma2/(ma2+Q2))*a_bs;*/

	double jacob=w1*w2*w3;

	e.weight=value*jacob/12/cm2*millibarn/GeV2;
	double Q2min=mlep*mlep*y*(1-y);                       

//According to original Berger-Sehgal the second set holds
  if(cc)
  {
		//double Corr = pow2(1.0 - 0.5 * Q2min / (Q2 + mpi2))+ y/4 * Q2min * (Q2 - Q2min) / pow2(Q2 + mpi2);
		//e.weight*=2*Corr*cos2thetac;
		double Corr = pow2(ma2/(ma2+Q2) - 0.5 * Q2min / (Q2 + mpi2))+ y/4 * Q2min * (Q2 - Q2min) / pow2(Q2 + mpi2);
		e.weight*=2*Corr/pow2(ma2/(ma2+Q2))*cos2thetac;
  }

  if ((par.coh_mass_correction) && ((y > 1.0 - mlep / E) || (Q2 < Q2min)))
      e.weight = 0.0;
	
	//Re-weight for A!=12 according to R-S model
	if(A!=12)
	{
		e.weight*=A/12.0*pow((1-z1),1.0-pow(12.0/A,2./3.))*exp(-9*(pow(A,1./3.)-pow(12.0,1./3.))/16/Pi/r0/r0*sigma(incl,q0));
	}
	
	
	if(e.weight>0)
	{  e.out.push_back (lepton);
	   e.out.push_back (pion);
	 }
 }
