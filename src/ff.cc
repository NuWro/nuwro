#include "ff.h"
#include "jednostki.h"
#include "pdg.h"
#include <cmath>
#include <utility>
#include <iostream>
#include "params.h"

using namespace std;

double inline pow2(double x)
	{ return x*x;
	}


static const double M2=pow2(0.5*(PDG::mass_proton+PDG::mass_neutron));

//________________________________________________________
/// Form Factor Parameters 
static double MV2=0.71*GeV2;
static double MA_cc=1030*MeV;
static double MA_nc=1030*MeV;
static double MA_s=1030*MeV;
static double delta_s=0;
static double mu_p=2.793;  // moment dipolowy protonu
static double mu_n=-1.913; // moment dipolowy neutronu
static double piMass2=pow2(PDG::mass_pi);
static double gA=-1.2673;
static int axialFFset=0;
static int strangeFFset=0;
static int strange=0; 
//stange =0 nie strange FF
//stange =1 old implementation (recover old bahaviour) 
//stange =2 new implementation (uses strange axial mass != nc axial mass) 

static int strangeEM=0;

struct FF
{
	double Q2;
	double GEp,GEn;
	double GMp,GMn;
	
	FF(): Q2(0), GEp(0), GEn(0), GMp(0), GMn(0)
	{
	}
	
    inline pair <double,double>  f12(int kind);    
};

//_________________________________________________________
/// Pointer to current vector form factors model
static FF (*FFfromq2)(const double)=0;

// _________________________________________________________
double axialcorr(int axialFF,double q2);

/*
/// strange form factors /??
struct FFs
{
	double F1s;
 	double F2s;
 	double FAs;
};

typedef FFs (*FFspointer)(const double);


/// strange correction

FFs strange_cor(double q2)//,double & f1,double & f2,double & fa) 
{
	double F1s_0 = 0.53;
	double F2s_0 = -0.40;
//?	double delta_s = -0.21;
	double Ma2=MA_s*MA_s;
    double tau=-q2/M2;
    double mian1=pow2(1-q2/MV2);
    FFs f;
	f.F1s= F1s_0/(1+tau)/mian1;
	f.F2s= F2s_0/(1+tau)/mian1;
	f.FAs= delta_s/mian1;
	return f;
}
*/


/// Functions calculating form factors
FF DipoleFF(const double q2);  // 1. dipole electric form factor G_E^V
FF bba03_FF(const double q2);  // 2. hep-ex/0308005 BBA-2003 for Q2<6 GeV
FF bbba05_FF(const double q2); // 3. arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV	
FF JLab_FF(const double q2);   // 4. PHYSICAL REVIEW C, VOLUME 65, 051001(R)
                               // PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
FF kg_FF(const double q2) ;    // 5. K. Graczyk ...
FF npar_FF(const double q2);   // 6. nowa (1990:) parametryzacja JS z qelcc

//double axialcorr(int axialFF,double q2);
////////////////////////////////////////////////////////////////////////
//                          IMPLEMENTATION
////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////
///Calculate F1,F2
pair<double,double> FF::f12(int kind)
{
  double Ge,Gm,f1,f2,F1s=0,F2s=0;
  switch(kind)
  { case 0: //cc
	  Ge=GEp-GEn;
	  Gm=GMp-GMn;
	  break;
	 case 1: // nc proton
	  Ge=0.5*(GEp-GEn)-2*sin2thetaW*GEp;//-0.5*GEs;
	  Gm=0.5*(GMp-GMn)-2*sin2thetaW*GMp;//-0.5*GMs;
	  break;
	 case 2: //nc neutron
	  Ge=0.5*(GEn-GEp)-2*sin2thetaW*GEn;//-0.5*GEs;
	  Gm=0.5*(GMn-GMp)-2*sin2thetaW*GMn;//-0.5*GMs;
	  break; 
	}
	
    const double tau = Q2/(4*M2);
    
	f1= (Ge + tau*Gm)/(1 + tau);
	f2= (Gm - Ge)/(1 + tau) ;
    if(kind and strangeEM)// strangeness in F1, F2 (only for kind!=0 i.e. nc)
    switch(strangeEM)
    {
    case 1:{ 
			  double mian=(1+tau)*pow2(1+Q2/MV2);
			  f1-= 0.5*( 0.53/mian);  
			  f2-= 0.5*(-0.40/mian );
			  break;
		   }
	 case 2: // more complex version from KG
		 {
			double mian=(1+tau)*pow2(1+Q2/MV2);;
			double mus=-0.39;
			double f1s=0.49;
					
			f1-=0.5*mus/mian;
			f2-=0.5*f1s*Q2/GeV2/mian;	 	
			break;
		 }
	  default: break;// no strange correction
    }
	return pair<double,double>(f1,f2);
}


/////////////////////////////////////////////////////////////
FF DipoleFF(const double q2) //dipole electric form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	double a2 = a*a;

	FF ff;
	ff.Q2=-q2;
	ff.GEp= 1.0/a2;
	ff.GEn= 0;
	ff.GMp= mu_p/a2; 
	ff.GMn= mu_n/a2;
	
	return ff;
};

/////////////////////////////////////////////////////////////
FF bba03_FF(const double q2) 
{
	const double Q2=-q2/GeV2;
	const double tau= Q2/(4*M2);

    FF ff;
	ff.Q2=-q2;
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV2
    ff.GEp = 1.0/(1.0 + Q2*(3.253+Q2*(1.422+Q2*(0.08582+Q2*(0.3318+Q2*(-0.09371+Q2*0.01076))))));
	ff.GEn = -mu_n*0.942*tau/(1 + 4.61*tau)/pow2(1 - q2/MV2);
	ff.GMp=  mu_p/(1.0 + Q2*(3.104+Q2*(1.428+Q2*(0.1112+Q2*(-0.006981+Q2*(0.0003705+Q2*-0.7063e-5))))));
	ff.GMn = mu_n/(1.0 + Q2*(3.043+Q2*(0.8548+Q2*(0.6806+Q2*(-0.1287+Q2*0.008912)))));

    return ff;
}

///////////////////////////////////////////////////////////////
FF bbba05_FForig(const double q2) 
{
    const double tau = -q2/(4*M2);

	FF ff;
	ff.Q2=-q2;
	
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV	
	ff.GEp =     ( 1.0-tau*0.0578)/( 1.0 + tau*(11.1+tau*(13.6+tau*33.0)));
	ff.GEn = tau*( 1.25+tau*1.30 )/( 1.0 +tau*(-9.86+tau*(305.0+tau*(-758.0+tau*802.0))));
	ff.GMp = mu_p*( 1.0+tau*0.15 )/( 1.0 + tau*(11.1+tau*(19.6+tau*7.54)));
	ff.GMn = mu_n*( 1.0+tau*1.81 )/( 1.0 + tau*(14.1+tau*(20.7+tau*68.7)));
/*
double Gep=(1-0.0577*tau)/(1+11.2*tau+13.6*tau2+33*tau3);
double Gen=(1.38*tau-0.214*tau2)/(1+8.51*tau+59.9*tau2+13.6*tau3+2.57*tau4);
double Gmp=2.79*(1+0.15*tau)/(1+11.1*tau+19.6*tau2+7.54*tau3);
double Gmn=-1.91*(1+1.82*tau)/(1+14.1*tau+20.7*tau2+69.7*tau3);
	*/
	return ff;
}

FF bbba05_FF(const double q2) 
{
	FF ff;
	ff.Q2=-q2;
	double Q2=-q2;
	double tau=-q2/4.0/M2;
	ff.GEp =            (1.0-tau*0.0578)/( 1.0 + tau*(11.1+tau*(13.6+tau*33.0)));//wersja Artura
	ff.GEn =         tau*(1.25+tau*1.30)/( 1.0 + tau*(-9.86+tau*(305.0+tau*(-758.0+tau*802.0))));//wersja Artura
	ff.GMp = 2.792847351*( 1.0+tau*0.15)/( 1.0 + tau*(11.1+tau*(19.6+tau*7.54)));
	ff.GMn = -1.91304273*( 1.0+tau*1.81)/( 1.0 + tau*(14.1+tau*(20.7+tau*68.7)));
	return ff;
}


/////////////////////////////////////////////////////////////
FF JLab_FF(const double q2) 
{
	const double Q2=-q2/GeV2;
	const double Q =sqrt(Q2); 
	const double tau=-q2/4/M2;
	
	FF ff;
	ff.Q2=-q2;
	//PHYSICAL REVIEW C, VOLUME 65, 051001(R)
	ff.GEp = (1.0 - 0.13*(Q2 - 0.04))/(1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2);	
	ff.GMp = mu_p/(1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2);

	//PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
	ff.GEn =  -1.25*mu_n*tau/(1 + 18.3*tau)/(1 + Q2*1.0e+6/MV2)/(1 + Q2*1.0e+6/MV2);
	ff.GMn = mu_n/(1.0 - (1.74 + 7.63*Q2)*Q + (9.29 + 4.63*Q2)*Q2);

    return ff;	
};

/////////////////////////////////////////////////////////////
FF npar_FF(const double q2) 
{
	const double Q2=-q2/GeV2;
	const double Q=sqrt(Q2); 
	const double tau=-q2/4/M2;
	
	FF ff;
	ff.Q2=-q2;
	
	double mu_p=2.793;  // moment dipolowy protonu
	double mu_n=-1.913; // moment dipolowy neutronu

	ff.GEp = 1/(1+Q*(0.62+Q*(0.68+Q*(2.80+Q*0.83))));
	ff.GMp = mu_p/(1+Q*(0.35+Q*(2.44+Q*(0.50+Q*(1.04+Q*0.34)))));
	ff.GMn = mu_n/(1+Q*(-0.74+Q*(9.29+Q*(-7.63+Q*4.63))));

	const double a=1.25;
	const double b=10.4;
	double GD=1/pow2(1-q2/(0.71*GeV2));
	
	ff.GEn = -a *mu_n*tau*GD/(1+b*tau);
	 
	bool stange=false;
	if (stange)
	  ff.GEp=ff.GMp/mu_p *(1-0.14*(-q2/GeV2-0.3)); 

    return ff;	
};

/*
/////////////////////////////////////////////////////////////
double newFA(const double q2, double MA) //axial form factor
{ 
	//double MA=p.sf_axial_mass;
	double MA2=MA*MA;
	double a = 1.0/(1.0 - q2/MA2);
	return gA*a*a; 
}

/////////////////////////////////////////////////////////////
double newFP(const double q2,double MA) //pseudoscalar form factor
{
	return 2.0*M2*newFA(q2,MA)/(piMass2 - q2); 
}
*/

/////////////////////////////////////////////////////////////
static double funkcja1( double Q2, double *w )
{

double q2 = Q2/GeV/GeV;

return w[4]/(1 + exp( - (q2*w[0] + w[1]) ) )+ w[5]/(1 + exp( - (q2*w[2] + w[3]) ))  + w[6];

}

/////////////////////////////////////////////////////////////
static double funkcja2( double Q2, double *w )
{
double q2 = Q2/GeV/GeV;

return w[6]/(1 + exp( - (q2*w[0] + w[1]) ) ) 
     + w[7]/(1 + exp( - (q2*w[2] + w[3]) )) 
     + w[8]/(1 + exp( - (q2*w[4] + w[5]) )) 
     + w[9];

}

/////////////////////////////////////////////////////////////
FF kg_FF(const double q2)
{ 
	static double tab_gen[7]  = { 10.19704, 2.36812, -1.144266,-4.274101,0.8149924,2.985524,-0.7864434};
	static double tab_gmn[10] = { 3.19646, 2.565681, 6.441526, -2.004055, -0.2972361, 3.606737, -3.135199, 0.299523, 1.261638, 2.64747};
	static double tab_gep[10] = { 3.930227, 0.1108384, -5.325479, -2.846154, -0.2071328, 0.8742101, 0.4283194, 2.568322, 2.577635, -1.185632};
	static double tab_gmp[10] = {-2.862682, -1.560675, 2.321148, 0.1283189, -0.2803566, 2.794296, 1.726774, 0.861083, 0.4184286, -0.1526676};
	static double tab_axial[10] = {-26.10885 ,1.823041, -8.391283,-7.737312, 15.27646, 0.3992788,-1.350184,-0.2021121,-2.870517, 3.879841};

	static const double magneton_N_nnff       = 1;   // tego nie jestem pewien !!!!
	static const double magneton_proton_nnff  = 2.793*magneton_N_nnff;
	static const double magneton_neutron_nnff =-1.913*magneton_N_nnff;

	static double gA_nnff = -1.267;
	static double MA_nnff =  1.015*GeV;     // 1.03*GeV;	


	double Q2=-q2; 
	double y=pow2(1.0 +Q2/0.71/GeV/GeV);

	FF ff;
	ff.Q2=Q2;
	ff.GEp=funkcja2( Q2, tab_gep)/y;
	ff.GEn=funkcja1( Q2, tab_gen);
	ff.GMp=funkcja2( Q2, tab_gmp)*magneton_proton_nnff/y;
	ff.GMn=funkcja2( Q2, tab_gmn)*magneton_neutron_nnff/y;
	return ff;
}


///////////////////////////////////////////////////////////////
// Calculate vector form factors
pair <double,double> f12(double q2, int kind)
{
 FF ff=FFfromq2(q2);
 ff.Q2=-q2;
//  cout<<ff.GEp<< ' '<<ff.GMp<<' '<<ff.GEn<<' '<<ff.GMn<<endl;
 return ff.f12(kind);	
}

///////////////////////////////////////////////////////////////
// Calculate the axial form factors
pair<double,double> fap(double q2,int kind)
{
	double ksi=3.706;

	static const double M12=(PDG::mass_proton+PDG::mass_neutron)/2;
	static double MM=M12*M12;

	double Ga, Fpa, Gas,Fpas;
	double Fa=0,Fp=0;
	switch(kind)
	{case 0: // cc
       Fa = -1.267 / pow2 (1 - q2 / MA_cc / MA_cc);
	   Fa *= axialcorr(axialFFset,q2);
       Fp = 2*MM*Fa/(piMass2-q2);
       break;
	 case 1: //nc proton
       Fa=0.5* -1.267/pow2(1-q2/MA_nc/MA_nc);
       //Fp=2.0*M2*Fa/(piMass2 - q2) ;

      break;
	 case 2: //nc neutron
      Fa=0.5* 1.267/pow2(1-q2/MA_nc/MA_nc);
      //Fp=2.0*M2*Fa/(piMass2 - q2) ;
      break;
    } 
    if(kind and strange)
    {
    	switch(strange)
    	{case 1:
    	  {
    	   Fa-= -0.5* delta_s/pow2(1-q2/MA_nc/MA_nc);
	       }  
    	   break;
    	 case 2: // new implementation
		  {  
		   Fa-= -0.5* delta_s/pow2(1-q2/MA_s/MA_s);
		  }
    	 default:break; //no strange content
	    }
	}
    
    return pair<double,double>(Fa,Fp); 
}

///////////////////////////////////////////////////////////////
double axialcorr(int axialFF,double q2)
{   
	//////////////////////////////
	double min;//maximal reduction
	double max; //maximal enhancement
	double szer;//reduction in Q2
	double dlug;//enhancement range
		
	switch(axialFF)
	{
	case 2:	min=0.9;//maximal reduction
		    max=1.1; //maximal enhancement
		    szer=0.2;//reduction in Q2
		    dlug=2.0;//enhancement range
		    break;

	case 3: min=0.8;//maximal reduction
		    max=1.2; //maximal enhancement
		    szer=0.2;//reduction in Q2
		    dlug=2.0;//enhancement range
		    
	default:return 1;
	}
	
	//parabolic approximation
	double x=-q2/GeV2;
	if (x<szer) return 1.0 + 4.0/szer/szer * (1.0 - min) * x * ( x - szer );
	  else
	if (x<dlug) return 1.0 - 4.0 * (max - 1.0)/pow2(dlug-szer) * ( x - szer ) * ( x - dlug );
      else        
                return 1.0;		
}


///////////////////////////////////////////////////////////////
/// Form Factor Configuration
void ff_configure(params& p)
{
  switch(p.qel_vector_ff_set)	
	 {
		case 1: FFfromq2=DipoleFF; break; 	
		case 2: FFfromq2=bbba05_FF;break;
		case 3: FFfromq2=bba03_FF; break;
		case 4: FFfromq2=JLab_FF;  break;
		case 5: FFfromq2=kg_FF;    break;
		case 6: FFfromq2=npar_FF;  break;
		default: throw("bad ffset");
	 };
  axialFFset=p.qel_axial_ff_set;
  strange=p.qel_strange;
  strangeEM=p.qel_strangeEM;
  delta_s=p.delta_s;
  
  MA_cc=p.qel_cc_axial_mass;
  MA_nc=p.qel_nc_axial_mass;
  MA_s=p.qel_s_axial_mass;
}
//________________________________________________________
