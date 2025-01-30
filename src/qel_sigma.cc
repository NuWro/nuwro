#include "qel_sigma.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "pdg.h"
#include "ff.h"

using namespace std;

inline static double pow2 (double x)          {  return x * x;}

           
            
double qel_amp(double Enu,double q2,
           double m,double M,
           double F1,double F2,double Fp,double Fa);

double qel_amp_orig(double Enu,double q2,
           double m,double M,
           double F1,double F2,double Fp,double Fa);

double amp2(double s,double t,
           double m,double M,double Mp,double Mn,
           double F1,double F2,double FP,double FA);

double amp3(double s,double t,
           double m,double M,double Mp,double Mn,
           double F1,double F2,double FP,double FA);

///
///(semi)elastic neutrino nucleon scattering cross section (Llewelyn-Smith)
///
double qel_sigma ( double Enu, ///< neutrino energy in the target frame
					double q2, ///< 
					int kind,  ///< process type: 0 - cc, 1 - nc proton, 2 - nc neutron
					bool anty, ///< true for antyneutrinos
					double m,  ///< outgoing lepton mass
					double M   ///< nucleon (effective) mass
				  )
{ 
    const double static M12=(PDG::mass_proton+PDG::mass_neutron)/2;
    //M=M12;    
    double F1,F2,Fa,Fp; 

    // Calculate Form Factors 
    list(F1,F2)=f12(q2,kind);
    list(Fa,Fp)=fap(q2,kind);

	if(not anty)
	{
		Fa=-Fa;
		Fp=-Fp;
	}

	double s=pow2(Enu+M)-pow2(Enu);
	
    //double ABC1= qel_amp(Enu,q2,  m,M12,      F1,F2,Fp,Fa);
    double ABC2=    amp2(s  ,q2,  m,M12,M,M,  F1,F2,Fp,Fa);
    //double ABC3=    amp3(s  ,q2,  m,M12,M,M,  F1,F2,Fp,Fa);
    
    //double w1=(G*G*M*M/8/Pi/Enu/Enu)*ABC1 ;
    double w2=(G*G/M/M/32/8/Pi/Enu/Enu)*ABC2;
    //double w3=(G*G/M/M/32/8/Pi/Enu/Enu)*ABC3;
    
	//cout<<w1<<endl;
	//cout<<w2<<endl;
	//cout<<w3<<endl;
	//cout<<w1/w3<<' '<<w2/w3<<endl<<endl;
    
    switch(kind)
	{
		case 0: case 3: case 6: return  w2*cos2thetac; // cc
        default: return  w2;            // nc
	}
}

////////////////////////////////////////////////////////////////////////
double qel_amp(double Enu,double q2,
           double m,double M,
           double F1,double F2,double Fp,double Fa)
{
	double MM=M*M;
	double mm=m*m;
    double x = q2 / MM ;   
    double faa=Fa*Fa;
    double fpp=Fp*Fp;
    double f11=F1*F1;
    double f22=F2*F2;
    
    double A = 
			(+ (4 - x) * faa
			 - (4 + x) * f11
			 - x * f22 * (1 + x / 4)
			 - 4 * x * F1 * F2
			 - mm / MM * ( pow2(F1+F2) + pow2(Fa + 2*Fp) + (x-4)*fpp )
			)
			* (mm - q2)/(4 * MM) 
		   ;

    double B = -x * Fa * (F1 + F2);

    double C = (faa + f11  - x/4 * f22) / 4;

	///  s - u
    double su = 4 * M * Enu + q2 - mm;
    double suM2 = su / MM;
    double ABC = A + suM2*(B  + suM2*C) ;
    return ABC;
}

 
inline double Power(double x, double n)
{ 
	return pow(x,n);
}


////////////////////////////////////////////////////////////////////////
double amp3(double s,double t,
           double m,double M,double Mp,double Mn,
           double F1,double F2,double FP,double FA)
{ 
	FP*=2;// KG has different convention for Fp	
	double AMP =
	(2*(-3*Power(F2,2)*Power(m,4)*Power(Mn,2) + Power(FP,2)*Power(m,4)*Power(Mn,2) -
		   2*Power(F2,2)*Power(m,2)*Power(Mn,4) - 2*Power(F2,2)*Power(m,4)*Mn*Mp -
		   2*Power(FP,2)*Power(m,4)*Mn*Mp + Power(F2,2)*Power(m,4)*Power(Mp,2) +
		   Power(FP,2)*Power(m,4)*Power(Mp,2) + 2*Power(F2,2)*Power(m,2)*Power(Mp,4) +
		   4*Power(F2,2)*Power(m,2)*Power(Mn,2)*s -
		   4*Power(F2,2)*Power(m,2)*Power(Mp,2)*s - Power(F2,2)*Power(m,4)*t -
		   Power(FP,2)*Power(m,4)*t + Power(F2,2)*Power(m,2)*Power(Mn,2)*t -
		   Power(FP,2)*Power(m,2)*Power(Mn,2)*t - 2*Power(F2,2)*Power(Mn,4)*t -
		   2*Power(F2,2)*Power(m,2)*Mn*Mp*t + 2*Power(FP,2)*Power(m,2)*Mn*Mp*t -
		   3*Power(F2,2)*Power(m,2)*Power(Mp,2)*t -
		   Power(FP,2)*Power(m,2)*Power(Mp,2)*t - 2*Power(F2,2)*Power(Mp,4)*t +
		   4*Power(F2,2)*Power(m,2)*s*t + 4*Power(F2,2)*Power(Mn,2)*s*t +
		   4*Power(F2,2)*Power(Mp,2)*s*t - 4*Power(F2,2)*Power(s,2)*t +
		   Power(F2,2)*Power(m,2)*Power(t,2) + Power(FP,2)*Power(m,2)*Power(t,2) +
		   2*Power(F2,2)*Power(Mn,2)*Power(t,2) + 4*Power(F2,2)*Mn*Mp*Power(t,2) +
		   2*Power(F2,2)*Power(Mp,2)*Power(t,2) - 4*Power(F2,2)*s*Power(t,2) +
		   8*Power(FA,2)*Power(M,2)*
				(-2*Power(Mp,2)*s + 2*Power(s,2) +
				  Power(m,2)*(Power(Mn,2) + 2*Mn*Mp + Power(Mp,2) - 2*s - t) +
				  Power(Mn,2)*(2*Power(Mp,2) - 2*s - t) 
				  - 2*Mn*Mp*t - Power(Mp,2)*t +
				  2*s*t + Power(t,2)
				  ) 
			  + 8*Power(F1,2)*Power(M,2)*
						(-2*Power(Mp,2)*s + 2*Power(s,2) +
						  Power(m,2)*(Power(Mn,2) - 2*Mn*Mp + Power(Mp,2) - 2*s - t) +
						  Power(Mn,2)*(2*Power(Mp,2) - 2*s - t) + 2*Mn*Mp*t - Power(Mp,2)*t +
						  2*s*t + Power(t,2)
						 ) 
		  - 8*FA*M*
			(-(F2*(Mn + Mp)*((Power(Mn,2) + Power(Mp,2) - 2*s - t)*t +
				   Power(m,2)*(Power(Mn,2) - Power(Mp,2) + t))) +
			  FP*Power(m,2)*(Power(m,2)*Mn + Power(Mn,3) - Power(Mn,2)*Mp +
			  Mp*s - Mn*(s + t))
			 ) 
		 - 8*F1*M*
			(-2*FA*M*((Power(Mn,2) + Power(Mp,2) - 2*s - t)*t +
				 Power(m,2)*(Power(Mn,2) - Power(Mp,2) + t)) +
			  F2*(Power(m,4)*Mn + (Mn + Mp)*(Power(Mn,2) - 2*Mn*Mp +
					Power(Mp,2) - t)*t 
					+ Power(m,2)*(Mn*(Power(Mp,2) - s) + Mp*(-Power(Mp,2) + s+ t))
				  )
			 )
	  )
	)/ Power(M,2);
   return AMP;
}

////////////////////////////////////////////////////////////////////////
double amp2(double s,double t,
           double m,double M,double Mp,double Mn,
           double F1,double F2,double FP,double FA)
{ 
	FP*=2;// KG has different convention for Fp	

	double F1F1=F1*F1;
	double F2F2=F2*F2;
	double FPFP=FP*FP;
	double FAFA=FA*FA;
	double mm=m*m;
	double MM=M*M;
	double MnMn=Mn*Mn;
	double MpMp=Mp*Mp;
	double ss=s*s;
	double tt=t*t;

	double AMP =
	(2*(-3*F2F2*mm*mm*MnMn + FPFP*mm*mm*MnMn -
		   2*F2F2*mm*MnMn*MnMn - 2*F2F2*mm*mm*Mn*Mp -
		   2*FPFP*mm*mm*Mn*Mp + F2F2*mm*mm*MpMp +
		   FPFP*mm*mm*MpMp + 2*F2F2*mm*MpMp*MpMp +
		   4*F2F2*mm*MnMn*s -
		   4*F2F2*mm*MpMp*s - F2F2*mm*mm*t -
		   FPFP*mm*mm*t + F2F2*mm*MnMn*t -
		   FPFP*mm*MnMn*t - 2*F2F2*MnMn*MnMn*t -
		   2*F2F2*mm*Mn*Mp*t + 2*FPFP*mm*Mn*Mp*t -
		   3*F2F2*mm*MpMp*t -
		   FPFP*mm*MpMp*t - 2*F2F2*MpMp*MpMp*t +
		   4*F2F2*mm*s*t + 4*F2F2*MnMn*s*t +
		   4*F2F2*MpMp*s*t - 4*F2F2*ss*t +
		   F2F2*mm*tt + FPFP*mm*tt +
		   2*F2F2*MnMn*tt + 4*F2F2*Mn*Mp*tt +
		   2*F2F2*MpMp*tt - 4*F2F2*s*tt +
		   8*FAFA*MM*
				(-2*MpMp*s + 2*ss +
				  mm*(MnMn + 2*Mn*Mp + MpMp - 2*s - t) +
				  MnMn*(2*MpMp - 2*s - t) 
				  - 2*Mn*Mp*t - MpMp*t +
				  2*s*t + tt
				  ) 
			  + 8*F1F1*MM*
						(-2*MpMp*s + 2*ss +
						  mm*(MnMn - 2*Mn*Mp + MpMp - 2*s - t) +
						  MnMn*(2*MpMp - 2*s - t) + 2*Mn*Mp*t - MpMp*t +
						  2*s*t + tt
						 ) 
		  - 8*FA*M*
			(-(F2*(Mn + Mp)*((MnMn + MpMp - 2*s - t)*t +
				   mm*(MnMn - MpMp + t))) +
			  FP*mm*(mm*Mn + Mn*MnMn - MnMn*Mp + Mp*s - Mn*(s + t))
			 ) 
		 - 8*F1*M*
			(-2*FA*M*((MnMn + MpMp - 2*s - t)*t +
				 mm*(MnMn - MpMp + t)) +
			  F2*(mm*mm*Mn + (Mn + Mp)*(MnMn - 2*Mn*Mp +
					MpMp - t)*t 
					+ mm*(Mn*(MpMp - s) + Mp*(-MpMp + s+ t))
				  )
			 )
	  )
	)/MM;
   return AMP;
}
