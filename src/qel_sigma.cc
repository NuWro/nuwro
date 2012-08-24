#include "qel_sigma.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "pdg.h"
//#include "qel_amp.cc"
#include "ff.h"

using namespace std;

inline static double pow2 (double x)          {  return x * x;}

/// (semi)elastic neutrino nucleon scattering cross section (Llewelyn-Smith)
double qel_sigma ( double Enu, ///< neutrino energy in the target frame
					double q2, ///< 
					int kind,  ///< process type: 0 - cc, 1 - nc proton, 2 - nc neutron
					bool anty, ///< true for antyneutrinos
					double m,  ///< outgoing lepton mass
					double M   ///< nucleon mass
				  )
{ 
	double mm=m*m;
	double MM=M*M;
    double F1,F2,Fa,Fp; 
    
    list(F1,F2)=f12(q2,kind);// Calculate Form Factors 
    list(Fa,Fp)=fap(q2,kind);
    
//cout<<F1<< ' '<<F2<<' '<<Fa<<' '<<Fp<<" new"<<endl;
    
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
			 - mm / MM * ( pow2(F1 + F2) + pow2(Fa + 2*Fp) + (x-4)*fpp )
			)
			* (mm - q2)/(4 * MM) 
		   ;

    double B = -x * Fa * (F1 + F2);

    double C = (faa + f11  - x/4 * f22) / 4;

///  s - u
    double su = 4 * M * Enu + q2 - mm;
    double suM2 = su / MM;
    double ABC = A + suM2*( (anty? B : -B)   + suM2*C );
//    ABC=amp2(s,q2,m,M,M,M,F1,F2,Fa,Fp);
    if(kind) // nc proton or nc neutron
	  return  (G*G*MM/8/Pi/Enu/Enu)*ABC;
    else
      return  (G*G*MM/8/Pi/Enu/Enu)*ABC *cos2thetac;
 }
