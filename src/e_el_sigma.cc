#include "e_el_sigma.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "pdg.h"
#include "ff.h"

using namespace std;

inline static double pow2 (double x)          
{  
    return x * x;
}

static double alpha = 1.0/137;

/// elastic electron - nucleon scattering cross section 

double e_el_sigma ( double Ee, // electron energy in the target frame
                    double q2, // 
                    int kind,  // process type: 10 - ep, 11 - en elastic scattering
                    bool anty, // true for positrons
                    double m,  // lepton mass
                    double M,   // nucleon (effective) mass
                    int what  // what to return
                  )
{ 
   /*
    * elastic electron - nucleon scattering cross section 
    */

    double F1,F2;
    list(F1,F2)=f12(q2,kind);

	if(not anty)
    {
        /// In the first order there is no difference between electron-nucleon and electron-proton scattering
    }
	    
//	double s = pow2(Ee+M)-pow2(Ee);
	double s = m*m + 2*Ee*M + M*M;   // 
   
    double t = q2;	
    double m2=m*m;
    double m4=m2*m2;
    double M2=M*M;  
    double ABC2  =    
    4*(
      4*F1*F2*(2*m2 + t)*t
    + 2*F1*F1*(2*pow2(m2+M2-s) + (2*s+t)*t) 
    +   F2*F2*(-pow2(m2-M2) + (2*M2-s)*(s + t) +m2*(2*s + t))*t/M2
     ); 
     // nie podzielone przez 4, ale usrednione po spinach
    
//    cout<<"ABC3="<<ABC2<<endl;
    double Eprim = Ee + t/2/M;  /// energy of outgoing electron

    double w2= (Pi*alpha*alpha/16/Ee/Ee/M/M) * (ABC2/t/t);    
    
    /// Powyższe różni się o czynnik (16*Pi*Pi*Ee)/Eprim od cross section dsigma/dQ2
    /// A powinien różnić się o 2*Eprim*Eprim;
    /// 4*Pi to powierzchnia sfery 
    /// pozostaje 4*Pi*Ee/Eprim
//    w2*=Eprim*Eprim/Pi;           /// uncomment to obtain dsigma/dTheta_Lab
//    w2*=Eprim*Eprim;              /// uncomment to obtain dsigma/dTheta_Lab
//    cout<<Ee<<endl;
    double w2dq2= alpha*alpha*( Pi/Ee/Eprim )    *( Eprim*Eprim/64/Pi/Pi/Ee/Ee/M/M ) * (ABC2/4/t/t);  /// cross section dsigma/dQ2
    double w2q2=  alpha*alpha/64/Pi * (Eprim/Ee/Ee/Ee/M/M)  * (ABC2/4/t/t);  /// cross section dsigma/dQ2 simplified version but

    /* 
     * Na razie wstawilem Jakobian przejscia dla w przypadku gdy zaniedba sie mase leptonu (KG)
     */
    // from the Handout 
    // dsigma_dq2 = 1/(64*Pi*s*pstar2) * Mfi2 
    // and s*pstar2= (Ee*M)^2 in lab frame and m=0
    // so Mfi may be identified as   alpha*alpha * (Eprim/Ee)  * (ABC2/4/t/t);
    
    double dsigma_dQ2=  alpha*alpha/64/Pi * (Eprim/Ee/Ee/Ee/M/M)  * (ABC2/4/t/t);
    double dsigma_dOmega=2*Eprim*Eprim* dsigma_dQ2;
    
//    cout<< w2q2<<' '<<dsigma_dQ2<<' '<<w2dq2<<endl;
    //return dsigma_dOmega;
	switch(what)
    {
        case 0: return w2;  
        case 1: return dsigma_dQ2;
        case 2: return w2q2;
        case 3: return w2dq2;
        case 4: return dsigma_dOmega;
        default: return w2;
    }
}


double dsigma_dOmega_ksiazka(double En, // 
                             double q2, // 	
                             int kind,  // process type: 10 - ep, 11 - en elastic scattering
                             bool anty, // true for positrons
                             double m,  // lepton mass
                             double M   // nucleon (effective) mass
                            )
{ 
    double t = q2;

    double F1,F2;
    
    list(F1,F2)=f12(q2,kind); // calculate form factors

    if(not anty)
    {
        /// In the first order there is no difference between electron-nucleon and electron-proton scattering
    }

    double par[1]={0.84};
    double s = 2.*En*M + M*M;
    double Eprim = En + t/2./M;
    double sin2theta2 = -t/4./En/Eprim;
    double cos2theta2 =  1. - sin2theta2; 
    double Mott = alpha*alpha*(Eprim/En)*cos2theta2
                  /pow2(2.*En*sin2theta2);
    
//    cout<<" mottks="<<Mott<<endl;
    //~ cout<<"E="<<En<<"  E'="<<Eprim<<" sinth2="<<sin2theta2<<endl;          

    //cerr << ope << "\t" << Eprim << "\t" << t << endl;

    double prze = F1*F1  - (t/4./M/M)*( F2*F2 + 2.0*pow2(F1 + F2)*sin2theta2/cos2theta2 );
    return prze*Mott;

};
