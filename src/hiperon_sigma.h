#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include "jednostki.h"
#include "pdg.h"

using namespace std;

inline double hiperon_sigma(double E,double Q2,double m, int N1pdg, int ffset=0)
{
    /*
     * hiperon differential cross section ds/dq2
     *
     * ffset=0 //
     * ffset=1 // BBBA05
     *
     */
    double m2=m*m;
    double g_a=1.267;
    double pi=M_PI;

    double MV=0.84*GeV;
    double Mv2=0.71*GeV2;

    double Ma=1.03*GeV;
    double Ma2=Ma*Ma;


    double M_out=0;
    double M_in=0;

    switch(N1pdg)
    {
        case PDG::pdg_Lambda:
            M_out =PDG::mass_Lambda;
            M_in = PDG::mass_proton;
            break;
        case PDG::pdg_Sigma:
            M_out =PDG::mass_Sigma;
            M_in = PDG::mass_proton;
            break;
        case PDG::pdg_SigmaM:
            M_out =PDG::mass_SigmaM;
            M_in = PDG::mass_neutron;
            break;
        default:
            cerr<<"wrong hiperon pdg code "<<N1pdg<<". (allowed values: 1,2,3)."<<endl;
            exit(13);
    }
//    cout<<M_out<<" "<<M_in<<" "<<m<<endl;
    double Ge=1.0/pow(1+Q2/Mv2,2.0);
    double Gm=4.71*Ge;
    double Gep=0,Gmp=0,Gen=0,Gmn=0,G_D=0;

    //nowe; BBBA05
    double M12=(PDG::mass_neutron+PDG::mass_proton)/2;
    double M2=M12*M12;


    double tau=Q2/4.0/M2;
//    cout<<Q2<<"->"<<M2<<"->"<<tau<<endl;
    if(ffset==1)
    {
        double tau2=tau*tau;
        double tau3=tau*tau2;
        double tau4=tau2*tau2;

        Gep=(1-0.0577*tau)/(1+11.2*tau+13.6*tau2+33*tau3);
        Gmp=2.79*(1+0.15*tau)/(1+11.1*tau+19.6*tau2+7.54*tau3);
        Gen=(1.38*tau-0.214*tau2)/(1+8.51*tau+59.9*tau2+13.6*tau3+2.57*tau4);
        Gmn=-1.91*(1+1.82*tau)/(1+14.1*tau+20.7*tau2+69.7*tau3);
    }
    if(ffset==0)
    {
        G_D = 1/(1+Q2/Mv2)/(1+Q2/Mv2);
        Gep = G_D;
        Gmp = 2.79*G_D;
        Gen = 0;
        Gmn = -1.91*G_D;
    }

    //~ double Ge_B=Gep-Gen;
    //~ double Gm_B=Gmp-Gmn;
    
//    cout<<"GEM="<<Ge<<" "<<Gm<<" "<<tau<<" "<<Gmp<<" "<<Gmn<<" "<<Gep<<" "<<Gen<<endl;

    double Fa=-g_a/pow(1+Q2/Ma2,2.0);
    //~ double Fp=0;
    double F2=(Gm-Ge)/(1+tau);
    double F2p=(Gmp-Gep)/(1+tau);
    double F2n=(Gmn-Gen)/(1+tau);

    double F1=(Ge+Gm*tau)/(1+tau);
    double F1p=(Gep+Gmp*tau)/(1+tau);
    double F1n=(Gen+Gmn*tau)/(1+tau);

//    cout <<F1<<" "<<F2<<endl;

    //~ double M1=M;
    M12 = (M_in+M_out)/2.;
    double M_minus = M_out-M_in;
    double M_plus = M_out+M_in;

    if(N1pdg==PDG::pdg_Lambda)
    {
        F1 = F1p;
        F2 = F2p;
    }
    else
    {
        F1 = (F1p + 2*F1n)/sqrt(3);
        F2 = (F2p + F2n)/sqrt(3);
    }
//    cout <<F1<<" "<<F2<<endl;

    //double xfd =1.;
    double F=0.463;
    double D=0.804;
    //~ double DF=D+F;
    //double xfd =0.2;
    //double xfd =F/(D+F);
    double xfd =F/(D);

    //Fa = (1+2*xfd)/3.*Fa;
    Fa = (1+2*xfd)/3.*Fa;


    double omega_1 = (M_minus*M_minus+Q2)/4/M12/M12*(F1+F2)*(F1+F2) + (M_plus*M_plus + Q2)/4/M12/M12*Fa*Fa;
    double omega_2 = Fa*Fa + (F1+F2- M_plus/2/M12*F2)*(F1+F2- M_plus/2/M12*F2) + Q2/M12/M12*(F2/2)*(F2/2);
    double omega_3 = 2*Fa*(F1+F2);

//    cout <<omega_1<<" "<<omega_2<<" "<<omega_3<<endl;

    double s_u = 4*M12*E-Q2-m2;
    double nu = (M_out*M_out -M_in*M_in+Q2)/2/M_in;
    double dif_1 = 8*M12*M12*Q2*omega_1;
    double dif_2 = (-4*(M12*M12/M_in/M_in*nu*nu+Q2) + s_u*s_u)*omega_2;
    double dif_3 = (2*s_u*Q2)*omega_3;

    double dif_0 = (G*G*(1-cos2thetac)/32/M12/M12/M_PI/E/E)*(dif_1+dif_2+dif_3);

//    cout<<dif_0<<endl;

    if(N1pdg==PDG::pdg_Sigma)
        dif_0/=2;
    
    if(dif_0==dif_0) //never return nan
        return  dif_0;
    else
        return 0;
}
