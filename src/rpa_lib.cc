#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <cassert>

#include "jednostki.h"
#include "rpa_lib.h"
#include "ff.h"

namespace rpa
{
	const double  M1=938.272*MeV; //masa protonu
	const double  M2=939.565*MeV; //masa neutronu
	const double  M12=(M1+M2)/2;  //¶rednia masa nukleonu 
	const double  m_e=0.5109989*MeV;      // masa elektronu
	const double  m_mu=105.658357*MeV;    // masa muonu           
	const double  m_tau=1777*MeV;         // masa tau 
	const double  m_pi = 139.57018*MeV;   // masa mezonu pi
	const double  m_r  = 769.3 *MeV;      // masa mezonu rho  
	using namespace std;

	typedef complex<double> zesp;
	
	
	double M=M12;             
	double M12_2=M12*M12;
	//double magneton=4.71;      
	double m_pi2= m_pi*m_pi;
	double m_r2 =m_r*m_r;
	double gr =sqrt(1.64*4*Pi); 
	double gr2=gr*gr;
	double fr = 6.1*gr;          
	double fr2 =fr*fr;
	double fpi = sqrt(4*Pi*0.075);
	double fpi2= fpi*fpi;
	double g_prim = 0.7;           
	bool ratio=false;      // used by ratio_rpa_fg;
	bool use_rpa = true;   
	bool use_mf = true;
	double En;             // changed by configure  
	double m=m_mu;         // changed by configure
	int znak = 1;   	   // changed by configure: 1 - particle, -1 - antiparticle
	double kf;             // changed by configure 
	double mf=M;          // changed by configure
	double sqrt_2=sqrt(2);

	int lambda_l =1;     
	 
	int l(1), t(1);
	 
	double max3(double a, double b, double c)
	{  
		return max(a,max(b,c));
	}

				
	double stala = 1.14e-11;   /// tak bylo

	//double stala = G*0.9737;

	double sigma_qv_q0(double qv,double q0)
	{
	//    cout<<"En="<<En/MeV<<endl;
	//    cout<<"qv="<<qv/MeV<<endl;
	//    cout<<"q0="<<q0/MeV<<endl;
		double q02=q0*q0;
		double q04=q02*q02;
		double qv2=qv*qv;
		double qv4=qv2*qv2;
		double qv3=qv2*qv;
		double q2 = q02-qv2;
		double q4 =q2*q2;
		double Mef=use_mf? mf :M;
		double Mef2=Mef*Mef;              		
		double Ef= sqrt(kf*kf + Mef2);   
		double Ef2=kf*kf + Mef2;
		double mm2=m*m;

		double alfa=log(fabs(q4 - 4*pow( q0*Ef - qv*kf, 2 ) )/fabs(q4 - 4*pow(q0*Ef + qv*kf,  2)) );
		double beta=log(fabs((pow(q2 + 2*qv*kf, 2) - 4*q02*Ef2))/fabs(pow(q2 - 2*qv*kf, 2) - 4*q02*Ef2 ) ) ;
		double stal=sqrt(q2*(q2 -4*Mef2));
		double czynnik1=fabs(q4*Ef - 4*q0*Mef2*(q0*Ef-qv*kf) - kf*q2*stal)/fabs(q4*Ef - 4*q0*Mef2*(q0*Ef-qv*kf) + kf*q2*stal);
		double czynnik2=fabs(q4*Ef - 4*q0*Mef2*(q0*Ef+qv*kf) - kf*q2*stal)/fabs(q4*Ef - 4*q0*Mef2*(q0*Ef+qv*kf) + kf*q2*stal);

		double LAMBDA =(log(czynnik1 ) + log(czynnik2)  )/2/stal; 



		double nawias_ReHs = kf*Ef - (3*Mef2 - q2/2)*log((kf+Ef)/Mef) - q0*(4*Mef2-q2)*alfa/8/qv + 
									Ef*(4*Mef2-q2)*beta/4/qv +(4*Mef2-q2)*(4*Mef2-q2)*LAMBDA/4; 
		double ReHs= lambda_l*nawias_ReHs /2/Pi/Pi;
		double nawias_ReHl= 2*kf*Ef/3 - (qv2/6)*log((kf+Ef)/Mef) + (q0/(4*qv))*(Ef2 + (q02-3*qv2)/12 )*alfa 
			- Ef*(3*q2 + 4*Ef2)*beta/(24*qv) + qv2*(2*Mef2+q2)*(4*Mef2-q2)*LAMBDA/12/q2 ;
		double ReHl=nawias_ReHl*lambda_l*q2/Pi/Pi/qv2; 
		double nawias_ReHt= kf*Ef*(1+ 2*q02/qv2 )/3 + (q2/3)*log((kf+Ef)/Mef)
			+ (q0*q2/(4*qv3))*((q02 + 3*qv2)/12 + qv2*Mef2/q2+ Ef2)*alfa
			+ (Ef/qv)*((qv4-q04)/(8*qv2) - Mef2/2 - q2*Ef2/(6*qv2))*beta
			-  (2*Mef2+q2)*(4*Mef2 - q2)*LAMBDA/6; 

		double ReHt= 2*nawias_ReHt*lambda_l/2/Pi/Pi; 
		double ReH0=(qv*kf + q0*Ef*alfa/2 -(q2/8 + Ef2/2 + qv2*Mef2/(2*q2))*beta)* lambda_l*Mef/(2*Pi*Pi*qv);


		double E_=min(Ef, max3( Mef , Ef - q0, 0.5*(- q0+qv * sqrt(1-4*Mef2/q2 ))));                     
		double E1= Ef - E_; 
		double E2=(Ef*Ef -E_*E_)/2;
		double E3=(Ef*Ef*Ef -E_*E_*E_)/3; 
				 

		double Re_H_a=( - 2*ReHs - ReHl + ReHt)/3;
		double Re_H_vv_l=ReHl;
		double Re_H_vv_t=ReHt;
		double Re_H_va=q2*ReH0/2/qv/qv/Mef;
		double Re_H_tt_l=q2*( -ReHl + ReHs - 2*ReHt)/ 12 / M/M;
		double nawias1_Re_H_tt_t = ReHl*(q2 - 8*Mef2) + ReHt*(2*Mef2 - q2) + 2*ReHs*(q2-2*Mef2);
		double nawias2_Re_H_tt_t = q2*ReH0 - 2*qv*qv*Mef*Re_H_va;
		double Re_H_tt_t=q2*(nawias1_Re_H_tt_t/24/M/M/Mef2 + nawias2_Re_H_tt_t/4/q0/M/M/Mef);
		double Re_H_vt_l=-q2*Re_H_a/4/M/Mef;
		double Re_H_vt_t=q2*Re_H_a/2/M/Mef ;

		double Im_H_vv_l=(q2/(2*Pi*qv3))*( E3 + q0*E2 + q2*E1/4 ) ; 
		double Im_H_vv_t=2*(q2/(4*Pi*qv3))*(E3+q0*E2+(qv2*Mef2/q2+0.25*(q02+qv2))*E1) ;
		double Im_H_tt_l=- ( q2/( 8*Pi*qv3*M*M))*( q2*E3 + q0*q2*E2+( qv2*Mef2 +0.25*q2*q02 )*E1 );
		double Im_H_tt_t=( q2/(8*Pi*qv3*M*M) )*(  (Mef2*qv2- 0.25*q4)*E1- q0*q2*E2 -q2*E3 ) ;
		double Im_H_vt_l=- q2*E1*Mef /(8*Pi*qv*M); 
		double Im_H_vt_t= -2*Im_H_vt_l;                             // dwojka
		double Im_H_va=qv*(q2/(8*Pi*qv3))*(2*E2 + q0*E1);
		double Im_H_a=(Mef2 * E1 )/(2*Pi*qv);

		zesp H_vv_l=zesp(Re_H_vv_l,Im_H_vv_l); 
		zesp H_vv_t=-0.5* zesp(Re_H_vv_t,Im_H_vv_t); 
		zesp H_tt_l=zesp(Re_H_tt_l,Im_H_tt_l); 
		zesp H_tt_t=-0.5*zesp(Re_H_tt_t,Im_H_tt_t); 
		zesp H_vt_l=      zesp(Re_H_vt_l,Im_H_vt_l); 
		zesp H_vt_t= -0.5*zesp(Re_H_vt_t,Im_H_vt_t); 
		zesp H_va  =  zesp(qv*Re_H_va, Im_H_va) ;
		zesp H_a   =  zesp(Re_H_a,Im_H_a); 

		zesp H_rr_l=(0.5*gr2*H_vv_l + fr*gr*H_vt_l + 0.5*fr2*H_tt_l);
		zesp H_rr_t=(0.5*gr2*H_vv_t + fr*gr*H_vt_t + 0.5*fr2*H_tt_t);
		zesp H_rp_va=-fpi*(gr + fr *Mef/M)*H_va/m_pi;
		zesp H_pp_l=2.*fpi2*H_vv_l/m_pi2;
		zesp H_pp_t=2.*fpi2*H_vv_t/m_pi2;
		zesp H_pp_a=2.*fpi2*H_a/m_pi2;

		double F_1_,F_2_,G_A_,F_P_;
		list(F_1_,F_2_)=f12(q2,0); // 0-cc, 1-nc proton, 2-nc neutron  
		list(G_A_,F_P_)=fap(q2,0);

		double V_l=- q2/(q2 - m_pi2);
		double V_t=V_l;                                        	        
		double V_a=-V_l - g_prim;
		double R_t= -q2/(q2 - m_r2)/m_r2;              
		double R_a=1/m_r2;
		double R_l= R_t;      

		double R_ta=R_t+R_a;
		double R_la=R_l+R_a;
		double V_ta=V_t+V_a;
		double V_la=V_l+V_a;

		zesp H_r_l= gr*(F_1_*H_vv_l + F_2_*H_vt_l)/sqrt_2+
				  fr*(F_1_*H_vt_l + F_2_*H_tt_l)/sqrt_2;

		zesp H_r_t=gr*(F_1_*H_vv_t + F_2_*H_vt_t)/sqrt_2 + fr*(F_1_*H_vt_t + F_2_*H_tt_t)/sqrt_2;

		zesp H_r_va=G_A_*(gr+fr*Mef/M)*H_va/sqrt_2;
		  
		zesp H_p_l= -sqrt_2*fpi*G_A_*H_vv_l/m_pi;
		  
		zesp H_p_t= -sqrt_2*fpi*G_A_*H_vv_t/m_pi;

		zesp H_p_va= -sqrt_2*fpi*(F_1_+F_2_*Mef/M)*H_va/m_pi;
	 
		zesp H_p_a=-sqrt_2*fpi*H_a*G_A_/m_pi;
	  
		zesp delta_a_1=R_a;
		  
		zesp delta_a_4= V_a/(1.-V_a*H_pp_a);
		zesp delta_l_1=(R_l+R_la*H_rr_l*R_a)/(1. - R_la*H_rr_l );
		  
		zesp licznik_delta_l_4 =V_l+V_la*H_pp_l*V_a;
		zesp mianownik_delta_l_4 = (1.0 - V_a*H_pp_a )*(1.0 - V_la*(H_pp_l+H_pp_a));
		zesp delta_l_4=licznik_delta_l_4/mianownik_delta_l_4;

		 
		zesp licznik_delta_t_1=(1.- V_ta*(H_pp_a+H_pp_t)) *( R_t + R_a* R_ta*H_rr_t  )
								+R_a*R_ta*V_ta*pow(H_rp_va,2);
		zesp mianownik_delta_t_1 = (1. - V_ta*(H_pp_a + H_pp_t))*(1. - R_ta*H_rr_t) - R_ta*V_ta*pow(H_rp_va,2);
		zesp delta_t_1 = licznik_delta_t_1/mianownik_delta_t_1;
		 
		zesp licznik_delta_t_4 =(1. - R_ta*H_rr_t)*(V_t +delta_a_4*(V_ta*H_pp_t+V_t*H_pp_a))
											 +delta_a_4*R_ta*V_ta*pow(H_rp_va,2);
		zesp mianownik_delta_t_4= (1. - V_ta*(H_pp_a + H_pp_t ) )*(1. - R_ta*H_rr_t)- R_ta*V_ta*pow(H_rp_va,2);
		zesp delta_t_4 =licznik_delta_t_4/mianownik_delta_t_4;
		 

		 
		zesp delta_va_2=R_ta*H_rp_va*(delta_t_4+delta_a_4)/(1. - R_ta*H_rr_t );
		 
		zesp delta_va_3=V_ta*H_rp_va*( R_a + delta_t_1)/(1. - V_ta*(H_pp_a + H_pp_t) );


		double poprawka_a=imag(delta_a_4*H_p_a*H_p_a);

		double poprawka_l= imag( ( delta_a_4+delta_l_4)*(pow(H_p_l ,2) + 2.*H_p_l*H_p_a ))
					   +imag( pow(H_p_a,2 )*delta_l_4 + pow( H_r_l,2 )*( delta_a_1+delta_l_1));

		double poprawka_t =  imag((pow(H_p_t,2) + 2.*H_p_t*H_p_a +pow(H_p_va,2))*
						(delta_a_4 +delta_t_4) + delta_t_4*pow(H_p_a,2))+
							imag((pow(H_r_t,2) + pow(H_r_va,2))*(delta_t_1+delta_a_1))+				    
							imag( (H_r_va*H_p_t+H_r_va*H_p_a+H_r_t*H_p_va)*(delta_va_2+delta_va_3));
			   
		double poprawka_va = imag((H_r_t*H_p_t + H_r_t*H_p_a+H_r_va*H_p_va)*(delta_va_2 + delta_va_3))   
							+2.*imag( H_p_va*( H_p_t+H_p_a )*( delta_t_4 + delta_a_4)
						  +H_r_va*H_r_t*(delta_t_1 + delta_a_1));
		double R_L[2],R_T[2],R_A[2],R_VA[2];  
		
		R_L[0]=( F_1_*F_1_ + G_A_*G_A_ )*Im_H_vv_l+ F_2_*F_2_*Im_H_tt_l+ 2*F_1_*F_2_*Im_H_vt_l;
		
		R_L[1]=R_L[0] + poprawka_l;

		R_T[0]=Im_H_vv_t*(G_A_*G_A_+ F_1_*F_1_)/2 + Im_H_tt_t*F_2_*F_2_/2 + F_1_*F_2_*Im_H_vt_t;

		R_T[1]=R_T[0] - poprawka_t; 

		R_A[0]=G_A_*G_A_*Im_H_a;

		R_A[1]=R_A[0]  + poprawka_a;

		R_VA[0]=(F_1_+Mef/M*F_2_)*G_A_*Im_H_va;
		
		R_VA[1]=R_VA[0] +poprawka_va;


		double L_L=-(16*En*(En-q0)-4*(mm2-q2))*q2/qv2 - 4*mm2*q0*(4*En - q0 + q0*mm2/q2 )/qv2;

		double L_T=-(16*En*(En-q0)-4*(mm2-q2))*q2/qv2 
					- 4*mm2*(4*En*q0- (q0*q0-qv*qv) + mm2)/qv2 - 8*(q2-mm2);

		double L_A=8*(q2-mm2);

		double L_VA=-16*(q2*(2*En-q0) + q0*mm2)/qv;


		double amplituda[2];
		amplituda[0]=    L_L*R_L[0] +(L_A+L_T)*R_A[0] + L_T*(R_T[0] -R_A[0]) -znak*L_VA*R_VA[0];   
		amplituda[1]=    L_L*R_L[1] +(L_A+L_T)*R_A[1] + L_T*(R_T[1] -R_A[1]) -znak*L_VA*R_VA[1];   
	  
	/*	  if(amplituda > 0)
		  cerr<<"amplituda zla( E="<< En/GeV 
			  <<", q0 = "<< q0/GeV<<", Q2 = " 
			  << (qv*qv-q0*q0)/GeV/GeV 
			  << ", (En-q0)*(En-q0)-mm2 = " << 
		  ((En-q0)*(En-q0)-mm2)/GeV/GeV<<")"<<endl;
	*/	 if(!ratio)
		  return -amplituda[use_rpa]*stala*stala * qv /(16 * Pi * Pi * (kf*kf*kf/3/Pi/Pi)  * En * En)/cm2;
		 else
		  return amplituda[0] ? min(amplituda[1]/amplituda[0],10.0) :1;
	}

	void configure(double E, int nu_pdg, double kf0, double mf0)
	{
		En = E;  
		kf=kf0;
		mf=mf0;
		
		znak= nu_pdg>0 ? 1 :-1;
		
		switch(nu_pdg)
		{ 
			case 12:case -12: m=m_e;break;
		  	case 14:case -14: m=m_mu;break;
		  	case 16:case -16: m=m_tau;break;
		  	default: m=0;break;
	    }
	}

	double ratio_rpa_fg(int qel_rpa, double q0, double qv)
	{

	    switch(qel_rpa)
	    {	
			case 0:return 1;
			case 1:
			{
				use_rpa = use_mf=false;
				ratio=true;
				double res=sigma_qv_q0(qv,q0);
				ratio=false;
				return res;
			}
			case 2:
			{
				ratio= use_rpa = use_mf=false;
				double sigma_fg=sigma_qv_q0(qv,q0);
				//assert(sigma_fg>=0);
				use_rpa = use_mf=true ; 
				double sigma_rpa=sigma_qv_q0(qv,q0);
				//assert(sigma_rpa>=0);
				return sigma_fg>0 && sigma_rpa>=0 ? min(max(sigma_rpa/sigma_fg,0.),10.) : 1;
			}
			default: return 1;
		}
	}


}
