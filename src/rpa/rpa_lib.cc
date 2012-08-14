#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>

#include "jednostki.h"
#include "form_faktory_rozdzielnik.h"
#include "calg.h"
#include "density_of_nucleus.h"
#include "rpa_lib.h"
#include "../calg5.h"
#include "../util2.h"

#ifndef RPA_MAIN
	#include "../ff.h"
#endif

namespace rpa
{
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
	double En;             // changed by configure  
	int kFF=0;             // changed by configure
	double kf;             // changed by configure 
	int znak = 1;   	   // changed by configure 
						   // znak = 1 - particle 
						   // znak =-1 - antiparticle
	double m=m_mu;         // changed by configure
	double mm2 = m*m;      // changed by configure
	bool use_rpa = true;   // changed by configure
	double Mef;            // changed by configure
	double Mef2;           // changed by configure      
	double Ef;             // changed by configure
	double Ef2;            // changed by configure
	double sqrt_2=sqrt(2);

	int lambda_l =1;     
	 
	int l(1), t(1);
	
	 
	double Ming(double a,double b)
	{ 
		return (a<b ? a : b); 
	}

	double Maxg(double a, double b, double c)
	{  
		double x = (a>b ? a : b); 
		return (x > c ? x : c) ;
	}

				
	double qMin(double q0)
	{      
		if( ((En-q0)*(En-q0)-mm2) < 0 ) 
			return 0;  
		else
			return sqrt(2*En*En-2*En*q0+q0*q0-mm2-2*En*sqrt((En-q0)*(En-q0)-mm2 ));
			
	}

	double qMax(double q0)
	{
		if((En-q0)*(En-q0)-mm2 < 0 ) 
			return 0;	  
		else
			return sqrt(2*En*En - 2*En*q0 + q0*q0 - mm2 + 2*En*sqrt((En-q0)*(En-q0)-mm2 ));
	}

	double qMin2(double q0)
	{
		return 2*En*En - 2*En*q0 + q0*q0 - mm2 - 2*En*sqrt((En-q0)*(En-q0)-mm2 );
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


		double E_=Ming(Ef, Maxg( Mef , Ef - q0, 0.5*(- q0+qv * sqrt(1-4*Mef2/q2 ))));                     
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

		//double MA= 1.03*GeV;
		//double MV2=0.71*GeV*GeV;
		//double G_E= 1/( pow((1-q2/MV2),2) ); 
		//double G_M=G_E*magneton ;       
#ifdef RPA_MAIN		
		double F_1_= F_1(-q2,kFF);//(q2*G_M - 4*M12_2*G_E )/(q2 - 4*M12_2) ;
		double F_2_= F_2(-q2,kFF);//4*M12_2*( G_M - G_E )/(4*M12_2 - q2);    
		double G_A_= G_A(-q2,kFF); //-1.26/(pow(1-q2/(MA*MA),2));                 
#else
		double F_1_,F_2_,G_A_,F_P_;
		list(F_1_,F_2_)=f12(q2,0); // 0-cc, 1-nc proton, 2-nc neutron  
		list(G_A_,F_P_)=fap(q2,0);
#endif

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

	double sigma_q0( double q0)
	{ 

		if( (En-q0)*(En-q0)-mm2 < 0 )
		{ 
			return 0; cerr<<"amplituda nie liczona"<<endl;
		}
		else   
			return  calg5a(fix2(sigma_qv_q0,q0),qMin(q0),qMax(q0),100);
	}    

	struct config{
		int kFF;
		double kf;
		int use_Mf;
		int nu_pdg;
		int kNucleus;
		bool use_rpa;
		
		config(
			int	   FF0,
			double kf0,
			bool   use_Mf0,
			int    nu_pdg0,
			int    kNucleus0,
			bool   use_rpa0):
				kFF (FF0),
				kf(kf0),
				use_Mf(use_Mf0),
				nu_pdg(nu_pdg0),
				kNucleus (kNucleus0),
				use_rpa(use_rpa0)
		{
		}

	};


	config ust1a  (NNFF, 225*MeV, true,  14, Ar, true);
	config ust10a (NNFF, 225*MeV, true, -14, Ar, true);
	config ust1b  (NNFF, 225*MeV, false, 14, Ar, true);
	config ust10b (NNFF, 225*MeV, false,-14, Ar, true);
	config ust1c  (NNFF, 225*MeV, true,  14, Ar, false);

	config ust2a  (NNFF2, 225*MeV, true,  14, Ar, true);
	config ust20a (NNFF2, 225*MeV, true, -14, Ar, true);
	config ust2b  (NNFF2, 225*MeV, false, 14, Ar, true);
	config ust20b (NNFF2, 225*MeV, false,-14, Ar, true);
	config ust2c  (NNFF2, 225*MeV, true,  14, Ar, false);




	void configure(const config &ust)
	{
		kFF=ust.kFF;
		
		kf=ust.kf>0?ust.kf : mean_kf(ust.kNucleus);       

		if(ust.use_Mf)
		{
			Mef= ust.kf==225*MeV
			 		? 638 * MeV 
					: mean_Mef(ust.kNucleus);

		}  
		else 
			Mef=M;
			
		use_rpa =ust.use_rpa;
		   
		Mef2=Mef*Mef;              
		Ef= sqrt(kf*kf + Mef2);   
		Ef2=kf*kf + Mef2;

		znak=ust.nu_pdg>0 ?1 : -1;
		switch(ust.nu_pdg)
		{ 
			case 12:case -12: m=m_e;break;
		  	case 14:case -14: m=m_mu;break;
		  	case 16:case -16: m=m_tau;break;
		  	default: m=0;break;
	    }
		mm2=m*m;
	}


	void configure(double E, const int kNucleus, int nu_pdg, int use_Mf, int kFF0,double kf0, double mf0)
	{
		En = E;  
		kFF=kFF0;

		kf = kf0>0 ? kf0 : mean_kf(kNucleus);       

		Mef = use_Mf ? (mf0>0? mf0 :638 * MeV) : M; 
		Mef2=Mef*Mef;              		
		Ef= sqrt(kf*kf + Mef2);   
		Ef2=kf*kf + Mef2;
		
		znak=nu_pdg>0 ? 1 :-1;
		switch(nu_pdg)
		{ 
			case 12:case -12: m=m_e;break;
		  	case 14:case -14: m=m_mu;break;
		  	case 16:case -16: m=m_tau;break;
		  	default: m=0;break;
	    }
		mm2=m*m;
	}   

	double sigma_q0_q2(double q0,double q2)
	{
		double qv=sqrt(q2+q0*q0);
		return (sigma_qv_q0(qv,q0)/2/qv);
	}

	double sigma_q2(double q2, double dokl)
	{   
		double q01 =0;
		double q02 = En - (q2+mm2)/4/En - En*mm2/(q2+mm2);
		return calg5x(sigma_q0_q2,q2,min(q01,q02),max(q01,q02),dokl,30);
	}    


/*
	void plot_sigma_q2(double E,double dokl, int ilosc,const int kNucleus, int nu_pdg, int use_Mf, int kFF0=0,double kf0=0 )
	//void plot_sigma_q2(double E,double dokl, config &ust )
	{   
	 
		configure(E, kNucleus, nu_pdg, use_Mf, kFF0,kf0,0);    
		  
		string s2= znak>0?"_n":"_an";
		
		stringstream s;
		
		s<<"rpa_Q2_"<<nazwa_jadra(kNucleus)<<"E="<<En/GeV<<"_g="<<g_prim
		 <<"kfsr="<<kf<<"Mef="<<Mef<<s2<<form(kFF)<<".dat"<<flush;
		
		ofstream out(s.str().c_str());

		out << "# Q2[GeV2]     cross rpa    crosss mf  cratio rpa/mf " <<endl;

		double crossrpa1, crossmf1, ratio1;
		double crossrpa2, crossmf2, ratio2;
		double ratiorpa_ff, ratiomf_ff;
		
		double krok= 2*M*(En-m)/ilosc;
		
		for(double q2 =0.01; q2<2*M*(En-m); q2+=krok)
		{
			kFF = NNFF;	use_rpa  = true;	crossrpa1   = sigma_q2(q2,dokl)/(1e-38);
						use_rpa = false;   crossmf1  = sigma_q2(q2,dokl)/(1e-38);

			kFF = NNFF2; use_rpa  = true;	crossrpa2   = sigma_q2(q2,dokl)/(1e-38);
						use_rpa = false;   crossmf2  = sigma_q2(q2,dokl)/(1e-38);
		
			ratio1      = crossmf1 ? crossrpa1/crossmf1 : 0;
			ratio2      = crossmf2 ? crossrpa2/crossmf2 : 0;
			ratiorpa_ff = crossrpa2? crossrpa1/crossrpa2: 0;
			ratiomf_ff  = crossmf2 ? crossmf1 /crossmf2 : 0;


			out<< q2/GeV/GeV 
			   <<"\t"<<crossrpa1 << "\t" << crossmf1<< "\t"// << ratio1 
			   //<<"\t"<<crossrpa2 << "\t" << crossmf2<< "\t" << ratio2 
			   //<<"\t"<<ratiorpa_ff<< "\t" << ratiomf_ff
			   <<endl;

		}
	}   
		  
	void plot_sigma_q2alt(double E,double dokl, int ilosc, const config &ust )
	{   
		En=E;  
		configure(ust);

		string s1 = use_rpa? "rpa" : "FG";
		string s2= znak>0 ? "_n":"_an";
		 
		stringstream s;
		  s<< "rpa_Q2_"<< nazwa_jadra(ust.kNucleus)
		   << "E="	   << En/GeV
		   << "_g="	   << g_prim
		   << "kfsr="  << kf
		   << "Mef="   << Mef
		   << s2	   << form(kFF)  
		   <<"_"	   <<s1<<"_"
		   << ".dat"<<flush;
		 
		ofstream out(s.str().c_str());
		double krok= 2*M*(En-m)/ilosc;
		for(double q2 =0.01; q2<2*M*(En-m); q2+=krok)
	   out<< q2/GeV/GeV <<' ' <<GeV*GeV*0.5*liczba_atomowa(ust.kNucleus)*sigma_q2(q2,dokl)/(1e-38)<<endl;
		// out<< q2/GeV/GeV <<' ' <<0.5*liczba_atomowa(kNucleus )*sigma_q2(q2,dokl)/(1e-38)<<endl;

	} 
*/		  
		  
	void plot_sigma_q0(double E,double dokl, double krok,const int kNucleus, int nu_pdg, int use_Mf, int kFF0,double kf0=0 )
	{
		 configure(E, kNucleus, nu_pdg, use_Mf, kFF0,kf0);

		 string s1=znak==1?"n":"an";
		 
		 stringstream s;

		 s<<"rpa_q0_"<<nazwa_jadra(kNucleus)<<",E="<<En/GeV
		  <<",g="<<g_prim<<",kfsr="<<kf<<",Mef="<<Mef<<form(kFF)<<",rpa="<<use_rpa<<",nu="<<nu_pdg<<".dat"<<flush;

		ofstream out(s.str().c_str());
		double frac=0;
        switch(kNucleus)
        {
			case C:
			case O: frac= 0.5;break;
			case Ar:frac= 18./40.;break;
			case Fe:frac= 26./56.;break;
		}
		if(nu_pdg>0)
			frac=1-frac;
		for(double q0=5*MeV; q0 < En-m; q0+=krok)
			out<< q0/GeV <<' ' <<frac*sigma_q0(q0)*GeV/(1e-38)<<endl;
	}



	double sigma()
	{
		return calg5a(sigma_q0,0 ,En - m);	
	}
/*		  
	double sigma_Atom(double E,int kNucleus, int nu_pdg, double dokl, int use_Mf, int kFF0, double kf0)
	{

		configure(E, kNucleus, nu_pdg, use_Mf, kFF0,kf0);

		return 0.5*liczba_atomowa(kNucleus)*sigma();      
	}

		  
	void plot_sigma(double zakres,double krok, int kNucleus,int nu_pdg,double dokl,int use_Mf,int kFF0, double kf0=0 )
	{  
	  
		configure(m+50*MeV, kNucleus, nu_pdg, use_Mf, kFF0,kf0);

		string s2=znak>0? "_n":"_an";
		 
		stringstream s;
		s<<"rpa_cal_"<<nazwa_jadra(kNucleus)<<"_g="<<g_prim<<"kfsr="<<kf<<"Mef="
		 <<Mef<<s2<<form(kFF)<<".dat"<<flush;
		
		ofstream out(s.str().c_str());
		
		for(double En =  m+50*MeV; En<zakres ;En+=krok)
			out<< En/GeV  <<' '<<0.5*liczba_atomowa(kNucleus)*sigma()/(1e-38)<<endl;
	}
*/
	double ratio_rpa_fg(int qel_rpa, double q0, double qv)//, double kf0, int kFF0=0,int nu_pdg=1,int kNucleus=0, bool use_Mf0=0)
	{

//		configure(E, kNucleus, nu_pdg, use_Mf0, kFF0, kf0);
	    switch(qel_rpa)
	    {	
			case 0:return 1;
			case 1:
			{
				ratio=true;
				double res=sigma_qv_q0(qv,q0);
				ratio=false;
				return res;
			}
			case 2:
			{
				ratio=false;
				use_rpa = true;
				double sigma_rpa=sigma_qv_q0(qv,q0);
				use_rpa = false; 
				double sigma_fg=sigma_qv_q0(qv,q0);
				return sigma_fg ? sigma_rpa/sigma_fg : 1;
			}
			default: return 1;
		}
	}
}


#ifdef RPA_MAIN

using namespace rpa;

int main()
{
    double Enu=1*GeV;
    int kNucleus=Ar;
    double kf=mean_kf(kNucleus);
    int nu=14;
    
    cout<<"Nucleus= "<<nazwa_jadra(kNucleus)<<endl;
    cout<<"kf= "<<kf/MeV<<" MeV"<<endl;
    cout<<"E= "<<Enu/MeV<<" MeV"<<endl;
    cout<<"nupdg= "<<nu<<endl;
    
	use_rpa=false;
	plot_sigma_q0(Enu,1e-3, 1*MeV,kNucleus , nu , false, BBA, kf );// qel_rpa=0
	En=Enu;
	cout<<sigma()<<endl;
	use_rpa=true;                            
	plot_sigma_q0(Enu,1e-3, 1*MeV,kNucleus , nu , false, BBA, kf );//qel_rpa=1
	En=Enu;
	cout<<sigma()<<endl;
	plot_sigma_q0(Enu,1e-3, 1*MeV,kNucleus , nu , true , BBA, kf ); //qel_rpa=2
	En=Enu;
	cout<<sigma()<<endl;

	if(false)
	{  
		for(double q0=0*MeV;q0<rpa::En-rpa::m;q0+=10*MeV)
		{
			for(double qv=rpa::qMin(q0);qv<rpa::qMax(q0);qv+=10*MeV)
			{
				double x=rpa::ratio_rpa_fg(rpa::En, q0, qv);
				if(x==x)
					 cout<<x<<endl;
			}
		}
		                                                          
		 
	}
}

#endif
