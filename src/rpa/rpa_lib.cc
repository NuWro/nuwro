#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>


#include "./rpa/jednostki.h"
#include "./rpa/form_faktory_rozdzielnik.h"
#include "./rpa/calg.h"
#include "./rpa/density_of_nucleus.h"
#include "./rpa/tablica_Mef.h"
#include "./rpa/flux_class.h"

namespace rpa
{

using namespace std;

typedef complex<double> zesp;

double M=M12;             
double En;                 
double M12_2=M12*M12;
//double magneton=4.71;      
double m_pi2= m_pi*m_pi;
double m_r2 =m_r*m_r;
double kf;
double gr =sqrt(1.64*4*Pi); 
double gr2=gr*gr;
double fr = 6.1*gr;          
double fr2 =fr*fr;
double fpi = sqrt(4*Pi*0.075);
double fpi2= fpi*fpi;
double g_prim = 0.7;           
int FF;
double mm2 = m_mu*m_mu;
int znak = 1;  
int nucleus;
int neutrino=1;
int antyneutrino=-1;

bool RPA_OR_NOT = true;

double Mef;         
double Mef2;              
double iloraz;
double Ef;   
double Ef2;
double ro;
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

double Przekroj(double qv,double q0)
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
	zesp H_rp_va=-fpi*(gr + fr *iloraz)*H_va/m_pi;
	zesp H_pp_l=2.*fpi2*H_vv_l/m_pi2;
	zesp H_pp_t=2.*fpi2*H_vv_t/m_pi2;
	zesp H_pp_a=2.*fpi2*H_a/m_pi2;

	//double MA= 1.03*GeV;
	double G_A_= G_A(q0,qv,FF); //-1.26/(pow(1-q2/(MA*MA),2));                 
	//    double MV2=0.71*GeV*GeV;
	//    double G_E= 1/( pow((1-q2/MV2),2) ); 
	//    double G_M=G_E*magneton ;       
	double F_1_= F_1(q0,qv,FF);//(q2*G_M - 4*M12_2*G_E )/(q2 - 4*M12_2) ;
	double F_2_ = F_2(q0,qv,FF);//4*M12_2*( G_M - G_E )/(4*M12_2 - q2);    


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

    zesp H_r_va=G_A_*(gr+fr*iloraz)*H_va/sqrt_2;
  	  
    zesp H_p_l= -sqrt_2*fpi*G_A_*H_vv_l/m_pi;
	  
    zesp H_p_t= -sqrt_2*fpi*G_A_*H_vv_t/m_pi;

    zesp H_p_va= -sqrt_2*fpi*(F_1_+F_2_*iloraz)*H_va/m_pi;
 
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
	  
	double R_L=( F_1_*F_1_ + G_A_*G_A_ )*Im_H_vv_l+ F_2_*F_2_*Im_H_tt_l+
					2*F_1_*F_2_*Im_H_vt_l + RPA_OR_NOT*poprawka_l;

	double R_T=Im_H_vv_t*(G_A_*G_A_+ F_1_*F_1_)/2
			   + Im_H_tt_t*F_2_*F_2_/2 + F_1_*F_2_*Im_H_vt_t  -  RPA_OR_NOT*poprawka_t; 

	double R_A=G_A_*G_A_*Im_H_a  + RPA_OR_NOT*poprawka_a;

	double R_VA=(F_1_+iloraz*F_2_)*G_A_*Im_H_va  + RPA_OR_NOT*poprawka_va;


	double L_L=-(16*En*(En-q0)-4*(mm2-q2 ) )*q2/qv2 -4*mm2*q0*(4*En - q0 + q0*mm2/q2 )/qv2;

	double L_T=-(16*En*(En-q0)-4*(mm2 - q2))*q2/qv2
			- 4*mm2*(4*En*q0- (q0*q0-qv*qv) + mm2)/qv2 - 8*(q2-mm2);

	double L_A=8*(q2-mm2);

	double L_VA=-16*(q2*(2*En-q0) + q0*mm2)/qv;


	double amplituda=    L_L*R_L +(L_A+L_T)*R_A + L_T*(R_T -R_A) -znak*L_VA*R_VA;   
  
/*	  if(amplituda > 0)
	  cerr<<"amplituda zla( E="<< En/GeV 
	      <<", q0 = "<< q0/GeV<<", Q2 = " 
	      << (qv*qv-q0*q0)/GeV/GeV 
	      << ", (En-q0)*(En-q0)-mm2 = " << 
	  ((En-q0)*(En-q0)-mm2)/GeV/GeV<<")"<<endl;
*/	 
	  return -amplituda*stala*stala * qv /(16 * Pi * Pi * ro  * En * En)/cm2;
}

double przekroj( double q0)
{ 

	if( (En-q0)*(En-q0)-mm2 < 0 )
	{ 
		return 0; cerr<<"amplituda nie liczona"<<endl;
	}
	else   
		return  calg5x(Przekroj,q0,qMin(q0),qMax(q0),1e-40,50);
}    

struct USTAWIENIA{
	int form_factors;
	double fermimomentum;
	int WLACZNIK;
	int nu_or_anti;
	int nucleus_kind;
	bool rpa;
};

int effective_mass =1;
int no_effective_mass =0;

USTAWIENIA ust1a ={NNFF, 225*MeV, effective_mass, neutrino, Ar, true };
USTAWIENIA ust1b ={NNFF, 225*MeV, no_effective_mass, neutrino, Ar, true};
USTAWIENIA ust1c ={NNFF, 225*MeV, effective_mass, neutrino, Ar, false };

USTAWIENIA ust2a ={NNFF2, 225*MeV, effective_mass, neutrino, Ar, true };
USTAWIENIA ust2b ={NNFF2, 225*MeV, no_effective_mass, neutrino, Ar, true };
USTAWIENIA ust2c ={NNFF2, 225*MeV, effective_mass, neutrino, Ar, false};

USTAWIENIA ust10a ={NNFF, 225*MeV, effective_mass, antyneutrino, Ar, true };
USTAWIENIA ust10b ={NNFF, 225*MeV, no_effective_mass, antyneutrino, Ar, true };

USTAWIENIA ust20a ={NNFF2, 225*MeV, effective_mass, antyneutrino, Ar, true };
USTAWIENIA ust20b ={NNFF2, 225*MeV, no_effective_mass, antyneutrino,Ar, true };


void roz(USTAWIENIA &ust)
{
       
	FF=ust.form_factors;

	if(ust.fermimomentum==0) 
		kf=sredni_ped_fermiego(ust.nucleus_kind);       
	else 
		kf = ust.fermimomentum; 

	if(ust.WLACZNIK==1)
	{
		if(ust.fermimomentum==(225*MeV))
			Mef =638 * MeV;
		else 
			sredni_Mef(ust.nucleus_kind);//638 * MeV;        

		RPA_OR_NOT =ust. rpa;
	}  
	else 
		Mef=M;
	   
	Mef2=Mef*Mef;              
	iloraz=Mef/M;
	Ef= sqrt(kf*kf + Mef2);   
	Ef2=kf*kf + Mef2;
	ro=(kf*kf*kf)/(3*Pi*Pi);

	znak=ust.nu_or_anti;
}


void rozdzielnik(double E, const int rodzaj_jadra, int rodzaj, int wlacznik, int FF_form,double ped)
{
	En = E;  
	FF=FF_form;

	if(ped==0) 
		kf = sredni_ped_fermiego(rodzaj_jadra);       
	else 
		kf = ped; 
	if(wlacznik==1)
		Mef = 638 * MeV;//sredni_Mef(rodzaj_jadra);//638 * MeV;          
	else 
		Mef = M;

	Mef2=Mef*Mef;              
	iloraz=Mef/M;
	Ef= sqrt(kf*kf + Mef2);   
	Ef2=kf*kf + Mef2;
	ro=(kf*kf*kf)/(3*Pi*Pi);

	znak=rodzaj;

}   

double Przekroj_q2(double q0,double q2)
{
	double qv=sqrt(q2+q0*q0);
	return (Przekroj(qv,q0)/2/qv);
}

double przekroj_q2(double q2, double dokl)
{ 
	  
	double q01 =0;
	double q02 =	En - (q2+mm2)/4/En - En*mm2/(q2+mm2);
	double q0min=q01;
	double q0max=q02;
	if(q01>q02)
	{
		q0min=q02;
		q0max=q01;  
	}
	return  calg5x(Przekroj_q2,q2,q0min,q0max,dokl,30);
}    



void plot_przekroj_q2(double E,double dokl, int ilosc,const int rodzaj_jadra, int rodzaj, int wlacznik, int FF_form,double ped=0 )
//void plot_przekroj_q2(double E,double dokl, USTAWIENIA &ust )
{   
 
    rozdzielnik(E, rodzaj_jadra, rodzaj, wlacznik, FF_form,ped);    
	stringstream s;
      

    string s2;
     
    if(znak==1) 
    	s2 ="_n";
	else
		s2 ="_an";
     
//      s<<"rpa_Q2_"<<nazwa_jadra(ust.nucleus_kind)<<"E="<<En/GeV<<"_g="<<g_prim
//      <<"kfsr="<<kf<<"Mef="<<Mef<<s2<<form(FF)<<".dat"<<flush;
  
	s<<"rpa_Q2_"<<nazwa_jadra(rodzaj_jadra)<<"E="<<En/GeV<<"_g="<<g_prim
     <<"kfsr="<<kf<<"Mef="<<Mef<<s2<<form(FF)<<".dat"<<flush;
    
    double crossrpa1, crossmf1, ratio1;
    double crossrpa2, crossmf2, ratio2;
    double ratiorpa_ff, ratiomf_ff;
    ofstream out(s.str().c_str());
    
    out << " Q2[GeV2]     cross rpa    crosss mf  cratio rpa/mf " <<endl;
    
    double krok= 2*M*(En-m_mu)/ilosc;
    
    for(double q2 =0.01; q2<2*M*(En-m_mu); q2+=krok)
    {
		FF = NNFF;
		RPA_OR_NOT  = true;
		crossrpa1   = przekroj_q2(q2,dokl)/(1e-38);
		RPA_OR_NOT = false;
		crossmf1  = przekroj_q2(q2,dokl)/(1e-38);

		FF = NNFF2;
		RPA_OR_NOT  = true;
		crossrpa2   = przekroj_q2(q2,dokl)/(1e-38);
		
		RPA_OR_NOT = false;
		crossmf2  = przekroj_q2(q2,dokl)/(1e-38);
    
		if(crossmf1==0)
			ratio1 = 0;
		else 
			ratio1 = crossrpa1/crossmf1;

		if(crossmf2==0)
			ratio2 = 0;
		else 
			ratio2 = crossrpa2/crossmf2;

		if(crossrpa2==0)
			ratiorpa_ff =0;
		else 
			ratiorpa_ff = crossrpa1/crossrpa2;

		if(crossmf2==0)
			ratiomf_ff =0;
		else 
			ratiomf_ff = crossmf1/crossmf2;


		out<< q2/GeV/GeV 
		   <<"\t"<<crossrpa1 << "\t" << crossmf1<< "\t" << ratio1 
		   <<"\t"<<crossrpa2 << "\t" << crossmf2<< "\t" << ratio2 
		   <<"\t"<<ratiorpa_ff<< "\t" << ratiomf_ff
		   <<endl;

	}
}   
      
void plot_przekroj_q2alt(double E,double dokl, int ilosc, USTAWIENIA &ust )
{   
	En=E;  
	roz(ust);

	stringstream s;

	string s1;
	if(RPA_OR_NOT)
		s1="rpa";
	else  
		s1="FG";

	string s2;

	if(znak==1) 
		s2 ="_n";
	else 
		s2 ="_an";
 
     
     
      s<< "rpa_Q2_"
       << nazwa_jadra(ust.nucleus_kind)
       << "E="
       << En/GeV
       << "_g="
       << g_prim
       << "kfsr="
       << kf
       << "Mef="
       << Mef
       << s2
       << form(FF)
       <<"_"<<s1<<"_"
       << ".dat"<<flush;
  
     
    ofstream out(s.str().c_str());
    double krok= 2*M*(En-m_mu)/ilosc;
    for(double q2 =0.01; q2<2*M*(En-m_mu); q2+=krok)
   out<< q2/GeV/GeV <<' ' <<GeV*GeV*0.5*liczba_atomowa(ust.nucleus_kind)*przekroj_q2(q2,dokl)/(1e-38)<<endl;
	// out<< q2/GeV/GeV <<' ' <<0.5*liczba_atomowa(rodzaj_jadra )*przekroj_q2(q2,dokl)/(1e-38)<<endl;

} 
      
      
void plot_przekroj_density_q0(double E,double dokl, double krok,const int rodzaj_jadra, int rodzaj, int wlacznik, int FF_form,double ped=0 )
{
     rozdzielnik(E, rodzaj_jadra, rodzaj, wlacznik, FF_form,ped);

     stringstream s;
     //char* s1;
     string s1;
     
     if(znak==1)
	     s1 ="n";
     else 
     	s1 ="an";
     
     s<<"rpa_q0_"<<nazwa_jadra(rodzaj_jadra)<<"_E="<<En/GeV
      <<"g="<<g_prim<<"kfsr="<<kf<<"Mef="<<Mef<<form(FF)<<".dat"<<flush;

    ofstream out(s.str().c_str());

    
    for(double q0=5*MeV; q0 < En-m_mu; q0+=krok)
    	out<< q0/GeV <<' ' <</*0.5*liczba_atomowa(rodzaj_jadra)*/przekroj(q0)*GeV/(1e-38)<<endl;
}



double calkowity(double E, double dokl)
{

	En=E;
	return calg5(przekroj,0 ,En - m_mu, dokl, 20);	
}
	  
double calkowity_liczba(double E,int rodzaj_jadra, int rodzaj, double dokl, int wlacznik, int FF_form, double ped)
{

	rozdzielnik(E, rodzaj_jadra, rodzaj, wlacznik, FF_form,ped);

	return 0.5*liczba_atomowa(rodzaj_jadra)*calkowity(En,dokl);      
}

	  
void plot_calkowity_density(double zakres,double krok, int rodzaj_jadra,int rodzaj,double dokl,int wlacznik,int FF_form, double ped=0 )
{  
	  /*
      FF=FF_form;
      
      if(ped==0) kf=sredni_ped_fermiego(rodzaj_jadra);       
	  else kf =ped; 
	  
	  if(wlacznik)
	   Mef =638*MeV;//sredni_Mef(rodzaj_jadra);          
           else Mef= M;
	   
	   Mef2=Mef*Mef;              
           iloraz=Mef/M;
           Ef= sqrt(kf*kf + Mef2);   
           Ef2=kf*kf + Mef2;
           ro=(kf*kf*kf)/(3*Pi*Pi);
        
      znak = rodzaj;
      */
  
	rozdzielnik(m_mu+50*MeV, rodzaj_jadra, rodzaj, wlacznik, FF_form,ped);

	stringstream s;

	//char* s2;
	string s2;
 
	if(znak==1) 
    	s2 ="_n";
    else 
     	s2 ="_an";
     
	s<<"rpa_cal_"<<nazwa_jadra(rodzaj_jadra)<<"_g="<<g_prim<<"kfsr="<<kf<<"Mef="
	 <<Mef<<s2<<form(FF)<<".dat"<<flush;
    
    ofstream out(s.str().c_str());
    
    for(double En =  m_mu+50*MeV; En<zakres ;En+=krok)
    	out<< En/GeV  <<' '<</*0.5*liczba_atomowa(rodzaj_jadra)*/calkowity(En,dokl)/(1e-38)<<endl;
}

double ratio_rpa_fg(double E, double q0, double qv, double kf)
{

    int wlacznik=0;
    int rodzaj_jadra=0;
    int FF_form=0;
    int rodzaj=-1;
    rozdzielnik(E, rodzaj_jadra, rodzaj, wlacznik, FF_form, kf);

	RPA_OR_NOT = true;
	double a=Przekroj(qv,q0);
	
	RPA_OR_NOT = false;
	double b=Przekroj(qv,q0);
        
    //cout<<"a="<<a<<"\nb="<<b<<endl;
	return a/b;
}


}


#ifdef RPA_MAIN

using namespace rpa;

int main()
{
	cout<<"as"<<endl;
	//cout<<" Tlen RPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 0, O,Dipol,0)<<endl;
	//cout<<" Tlen MFRPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 1, O,Dipol,0)<<endl;
	//cout<<" Argon RPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 0, Ar,Dipol,0)<<endl;
	//cout<<" Argon MFRPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 1, Ar,Dipol,0)<<endl;
	//cout<<" Iron RPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 0, Fe,Dipol,0)<<endl;
	//cout<<" Iron MFRPA "<<calkowity_liczba(1*GeV, 1e-39, neutrino, 1, Fe,Dipol,0)<<endl;

	//cout << "Mef_sr Carbon =" << sredni_Mef(C) <<endl;

	//cout << Ar << "\t" << O << "\t" << endl;

	//(double E,double dokl, double krok,int rodzaj_jadra, int rodzaj, int wlacznik, int FF_form,double ped=0 )
	// plot_przekroj_density_q0(0.7*GeV,1e-3, 1*MeV,  Ar, neutrino, 0, BBA,225*MeV );


	//plot_przekroj_q2(0.5*GeV,1e-3, 500,Ar ,neutrino, 1,NNFF,225*MeV );
	//plot_przekroj_q2(0.5*GeV,1e-3, 500,Ar ,neutrino, 0,NNFF,225*MeV );

	//plot_przekroj_density_q0(0.5*GeV,1e-3, 1*MeV,Ar ,neutrino, 1,NNFF,225*MeV );
	//plot_przekroj_density_q0(0.5*GeV,1e-3, 1*MeV,Ar ,neutrino, 0,NNFF,225*MeV );



	for(double E=0.2*GeV; E < 2*GeV; E+=20*MeV)
	{  
	
	  for(double q0=55*MeV;q0<E;q0+=10*MeV)
	  {

	   double qv=2*q0;
	   double kf=20*MeV;
	   double x=ratio_rpa_fg(E, q0, qv, kf);
	   if(x==x)
  	     cout<< ratio_rpa_fg(E, q0, qv, kf)<<endl;
  	  }
//	plot_przekroj_q2alt(E, 1e-3, 500 , ust1c);
		//plot_przekroj_q2alt(E, 1e-3, 500 , ust2c);

		/*
		plot_przekroj_q2alt(E, 1e-3, 500 , ust1a);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust1b);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust2a);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust2b);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust10a);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust10b);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust20a);
		plot_przekroj_q2alt(E, 1e-3, 500 , ust20b);
		*/
	}

	// plot_przekroj_density_q0(1*GeV,1e-38, 1*MeV,C,neutrino, 1,BBA,225*MeV );
	//plot_przekroj_density_q0(1*GeV,1e-38, 10*MeV,Ar,neutrino, 1,0 );
	//plot_przekroj_density_q0(1*GeV,1e-38, 10*MeV,Fe,neutrino, 1,0 );
	//plot_przekroj_q2(1*GeV,1e-39, C,neutrino,0, 225*MeV ); 
	//plot_przekroj_q2(0.7*GeV,1e-39, C,neutrino,1, 225*MeV );
	//plot_przekroj_q2(1*GeV,1e-39, C,neutrino,0, 225*MeV ); 
	//plot_przekroj_q2(1*GeV,1e-39, C,neutrino,1, 225*MeV );
	//plot_przekroj_q2(1.2*GeV,1e-39, C,neutrino,0, 225*MeV );
	//plot_przekroj_q2(0.7*GeV,1e-3, C,neutrino,1, 225*MeV );
	//plot_przekroj_q2(1.5*GeV,1e-39, C,neutrino,0, 225*MeV );
	//plot_przekroj_q2(1.5*GeV,1e-39, C,neutrino,1, 225*MeV );
	//plot_calkowity_density(10*GeV,50*MeV, C,neutrino,1e-39,0,BBA,225*MeV);
	//plot_calkowity_density(10*GeV,50*MeV, O,neutrino,1e-38,1,0);
	//plot_calkowity_density(10*GeV,50*MeV, Fe,neutrino,1e-38,1,0);
	//plot_calkowity_density(10*GeV,50*MeV, C,neutrino,1e-40,1,225);

	//plot_calkowity_density(5*GeV,100*MeV, 1e-4,ust1a);
	//plot_calkowity_density(5*GeV,100*MeV, 1e-4,ust1b);

	//plot_calkowity_density(5*GeV,100*MeV, 1e-4,ust2a);
	//plot_calkowity_density(5*GeV,100*MeV, 1e-4,ust2b);
		 
}


#endif
