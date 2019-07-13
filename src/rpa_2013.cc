#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>

#include "rpa/jednostki.h"
#include "rpa/form_faktory_rozdzielnik.h"
#include "calg5.h"
#include "util2.h"
#include "ff.h"

using namespace std;

typedef complex<double> zesp;

static double  M=M12;
static double  m_zapach= m_mu;
static double  mm2=m_zapach*m_zapach;
static double  En;
//static double  M12_2=M12*M12;
static double  m_pi2= m_pi*m_pi;
static double  m_r2 =m_r*m_r;
static double  kf=225*MeV;
static double  gr = sqrt(1.64*4*Pi);
static double  gr2=gr*gr;
static double  fr = 6.1*gr;
static double  fr2 =fr*fr;
static double  fpi = sqrt(4*Pi*0.075);
static double  fpi2= fpi*fpi;
static double  g_prim = 0.7;
static int znak = 1;
//~ static int neutrino=1;
//~ static int antyneutrino=-1;

static double  Mef = M12;				 //638 * MeV;
static double  Mef2=Mef*Mef;
static double  Ef= sqrt(kf*kf + Mef2);
static double  Ef2=kf*kf + Mef2;
static double  ro=(kf*kf*kf)/(3*Pi*Pi);
static double  sqrt_2=sqrt(2);

//~ static int FF= Dipol;
static int FF= BBA;
static int rpa_switch;
static int RPA=1;
static int FM=0;
static int l(1);
static int t(1);
static int va(1);
static int a(1);
static double stala = 1.14e-11;


bool nuwro_ff=true;

//--------------------------------------------------------------------------
// Wprowadzamy funkcje okreslajace czesci rzeczywiste
//--------------------------------------------------------------------------

static int lambda_l =1;

////// ----------- Pomocnicze  Czesci rzeczywiste ----------////////

static int fey =0;						 //Czesci Faynmanowskie sa rowne zero dla fey = 0;

// -------------------------------------------------------------------------



//      -----------------------------------------------------------------------------

double qMin(double q0)
{
	if( ((En-q0)*(En-q0)-mm2) < 0 ) return 0;
	else
		return sqrt(2*En*En- 2*En*q0 + q0*q0- mm2- 2*En*sqrt((En-q0)*(En-q0) - mm2 ));
}


double qMax(double q0)
{
	if((En-q0)*(En-q0)-mm2 < 0 ) return 0;
	else
		return sqrt(2*En*En - 2*En*q0 + q0*q0 - mm2 + 2*En*sqrt((En-q0)*(En-q0)-mm2 ));
}


double sigma_qv_q0_2013(double qv,double q0)
{
	double q2=q0*q0-qv*qv;
	double qv2= qv*qv;
	double q02= q0*q0;

	/// Czynniki postaci
	double F1,F2,GA,Fp;
//	{  	// Form Factor debug
//		double F1n,F2n,GAn,Fpn;
//		list(F1n,F2n)=f12(q2,0);
//		list(GAn,Fpn)=fap(q2,0);
//		Fpn/=M12;
//		F1 = F_1(q0,qv,FF);
//		F2 = F_2(q0,qv,FF);
//		GA = G_A(q0,qv,FF);
//		Fp = F_p(q0,qv,FF);
//		printf(" %15.8f %15.8f %15.8f %15.8f \n",F1,F2,GA,Fp);
//		printf(" %15.8f %15.8f %15.8f %15.8f \n\n",F1n,F2n,GAn,Fpn);
//	}		
	if(nuwro_ff)
	{  
		list(F1,F2)=f12(q2,0);
		list(GA,Fp)=fap(q2,0);
		Fp/=M12;
	}	
	else
	{
		F1 = F_1(q0,qv,FF);
		F2 = F_2(q0,qv,FF);
		GA = G_A(q0,qv,FF);
		Fp = F_p(q0,qv,FF);
	}
	/// tensor leptonowy

	double L_L_ = -(16*En*(En-q0)-4*(mm2-(q0*q0-qv*qv) ) )*(q0*q0-qv*qv)/(qv*qv)
		-4*mm2*q0*(4*En - q0 + q0*mm2/(q0*q0-qv*qv) )/(qv*qv);

	double L_T_ = -(16*En*(En-q0)-4*(mm2-(q0*q0-qv*qv)))*(q0*q0-qv*qv)/(qv*qv)
		- 4*mm2*(4*En*q0- (q0*q0-qv*qv) + mm2)/(qv*qv) -8*((q0*q0-qv*qv)-mm2);

	double L_A_   = 8*((q0*q0-qv*qv)-mm2);

	double L_VA_ = -16*((q0*q0-qv*qv)*(2*En-q0) + q0*mm2)/qv;

	//// tensory

	//--------------------------------------------------------------------------
	// Wprowadzamy funkcje okreslajace czesci rzeczywiste
	//--------------------------------------------------------------------------

	double licznikalfa = fabs(q2*q2 - 4*pow( q0*Ef - qv*kf, 2 ) );
	double mianownikalfa =fabs(q2*q2 - 4*pow(q0*Ef + qv*kf,  2));
	double alfa_= log(licznikalfa/mianownikalfa );

	double licznikbeta =fabs((pow(q2 + 2*qv*kf, 2) - 4*q02*Ef2));
	double mianownikbeta =fabs(pow(q2 - 2*qv*kf, 2) - 4*q02*Ef2 );
	double beta_ = log(licznikbeta/mianownikbeta ) ;

	double stalLAMBDA=sqrt(q2*(q2 -4*Mef2));
	double licznik1LAMBDA =      fabs(q2*q2*Ef - 4*q0*Mef2*(q0*Ef-qv*kf) - kf*q2*stalLAMBDA);
	double Mianownik1LAMBDA=fabs(q2*q2*Ef - 4*q0*Mef2*(q0*Ef-qv*kf) + kf*q2*stalLAMBDA);
	double licznik2LAMBDA =      fabs(q2*q2*Ef - 4*q0*Mef2*(q0*Ef+qv*kf) - kf*q2*stalLAMBDA);
	double Mianownik2LAMBDA=fabs(q2*q2*Ef - 4*q0*Mef2*(q0*Ef+qv*kf) + kf*q2*stalLAMBDA);

	double LAMBDA_=0;

	if(q2==0) { LAMBDA_ = 0;}
	else
	{
		LAMBDA_= (log(licznik1LAMBDA/Mianownik1LAMBDA )
			+         log(licznik2LAMBDA/Mianownik2LAMBDA)  )/2./stalLAMBDA;
	}

	////// ----------- Pomocnicze  Czesci rzeczywiste ----------////////

	double nawiasReHs = kf*Ef - (3*Mef2 - q2/2)*log((kf+Ef)/Mef) -
		q0*(4*Mef2-q2)*alfa_ /8/qv +
		Ef*(4*Mef2-q2)*beta_  /4/qv +
		(4*Mef2-q2)*(4*Mef2-q2)*LAMBDA_/4;

	double ReHs_= lambda_l*nawiasReHs /(2*Pi*Pi);

	double nawiasReHl= 2*kf*Ef/3 - (qv2/6)*log((kf+Ef)/Mef) +
		(q0/(4*qv))*(Ef2 + (q02-3*qv2)/12 )*alfa_
		- Ef*(3*q2 + 4*Ef2)*beta_/(24*qv)
		+ qv2*(2*Mef2+q2)*(4*Mef2-q2)*LAMBDA_/(12*q2) ;

	double ReHl_= nawiasReHl*lambda_l*q2/(Pi*Pi*qv2);

	double nawiasReHt= kf*Ef*(1+ 2*q02/qv2 )/3 + (q2/3)*log((kf+Ef)/Mef)
		+ (q0*q2/(4*qv2*qv))*((q02 + 3*qv2)/12 + qv2*Mef2/q2+ Ef2)*alfa_
		+ (Ef/qv)*((qv2*qv2-q02*q02)/(8*qv2) - Mef2/2 - q2*Ef2/(6*qv2))*beta_
		-  (2*Mef2+q2)*(4*Mef2 - q2)*LAMBDA_/6;

								 // tutaj przemnozylem przez 2
	double  ReHt_= 2*nawiasReHt*lambda_l/(2*Pi*Pi);

	double nawiasReH0 = qv*kf + q0*Ef*alfa_/2
		-(q2/8 + Ef2/2 + qv2*Mef2/(2*q2))*beta_;

	double ReH0_= nawiasReH0 * lambda_l*Mef/(2*Pi*Pi*qv);

	// -------------------------------------------------------------------------

	double E_max_ = max( Mef , max(Ef - q0, 0.5*(- q0+qv * sqrt(1-4*Mef2/q2 ))));

	double E__ = min(Ef, E_max_);

	double E1_ = Ef - E__;

	double E2_ = (Ef*Ef -E__*E__)/2;

	double E3_=  (Ef*Ef*Ef -E__*E__*E__)/3;

	//              Czesci Rzeczywiste tensorw Hadronowych

	//                        ------------------

	double Re_H_a_ =  ( - 2*ReHs_ - ReHl_ + ReHt_)/3;

	double Re_H_vv_l_= ReHl_;

	double Re_H_vv_t_= ReHt_;

	double Re_H_va_= (q0*q0-qv*qv)*ReH0_/2/qv/qv/Mef;

	double Re_H_tt_l_=  (q0*q0-qv*qv)*(-ReHl_ + ReHs_ - 2*ReHt_)/ (12 * M*M)  ;

	double nawias1_Re_H_tt_t = ReHl_*(q2 - 8*Mef2) + ReHt_*(2*Mef2 - q2)
		+ 2*ReHs_*(q2-2*Mef2);
	double nawias2_Re_H_tt_t = q2*ReH0_ - 2*qv*qv*Mef*Re_H_va_;
	double Re_H_tt_t_ =  q2*(nawias1_Re_H_tt_t/24/M/M/Mef2 + nawias2_Re_H_tt_t/4/q0/M/M/Mef);

	double Re_H_vt_l_ =-q2*Re_H_a_/(4*M*Mef);

	double Re_H_vt_t_= (q0*q0-qv*qv)*Re_H_a_/(2*M*Mef) ;

	//~ double Re_H_ppf_l_= ReHs_ +0.5*Re_H_a_;
//~ 
	//~ double Re_H_ppf_t_= ReHs_ +0.5*Re_H_a_;
//~ 
	//~ double Re_H_ppf_a_= ReHs_ +0.5*Re_H_a_;
//~ 
	//~ double Re_H_apf_l_=  0.5*Re_H_a_/Mef;
//~ 
	//~ double Re_H_apf_t_= 0.5*Re_H_a_/Mef;
//~ 
	//~ double Re_H_apf_a_= 0.5*Re_H_a_/Mef;

	//    ----------------------------------------------------------------------

	//                   Czesci Urojone tensorow Hadronowych

	//                   -----------------------------------

	double Im_H_vv_l_= ((q0*q0-qv*qv)/(2*Pi*qv*qv*qv))*
		( E3_ + q0*E2_ + (q0*q0-qv*qv)*E1_/4 ) ;

								 // dwojka
	double Im_H_vv_t_= 2*((q0*q0-qv*qv)/(4*Pi*qv*qv*qv))*
		(E3_+q0*E2_+(qv*qv*Mef2/(q0*q0-qv*qv)+0.25*(q0*q0+qv*qv))*E1_) ;

	double Im_H_tt_l_=  - ( ( q0*q0-qv*qv )/( 8*Pi*qv*qv*qv*M*M))*
		( (q0*q0-qv*qv)*E3_ + q0*(q0*q0-qv*qv)*E2_+
		( qv*qv*Mef2 +0.25*(q0*q0-qv*qv)*q0*q0 )*E1_ );

								 // bylo 16 zamiast 8;
	double Im_H_tt_t_ =( (q0*q0-qv*qv)/(8*Pi*qv*qv*qv*M*M) )*

		( ( Mef2*qv*qv - 0.25*(q0*q0-qv*qv)*(q0*q0-qv*qv) )*E1_

		- q0*(q0*q0-qv*qv)*E2_ -(q0*q0-qv*qv)*E3_ ) ;

	double Im_H_vt_l_= -( (q0*q0-qv*qv)*E1_ *Mef )/(8*Pi*qv*M);

								 // dwojka
	double Im_H_vt_t_= -2*Im_H_vt_l_;

	double Im_H_va_ = qv*((q0*q0-qv*qv)/(8*Pi*qv*qv*qv))*(2*E2_ + q0*E1_);

	double Im_H_a_= (Mef2 * E1_ )/(2*Pi*qv);

	double Im_H_pp_ = (q0*q0-qv*qv)*E1_/8/Pi/qv ;

	double Im_H_ap_ =  Mef*E1_/4/Pi/qv;

	// --------------------------------------------------------------------------
	//           Wprowadzamy kompletne postacie Tensorw Hadronowych
	//           ----------------------------------------------------------------

	zesp H_vv_l_=  zesp(Re_H_vv_l_,Im_H_vv_l_);

	zesp H_vv_t_=  -0.5* zesp(Re_H_vv_t_,Im_H_vv_t_);

	zesp H_tt_l_=   zesp(Re_H_tt_l_,Im_H_tt_l_);

	zesp H_tt_t_ = -0.5*zesp(Re_H_tt_t_,Im_H_tt_t_);

	zesp H_vt_l_= zesp(Re_H_vt_l_,Im_H_vt_l_) ;

	zesp H_vt_t_= -0.5*zesp(Re_H_vt_t_,Im_H_vt_t_);

								 // czesc urojona pomnozylem przez
	zesp H_va_=   zesp(qv*Re_H_va_, Im_H_va_) ;
								 // qv wczaesniej w definicji Im

	zesp H_a_=  zesp(Re_H_a_,Im_H_a_);

	zesp H_apf_l_=-(q0*q0-qv*qv)*0.5*zesp(Re_H_a_,Im_H_a_)/Mef;

	zesp H_apf_t_= -(q0*q0-qv*qv)*0.5*zesp(Re_H_a_,Im_H_a_)/Mef;

	zesp H_apf_a_=(q0*q0-qv*qv)*0.5*zesp(Re_H_a_,Im_H_a_)/Mef;

	//    -----------------------------------------------------------------------

	zesp H_rr_l_ = (gr2*H_vv_l_*2. + fr*gr*H_vt_l_ + fr2*H_tt_l_/2.);

	zesp H_rr_t_ = (gr2*H_vv_t_/2. + fr*gr*H_vt_t_ + fr2*H_tt_t_/2.);

	zesp H_rp_va_ =  -fpi*(gr + fr *Mef/M)*H_va_/m_pi;

	zesp H_pp_l_ =  2.*fpi2*H_vv_l_/m_pi2;

	zesp H_pp_t_ = 2.*fpi2*H_vv_t_/m_pi2;

	zesp H_pp_a_ = 2.*fpi2*H_a_/m_pi2;

	//-----------------------------// PROPAGATORY SWOBODNE  //-----------------------------------------------//

	double V_l_ = - q2/(q2 - m_pi2);;

								 // zmiana znaku dla elemenu macierzowego
	double V_t_ = -q2/(q2 - m_pi2);

	double V_a_ = q2/(q2 - m_pi2) - g_prim;

	double R_l_ = -q2/(q2 - m_r2)/m_r2;

								 // zmiana znaku dla elemenu macierzowego
	double R_t_ = -q2/(q2 - m_r2)/m_r2;

	double R_a_ = 1/m_r2;

	//--------------------------------------// Tensory Pi  i  Rho  //-------------------------------------------//

	zesp H_r_l_= gr*(F1*H_vv_l_ + F2*H_vt_l_)/sqrt(2)+
		fr*(F1*H_vt_l_ + F2*H_tt_l_)/sqrt(2);

	zesp H_r_t_ =  gr*(F1*H_vv_t_ + F2*H_vt_t_)/sqrt(2)
		+ fr*(F1*H_vt_t_ + F2*H_tt_t_)/sqrt(2);

	zesp H_r_va_ = GA*(gr+fr*Mef/M)*H_va_/sqrt(2);

	zesp H_p_l_ = -sqrt(2)*fpi*(GA*H_vv_l_+Fp*H_apf_l_)/m_pi;

	zesp H_p_t_ =  -sqrt(2)*fpi*( GA*H_vv_t_+Fp*H_apf_t_ )/m_pi;

	zesp H_p_a_ = -sqrt(2)*fpi*( H_a_*GA+Fp*H_apf_a_ )/m_pi;

	zesp H_p_va_ = -sqrt(2)*fpi*(F1+F2*Mef/M)*H_va_/m_pi;

	////   delty

	zesp delta_a_1_ = R_a_;

	zesp delta_a_4_ = V_a_/(1.-V_a_*H_pp_a_);

	zesp delta_l_1_ = (R_l_+(R_l_+R_a_)*H_rr_l_*R_a_)/
		(1. - (R_l_+R_a_)*H_rr_l_ );

	zesp licznik_delta_l_4 =V_l_+( V_l_+V_a_ )*H_pp_l_*V_a_;
	zesp mianownik_delta_l_4 = (1. - V_a_*H_pp_a_ )*(1. - (V_l_+V_a_)*(H_pp_l_+H_pp_a_));
	zesp delta_l_4_ = licznik_delta_l_4/mianownik_delta_l_4;

	zesp R_ta =R_t_+R_a_;		 ///
	zesp V_ta =V_t_+V_a_;		 ///

	zesp licznik_delta_t_1 =(1.- V_ta*(H_pp_a_+H_pp_t_)) *( R_t_ + R_a_* R_ta*H_rr_t_  )
		+ R_a_*R_ta*V_ta*pow(H_rp_va_,2);

	zesp mianownik_delta_t_1 = (1. - V_ta*(H_pp_a_ + H_pp_t_ )) *(1. - R_ta*H_rr_t_)
		- R_ta*V_ta*pow(H_rp_va_,2);

	zesp delta_t_1_ = licznik_delta_t_1/mianownik_delta_t_1;

	zesp licznik_delta_t_4 = (1. - R_ta*H_rr_t_ )*
		( V_t_ +delta_a_4_*(V_ta*H_pp_t_ + V_t_*H_pp_a_))
		+ delta_a_4_*R_ta*V_ta*pow(H_rp_va_,2);

	zesp mianownik_delta_t_4= (1. - V_ta*(H_pp_a_ + H_pp_t_ ) )*(1. - R_ta*H_rr_t_)
		- R_ta*V_ta*pow(H_rp_va_,2);

	zesp delta_t_4_ = licznik_delta_t_4/mianownik_delta_t_4;

	zesp delta_va_2_ = R_ta*H_rp_va_*(delta_t_4_+delta_a_4_)/
		(1. - R_ta*H_rr_t_ );

	zesp delta_va_3_ =  V_ta*H_rp_va_*( R_a_ + delta_t_1_)/(1. - V_ta*(H_pp_a_ + H_pp_t_) );

	//// Poprawki RPA
	double poprawka_l_ =  imag( ( delta_a_4_+delta_l_4_)*(H_p_l_*H_p_l_  + 2.*H_p_l_*H_p_a_ ))
		+imag( H_p_a_*H_p_a_*delta_l_4_ + H_r_l_*H_r_l_*( delta_a_1_+delta_l_1_));

	double poprawka_t_ = imag((pow(H_p_t_,2) + 2.*H_p_t_*H_p_a_ +pow(H_p_va_,2))*
		(delta_a_4_ +delta_t_4_) + delta_t_4_*pow(H_p_a_,2))
		+
		imag((pow(H_r_t_,2) + pow(H_r_va_,2))*(delta_t_1_+delta_a_1_))
		+
		imag( (H_r_va_*H_p_t_+H_r_va_*H_p_a_+H_r_t_*H_p_va_)*
		(delta_va_2_+delta_va_3_));

	double poprawka_a_ = imag(delta_a_4_*H_p_a_*H_p_a_);

	double poprawka_va_ = imag((H_r_t_*H_p_t_ + H_r_t_*H_p_a_+H_r_va_*H_p_va_)*
		(delta_va_2_ + delta_va_3_)) +
		2.*imag( H_p_va_*( H_p_t_+H_p_a_ )*( delta_t_4_ + delta_a_4_)
		+H_r_va_*H_r_t_*(delta_t_1_ + delta_a_1_));

	//poprawka_va = imag((H_r_t*H_p_t + H_r_t*H_p_a+H_r_va*H_p_va)*(delta_va_2 + delta_va_3))
	//                        +2.*imag( H_p_va*( H_p_t+H_p_a )*( delta_t_4 + delta_a_4)
	//		              +H_r_va*H_r_t*(delta_t_1 + delta_a_1));
	//----------------------------------------------------------------------------

	/// Tensory odpowiedzi
	double R_L_[2];
	double R_T_[2];
	double R_A_[2];
	double R_VA_[2];
    
	R_L_[0] = ( F1*F1 + GA*GA )*Im_H_vv_l_+
		F2*F2*Im_H_tt_l_+
		2*F1*F2*Im_H_vt_l_
		-q2*Fp*Fp*Im_H_pp_
		-2*q2*GA*Fp*Im_H_ap_;
	R_L_[1] = R_L_[0] + poprawka_l_;

	R_T_[0] = Im_H_vv_t_*(GA*GA+ F1*F1)/2
		+ Im_H_tt_t_*F2*F2/2
		+ F1*F2*Im_H_vt_t_
		+ q2*Fp*Fp*Im_H_pp_
		+ 2*q2*GA*Fp*Im_H_ap_;
	R_T_[1] = R_T_[0] - poprawka_t_;

	R_A_[0] = pow(GA,2)*Im_H_a_ + q2*pow(Fp,2)*Im_H_pp_
		+ 2*q2*GA*Fp*Im_H_ap_;
	R_A_[1] = R_A_[0] +  poprawka_a_;

	R_VA_[0] = ( F1 +  Mef*F2/M )*GA*Im_H_va_;
								 //wydaje sie ze musi tu byc 0.5 bo tensor leptonowy jest tak unormowany
	R_VA_[1] = R_VA_[0] + 0.5*poprawka_va_;

	double Ampl[2];
	Ampl[0] = L_L_*R_L_[0] + L_T_*R_T_[0] + L_A_*R_A_[0] - znak*L_VA_*R_VA_[0];
	Ampl[1] = L_L_*R_L_[1] + L_T_*R_T_[1] + L_A_*R_A_[1] - znak*L_VA_*R_VA_[1];

	//~ if(Amplituda > 0)
		//~ cerr<<"amplituda zla( q0="
			//~ << q0/GeV<<" qv="<<qv/GeV<<','<<-Amplituda*(stala*stala * qv )/(16 * Pi * Pi * ro  * En * En)*GeV/cm2<<')'<<endl;
	switch(rpa_switch)
	{
		case 0: return max(-Ampl[0]*(stala*stala * qv )/(16 * Pi * Pi * ro  * En * En)*GeV/cm2, 0.);
		case 1: return max(-Ampl[1]*(stala*stala * qv )/(16 * Pi * Pi * ro  * En * En)*GeV/cm2, 0.);
		case 2: if(Ampl[0]>=0) return 1; 
			if(Ampl[1]>=0) return 0;
			return  min(Ampl[1]/Ampl[0],10.0);
		default: return 1;
	}
}


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


	double E_=min(Ef, max( Mef , max(Ef - q0, 0.5*(- q0+qv * sqrt(1-4*Mef2/q2 )))));                     
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

	double F_1_, F_2_, G_A_,Fp;

	if(nuwro_ff)
	{  
		list(F_1_,F_2_)=f12(q2,0);
		list(G_A_,Fp)=fap(q2,0);
		Fp/=M12;
	}	
	else
	{
		F_1_ = F_1(q0,qv,FF);
		F_2_ = F_2(q0,qv,FF);
		G_A_ = G_A(q0,qv,FF);
		Fp = F_p(q0,qv,FF);
	}

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
	
	R_VA[1]=R_VA[0] +poprawka_va/2.0;   /// UWAGA TUTAJ BYL BLAD, TRZEBA PODZIELIC PRZEZ 2!

	
	double L_L=-(16*En*(En-q0)-4*(mm2-q2))*q2/qv2 - 4*mm2*q0*(4*En - q0 + q0*mm2/q2 )/qv2;

	double L_T=-(16*En*(En-q0)-4*(mm2-q2))*q2/qv2 
				- 4*mm2*(4*En*q0- (q0*q0-qv*qv) + mm2)/qv2 - 8*(q2-mm2);

	double L_A=8*(q2-mm2);

	double L_VA=-16*(q2*(2*En-q0) + q0*mm2)/qv;


	double Ampl[2];
	Ampl[0]=    L_L*R_L[0] +(L_A+L_T)*R_A[0] + L_T*(R_T[0] -R_A[0]) -znak*L_VA*R_VA[0];   
	Ampl[1]=    L_L*R_L[1] +(L_A+L_T)*R_A[1] + L_T*(R_T[1] -R_A[1]) -znak*L_VA*R_VA[1];   
  
/*	  if(amplituda[0] > 0 || amplituda[1] > 0)
	  cerr<<"amplituda zla( E="<< En/GeV 
		  <<", q0 = "<< q0/GeV<<", Q2 = " 
		  << (qv*qv-q0*q0)/GeV/GeV 
		  << ", (En-q0)*(En-q0)-mm2 = " << 
	  ((En-q0)*(En-q0)-mm2)/GeV/GeV<<") mm2="<<mm2<<endl;
*/
	switch(rpa_switch)
	{
		case 0: return max(-Ampl[0]*(stala*stala * qv )/(16 * Pi * Pi * ro  * En * En)*GeV/cm2, 0.);
		case 1: return max(-Ampl[1]*(stala*stala * qv )/(16 * Pi * Pi * ro  * En * En)*GeV/cm2, 0.);
		case 2: if(Ampl[0]>=0) return 1; 
			if(Ampl[1]>=0) return 0;
			return  min(Ampl[1]/Ampl[0],10.0);
		default: return 1;
	}
}

double przekroj(double q0)
{

	if( (En-q0)*(En-q0)-mm2 < 0 )
		{ return 0; cerr<<"amplituda nie liczona"<<endl;}

		else   return  calg5a(fix2(sigma_qv_q0_2013,q0),qMin(q0),qMax(q0),1000);
}

double sigma_q0( double q0)
{ 

	if( (En-q0)*(En-q0)-mm2 < 0 )
	{ 
		return 0; cerr<<"amplituda nie liczona"<<endl;
	}
	else   
		return  calg5a(fix2(sigma_qv_q0,q0),qMin(q0),qMax(q0),1000);
}    


void plot_przekroj_q0(double E, int ilosc_punktow, int ZNAK, bool nowy)
{

	znak=ZNAK;
	En =E;

	
	int tmp=rpa_switch;
	
	stringstream s[3];
	s[0]<<"rpa=0,E="<<En/MeV<<",g="<<g_prim<<",kf="<<kf<<",Mef="<<Mef<<",znak="<<znak<<".v="<<nowy<<".dat"<<flush;
	s[1]<<"rpa=1,E="<<En/MeV<<",g="<<g_prim<<",kf="<<kf<<",Mef="<<Mef<<",znak="<<znak<<".v="<<nowy<<".dat"<<flush;
	s[2]<<"rpa=r,E="<<En/MeV<<",g="<<g_prim<<",kf="<<kf<<",Mef="<<Mef<<",znak="<<znak<<".v="<<nowy<<".dat"<<flush;
	cout<<s[0].str()<<endl;
	ofstream out0(s[0].str().c_str());
	ofstream out1(s[1].str().c_str());
	ofstream out2(s[2].str().c_str());

	double cross;
	double crossFG;
	
	double(*f)(double)= nowy ? przekroj : sigma_q0;

	double wsp=nowy ?1 :500;
	
	for(double q0=0.001; q0 < 0.06*GeV; q0+=0.6*MeV)
	{
		
		rpa_switch=0;
		cross   = f(q0)/(1e-38);
		rpa_switch=1;
		crossFG = f(q0)/(1e-38);

		out0<< q0/GeV <<"\t" << crossFG*wsp <<endl;
		out1<< q0/GeV <<"\t" << cross*wsp <<endl;
		out2<< q0/GeV <<"\t" << cross/crossFG<<endl;
	}

	double krok=fabs(En - m_zapach)/ilosc_punktow;
	for(double q0=0.06*GeV; q0 < En-m_zapach+krok; q0+=10*MeV)
	{

		rpa_switch=true;
		cross   = f(q0)/(1e-38);
		rpa_switch=false;
		crossFG = f(q0)/(1e-38);

		out0<< q0/GeV <<"\t" << crossFG*wsp <<endl;
		out1<< q0/GeV <<"\t" << cross*wsp <<endl;
		out2<< q0/GeV <<"\t" << cross/crossFG<<endl;

	}
	rpa_switch=tmp;
}

double ratio_rpa(double qv,double q0,double E, int pdg,double m_lepton, double m_eff, double kf0,  bool nowy=true)
{	
	znak=(pdg>0)-(pdg<0);
	En =E;
	m_zapach= m_lepton;
	mm2=m_zapach*m_zapach;
	Mef = m_eff;
	Mef2=Mef*Mef;
	Ef= sqrt(kf*kf + Mef2);
	Ef2=kf*kf + Mef2;
	ro=(kf*kf*kf)/(3*Pi*Pi);
	
	double(*f)(double,double)= nowy ? sigma_qv_q0_2013 : sigma_qv_q0;

	
/*	{ // DEBUG RPA
		rpa_switch=0; double crossFG = f(qv,q0)/(1e-38);
		rpa_switch=1; double crossRPA= f(qv,q0)/(1e-38);
		rpa_switch=2; double ratio   = f(qv,q0);
		cout<<ratio<< ' '<<crossRPA/crossFG<<"    "<<crossFG<<"  "<<crossRPA<<endl;
		if(crossFG<=0)	
			return 1;
		if(crossRPA<=0)
			return 0;
		return min(10.,crossRPA/crossFG);	
	}
*/	
	rpa_switch=2; return f(qv,q0);
}
	

void main1()
{
	for(double E=1000*MeV;E<1100*MeV;E+=500*MeV)
	for(int znak=-1;znak<2;znak+=2)
	for(int v=0;v<2;v++)
	{	
		plot_przekroj_q0(E,  500, znak, v);
	}
}
