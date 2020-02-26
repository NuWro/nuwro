#ifndef _FormFaktoryHorowitz_h_
#define _FormFaktoryHorowitz_h_
#include <iostream>
#include <cmath>
#include "jednostki.h"

namespace rpa
{
const char FormFaktor = 'H';
double G_A(double q0,double qv)                                         
       { 
       double MA_H= 1.3*GeV;                                                               
       return  -1.26/(pow(1-(q0*q0-qv*qv)/(MA_H*MA_H),2));                 
       }                                                                
double lp = 1.793;

	
double G_(double q0, double qv)
       {
       double tau = (qv*qv-q0*q0)/4/M/M;
       return 1/pow(1+4.97*tau,2);
       }

double F_1p(double q0, double qv)
       {
       double tau = (qv*qv-q0*q0)/(4*M12_2);
       return (1+tau*(1+lp))*G_(q0,qv)/(1+tau);
       }

double F_2p(double q0, double qv)
       {
       double tau = (qv*qv-q0*q0)/(4*M12_2);
       return lp*G_(q0,qv)/(1+tau);
       }

double ln= -1.913;
       
double F_1n(double q0, double qv)
       {
       double tau = (qv*qv-q0*q0)/(4*M12_2);
       double ni  = 1/(1+5.6*tau);
       return tau*ln*(1-ni)*G_(q0,qv)/(1 + tau);
       }
double F_2n(double q0, double qv)
       {
       double tau = (qv*qv-q0*q0)/(4*M12_2);
       double ni  = 1/(1+5.6*tau);
       return ln*(1 + tau*ni)*G_(q0,qv)/(1 + tau);
       }
double F_1(double q0, double qv)
        {
	return   F_1p(q0,qv) - F_1n(q0,qv) ;
	}

double F_2(double q0, double qv)
        {
	return   F_2p(q0,qv) - F_2n(q0,qv);
	}
}
#endif
