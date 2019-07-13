#ifndef _czesci_rzeczywist_h
#define _czesci_rzeczywist_h
#include <iostream>
#include <cmath>
#include "jednostki.h"
#include "alfa_beta_LAM.h"

namespace rpa
{
double Re_H_a(double q0, double qv)
{
	//return -LinRelRe( q0,  qv);
	return ( - 2*ReHs(q0,qv) - ReHl(q0,qv) + ReHt(q0,qv))/3;
}

double Re_H_vv_l(double q0, double qv)
{
	return ReHl(q0,qv);
}

double Re_H_vv_t(double q0, double qv)
{
	return  ReHt(q0,qv);
}

double Re_H_va(double q0, double qv)
{
	return (q0*q0-qv*qv)*ReH0(q0, qv)/2/qv/qv/Mef;
}	     
 
double Re_H_tt_l(double q0, double qv)
{
	double nawias =  -ReHl(q0,qv) + ReHs(q0,qv) - 2*ReHt(q0,qv);
	return    (q0*q0-qv*qv)*nawias/ 12 / M/M  ;
}

double Re_H_tt_t(double q0, double qv)
{
	double q2 = q0*q0-qv*qv;
	double nawias1 = ReHl(q0,qv)*(q2 - 8*Mef2) + ReHt(q0,qv)*(2*Mef2 - q2) + 2*ReHs(q0,qv)*(q2-2*Mef2);
	double nawias2 = q2*ReH0(q0,qv) - 2*qv*qv*Mef*Re_H_va(q0,qv);  
	return  q2*(nawias1/24/M/M/Mef2 + nawias2/4/q0/M/M/Mef);
}

double Re_H_vt_l(double q0, double qv)
{
	return -(q0*q0-qv*qv)*Re_H_a(q0,qv)/(4*M*Mef);
}

double Re_H_vt_t(double q0, double qv)
{
	return (q0*q0-qv*qv)*Re_H_a(q0,qv)/2/M/Mef ;
}

}
#endif

