// Forma Faktory wziete z pracy: A. Bodek, H. Budd, and J. Arrington hep-ex/0308005
// Zwane pozniej BBA - 2003 
#ifndef _form_faktory_rozdzielnik_h_
#define _form_faktory_rozdzielnik_h_
#include "jednostki.h"
#include "form_faktory_BBA_2003_long.h"
#include "form_faktory_Dipol_long.h"
#include "form_faktory_bound_long.h"
#include <math.h>

namespace rpa
{
const int BBA   = 0;
const int Dipol = 1;
const int Bound = 2;

char* form(int switch_FF) 
{
	switch(switch_FF)
	{
		 case BBA     : return "BBA"; break;
		 case Dipol   : return "Dipol"; break;
		 case Bound   : return   "Bound"; break;
		 default        : return "Zero"; break; 
	}
}

//char* FF_char;
long double F_1(long double q0, long double qv, int switch_FF)
{
	switch(switch_FF)
	{
		 case BBA     : return F_1_BBA_2003(qv*qv-q0*q0); break;
		 case Dipol   : return F_1_Dipol(qv*qv-q0*q0); break;
		 case Bound : return F_1_bound(qv*qv-q0*q0); break;
		 default        : return 0; break; 
	}
	//if(switch_FF==0) { return  F_1_BBA_2003(qv*qv-q0*q0); }
	//else {  return F_1_Dipol(qv*qv-q0*q0);}
}

long double F_2(long double q0, long double qv, int switch_FF)
{
	switch(switch_FF)
	{
		 case BBA         : return F_2_BBA_2003(qv*qv-q0*q0); break;
		 case Dipol       : return F_2_Dipol(qv*qv-q0*q0); break;
		 case Bound     : return F_2_bound(qv*qv-q0*q0); break;
		 default            : return 0; break;
	 }
	//if(switch_FF==0) { return F_2_BBA_2003(qv*qv-q0*q0);}
	//else { return F_2_Dipol(qv*qv-q0*q0);}
}

long double F_p(long double q0, long double qv, int switch_FF)
{
	 switch(switch_FF)
	{
		 case BBA     : return F_p_BBA_2003(qv*qv-q0*q0); break;
		 case Dipol   : return F_p_Dipol(qv*qv-q0*q0); break;
		 case Bound : return F_p_bound(qv*qv-q0*q0); break;
		 default        :return 0; break;
	}
	//if(switch_FF==0) { return F_p_BBA_2003(qv*qv-q0*q0);}
	//else {  return F_p_Dipol(qv*qv-q0*q0);}
}

long double G_A(long double q0, long double qv, int switch_FF)
{
	switch(switch_FF)
	{
		 case BBA     : return G_A_BBA_2003(qv*qv-q0*q0); break;
		 case Dipol   : return G_A_Dipol(qv*qv-q0*q0); break;
		 case Bound : return G_A_bound(qv*qv-q0*q0); break;
		 default        : return 0; break;
	}
	//if(switch_FF==0) { return G_A_BBA_2003(qv*qv-q0*q0);}
	//else {  return G_A_Dipol(qv*qv-q0*q0);}
}
}
#endif














