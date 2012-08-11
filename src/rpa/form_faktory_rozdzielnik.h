// Forma Faktory wziete z pracy: A. Bodek, H. Budd, and J. Arrington hep-ex/0308005 
// Zwane pozniej BBA - 2003 
#ifndef _form_faktory_rozdzielnik_h_
#define _form_faktory_rozdzielnik_h_
#include "./jednostki.h"
#include "./form_faktory_BBA_2003.h"
#include "./form_faktory_Dipol.h"
#include "./form_faktory_bound.h"
#include "./nnff.h"

#include <cmath>

using namespace std;

namespace rpa
{
	
	const int BBA   = 0;
	const int Dipol = 1;
	const int Bound = 2;
	const int NNFF  = 3;
	const int NNFF2  = 4;

	const char* form(int kFF) 
	{
		cout << kFF <<endl;
		
		switch(kFF)
		{
			 case BBA     :   return "BBA";   break;
			 case Dipol   :   return "Dipol"; break;
			 case Bound   :   return "Bound"; break;
			 case NNFF    :   return "NNFF";  break;
			 case NNFF2   :   return "NNFF2";  break;
			 
			 default        : return "Zero";  break; 
		}
	}

	double F_1(double Q2, int kFF)
	{
		switch(kFF)
		{
			 case BBA     : return F_1_BBA_2003(Q2); break;
			 case Dipol   : return F_1_Dipol(Q2); break;
			 case Bound   : return F_1_bound(Q2); break;
			 case NNFF    : return F_1_nnff(Q2); break;
			 case NNFF2    : return F_1_nnff(Q2); break;
			 
			 default        : return 0; break; 
		}
	}

	double F_2(double Q2, int kFF)
	{
		switch(kFF)
		{
			 case BBA         : return F_2_BBA_2003(Q2); break;
			 case Dipol       : return F_2_Dipol(Q2);    break;
			 case Bound       : return F_2_bound(Q2);    break;
			 case NNFF        : return F_2_nnff(Q2);     break;
			 case NNFF2        : return F_2_nnff(Q2);    break;

			 default            : return 0; break;
		 }
	}

	double F_p(double Q2, int kFF)
	{
		switch(kFF)
		{
			case BBA     : return F_p_BBA_2003(Q2); break;
			case Dipol   : return F_p_Dipol(Q2);    break;
			case Bound   : return F_p_bound(Q2);    break;
			case NNFF    : return F_p_nnff(Q2);     break;
			case NNFF2   : return F_p_nnff2(Q2);    break;

			default        :return 0; break;
		}
	}

	double G_A(double Q2, int kFF)
	{
		switch(kFF)
		{
			case BBA     : return G_A_BBA_2003(Q2); break;
			case Dipol   : return G_A_Dipol(Q2);    break;
			case Bound   : return G_A_bound(Q2);    break;
			case NNFF    : return G_A_nnff(Q2);     break;
			case NNFF2    : return G_A_nnff2(Q2);     break;

			default        : return 0; break;
		}
	}
}
#endif
