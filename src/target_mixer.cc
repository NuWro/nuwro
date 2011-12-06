#include "params.h"
#include "jednostki.h"

////////////////////////////////////////////////
double get_Eb(int p, int n)
{
	switch(p*1000+n)
	{
		case 1000: return 2;
		case 1001: return 2;
		
		default: return 8*MeV;
	}
	
}
/////////////////////////////////////////
double get_kf(int p, int n)
{
	switch(p*1000+n)
	{
		case 1000: return 0;
		case 1001: return 0;
		
		default: return 222*MeV;
	}
	
}

/*
# Models for the description of nucleus as a target
# 0 is free target; 
# 1 is Fermi gas; 
# 2 is local Fermi gas; 
# 3 is Bodek-Ritchie; 
# 4 is "effective" spectral function (carbon or oxygen); 
# 5 is deuterium; 
# 6 is deuterium with constant binding energy nucleus_E_b (for tests only!)
*/
///////////////////////////////////////
double get_model(int p, int n)
{
	switch(p*1000+n)
	{
		case 1000: case 1: return 0;
		case 1001: return 5;
		
		default: return 2;
	}	
}

