#include "DM.h"
#include <iostream>
//Dirac matrices in the Dirac representation!
DM Gamma(unsigned short i)
{
    switch(i)
    { 	
		case 0: //gamma0
			return DM(
				  1,0,0,0,
				  0,1,0,0,
				  0,0,-1,0,
				  0,0,0,-1);
		case 1: //gamma1
			return DM(
				  0,0,0,1,
				  0,0,1,0,
				  0,-1,0,0,
				  -1,0,0,0);
		case 2: //gamma2
			return DM(
				  0,0,0,comp(0,-1),
				  0,0,comp(0,1),0,
				  0,comp(0,1),0,0,
				  comp(0,-1),0,0,0);
		case 3: //gamma3
			return DM(
				  0,0,1,0,
				  0,0,0,-1,
				  -1,0,0,0,
				  0,1,0,0);
		case 5: //gamma5
			return DM(
				  0,0,1,0,
				  0,0,0,1,
				  1,0,0,0,
				  0,1,0,0);
		default: throw "bad argument to Gamma";
	}
}


//Hermitian adjounsigned short
void DM::hermit()
{
	transp();
	conj();
}

//complex conjuation
void DM::conj()
{
	for(comp *a=&a0;a<=&d3;a++)
		*a=std::conj(*a);
}

//transposition
void DM::transp()
{
	swap(a1,b0);swap(a2,c0);swap(a3,d0);
	            swap(b2,c1);swap(b3,d1);
	                        swap(c3,d2);
}

