#ifndef LEPTONIC_H
#define LEPTONIC_H
/*****************************************************************************

*                            THE LEPTON TENSOR                               *

*****************************************************************************/

/*****************************************************************************
The function "Leptonic_Tensor" builds the leptonic tensor for neutrino reac-
tions:

L_[r,s] = (1/(E_nu*Ml))*(knu_[r]*kl_[s] + knu_[s]*kl_[r] - knu dot kl g_[r,s] 
          - i*e_[a,r,b,s]*knu^[a]*kl^[b]).

First, the four vectors of neutrino and outgoing lepton are calculated, in con-
travariant notation (upper indices).  Next, the leptonic tensor is stored in a 
complex<double> array, where a distinction is made between the symmetric and 
asymmetric part of the tensor.  This allows also to build the electron scatte-
ring tensor:

 
*****************************************************************************/

void Leptonic_Tensor( const int &process, const double &El, const double &leptonmass, 
		      const double &radThetal, const double &w, const double &q, 
		      double Lepton_S[][4], double Lepton_A[][4] );


void Leptonic_Tensor_fourvector(const int &process, double kl[4], double kl_inc[4], 
		      double Lepton_S[][4], double Lepton_A[][4] );

#endif
