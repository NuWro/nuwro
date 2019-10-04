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

#include <stdlib.h>
#include <math.h>
void Leptonic_Tensor( const int &process, const double &El, const double &leptonmass, 
		      const double &radThetal, const double &w, const double &q, 
		      double Lepton_S[][4], double Lepton_A[][4] )
{

  // the magnitude of the three-momenta. for the neutrino: E_nu = k_nu and within ERL,
  // E_l = k_l

  double k_l = sqrt( pow(El,2) - pow(leptonmass,2) );
  double k_l_inc = w + El;

  // the four vectors, in CONTRAVARIANT notation

  double kl_inc[4], kl[4];

  kl_inc[0] = k_l_inc;
  kl_inc[1] = 0.;
  kl_inc[2] = 0.;
  kl_inc[3] = k_l_inc;

  kl[0] = El;
  kl[1] = sin(radThetal)*k_l;
  kl[2] = 0.;
  kl[3] = cos(radThetal)*k_l;

  double kl_inc_dot_kl = kl_inc[0]*kl[0] - kl_inc[1]*kl[1] - kl_inc[2]*kl[2] - kl_inc[3]*kl[3];


  double Norm_factor;
  
  if( process == 0 ){
    Norm_factor = 1./(2.*kl_inc[0]*kl[0]);  // normally: also divide by m_e*m_e -> see General_FF, where m_e*m_e is left out in numerator
  }else{
    Norm_factor = 2./(kl_inc[0]*kl[0]);  // see previous remark about the appearance of a factor 2.
  }

  
  Lepton_S[0][0] = Norm_factor*(2.*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  Lepton_S[0][1] = Norm_factor*(kl_inc[0]*kl[1] + -kl[0]*kl_inc[1]);  // mind the sign difference for the spatial compts!
  Lepton_S[0][2] = Norm_factor*(kl_inc[0]*kl[2] + -kl[0]*kl_inc[2]);
  Lepton_S[0][3] = Norm_factor*(kl_inc[0]*kl[3] + -kl[0]*kl_inc[3]);

  Lepton_S[1][0] = Lepton_S[0][1];
  Lepton_S[1][1] = Norm_factor*(2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  Lepton_S[1][2] = Norm_factor*(kl_inc[1]*kl[2] + kl[1]*kl_inc[2]);
  Lepton_S[1][3] = Norm_factor*(kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);

  Lepton_S[2][0] = Lepton_S[0][2];
  Lepton_S[2][1] = Lepton_S[1][2];
  Lepton_S[2][2] = Norm_factor*(2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  Lepton_S[2][3] = Norm_factor*(kl_inc[2]*kl[3] + kl[2]*kl_inc[3]);

  Lepton_S[3][0] = Lepton_S[0][3];
  Lepton_S[3][1] = Lepton_S[1][3];
  Lepton_S[3][2] = Lepton_S[2][3];
  Lepton_S[3][3] = Norm_factor*(2*kl_inc[3]*kl[3] + kl_inc_dot_kl);

  
  if( process == 0 ){
    
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
	
	Lepton_A[i][j] = 0.;
	
      }
    }
    
  }
  else{
    
    Lepton_A[0][0] = 0.;
    Lepton_A[0][1] = Norm_factor*(kl_inc[2]*kl[3] - kl[2]*kl_inc[3]);
    Lepton_A[0][2] = Norm_factor*(-kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);
    Lepton_A[0][3] = Norm_factor*(kl_inc[1]*kl[2] - kl[1]*kl_inc[2]);
    
    Lepton_A[1][0] = -Lepton_A[0][1];
    Lepton_A[1][1] = 0.;
    Lepton_A[1][2] = Norm_factor*(kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
    Lepton_A[1][3] = Norm_factor*(-kl_inc[0]*kl[2] + kl[0]*kl_inc[2]);
    
    Lepton_A[2][0] = -Lepton_A[0][2];
    Lepton_A[2][1] = -Lepton_A[1][2];
    Lepton_A[2][2] = 0.;
    Lepton_A[2][3] = Norm_factor*(kl_inc[0]*kl[1] - kl[0]*kl_inc[1]);
    
    Lepton_A[3][0] = -Lepton_A[0][3];
    Lepton_A[3][1] = -Lepton_A[1][3];
    Lepton_A[3][2] = -Lepton_A[2][3];
    Lepton_A[3][3] = 0.;

  }

  
}


void Leptonic_Tensor_fourvector(const int &process, double kl[4], double kl_inc[4], 
		      double Lepton_S[][4], double Lepton_A[][4] )
{
  //Identical to the function above, but takes the four vectors as arguments, in this way we can change reference system easily
  // the magnitude of the three-momenta. for the neutrino: E_nu = k_nu and within ERL,
  // E_l = k_l

  double kl_inc_dot_kl = kl_inc[0]*kl[0] - kl_inc[1]*kl[1] - kl_inc[2]*kl[2] - kl_inc[3]*kl[3];

  double Norm_factor;
  
  if( process == 0 ){
    Norm_factor = 1./4;  // normally: also divide by m_e*m_e -> see General_FF, where m_e*m_e is left out in numerator
  }else{
    Norm_factor = 1.;  // see previous remark about the appearance of a factor 2.
  }

  //NOTE: set norm factor to one, to restore lorentz invariance of matrix element defined as Leptonic*Hadron contraction in Calculations.cpp, 

  
  Lepton_S[0][0] = Norm_factor*(2.*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  Lepton_S[0][1] = Norm_factor*(-kl_inc[0]*kl[1] + -kl[0]*kl_inc[1]);  // mind the sign difference for the spatial compts!
  Lepton_S[0][2] = Norm_factor*(-kl_inc[0]*kl[2] + -kl[0]*kl_inc[2]);
  Lepton_S[0][3] = Norm_factor*(-kl_inc[0]*kl[3] + -kl[0]*kl_inc[3]);

  Lepton_S[1][0] = Lepton_S[0][1];
  Lepton_S[1][1] = Norm_factor*(2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  Lepton_S[1][2] = Norm_factor*(kl_inc[1]*kl[2] + kl[1]*kl_inc[2]);
  Lepton_S[1][3] = Norm_factor*(kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);

  Lepton_S[2][0] = Lepton_S[0][2];
  Lepton_S[2][1] = Lepton_S[1][2];
  Lepton_S[2][2] = Norm_factor*(2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  Lepton_S[2][3] = Norm_factor*(kl_inc[2]*kl[3] + kl[2]*kl_inc[3]);

  Lepton_S[3][0] = Lepton_S[0][3];
  Lepton_S[3][1] = Lepton_S[1][3];
  Lepton_S[3][2] = Lepton_S[2][3];
  Lepton_S[3][3] = Norm_factor*(2*kl_inc[3]*kl[3] + kl_inc_dot_kl);

  
  if( process == 0 ){
    
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
	
	Lepton_A[i][j] = 0.;
	
      }
    }
    
  }
  else{
    Lepton_A[0][0] = 0.;
    Lepton_A[0][1] = Norm_factor*(kl_inc[2]*kl[3] - kl[2]*kl_inc[3]);
    Lepton_A[0][2] = Norm_factor*(-kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);
    Lepton_A[0][3] = Norm_factor*(kl_inc[1]*kl[2] - kl[1]*kl_inc[2]);
    
    Lepton_A[1][0] = -Lepton_A[0][1];
    Lepton_A[1][1] = 0.;
    Lepton_A[1][2] = Norm_factor*(kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
    Lepton_A[1][3] = Norm_factor*(-kl_inc[0]*kl[2] + kl[0]*kl_inc[2]);
    
    Lepton_A[2][0] = -Lepton_A[0][2];
    Lepton_A[2][1] = -Lepton_A[1][2];
    Lepton_A[2][2] = 0.;
    Lepton_A[2][3] = Norm_factor*(kl_inc[0]*kl[1] - kl[0]*kl_inc[1]);
    
    Lepton_A[3][0] = -Lepton_A[0][3];
    Lepton_A[3][1] = -Lepton_A[1][3];
    Lepton_A[3][2] = -Lepton_A[2][3];
    Lepton_A[3][3] = 0.;

  }

  
}

