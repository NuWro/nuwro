#include <iostream>
#include <iomanip>  // this include file enables the use of I/O manipulators for controlling the format of the output stream
#include <fstream>
#include <complex>
#include <cstdarg>  // this include file enables working with variable length parameter lists

using namespace std;

#include "diracMatrices.h"
#include "Matrix.h"
#include "Calculations.h"
#include "Constants.h"
#include "Leptonic.h"
#include "Operators.h"
#include "ABCDE.h"
// 


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


// // // // 
//

double CS_interpolate(double QQ,double W,CS_tabulated* CS_table, Lepton_kin* Lep, Reaction_parameters* Reac)
{

	//We build the cross section from the table. The lepton kinematics are determined by the Leptonic tensor, and we take the lazy route and just build the whole lepton tensor. One can in principle just calculate the 5 terms directly from the incoming, outgoing energies and scattering angle in the CMS. This has been done in Nakamura/Sato paper, where the different terms are given in terms of omega, q etc.

	double helic = Reac->Helicity*1.; // -1 (neut), 1 (antineut)
	double Lepton_S[4][4];
	double Lepton_A[4][4];

	//Compute symmetric and antisymmetric part of Lepton tensor based on the Fourmomenta kl and kl_inc (in CMS!). The Lepton Tensor contains different normalization for electron versus neutrino interactions ( 2 and 1/2 respect)	
	Leptonic_Tensor_fourvector( Reac->process, Lep->kl, Lep->kl_inc ,Lepton_S, Lepton_A );

	complex<double> Lepton[4][4];
	
	for(int i=0; i<4; i++){
	  for(int j=0; j<4; j++){

	    if( Reac->process == 0 ){
	      Lepton[i][j] = Lepton_S[i][j];
	    }else{      
	      Lepton[i][j] = Lepton_S[i][j] -I*helic* Lepton_A[i][j];
	    }

	  }
	}

        double	Lepton_factors[5];
	
	Lepton_factors[0] = real(Lepton[0][0]);
	Lepton_factors[1] = real(2.*Lepton[3][0]);
	Lepton_factors[2] = real(Lepton[3][3]);
	Lepton_factors[3] = real(0.5*(Lepton[1][1]+Lepton[2][2]));
	cout << Lepton_factors[3] << endl;
	Lepton_factors[4] = real(2.*I*Lepton[1][2]); 


	double W_step = CS_table->W_step;
	double QQ_step = CS_table->QQ_step;

	double *QQ_vec = CS_table->QQ_vec;
	double *W_vec = CS_table->W_vec;

	//interpolate the table on a grid of constant stepsize:

    if( W < W_vec[0]){W = W_vec[0];}
    else if(W >= CS_table->W_fin){W=CS_table->W_fin-0.00000001;} //should never happen, if it does its problems


    if( QQ < QQ_vec[0]){QQ = QQ_vec[0];}
    else if(QQ >= CS_table->QQ_fin){QQ=CS_table->QQ_fin-0.00000001;} //should never happen, if it does its problems

    int index_W = (int)((W-W_vec[0])/W_step);
    int index_QQ = (int)((QQ-QQ_vec[0])/QQ_step);


    double x2x, y2y, yy1, xx1, x1,x2,y1,y2;
    x1 = W_vec[0] + index_W*W_step;
    x2 = x1 + W_step;
    y1 = QQ_vec[0] + index_QQ*QQ_step;
    y2 = y1 + QQ_step;
   
    x2x = x2 - W;
    xx1 = W - x1;
    y2y = y2 - QQ;
    yy1 = QQ - y1;
    
    double den = W_step * QQ_step;
    x2x = x2x/den;
    xx1 = xx1/den;
    double x2x_y2y = x2x * y2y;
    double xx1_y2y = xx1 * y2y;
    double x2x_yy1 = x2x * yy1;
    double xx1_yy1 = xx1 * yy1;
    
    double q11,q12,q22,q21;

    double CS = 0;
    for (int i = 0 ; i < 5 ;  i++)
    {
	    q11 = CS_table->table[index_QQ][index_W][i];
	    q12 = CS_table->table[index_QQ+1][index_W][i];
	    q21 = CS_table->table[index_QQ][index_W +1][i];
	    q22 = CS_table->table[index_QQ+1][index_W+1][i];

    //sum the five different inclusive responses and multiply by lepton factors
    CS += Lepton_factors[i]*(
	q11 * x2x_y2y +
	q21 * xx1_y2y +
	q12 * x2x_yy1 +
	q22 * xx1_yy1 );

    }

    return CS;
}

double Do_Fivefold_Calc( Reaction_parameters *Reac, Lepton_kin* Lep, Hadron_prime * Had)
{ 
  
// // // // NOTATION // // // // 
// process = 1   CC interaction
//   HELICITY == -1 --> W^+ induced 1-pion production 
// 	nucleon = 1 --> proton initial state
// 		decay = 1 or 2 --> p + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
//   
//   HELICITY == 1 --> W^- induced 1-pion production 
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 or 2 --> n + pi^-
// 	nucleon = 1 --> proton initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^+
// //    
// process = 0   EM interaction --> photon induced 1-pion production 
// 	nucleon = 1 --> proton initial state 
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// // // // // // // // // // //
      
  
	/******************************************
	*       Construct the lepton tensor       *        
	******************************************/


	double helic = Reac->Helicity*1.; // -1 (neut), 1 (antineut)
	double Lepton_S[4][4];
	double Lepton_A[4][4];

	//Compute symmetric and antisymmetric part of Lepton tensor based on the Fourmomenta kl and kl_inc. The Lepton Tensor contains different normalization for electron versus neutrino interactions ( 2 and 1/2 respect)	
	Leptonic_Tensor_fourvector( Reac->process, Lep->kl, Lep->kl_inc ,Lepton_S, Lepton_A );

	complex<double> Lepton[4][4];
	
	for(int i=0; i<4; i++){
	  for(int j=0; j<4; j++){

	    if( Reac->process == 0 ){
	      Lepton[i][j] = Lepton_S[i][j];
	    }else{      
	      Lepton[i][j] = Lepton_S[i][j] -I*helic* Lepton_A[i][j];
	    }

	  }
	}
    
  
  /******************************************
   *      Construct the Operators      *
   ******************************************/

   Matrix Operator[4], Op_delta[4], Op[4];
   GetOperator(Had, Reac, Operator, Op_delta); //The type of operator depends on the flag iModel which is set in constants.h!


   //Add the Delta (contains Delta + Crossed delta) to the Operator, contains all other resonances and ChPT background
   for (int i = 0 ; i < 4 ; i++)
   {
	
   	for (int j = 0 ; j < 4 ; j++)
	{
		
   		for (int k = 0 ; k < 4 ; k++)
		{
			Op[i].M[j][k] = Op_delta[i].M[j][k] + Operator[i].M[j][k] ; 
		}
	}
   }
	
	/**************************
	  BUILD THE HADRON TENSOR
	**************************/

	//Building the matrices (kslash + M)*Gamma
	
	complex<double> block1[4][4], block2[4][4];
	
	double Mfac;
	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j = 0 ; j < 4 ; j++)
		{
			Mfac = real(MN*Id.M[i][j]);
			block1[i][j] = (1./(4.*MN))*(Had->kiSlash.M[i][j] + Mfac)*Gamma[0].M[j][j];
			block2[i][j] = (1./(2.*MN))*(Had->kNSlash.M[i][j] + Mfac)*Gamma[0].M[i][i];
		}
	}
	//Block 1 and 2 contain the factors 1/2 and 1/4 from summation over helicities already, they are not added later.




	//Conjugating the Operator, and then computing the "Sandwich" which is (ki_slash + M)*Op^+*(kN_slash +M)
	Matrix Conj_Op[4];
	for (int i = 0 ; i < 4 ; i++)
	{
		Conj_Op[i] = ConjMatrix(Op[i]);
	}

	Matrix Sandwich_Op[4];
	for (int i = 0 ; i < 4 ; i++)
	{
		for (int m = 0 ; m < 4 ; m++)
		{
			for (int j = 0 ; j < 4 ; j++)
			{
				Sandwich_Op[i].M[j][m] = 0.0;	
				for (int k = 0 ; k < 4 ; k++)
				{
					for (int l = 0 ; l < 4 ; l++)
					{
						Sandwich_Op[i].M[j][m] += block1[j][l]*Conj_Op[i].M[l][k]*block2[k][m];
					}
				}
				
			}
		}
	}
	
	
	complex<double> Hadron[4][4] = {}; //Hadron tensor is Trace(Sandwich[i]*Op[j])

	for(int i=0; i<4; i++){
	  for(int j=0; j<4; j++){
	    Trace_product(Sandwich_Op[i], Op[j], Hadron[i][j]); // The trace op Sandwich_Op*Op is stored in Hadron[i][j]
	  }
	}
	//Calculate the factors ABCDE, defined as in J.E. Sobczycks paper (angular distributions ...) but without any prefactors except for the additional normalization in Lepton tensor an 1/(8MN^2) from hadron
   	ABCDE(Hadron, Lepton, Reac->ABCDE, Reac->SL_factors, Reac->e_responses, Lep, Lep);

	//The SL_factors are five terms one needs to calculate the double differential cross section in an energy-independent way. One can combine the relevant leptonic expressions with them later as detailed in SL paper or just from the definition of the structure function A to actually get the cross section as function of energy.
	

    return 0;
    
  
}


