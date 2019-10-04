#include "Operators.h"
#include "diracMatrices.h"
#include "Delta_free.h"
#include "D13.h"
#include "S11.h"
#include "P11.h"
#include "NON-Resonant.h"
#include "CT_PF_PP.h"
#include "Regge_CPT.h"


// // // // // // // // // Operators which do not depend on "r" // // // // // // // // // 
void GetOperator(Hadron_prime* kin, Reaction_parameters* par, Matrix Operator[], Matrix Op_delta[])
{
    
      
// // // // // // // // LEM computed
      
    Matrix ChPT[4], Resonances[4];
    if( iModel == 2){
        
// // // LEM operators: ChpT and Resonances

      if (kin->W < W_upper)
      {
	      Get_Res_ChPT(kin, par, ChPT, Resonances, Op_delta);
      }

// // //  // // // ReChi-OPERATOR
      Matrix ReChi[4];
      if (kin->W > W_lower){
      		Get_ReChi(kin, par, ReChi);
      }
// // // // ReChi computed /// /// // / // 



      // // Hyb operator

//       double W0 = 1.7E3; // MeV, center of the transition
//       double leffe = 0.1E3; // MeV, width of the transition (divided by 2)
      double Jupiler = Pi/2.*( 1.- 1./( 1.+exp( (kin->W-W0)/leffe ) ) );
      double cos2Jup = pow(cos(Jupiler),2);
      double sin2Jup = pow(sin(Jupiler),2);

//       for( int i=0; i<4; i++ ){
// 
//      ChPT_ReChi[i] = (cos2Jup*ChPT[i]) + (sin2Jup*ReChi[i]);
// 
//       }
      // // // // // // // //      

      // // // Defining the current operator (without the r-dependent DeltaPole) The delta is already set at this point (in get_resonances)
      for( int i=0; i<4; i++ ){
	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
		Operator[i].M[l][k] = Resonances[i].M[l][k] + cos2Jup*ChPT[i].M[l][k] + sin2Jup*ReChi[i].M[l][k]; // Hybrid model
//		Operator[i].M[l][k] = Resonances[i].M[l][k];
//		Operator[i].M[l][k] = ReChi[i].M[l][k];
	    }
	  }

      }
      
    }else{ // if( iModel == 0 || iModel == 1 ){

        Get_Res_ChPT(kin, par, ChPT, Resonances, Op_delta);

      // // // Defining the current operator (without the r-dependent DeltaPole) The delta is already set at this point (in get_resonances)
      for( int i=0; i<4; i++ ){
	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
		Operator[i].M[l][k] = Resonances[i].M[l][k]+ ChPT[i].M[l][k]; //Low-energy model w/o ff
//		Operator[i].M[l][k] =  ChPT[i].M[l][k]; //Low-energy model w/o ff
	    }
	  }

      }
        
    } // iModel == 0


}

void Get_ReChi(Hadron_prime* kin, Reaction_parameters* par, Matrix ReChi[])
{

      //Set parameters:
      int process = par->process;
      int nucleon = par->nucleon;
      int decay = par->decay;
      int Helicity = par->Helicity;

      if( -(kin->t) > kin->s ){ return; }

        Matrix PF_regge[4], CT_regge[4], NP_regge[4], CNP_regge[4]; // Making a lot of matrix here ...

        Regge( process, nucleon, decay, Helicity, kin->Qsq, kin->ki, kin->kN, kin->Q, kin->kpi, kin->QSlash, kin->kpiSlash, kin->u, kin->s, kin->t, kin->sSlash, kin->uSlash, kin->tSlash, PF_regge, CT_regge, NP_regge, CNP_regge );

	//We sum the matrices using the acces to the private member, should be slightly faster
        for( int i=0; i<4; i++ ){
	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
          	ReChi[i].M[l][k] = PF_regge[i].M[l][k] + CT_regge[i].M[l][k] + NP_regge[i].M[l][k] + CNP_regge[i].M[l][k];
	    }
	  }
	  
//             cout << "ReChi= " << ReChi[i] << endl;
        
        }

}


void Get_Res_ChPT(Hadron_prime* kin, Reaction_parameters* par, Matrix ChPT[], Matrix Resonances[], Matrix Op_delta[])
{
      //Set parameters:
      int process = par->process;
      int nucleon = par->nucleon;
      int decay = par->decay;
      int Helicity = par->Helicity;
      int medmod = par->medmod;
      medmod = 0;

      //Some kinematic invariants:
      double s = kin->s;
      double t = kin->t;
      double u = kin->u;
      double Qsq = kin->Qsq;

      int cross;
      // ///////////////////////////////////
      // // Non-resonant contributions // // 
      /////////////////////////////////////

      // // NUCLEON POLE (non-resonant)
      cross = 0;
      Matrix Op_NP[4];
        NP_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, Op_NP );

      // // CROSS NUCLEON POLE (Non-resonant)
      cross = 1;
      Matrix Op_NP_cross[4];
        CNP_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->uMan, kin->kpiSlash, Op_NP_cross );

      Matrix Op_CT[4], Op_PF[4], Op_PP[4];
        CT_PF_PP( process, nucleon, decay, Helicity, kin->Q, kin->kpi, t, Qsq, kin->QSlash, Op_CT, Op_PF, Op_PP );

      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // HIGHER RESONANCES       
      cross = 0;
      Matrix Op_delta_ffNO[4];
        DP_current( medmod, Pascalutsa, process, nucleon, decay, Helicity, cross, Qsq, s, u, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_delta_ffNO );
        double Deltaff;             //Delta cutoff form factor 
        Delta_ff(s, u, Deltaff);

// // D13 POLE
      cross = 0;
      Matrix Op_D13[4];
        D13P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_D13 );
        double D13ff;
        D13_ff( s, u, D13ff);

      // // // P11 POLE
      cross = 0;
      Matrix Op_P11[4];
        P11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_P11 );
        double P11ff;
        P11_ff( s, u, P11ff);

      // // // // S11 POLE
      cross = 0;
      Matrix Op_S11[4];
        S11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_S11 );
        double S11ff;
        S11_ff( s, u, S11ff);

      // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // HIGHER RESONANCES (u-channel)
      //       
      // // CROSS D13 POLE
      cross = 1;
      Matrix Op_D13_cross[4];
        CD13P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->minusQ, kin->uMan, kin->kN, kin->kpi, kin->uSlash, Op_D13_cross );


      // // // // CROSS P11 POLE
      cross = 1;
      Matrix Op_P11_cross[4];
        CP11P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->uMan, kin->kpiSlash, kin->QSlash,kin->uSlash,Op_P11_cross );

      // // // // CROSS S11 POLE
      cross = 1;
      Matrix Op_S11_cross[4];
        CS11P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->uMan, kin->kpiSlash, kin->QSlash, kin->uSlash, Op_S11_cross );

// // // CROSS DELTA POLE        
       cross = 1;
      Matrix Op_delta_cross[4];
        CDP_current( medmod, Pascalutsa, process, nucleon, decay, Helicity, cross, Qsq, s, u, kin->minusQ, kin->uMan, kin->kN, kin->kpi, kin->uSlash, Op_delta_cross );

//       
      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // // ChPT-background and resonances operator
	// Summation through acces to matrix M inside the matrix class
        
    if( iModel == 0 ){ //LEM w/o resonance form factors
      for( int i=0; i<4; i++ ){
	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
		ChPT[i].M[l][k] = Op_NP[i].M[l][k] + Op_NP_cross[i].M[l][k] + Op_PP[i].M[l][k] + Op_CT[i].M[l][k] + Op_PF[i].M[l][k];

		Resonances[i].M[l][k] = Op_D13[i].M[l][k] + Op_D13_cross[i].M[l][k] +
                			Op_P11[i].M[l][k] + Op_P11_cross[i].M[l][k] +
                                        Op_S11[i].M[l][k] + Op_S11_cross[i].M[l][k]; 

		Op_delta[i].M[l][k] = Op_delta_ffNO[i].M[l][k]+Op_delta_cross[i].M[l][k];
	    } //k
          } //l
      } //i

    }else{ //LEM w ff and Hybrid
        
      for( int i=0; i<4; i++ ){



	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
		ChPT[i].M[l][k] = Op_NP[i].M[l][k] + Op_NP_cross[i].M[l][k] + Op_PP[i].M[l][k] + Op_CT[i].M[l][k] + Op_PF[i].M[l][k];
		
		Resonances[i].M[l][k] = D13ff*(Op_D13[i].M[l][k] + Op_D13_cross[i].M[l][k]) + 
					S11ff*(Op_S11[i].M[l][k]+ Op_S11_cross[i].M[l][k]) +
					P11ff*(Op_P11[i].M[l][k] + Op_P11_cross[i].M[l][k]);
					//Deltaff*Op_delta_cross[i].M[l][k];


		Op_delta[i].M[l][k] = Deltaff*( Op_delta_ffNO[i].M[l][k] + Op_delta_cross[i].M[l][k]);
	    } //k
          } //l
      } //i
        
    } //else (LEM and Hybrid)
    
    
}
