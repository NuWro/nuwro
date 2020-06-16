#include "Operators.h"
#include "diracMatrices.h"
#include "Delta_free.h"
#include "D13.h"
#include "S11.h"
#include "P11.h"
#include "NON-Resonant.h"
#include "CT_PF_PP.h"
//#include "Regge_CPT.h"


// // // // // // // // // Operators which do not depend on "r" // // // // // // // // // 
void GetOperator(Hadron_prime* kin, Reaction_parameters* par,  Matrix Operator[])
{
    //This function is used to select which model with which parameters are computed. 
    //Here it is reduced to a single function call: the code only computes LEM w/o FF, keep it to expand implementation later
      
        Get_Res_ChPT(kin, par, Operator);

}


void Get_Res_ChPT(Hadron_prime* kin, Reaction_parameters* par, Matrix Operator[])
{
      //Set parameters:
      int process = par->process;
      int nucleon = par->nucleon;
      int decay = par->decay;
      int Helicity = par->Helicity;
      int medmod = par->medmod; //OSMM Only for interaction with the nucleus, here medmod set to zero to make sure
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
        NP_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->QSlash, kin->sSlash, kin->kpiSlash, Op_NP );


      // // CROSS NUCLEON POLE (Non-resonant)
      cross = 1;
      Matrix Op_NP_cross[4];
        CNP_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->QSlash, kin->uSlash, kin->kpiSlash, Op_NP_cross );

       Matrix Op_CT[4], Op_PF[4], Op_PP[4];
        CT_PF_PP( process, nucleon, decay, Helicity, kin->Q, kin->kpi, t, Qsq, kin->QSlash, Op_CT, Op_PF, Op_PP );

      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // HIGHER RESONANCES       
      cross = 0;
      Matrix Op_delta[4];
      int Pascalutsa = 0; // Pascalutsa = 1 for Pascalutsa decay vertex, needs to be checked -> use 0
      DP_current( medmod, Pascalutsa, process, nucleon, decay, Helicity, cross, Qsq, s, u, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_delta );

// // D13 POLE
      cross = 0;
      Matrix Op_D13[4];
        D13P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_D13 );

      // // // P11 POLE
      cross = 0;
      Matrix Op_P11[4];
        P11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_P11 );

      // // // // S11 POLE
      cross = 0;
      Matrix Op_S11[4];
        S11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_S11 );

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

      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // // ChPT-background and resonances operator
	  // Summation through acces to matrix M inside the matrix class
        
      for( int i=0; i<4; i++ ){
          for (int l = 0 ; l < 4 ; l++)
          {
            for (int k = 0 ; k < 4 ; k++)
            {
            Operator[i].M[l][k] = Op_NP[i].M[l][k] + Op_NP_cross[i].M[l][k] + Op_PP[i].M[l][k] + Op_CT[i].M[l][k] + Op_PF[i].M[l][k] + //Non-resonant
                                 (Op_D13[i].M[l][k] + Op_D13_cross[i].M[l][k]) + //D13
                                 (Op_P11[i].M[l][k] + Op_P11_cross[i].M[l][k]) + //P11
                                 (Op_S11[i].M[l][k] + Op_S11_cross[i].M[l][k]) + //S11
                                  (Op_delta[i].M[l][k]+Op_delta_cross[i].M[l][k]); //Delta
            } //k
          } //l
      } //i

    
    
}
