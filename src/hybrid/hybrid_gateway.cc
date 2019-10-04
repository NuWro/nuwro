#include <complex>
#include <cstdarg>

using namespace std;

#include "Matrix.h"
#include "Calculations.h"
#include "Constants.h"
#include "diracMatrices.h"

// // // for the interpolations // // //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// // // // // // // // // // // // // // 

#include <cstdlib>
#include <cassert>
#include <cstdio>


int hybrid_ABCDE(double El_inc, double Q2, double W, double *costheta_pi_in, int N, int *params, double (*strucfuncs)[5])
{
  //Making the structs because it is what we use in the code, these will be removed at some point, because it may not be efficient
  Hadron_prime Hadronbag;
  Lepton_kin Leptonbag;

  //We need to set the pointers to something usefull that will not go out of scope at this point
  double kl[4],kl_inc[4], Q[4], Pa[4];
  Leptonbag.kl=kl, Leptonbag.kl_inc=kl_inc, Leptonbag.Q= Q, Leptonbag.Pa =Pa;

  //Same for Hadron Kin
  double xQ[4], xminusQ[4];
  double kpi[4], kN[4], ki[4];         //4vectors in the kN//z system
  Hadronbag.Q= xQ, Hadronbag.minusQ=xminusQ, Hadronbag.kpi = kpi, Hadronbag.kN=kN, Hadronbag.ki=ki;

  //Mandelstam variables
  double sMan[4], tMan[4], uMan[4];
  Hadronbag.sMan=sMan;
  Hadronbag.tMan=tMan;
  Hadronbag.uMan=uMan;

  Reaction_parameters Reacbag;
  double ABCDE[5], SL_factors[5], e_responses[5];
  //Fill Reac
  Reacbag.process = params[0], Reacbag.decay= params[1], Reacbag.Helicity = params[2] , Reacbag.nucleon = params[3];
  Reacbag.ABCDE = ABCDE, Reacbag.SL_factors = SL_factors, Reacbag.e_responses = e_responses;

  //Make the pointers:
  Hadron_prime* Had = &Hadronbag;
  Lepton_kin* Lep = &Leptonbag;
  Reaction_parameters* Reac = &Reacbag;

  //////////////////////////////////////////////
  ///////Lepton kinematics:
  //////////////////////////////////////////////////

  double leptonmass = 0.00;  //We dont neglect it, because everything is numerical anyway, woohoo!
  if( Reac->process == 1 ){leptonmass = muonmass; }

  double W2 = W*W;
  double QQ = Q2;

  double Enu = El_inc;
  double w = (W*W + QQ - MN*MN)/(2*MN); //This omega is  in LAB!
  if (w <= 0 || w > Enu){return -1;}
  if (pow((Enu - w),2) - leptonmass*leptonmass <= 0){return -1;}
  double k_l = sqrt(pow((Enu - w),2) - leptonmass*leptonmass);
  double E_l = Enu - w;
  if (fabs((QQ+w*w - pow(Enu,2) - pow(k_l,2))/(2.*Enu*k_l)) > 1){return -1;}
  double radThetal = acos((QQ+w*w - pow(Enu,2) - pow(k_l,2))/(-2.*Enu*k_l) );
  double q = sqrt(QQ + w*w);

  //define the boost parameters
  double v = q/(w+MN);
  double gamma = (w + MN)/W;

  //define CMS angles for the leptons
  double costheta = (Enu - k_l*cos(radThetal))/q;
  double sintheta = k_l*sin(radThetal)/q;

  //The leptonic side is now given as:
  Lep->k_l_inc = gamma*(El_inc - v*Enu*costheta);

  Lep->kl[0] = gamma*(E_l - v*(Enu*costheta - q));
  Lep->kl[1] = sintheta*Enu;
  Lep->kl[2] = 0.;
  Lep->kl[3] = gamma*(-1*v*E_l + (Enu*costheta - q));
  Lep->k_l = sqrt(pow(Lep->kl[0],2) - leptonmass*leptonmass);

  Lep->kl_inc[0] = Lep->k_l_inc; //Neglect incoming lepton rest mass
  Lep->kl_inc[1] = Lep->kl[1];
  Lep->kl_inc[2] = 0;
  Lep->kl_inc[3] = gamma*(-1*v*Enu+ Enu*costheta);

  Lep->Q[0] = Lep->kl_inc[0] - Lep->kl[0];
  Lep->Q[1] = Lep->kl_inc[1]-Lep->kl[1];
  Lep->Q[2] = 0;
  Lep->Q[3] = Lep->kl_inc[3] - Lep->kl[3];

  Lep->Qsq = QQ;
  Lep->leptonmass = leptonmass;
  Lep->El = Lep->kl[0];
  Lep->w = Lep->Q[0];
  Lep->q = Lep->Q[3];

  /////////////////////////////////////////////////////
  ///////HADRON kinematics (the part invariant under rotations)
  ////////////////////////////////////////////////////////////

  double E_pi = (W*W + pow(Mpi,2) - pow(MN,2))/(2*W);
  double E_N = W - E_pi;
  double p_N = MN*q/W;
  Lep->Pa[0] = W - Lep->Q[0];
  Lep->Pa[1] = 0;
  Lep->Pa[2] = 0;
  Lep->Pa[3] = -p_N;

  Had->Q[0] = Lep->Q[0];
  Had->Q[1] = Lep->Q[1];
  Had->Q[2] = Lep->Q[2];
  Had->Q[3] = Lep->Q[3];

  Had->minusQ[0] = -Lep->Q[0];
  Had->minusQ[1] = -Lep->Q[1];
  Had->minusQ[2] = -Lep->Q[2];
  Had->minusQ[3] = -Lep->Q[3];

  //Now the pion and final nucleon, the kinematics are like decay at rest, so momenta are opposite, and determined by energy balance
  double Tpi = E_pi - Mpi;
  if (Tpi < 0){return -1;}

  Had->ki[0] = W - Lep->Q[0];
  Had->ki[1] = 0;
  Had->ki[2] = 0;
  Had->ki[3] = -p_N;

  Had->kiSlash =  Had->ki[0]*Gamma[0]-Had->ki[1]*Gamma[1]-Had->ki[2]*Gamma[2]-Had->ki[3]*Gamma[3]; //we comment it out, because spatial components of ki are zero
  Had->QSlash = Had->Q[0]*Gamma[0]-Had->Q[1]*Gamma[1]-Had->Q[2]*Gamma[2]-Had->Q[3]*Gamma[3];

  //Mandelstam s is trivial in CMS:
  Had->sMan[0] = W;
  Had->sMan[1] = 0;
  Had->sMan[2] = 0;
  Had->sMan[3] = 0;
  Had->sSlash = Had->sMan[0]*Gamma[0] - Had->sMan[1]*Gamma[1]- Had->sMan[2]*Gamma[2] - Had->sMan[3]*Gamma[3]; 

  Had->Qsq = Lep->Qsq;

  Had->k_pi = sqrt(Tpi*(Tpi+2*Mpi));
  Had->kpi[0] = E_pi;
  Had->kpi[2] = 0;

  Had->kN[0] = W - E_pi;
  if (Had->kN[0] < MN){return -1;}
  Had->kN[2] = 0;

  Had->s = W*W;

  for (int i = 0 ; i < N ; i++)
  {
    double CT = costheta_pi_in[i];
    double radThetapi = acos(CT);

    Had->kpi[1] = sin(radThetapi)*Had->k_pi;
    Had->kpi[3] = cos(radThetapi)*Had->k_pi;

    Had->kN[1] = -Had->kpi[1];
    Had->kN[3] = -Had->kpi[3];

    double prefactor = 0.5;

    //Mandelstam
    for(int i=0; i<4; i++)
    {
      Had->uMan[i] = Had->kN[i]  - Had->Q[i];
      Had->tMan[i] = Had->Q[i]   - Had->kpi[i];
    }
    Had->u = pow(Had->uMan[0],2) - pow(Had->uMan[1],2) - pow(Had->uMan[2],2) - pow(Had->uMan[3],2);
    Had->t = pow(Had->tMan[0],2) - pow(Had->tMan[1],2) - pow(Had->tMan[2],2) - pow(Had->tMan[3],2);

    //Slashes (to be optimized!)
    Had->kNSlash = Had->kN[0]*Gamma[0]-Had->kN[1]*Gamma[1]-Had->kN[2]*Gamma[2]-Had->kN[3]*Gamma[3];
    Had->kpiSlash = Had->kpi[0]*Gamma[0]-Had->kpi[1]*Gamma[1]-Had->kpi[2]*Gamma[2]-Had->kpi[3]*Gamma[3];
    Had->tSlash = Had->tMan[0]*Gamma[0] - Had->tMan[1]*Gamma[1]- Had->tMan[2]*Gamma[2] - Had->tMan[3]*Gamma[3];
    Had->uSlash = Had->uMan[0]*Gamma[0] - Had->uMan[1]*Gamma[1]- Had->uMan[2]*Gamma[2] - Had->uMan[3]*Gamma[3];

    Do_Fivefold_Calc(Reac, Lep, Had); //Same as in function DDCS
    strucfuncs[i][0] = Reac->ABCDE[0]*prefactor; 
    strucfuncs[i][1] = Reac->ABCDE[1]*prefactor; 
    strucfuncs[i][2] = Reac->ABCDE[2]*prefactor; 
    strucfuncs[i][3] = Reac->ABCDE[3]*prefactor; 
    strucfuncs[i][4] = Reac->ABCDE[4]*prefactor; 
  }

  return 0;
}