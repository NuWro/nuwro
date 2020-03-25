#include <math.h> 
#include <iostream>
#include <fstream>
#include "hyperon_hadron_tensor.h" 
#include "ff.h"
#include "jednostki.h"
#include "pdg.h"
#include "vect.h"

/*

CALCULATION OF CONTRACTION OF HADRONIC AND LEPTONIC 
TENSORS IN HYPERON PRODUCTION

Author: Christopher Thorpe

Argument list:

Q2 = Q^2 
E_nu =  neutrino energy
h = selects channel
k1 = neutrino 4 momentum
p1 = nucleon 4 momentum
k2 = lepton 4 momentum
p2 = hyperon 4 momentum  

 */

//DO NOT USE
double Hyperon_Interaction(double Q2, double E_nu, int h, vect k1, vect p1, vect k2, vect p2,bool anti){

  //assume input particles are on shell
  //set Hyperon mass
  double My = sqrt(p2*p2);
  //set nucleon mass
  double Mp = sqrt(p1*p1);

  double M = Mp + My;

   // NOTATION:
   // p1 = nucleon 4 momentum
   // p2 = hyperon 4 momentum
   // k1 = neutrino 4 momentum
   // k2 = lepton 4 momentum

   // q = k1-k2  4 mometum transfer

  //invariant 4 momentum products;

  double p1k1;
  double p2k2;
  double p1k2; 
  double p2k1; 
  double p1p2;
  double k1k2;
  double q2; //q^2 = -Q^2
  double p1q;
  double p2q;
  double k1q;
  double k2q;

  //contractions with Lepton tensor

  double p1mp2n; //L^{mu nu}*(p_mu p'_nu) 
  double p2mp1n; //L^{mu nu}*(p'_mu p_nu)  
  double qmqn; //L^{mu nu}*(q_mu q_nu) 
  double p1mqn; //L^{mu nu}*(p_mu q_nu) 
  double p2mqn; //L^{mu nu}*(p'_mu q_nu)
  double qmp1n; //L^{mu nu}*(q_mu p_nu) 
  double qmp2n; //L^{mu nu)*(q_mu p'_nu_) 
  double p1mp1n; //L^{mu nu)*(p_mu p_nu_) 
  double p2mp2n; //L^{mu nu)*(p'_mu p'_nu_) 

  double gmn;  //L^{mu nu}*g_{mu nu}

  double Mat; //Contraction of hadronic and leptonic tensors

  //form factors
  double f1,f2,f3,f4,g1,Rg2,Ig2,g3;

  //calculate products of 4 momenta

  p1k1 = p1*k1;
  p2k2 = p2*k2;
  p1k2 = p1*k2; 
  p2k1 = p2*k1; 
  p1p2 = p1*p2;
  k1k2 = k1*k2;

  /*  
    std::cout << "Q2: " << Q2 << std::endl;
    std::cout << "p1k1: " <<  p1k1 << std::endl;
    std::cout << "p2k2: " << p2k2 << std::endl;
    std::cout << "p1k2: " << p1k2 << std::endl;
    std::cout << "p2k1: " << p2k1 << std::endl;
    std::cout <<"p1p2: " << p1p2 << std::endl;
    std::cout <<"k1k2: " << k1k2 << std::endl;
    std::cout << std::endl;
  */ 

  p1q = p1k1 - p1k2;
  p2q = p2k1 - p2k2;
  k1q = p2k1 - p1k1;
  k2q = p2k2 - p1k2;


  double s = (p1+k1)*(p1+k1);
  double t = (k1-k2)*(k1-k2);
  double ml = sqrt(k2*k2);


  //calculate contractions with Lepton tensor

  q2 = -Q2;
  p1mp2n = 8*(p1k1*p2k2+p1k2*p2k1-p1p2*k1k2);
  p2mp1n = 8*(p1k1*p2k2+p1k2*p2k1-p1p2*k1k2);
  qmqn = 16*k1q*k2q - 8*q2*k1k2;
  p1mqn = 8*(p1k1*k2q + p1k2*k1q - k1k2*p1q);
  p2mqn = 8*(p2k1*k2q + p2k2*k1q - k1k2*p2q);

  qmp1n = 8*(k1q*p1k2 + p1k1*k2q - k1k2*p1q);
  qmp2n = 8*(k1q*p2k2 + p2k1*k2q - k1k2*p2q);

  gmn = -16*k1k2;
  p1mp1n = 16*p1k1*p1k2 - 8*Mp*Mp*k1k2;
  p2mp2n = 16*p2k1*p2k2 - 8*My*My*k1k2;


  /*
    std::cout << "Q2: " << Q2 << std::endl;
    std::cout << "p1mp2n: " <<  p1mp2n << std::endl;
    std::cout << "p2mp1n: " << p2mp1n << std::endl;
    std::cout << "qmqn: " << qmqn << std::endl;
    std::cout << "p1mqn: " << p1mqn << std::endl;
    std::cout <<"p2mqn: " << p2mqn << std::endl;
    std::cout <<"qmp1n: " << qmp1n << std::endl;
    std::cout <<"qmp2n: " << qmp2n << std::endl;
   std::cout <<"gmn: " << gmn << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  */

  // h+11 = 12,13,14 for the hyperon channels

  //calculate form factors
  list(f1,f2)=f12(q2,h+11);
  list(g1,g3) = fap(q2,h+11);

  //SCC form factors
  //real and imaginary parts of g2 (axial SCC)
  list(Rg2,Ig2) = g2(q2,h+11);

  /*
   std::cout << std::endl;
   std::cout << f1 << std::endl;
   std::cout << f2 << std::endl;
  std::cout << g1 << std::endl;
  std::cout << g3 << std::endl;
  */

  //list of all of the variables needed by H tensor  

  double Hadron_Tensor_Variables[19];

 //contractions of all the different 4 momenta
 Hadron_Tensor_Variables[0] = p1k1;
 Hadron_Tensor_Variables[1] = p2k2;
 Hadron_Tensor_Variables[2] = p1k2;
 Hadron_Tensor_Variables[3] = p2k1;
 Hadron_Tensor_Variables[4] = p1p2;
 Hadron_Tensor_Variables[5] = k1k2;

 Hadron_Tensor_Variables[6] = q2;

 Hadron_Tensor_Variables[7] = p1q;
 Hadron_Tensor_Variables[8] = p2q;
 Hadron_Tensor_Variables[9] = k1q;
 Hadron_Tensor_Variables[10] = k2q;

 //contractions of different pairs of 4 momenta with the lepton tensor
 Hadron_Tensor_Variables[11] = p1mp2n;
 Hadron_Tensor_Variables[12] = p2mp1n;
 Hadron_Tensor_Variables[13] = qmqn;
 Hadron_Tensor_Variables[14] = p1mqn;
 Hadron_Tensor_Variables[15] = p2mqn;
 Hadron_Tensor_Variables[16] = qmp1n;
 Hadron_Tensor_Variables[17] = qmp2n;
 Hadron_Tensor_Variables[18] = gmn;

 //initialise the hadronic tensor
 Hyperon_Hadron_Tensor H(Hadron_Tensor_Variables,Mp,My);


 // select pieces of the hadronic tensor contracted with the leptonic tensor 
 // and multiply by the required form factor to calculate contractions

 //original

 //if antineutrino flip sign of vector/axial interference terms
 if(anti == true){
   g1 *= (-1);
   g3 *= (-1);
   Rg2 *= (-1);
   Ig2 *= (-1); 
 }

 Mat = H.F1F1()*f1*f1 + H.F2F2()*f2*f2/(M*M) + H.F1F2()*f1*f2/M 
   + H.G1G1()*g1*g1 + H.F1G1()*f1*g1 + H.F2G1()*f2*g1/M
   + 4*H.G3G3()*g3*g3/(M*M) + 2*H.G1G3()*g1*g3/M
   + H.MG2G2()*(Rg2*Rg2+Ig2*Ig2)/(M*M) + H.F1RG2()*f1*Rg2/M 
   + H.F2RG2()*f2*Rg2/(M*M) + H.G1RG2()*g1*Rg2/M +2*H.RG2G3()*Rg2*g3/(M*M);
 
 
  /*
 Mat = 64*f1*f1*(p2k1*p1k2 + p1k1*p2k2 - Mp*My*k1k2)

+  64*g1*g1*(p2k1*p1k2 + p1k1*p2k2 + Mp*My*k1k2) 

 - 128*f1*g1*(p1k1*p2k2 - p1k2*p2k1)

+ (32*f2*f2/(M*M))*(-2*Mp*My*k1q*k2q - Mp*My*q2*k1k2 - 2*k1q*k2q*p1p2
+ 2*p2k1*k2q*p1q + 2*p1k1*k2q*p2q + 2*k1q*p2k2*p1q + 2*k1q*p1k2*p2q
			- 2*q2*p2k1*p1k2 - 2*q2*p1k1*p2k2 + q2*p1p2*k1k2)

   + (64*f1*f2/M)*(-Mp*k1q*p2k2 - Mp*p2k1*k2q - Mp*k1k2*p2q 
		   + My*k1q*p1k2 + My*p1k1*k2q + My*k1k2*p1q)
 
 - (128*f2*g1/M)*(Mp*(p2k1*k2q - k1q*p2k2) + My*(p1k1*k2q - k1q*p1k2));
  */
 

 //if(Mat < 0){std::cout << Mat << std::endl; std::cin.get();}
 
  return Mat/32;
}

//implementation of Aligarh MU Model from https://arxiv.org/pdf/1806.08597.pdf

double Singh_Model(double Q2, double E_nu, int h, vect k1, vect p1, vect k2, vect p2,bool anti){

 //assume input particles are on shell
  //set Hyperon mass
  double My = sqrt(p2*p2);
  //set nucleon mass
  double Mp = sqrt(p1*p1);

  double M = Mp + My;

   // NOTATION:
   // p1 = nucleon 4 momentum
   // p2 = hyperon 4 momentum
   // k1 = neutrino 4 momentum
   // k2 = lepton 4 momentum

   // q = k1-k2 = p2-p1  4 mometum transfer
 
  double q2 = -Q2; //q^2 = -Q^2


  double Mat; //Contraction of hadronic and leptonic tensors

  //form factors
  double f1,f2,f3,f4,g1,Rg2,Ig2,g3;

// s and t are Mandelstam variables

  double s = (p1+k1)*(p1+k1);
  double t = (k1-k2)*(k1-k2);
  double ml = sqrt(k2*k2);
  double Delta = My - Mp;

  //also have M = Mp+My

// h+11 = 12,13,14 for the hyperon channels
 
  //calculate form factors
  list(f1,f2)=f12(q2,h+11);
  list(g1,g3) = fap(q2,h+11);

  //SCC form factors
  //real and imaginary parts of g2 (axial SCC)
  list(Rg2,Ig2) = g2(q2,h+11);

  //note sign convention for axial ffs in NuWro is different to Aligarh paper

  //if antineutrino flip signs of vector/axial interference terms
  if(anti == true)
    {
      g1 *= (-1);
      g3 *= (-1); 
      Rg2 *= (-1);
      Ig2 *= (-1); 
    }


Mat =  f1*f1*0.5*(2*(Mp*Mp -s)*(My*My-s) - t*(Delta*Delta-2*s) + t*t + ml*ml*(Delta*Delta - 2*s -t))

+  g1*g1*0.5*(2*(Mp*Mp -s)*(My*My-s) - t*(M*M-2*s) + t*t + ml*ml*(M*M - 2*s -t))

+ (f2*f2/(M*M))*0.25*(-2*t*(Mp*Mp*Mp*Mp-2*s*(Mp*Mp+My*My) + My*My*My*My+2*s*s) + 2*t*t*(M*M-2*s)
 			    + ml*ml*(2*Delta*M*(Mp*Mp+My*My-2*s)+t*((Mp-3*My)*(Mp+My)+4*s)+t*t)
 			    -ml*ml*ml*ml*((3*Mp-My)*M+t))

  + g1*f1*(-1)*(t*(Mp*Mp+My*My - 2*s - t)+ml*ml*(Mp*Mp-My*My+t))

+ (f1*f2/M)*(-1)*(t*M*(Delta*Delta-t) + ml*ml*( (-1)*Delta*(My*My-s)+My*t) +ml*ml*ml*ml*Mp)

  + (f2*g1/M)*(-M)*(t*(Mp*Mp+My*My-2*s-t)+ml*ml*(Mp*Mp-My*My+t)) 

 + (g3*g3/(M*M))*ml*ml*(ml*ml-t)*(Delta*Delta -t)

   + (g1*g3/M)*(-2)*(ml*ml)*(ml*ml*Mp + Mp*Mp*Mp - Mp*Mp*My - Mp*(s+t) + My*s)
  
  //second class current

 + ((Rg2*Rg2+Ig2*Ig2)/(M*M))*(0.25*(4*(Delta*Delta - t)*((Mp*Mp-s)*(My*My-s)+s*t) + ml*ml*(4*Delta*(Mp*Mp*Mp+Mp*Mp*My - Mp*(3*s+t)+My*s) + 2*Delta*Delta*(M*M - 2*s - t) - (4*s+t)*(Delta*Delta -t))  + 2*Delta*Delta*(-2*(Mp*Mp-s)*(My*My-s)-t*(M*M +2*s)+t*t) + ml*ml*ml*ml*(Delta*Delta + 4*Mp*Delta - t)))

  + (Rg2*f1/M)*(-Delta*(t*(Mp*Mp + My*My - 2*s-t) + ml*ml*(Mp*Mp - My*My + t)))

  + (Rg2*f2/(M*M))*(Delta*(-M)*(t*(Mp*Mp + My*My - 2*s-t)+ml*ml*(Mp*Mp-My*My+t)))

  + (Rg2*g1/M)*((Delta*(-t*(M*M)+t*t) + ml*ml*(Mp*Mp*Mp + Mp*Mp*My + Delta*(M*M - 2*s - t ) - 3*Mp*s - M*t + My*s) + ml*ml*ml*ml*Mp))

  + (Rg2*g3/(M*M))*(ml*ml*(-2*Delta*(ml*ml*Mp + Mp*Mp*Mp - Mp*Mp*My - Mp*(s+t) + My*s) - (Delta*Delta-t)*(ml*ml + 2*Mp*Mp - 2*s -t)));

							
//if sqaured matrix element is negative there is a problem!						 
//if(Mat < 0){std::cout << Q2 << "  "  << Mat << std::endl; std::cin.get();}


  return Mat;

}
