#ifndef _hyperon_hadron_tensor_h_
#define _hyperon_hadron_tensor_h_


#include <math.h>



class Hyperon_Hadron_Tensor
{
 private:


  // Mass parameters

  double Mp; //nucleon mass
  double My; //hyperon mass
  double M; //nucleon + hyperon mass

  //4 momenta contractions

  double p1k1;
  double p2k2;
  double p1k2;
  double p2k1;
  double p1p2;
  double k1k2;

  double q2;

  double p1q;
  double p2q;
  double k1q;
  double k2q;

  //contractions with the lepton tensor

  // eg: p1mp2n is p1^mu p2^nu L_{mu nu} //contracting the mu and nu 

  double p1mp2n;
  double p2mp1n;
  double qmqn;
  double p1mqn;
  double p2mqn;
  double qmp1n;
  double qmp2n;
  double gmn;


 public:

  //CONSTRUCTORS

  Hyperon_Hadron_Tensor();
  Hyperon_Hadron_Tensor(double initialise[19], double M1, double M2);


  //Used to check input quantities are correct (for debugging)
  void output();


  //ELEMENTS OF HADRON TENSOR

  double F1F1();
  double F2F2();
  double G1G1();
  double MG2G2();
  double G3G3();
  double F1F2();
  double F1G1();
  double F1RG2();
  double F2G1();
  double F2RG2();
  double G1RG2();
  double G1G3();
  double RG2G3();



};


//////////////////////////////////////
// CONSTRUCTORS
///////////////////////////////////////

//default constructor
Hyperon_Hadron_Tensor::Hyperon_Hadron_Tensor(){

}

//constructor with all kinematic variables 
Hyperon_Hadron_Tensor::Hyperon_Hadron_Tensor(double initialise[19],double M1,double M2){



  Mp = M1;
  My = M2;
  M = Mp + My;

  //keep things in the same order used in Hyperon_Interaction.cc

  //contractions of 4 momenta

  p1k1 = initialise[0];
  p2k2 = initialise[1];
  p1k2 = initialise[2];
  p2k1 = initialise[3];
  p1p2 = initialise[4];
  k1k2 = initialise[5];

  q2 = initialise[6];

  p1q = initialise[7];
  p2q = initialise[8];
  k1q = initialise[9];
  k2q = initialise[10];

  //contractions of leptonic tensor with differnt kinematic variables

  p1mp2n = initialise[11];
  p2mp1n = initialise[12];
  qmqn  = initialise[13];
  p1mqn = initialise[14];
  p2mqn = initialise[15];
  qmp1n = initialise[16];
  qmp2n = initialise[17];
  gmn = initialise[18];



}


//prints one of the inpu variables
//compare with the corresponding value in Hyperon_Interaction 
// (for debugging)

void Hyperon_Hadron_Tensor::output(){

  std::cout << k1k2 << std::endl;

}




/////////////////////////////////////////////////////////////////////
// Elements of hadronic tensor
/////////////////////////////////////////////////////////////////////

/*
 LABELLING SYSTEM:

Functions are labelled accoring to the combination of form factors 
they are multiplied by:

eg F1F1 is multiplied by f1f1 , G1RG2 is multiplied by g1*Re(g2) 

MG2G2 is multiplied by |g2|^2


*/


// contribution of f1 * f1

double Hyperon_Hadron_Tensor::F1F1(){

return  4*(p1mp2n + p2mp1n + gmn*Mp*My - gmn*p1p2);


}

// contribution from f2 * f2

double Hyperon_Hadron_Tensor::F2F2(){

  //original
  return -4*(2*gmn*p1q*p2q - gmn*p1p2*q2 + p1mp2n*q2 + p2mp1n*q2 + qmqn*p1p2 - p1mqn*p2q  - p2mqn*p1q  - qmp1n*p2q  - qmp2n*p1q - gmn*Mp*My*q2 + qmqn*Mp*My);

  


}

// contribution from g1 * g1

double Hyperon_Hadron_Tensor::G1G1(){

return  4*(p1mp2n + p2mp1n - gmn*Mp*My - gmn*p1p2);


}


//SECOND CLASS CURRENT
// contribution from |g2|^2

double Hyperon_Hadron_Tensor::MG2G2(){

return  -4*(2*gmn*p1q*p1q - gmn*p1p2*q2 + p1mp2n*q2+p2mp1n*q2+qmqn*p1p2 - p1mqn*p2q - p2mqn*p1q - qmp1n*p2q - qmp2n*p1q + gmn*Mp*My*q2 - qmqn*Mp*My);


}


// contribution from g3 * g3

double Hyperon_Hadron_Tensor::G3G3(){

return  4*(p1p2-Mp*My)*qmqn;


}

// contribution from f1 * f2

double Hyperon_Hadron_Tensor::F1F2(){

 
  //original
  return  4*(qmp1n*My - qmp2n*Mp -p2mqn*Mp + p1mqn*My +  2*gmn*p2q*Mp - 2*gmn*p1q*My);


}

// contribution from f1 * g1

double Hyperon_Hadron_Tensor::F1G1(){


  return 128*(p1k2*p2k1 - p1k1*p2k2);



}


//SECOND CLASS CURRENT
// contribution from f1 * R2(g2)

double Hyperon_Hadron_Tensor::F1RG2(){

return  128*(p1k2*p2k1-p1k1*p2k2)*(My-Mp);


}

// contribution from f2 * g1

double Hyperon_Hadron_Tensor::F2G1(){


  return  128*(p1k2*p2k1 - p1k1*p2k2)*(Mp+My);


}


//SECOND CLASS CURRENT
// contribution from f2 * Re(g2)

double Hyperon_Hadron_Tensor::F2RG2(){

return  128*(My*My - Mp*Mp)*(p1k2*p2k1 - p1k1*p2k2);


}

// SECOND CLASS CURRENT
// contribution from g1 * Re(g2)

double Hyperon_Hadron_Tensor::G1RG2(){

return  4*((-2)*Mp*gmn*p2q - 2*My*gmn*p1q + Mp*qmp2n + Mp*p2mqn + My*qmp1n + My*p1mqn);


}

// contribution from g1 * g3

double Hyperon_Hadron_Tensor::G1G3(){

return  4*(p1mqn*My - p2mqn*Mp + qmp1n*My - qmp2n*Mp);


}


// SECOND CLASS CURRENT
// contribution from Re(g2) * g3

double Hyperon_Hadron_Tensor::RG2G3(){

return  4*(p1mqn*p2q - p2mqn*p1q - qmp2n*p1q + qmp1n*(p2q));


}





#endif
