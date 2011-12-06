#ifndef FORM_FACTORS_H 
#define FORM_FACTORS_H

#include "GConstants.h"

double FA(const double q2til); //axial form factor
double FP(const double q2til); //pseudoscalar form factor



//-------------------------------------Dipole-------------------------------------------------------

double DipoleGE(const double q2til);		  //dipole electric form factor G_E^V
double DipoleGM(const double q2til);		  //dipole magnetic form factor G_E^V
double ProtonDipoleGE(const double q2til); //dipole proton electric form factor G_E^V
double ProtonDipoleGM(const double q2til); //dipole proton magnetic form factor G_E^V
double NeutronDipoleGE(const double q2til);//dipole neutron electric form factor G_E^V
double NeutronDipoleGM(const double q2til);//dipole neutron magnetic form factor G_E^V



//-------------------BBBA_05------------------------------------------------------------------------
//arXiv:hep-ex/0602017 BBBA05 for q2til<18 GeV

double BBBA05_GE(const double q2til);	   //electic form factor
double BBBA05_GM(const double q2til);	   //magnetic form factor
double ProtonBBBA05_GE(const double q2til); //proton electic form factor
double ProtonBBBA05_GM(const double q2til); //proton magnetic form factor
double NeutronBBBA05_GE(const double q2til);//neutron electic form factor
double NeutronBBBA05_GM(const double q2til);//neutron magnetic form factor



//--------------------BBA_03------------------------------------------------------------------------
//arXiv:hep-ex/0308005 BBA-2003 for q2til<6 GeV

double BBA03_GE(const double q2til);		  //electic form factor
double BBA03_GM(const double q2til);		  //magnetic form factor
double ProtonBBA03_GE(const double q2til); //proton electic form factor
double ProtonBBA03_GM(const double q2til); //proton magnetic form factor
double NeutronBBA03_GE(const double q2til);//neutron electic form factor
double NeutronBBA03_GM(const double q2til);//neutron magnetic form factor



//--------------------JLab---------------------------------------------------------------------------
//PHYSICAL REVIEW C, VOLUME 65, 051001(R)  //(proton)
//PHYSICAL REVIEW C, VOLUME 51, 409 (1995) //(neutron)

double JLabGE(const double q2til);		//electic form factor
double JLabGM(const double q2til);		//magnetic form factor
double ProtonJLabGE(const double q2til); //proton electic form factor
double ProtonJLabGM(const double q2til); //proton magnetic form factor
double NeutronJLabGE(const double q2til);//neutron electic form factor
double NeutronJLabGM(const double q2til);//neutron magnetic form factor


#endif
