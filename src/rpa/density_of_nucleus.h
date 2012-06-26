#ifndef _density_of_nucleus_h_
#define _density_of_nucleus_h_

//#include<iostream>
//#include<fstream>
#include<cmath>
#include   "jednostki.h"
//#include <strstream>
#include   "calg.h"
#include   "tablica_Mef.h"

using namespace std;

        //       Wspolczynniki okreslajace rozklad pedow w jadrze.
        //        Pochodza one z hep-ph/0005255 E. A. Paschos, L. Pasquali and J. Y. You 
	    //        orginalnie poraz pierwszy: S. L. Adler, S. Nussinov and E. A. Paschos, Phys. Rev. D 9, 2125 (1974)
       //-------    Tlen    -------//

namespace rpa
{
    const double a_O       =2.718*fm;
    const double C_O       =1.544;//*fm   --- Tutaj ma nie by fm a po prostu bez jednostek
    const double C1_O     =0;
    const double R_O       =1.833*fm;
    const double rho0_O  =0.141/fm/fm/fm;

         //-------    Argon   -----//

    const double a_Ar      =3.393*fm;
    const double C_Ar      =3.530*fm;  
    const double C1_Ar    =0.542*fm;  
    const double R_Ar      =4.830*fm;
    const double rho0_Ar =0.176/fm/fm/fm;

       //-------   Fe -zelazo ----//

   const double a_Fe      =3.801*fm;
   const double C_Fe      =4.111*fm; 
   const double C1_Fe    =0.558*fm; 
   const double R_Fe      =4.907*fm;
   const double rho0_Fe =0.163/fm/fm/fm;

   //-----------// Liczby porzadkowe pierwiastkow //-------------//  
   const int Ar =18;
   const int Fe =26;
   const int O =8;
   const int C = 6;
   
   //----------//Liczby atomowe pierwiastkow    //--------------//
   const int A_Ar=40;
   const int A_Fe=56;
   const int A_O =16;
   const int A_C =12;   
   /// Rozklad gestos ladunkow  jadra:
    /// W zaleznosci od promienia;
      
double density_O(double r)
{
	return rho0_O*exp(-pow(r/R_O,2))*(1+C_O*pow(r/R_O,2)+C1_O*pow(r/R_O,4));
}
  
double density_Ar(double r)
{
	return rho0_Ar/(1.+exp((r-C_Ar)/C1_Ar));
}
  
double density_Fe(double r)
{
	return rho0_Fe/(1.+exp((r-C_Fe)/C1_Fe));
}

double density(double r, int rodzaj_jadra)
{
	switch(rodzaj_jadra)
	{
		case O  : return rho0_O*exp(-pow(r/R_O,2))*(1+C_O*pow(r/R_O,2)+C1_O*pow(r/R_O,4));
		case Ar : return rho0_Ar/(1.+exp((r-C_Ar)/C1_Ar));
		case Fe : return rho0_Fe/(1+exp((r-C_Fe)/C1_Fe));
		default : return 0; 
	} 
}
  
  /// Pedy Fermiego zadane gestoscia jadra ///
  
double Fermi_momentum(double r,int nucleus)
{
	return cbrt(3*Pi*Pi*density(r,nucleus)/2 );
}

double density_calka(double r, int nucleus)
{
	return 4*Pi*r*r*density(r,nucleus); 	       
}

double density_r(double r, int nucleus)
{
	return r*density(r,nucleus);
}

    
double kf_density(double r, int nucleus)
{
	return Fermi_momentum(r,nucleus)*density_calka(r,nucleus);     
}      
      
double Mf_density(double r, int nucleus)
{
	return Masa_Efektywna(Fermi_momentum(r,nucleus))*density_calka(r,nucleus);     
}

double sredni_ped_fermiego(int nucleus)
{
   return calg5x_int(kf_density, nucleus, 0, 5*R_O,  0.01, 100)/calg5x_int(density_calka,nucleus, 0, 5*R_O, 0.01, 100);
}

const char* nazwa_jadra(int rodzaj_jadra)
//string nazwa_jadra(int rodzaj_jadra)
{
 	switch(rodzaj_jadra)
	{
		case O : return "Tlen"  ; break;
		case Ar: return "Argon" ; break;
		case Fe: return "Ferryt"; break;
		case C : return "Carbon"; break;
		default : return "_" ;    break;
	}
}

double sredni_Mef(int nucleus)
{
	return calg5x_int(Mf_density, nucleus, 0, 7*R_Fe,  0.01, 20)/calg5x_int(density_calka,nucleus, 0, 7*R_Fe, 0.01, 20);
}

double sredni_Mef_st(int nucleus)
{
	switch(nucleus)
	{
		case O : return 690.78*MeV  ; break;
		case Ar: return 631.37*MeV ; break;
		case Fe: return 634.84*MeV; break;
		//case C : return "Carbon"; break;
		default : return 1000 ;    break;
	}
}


int liczba_atomowa(int rodzaj_jadra)
{
	switch(rodzaj_jadra)
	{
		case O : return A_O  ; break;
		case Ar: return A_Ar ; break;
		case Fe: return A_Fe; break;
		case C : return A_C; break;
		default : return 0 ;    break;
	}
}       

}
#endif
