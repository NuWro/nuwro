// Jednostki fizyczne w uk³adzie hbar=c=MeV=1
// (Pozwale nie pisaæ hbar ani c we wzorach.
// aby otrzymaæ wynik np. w cm  pisz  cout<<wynik/cm <<"[cm]"<<endl;
// idea podpatrzona w geancie, ale tam przyjêto inne jednostki podstawowe.) 
//
// Plik sprawdzony: wszelkie zmiany proszê opatrywaæ komentarzem
// (C. Juszczak 2001)
#ifndef _jednostki_h_
#define _jednostki_h_
#include <cmath>
double  ksi=3.71;   
const double  MeV=1;
const double  eV=MeV/1E6;
const double  GeV=1000*MeV;
const double  GeV2=GeV*GeV;
const double  TeV=1000*GeV;
const double  hbar=1;
const double  c=1;
const double  sek=1/(6.58211889*1E-22*MeV);
const double  metr=sek/299792458;
const double  km=1000*metr;
const double  cm=metr/100;
const double  cm2=cm*cm;
const double  fm =1e-13*cm;
const double  barn=1E-28*metr*metr;
const double  mbarn=barn/1000;
const double  kg=eV/(c*c)/(1.782661731*1E-36);
const double  J=kg*metr*metr/(sek*sek);
const double  gram=kg/1000;
//const double  MA =1030*MeV; // Masa Aksialna
//const double  MA_2=MA*MA;   // Kwadrat masy aksialnej 
//const double  MV2= 710000*MeV*MeV; //Kwadrat masy wektorowej z Form Faktorow
const double  G=1.16639*1E-5/(GeV*GeV);
const double  M1=938.272*MeV; //masa protonu
const double  M2=939.565*MeV; //masa neutronu
const double  M12=(M1+M2)/2;  //¶rednia masa nukleonu 
const double  m_e=0.5109989*MeV;      // masa elektronu
const double  m_mu=105.658357*MeV;    // masa muonu           
const double  m_tau=1777*MeV;         // masa tau 
const double  m_pi = 139.57018*MeV;   // masa mezonu pi
const double  m_r  = 769.3 *MeV;      // masa mezonu rho  
const double  cos2thetac=0.9737*0.9737; //kwadrat cosinusa k±ta Cabbibo
const double  sin2thetac=1-cos2thetac;
const double  sin2thetaW=0.23113;
/// Dodane przez k.g.
const double magneton_N       = 1;   // tego nie jestem pewien !!!!
const double magneton_proton  = 2.793*magneton_N;
const double magneton_neutron =-1.913*magneton_N;
const double Pi=M_PI;

//~ const int O_16_n =8;
//~ const int O_16_p =8;
//~ const int Ar_40_n=22;
//~ const int Ar_40_p=18;
//~ const int Fe_56_n=30;
//~ const int Fe_56_p=26;
//~ const int C_12_n=6;
//~ const int C_12_p=6;
//~ 

#endif
