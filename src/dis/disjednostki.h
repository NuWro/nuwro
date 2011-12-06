// Jednostki fizyczne w uk�adzie hbar=c=MeV=1
// (Pozwale nie pisa� hbar ani c we wzorach.
// aby otrzyma� wynik np. w cm  pisz  cout<<wynik/cm <<"[cm]"<<endl;
// idea podpatrzona w geancie, ale tam przyj�to inne jednostki podstawowe.) 
//
// Plik sprawdzony: wszelkie zmiany prosz� opatrywa� komentarzem
// (C. Juszczak 2001)
#ifndef _jednostki_h_
#define _jednostki_h_
#include <cmath>
const double MeV = 1;
const double eV = MeV / 1E6;
const double GeV = 1000 * MeV;
const double GeV2 = GeV*GeV;
const double hbar = 1;
const double c = 1;
const double sek = 1 / (6.58211889 * 1E-22 * MeV);
const double metr = sek / 299792458;
const double km = 1000 * metr;
const double cm = metr / 100;
const double cm2 = cm * cm;
const double barn = 1E-28 * metr * metr;
const double mbarn = barn / 1000;
const double kg = eV / (c * c) / (1.782661731 * 1E-36);
const double J = kg * metr * metr / (sek * sek);
const double gram = kg / 1000;
//double  M_axial=1050*MeV;
const double G = 1.16639 * 1E-5 / (GeV * GeV);
const double Mp = 938.272 * MeV;	//masa protonu
const double Mn = 939.565 * MeV;	//masa neutronu
const double M12 = (Mp + Mn) / 2;	//�rednia masa nukleonu 
const double M2 = M12*M12;	        //kwadrat �rednej masy 
const double MW = 80.423*GeV;
const double MZ = 91.1876*GeV;
//double  MM=M*M; 
const double m_e = 0.5109989 * MeV;	//masa elektronu
const double m_mu = 105.658357 * MeV;	//masa muonu 
const double m_tau = 1777 * MeV;	//masa tau 
const double cos2thetac = 0.9737 * 0.9737;	//kwadrat cosinusa k�ta Cabbibo
const double Pi = M_PI;
//const double Pi = acos(0)*2;
#endif
