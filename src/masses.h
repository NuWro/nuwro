#ifndef _masses_h_
#define _masses_h_
#include "jednostki.h"
//double  M_axial=1050*MeV;
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
//const double Pi = M_PI;
#endif
