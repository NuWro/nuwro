#ifndef _stale_rpa_h_
#define _stale_rpa_h_
#include <math.h>
#include "jednostki.h"

namespace rpa
{
double M=M12;             

double mm2=m_mu*m_mu;    

//double MA_2=MA*MA;

double M12_2=M12*M12;

double magneton=4.71;      

double m_pi2= m_pi*m_pi;

double m_r2 =m_r*m_r;

double kf=225*MeV;       

double gr =sqrt(1.64*4*Pi); 

double gr2=1.64*4*Pi;

double fr = 6.1*gr;          

double fr2 =fr*fr;

double fpi = sqrt(4*Pi*0.075)*MeV;

//double fpi =  94*MeV;

double fpi2 = fpi*fpi;

double pion = fpi/m_pi; 

double pion2 = pion*pion;

double g_prim = 0.7;           

int znak = 1;  

const int neutrino=1;

const int antyneutrino=-1;

double Mef = 638 * MeV;          

double ilor = Mef/M;

double ilor2 = ilor*ilor;

double Mef2=Mef*Mef;              

double Ef= sqrt(kf*kf + Mef2);   

double Ef2=kf*kf + Mef2;

double ro=kf*kf*kf/3/Pi/Pi;
}
#endif
