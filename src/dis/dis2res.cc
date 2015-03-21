#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "pdg_name.h"
#include "jednostki.h"
#include "masses.h"
#include "parameters.h"
#include "dis_cr_sec.h"
#include "delta.h"
#include "charge.h"
#include "generatormt.h"


extern "C" void py2ent_(int *, const int *, const int *, double *);


double spp_p_cc[200][3];//single pion production produced number
double spp_n_cc[200][3];//single pion production produced number
double spp_p_nc[200][3];//single pion production produced number
double spp_n_nc[200][3];//single pion production produced number

double dis_p_cc[200];//dis event produced number; spp + >spp
double dis_n_cc[200];//dis event produced number; spp + >spp
double dis_p_nc[200];//dis event produced number; spp + >spp
double dis_n_nc[200];//dis event produced number; spp + >spp

double rsp_p_cc[200][3];//ratio od single pion = spp/dis. For W<1.2GeV equal 1
double rsp_p_nc[200][3];//ratio od single pion = spp/dis. For W<1.2GeV equal 1
double rsp_n_nc[200][3];//ratio od single pion = spp/dis. For W<1.2GeV equal 1
double rsp_n_cc[200][3];//ratio od single pion = spp/dis. For W<1.2GeV equal 1


double alfa0;
double W_min;
double W_max;

double cr_dis;
double cr_res;
bool disres;
//function combaning dis and delta contributions for small W
int p_plus=0;
int p_zero=0;

void funkcja_spp_0()
{
for(int  k=0; k<200; k++)
for(int  j=0; j<3; j++)
    {
       rsp_p_cc[k][j]=1;
       rsp_p_nc[k][j]=1;
       rsp_n_cc[k][j]=1;
       rsp_n_nc[k][j]=1;
   }
}


double funkcja_spp_3(double W, int nukleon_in, int pion, bool current)
{

for(int  k=0; k<200;k++)
    {
      if(W>1*GeV+10*MeV*k && W<1*GeV+10*MeV*(k+1))
       {
       if(pion== 211 && current == true && nukleon_in == proton){return rsp_p_cc[k][0];}
       if(pion== 111 && current == true && nukleon_in == proton){return rsp_p_cc[k][1];}
       if(pion==-211 && current == true && nukleon_in == proton){return rsp_p_cc[k][2];}

       if(pion== 211 && current == true && nukleon_in == neutron){return rsp_n_cc[k][0];}
       if(pion== 111 && current == true && nukleon_in == neutron){return rsp_n_cc[k][1];}
       if(pion==-211 && current == true && nukleon_in == neutron){return rsp_n_cc[k][2];}

       if(pion== 211 && current ==false && nukleon_in == proton){return rsp_p_nc[k][0];}
       if(pion== 111 && current ==false && nukleon_in == proton){return rsp_p_nc[k][1];}
       if(pion==-211 && current ==false && nukleon_in == proton){return rsp_p_nc[k][2];}

       if(pion== 211 && current ==false && nukleon_in == neutron){return rsp_n_nc[k][0];}
       if(pion== 111 && current ==false && nukleon_in == neutron){return rsp_n_nc[k][1];}
       if(pion==-211 && current ==false && nukleon_in == neutron){return rsp_n_nc[k][2];}

       }
       
    }
return 0;
}


double funkcja_spp_2(double W, int nukleon_in, bool current)
{
for(int  k=0; k<200;k++)
    {
      if(W>1*GeV+10*MeV*k && W<1*GeV+10*MeV*(k+1))
       {
       if(current == true){return (spp_p_cc[k][0]+spp_p_cc[k][1]+spp_p_cc[k][2])/dis_p_cc[k];}
       if(current ==false){return (spp_p_nc[k][0]+spp_p_nc[k][1]+spp_p_nc[k][2])/dis_p_nc[k];}
       if(current == true){return (spp_n_cc[k][0]+spp_n_cc[k][1]+spp_n_cc[k][2])/dis_n_cc[k];}
       if(current ==false){return (spp_n_nc[k][0]+spp_n_nc[k][1]+spp_n_nc[k][2])/dis_n_nc[k];}

       }
    }
return 0;
}


int meson_out_(double W, int lepton_in, int nukleon_in, bool current)
{

if(current == true && lepton_in > 0 && nukleon_in == neutron)
    {
    if(funkcja_spp_3(W, nukleon_in, piplus, current)>frandom()){return piplus;}
    else{return pizero;	}
    }
    
if(current == true && lepton_in > 0 && nukleon_in == proton){return piplus;}

if(current == true && lepton_in < 0 && nukleon_in == neutron){return piminus;}

if(current == true && lepton_in < 0 && nukleon_in == proton)
    {
    if(funkcja_spp_3(W, nukleon_in, pizero, current)>frandom()){return pizero;}
    else{return piminus;}
    }

//NC

if(current == false &&  nukleon_in == proton)
    {
    if(funkcja_spp_3(W, nukleon_in, piplus, current)>frandom()){return piplus;}
    else{return pizero;	}
    }

if(current == false &&  nukleon_in == neutron)
    {
    if(funkcja_spp_3(W, nukleon_in, pizero, current)>frandom()){return pizero;}
    else{return piminus;}
    }


return 0;
}


int nukleon_out_(double W, int lepton_in, int nukleon_in,int meson_out, bool current)
{
//CC
    if(current == true && lepton_in >0 && nukleon_in == neutron && meson_out == piplus){return neutron;}
    if(current == true && lepton_in >0 && nukleon_in == neutron && meson_out == pizero){return proton;}

    if(current == true && lepton_in >0 && nukleon_in == proton && meson_out == piplus){return proton;}

    if(current == true && lepton_in <0 && nukleon_in == proton && meson_out == piminus){return proton;}
    if(current == true && lepton_in <0 && nukleon_in == proton && meson_out == pizero){return neutron;}

    if(current == true && lepton_in <0 && nukleon_in == neutron && meson_out == piminus){return neutron;}

//NC

    if(current == false && nukleon_in == proton && meson_out == pizero){return proton;}
    if(current == false && nukleon_in == proton && meson_out == piplus){return neutron;}
    if(current == false && nukleon_in == neutron && meson_out == pizero){return neutron;}
    if(current == false && nukleon_in == neutron && meson_out == piminus){return proton;}
    
    return 0;    
}

double funkcja_dis (double W, double W_min, double W_max, double alfa)
{

    if ((W > M12 + m_pi) && W < W_min){  return alfa * (W - M12 - m_pi) / (W_min - M12 - m_pi);}
    if (W >= W_min && W < W_max){return alfa + (1 - alfa) * (W - W_min) / (W_max - W_min);}
    if (W >= W_max){return 1;}

    return 0;		
}
		

double kombinacja (int FFset, double delta_axial_mass, double delta_C5A, double E, double W, double nu, int lepton_in, int nukleon_in, int nukleon_out, int meson_out, bool current, double alfa, double W_min, double W_max)
{

  double cros = 0;
  double spp_f = funkcja_dis(W,W_min,W_max, alfa);
  double funkcja_spp_3_ = funkcja_spp_3(W, nukleon_in, meson_out, current);

  double cros_dis = cr_sec_dis (E, W, nu, lepton_in, nukleon_in, current);

  double cros_delta = cr_sec_delta (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, nukleon_out, meson_out, current);


if(funkcja_spp_3_==0)
    {
      cros = cros_dis * spp_f + cros_delta * (1-spp_f);
      cr_dis =cros_dis*spp_f;
      cr_res =cros_delta * (1-spp_f);
    }
else
    {
     cros = cros_dis * spp_f + cros_delta * (1-spp_f)/(funkcja_spp_3_);
     cr_dis =cros_dis*spp_f;
     cr_res =cros_delta * (1-spp_f)/(funkcja_spp_3_);

    }

//cout<<"tutuaj="<<cr_dis<<" "<<cr_res<<" "<<cros<<endl;

  return cros;
}

//nowa kombinacja dis i res z produkcja stanow koncowych

double cr_dis_res (int FFset, double delta_axial_mass, double delta_C5A, double E, double W, double nu, int lepton_in, int nukleon_in, int meson_out, int nukleon_out, bool current)
{

  int ip=0;
  double W1=W/GeV;
  W_min = 1300/MeV;
  W_max = 1600/MeV;
  double przekroj_=0.0;


//Neutrino + proton
//nu + p --> l- + p + pi+
//int nukleon_out=0;
//int meson_out=0;

  if (current == true && lepton_in > 0 && nukleon_in == proton && nukleon_out == proton && meson_out == piplus)
    {
      double  alfa = 0.;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, proton, piplus,  current, alfa, W_min, W_max);
      py2ent_(&ip, &proton, &piplus, &W1);
      return przekroj_;
    }

//neutrino + neutron
//nu + n --> l- + n + pi+
//nu + n --> l- + p + pi0


  if (current == true && lepton_in > 0 && nukleon_in == neutron && nukleon_out == neutron && meson_out == piplus)
    {
    double  alfa = 0.2;
    double przekroj_plus= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, neutron, piplus, current, alfa, W_min, W_max);
    py2ent_(&ip, &neutron, &piplus, &W1);
    return przekroj_plus;
    }

  if (current == true && lepton_in > 0 && nukleon_in == neutron && nukleon_out == proton && meson_out == pizero)
    {
    double  alfa = 0.3;
    double przekroj_zero= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, proton, pizero,  current, alfa, W_min, W_max);
    py2ent_(&ip, &proton, &pizero, &W1);
    return przekroj_zero;
    }


//anu

  if (current == true && lepton_in < 0 && nukleon_in == neutron && nukleon_out == neutron && meson_out == piminus)
    {
      double  alfa = 0.0;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, neutron, piminus,  current, alfa, W_min, W_max);
      py2ent_(&ip, &neutron, &piminus, &W1);
      return przekroj_;
    }



  if (current == true && lepton_in < 0 && nukleon_in == proton && nukleon_out == proton && meson_out == piminus)
    {
      double  alfa = 0.2;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, proton, piminus,  current, alfa, W_min, W_max);
      py2ent_(&ip, &proton, &piminus, &W1);
      return przekroj_;
    }

  if (current == true && lepton_in < 0 && nukleon_in == proton && nukleon_out == neutron && meson_out == pizero)
    {
      double  alfa = 0.3;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, neutron, pizero,  current, alfa, W_min, W_max);
      py2ent_(&ip, &neutron, &pizero, &W1);
      return przekroj_;
    }



//NC

  if (current == false &&  nukleon_in == proton && nukleon_out == neutron && meson_out == piplus)
    {
      double  alfa = 0.0;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, neutron, piplus,  current, alfa, W_min, W_max);
      py2ent_(&ip, &neutron, &piplus, &W1);
      return przekroj_;
    }

  if (current == false &&  nukleon_in == proton && nukleon_out == proton && meson_out == pizero)
    {
      double  alfa = 0.0;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, proton, pizero,  current, alfa, W_min, W_max);
      py2ent_(&ip, &proton, &pizero, &W1);
      return przekroj_;
    }

  if (current == false &&  nukleon_in == neutron && nukleon_out == neutron && meson_out == pizero)
    {
      double  alfa = 0.0;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, neutron, pizero,  current, alfa, W_min, W_max);
      py2ent_(&ip, &neutron, &pizero, &W1);
      return przekroj_;
    }

  if (current == false &&  nukleon_in == neutron && nukleon_out == proton && meson_out == piminus)
    {
      double  alfa = 0.0;
      przekroj_= kombinacja (FFset, delta_axial_mass, delta_C5A, E, W, nu, lepton_in, nukleon_in, proton, piminus,  current, alfa, W_min, W_max);
      py2ent_(&ip, &proton, &piminus, &W1);
      return przekroj_;
    }







cout<<"poza!!!!!!!!!!!!"<<endl;
cout<<"W="<<W/GeV<<endl;
cout <<"N0="<<nukleon_in<<endl;
cout <<"N1="<<nukleon_out<<endl;
cout <<"ME="<<meson_out<<endl;
return 0;
}


