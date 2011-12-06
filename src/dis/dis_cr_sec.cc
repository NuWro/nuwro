#include "dis_cc_proton.h"
#include "dis_cc_neutron.h"
#include "dis_nc.h"
#include "LeptonMass.h"
#include "pdg_name.h"
#include "parameters.h"
#include "dis_cr_sec.h"
#include <cstdlib>
#include "generatormt.h"
//////////////////////////////decomposition //////////////////////////

double hadr()
{

//double z = 0.25 + 0.75*frandom();
double z = frandom();

double a = 0.77;
double P_z = (1-a+3*a*kwad(1-z))/(1+2*a);

return P_z;

}

double prob_fey()
{
    double a = 0.77;
    double z = frandom();
    double s = z - a*z + 3*a*(z-z*z+z*z*z/3.);
    return s;
}
//cross section for nukleons with bodek corrections
//double cr_sec_dis_temp(double E, double W, double nu, int lepton, int nukleon, bool current)

double cr_sec_dis(double E, double W, double nu, int lepton, int nukleon, bool current)
{

double wynik;
double m;

if(current == true && lepton>0)
{

m = lepton_mass(abs(lepton),current);

    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_cc_nu_p (E, W, nu, m);
	break;
	case neutron:
	wynik =  cr_sec_cc_nu_n (E, W, nu, m);
	break;
    }
}

if(current == true && lepton<0)
{
m = lepton_mass(abs(lepton),current);

    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_cc_anu_p (E, W, nu, m);

	break;
	case neutron:
	wynik =  cr_sec_cc_anu_n (E, W, nu, m);
	
	break;
    }
}

if(current == false && lepton > 0)
{
 m = 0;
    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_nc_nu_p_grv (E, W, nu, m);
//	cout<<"current==false, nu, proton"<<endl;    
	break;
	case neutron:
	wynik =  cr_sec_nc_nu_n_grv (E, W, nu, m);
//	cout<<"current==false, nu, neutron"<<endl;    
	break;
    }
}

if(current == false && lepton<0)
{
 m = 0;
    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_nc_anu_p_grv (E, W, nu, m);
	break;
	case neutron:
	wynik =  cr_sec_nc_anu_n_grv (E, W, nu, m);
	break;
    }
}


return wynik/cm2*1e38 ;
//return wynik ;

}

///////////////////////////////only grv////////////////////////////////
//double cr_sec_dis(double E, double W, double nu, int lepton, int nukleon, bool current)
double cr_sec_dis_grv(double E, double W, double nu, int lepton, int nukleon, bool current)
{

double wynik;
double m;

if(current == true && lepton>0)
{

m = lepton_mass(abs(lepton),current);

    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_cc_nu_p_grv (E, W, nu, m);
	break;
	case neutron:
	wynik =  cr_sec_cc_nu_n_grv (E, W, nu, m);
	break;
    }
}

if(current == true && lepton<0)
{
m = lepton_mass(abs(lepton),current);

    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_cc_anu_p_grv (E, W, nu, m);

	break;
	case neutron:
	wynik =  cr_sec_cc_anu_n_grv (E, W, nu, m);
	
	break;
    }
}

if(current == false && lepton > 0)
{
 m = 0;
    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_nc_nu_p_grv (E, W, nu, m);
//	cout<<"current==false, nu, proton"<<endl;    
	break;
	case neutron:
	wynik =  cr_sec_nc_nu_n_grv (E, W, nu, m);
//	cout<<"current==false, nu, neutron"<<endl;    
	break;
    }
}

if(current == false && lepton<0)
{
 m = 0;
    switch(nukleon)
    {
	case proton:
	wynik =  cr_sec_nc_anu_p_grv (E, W, nu, m);
	break;
	case neutron:
	wynik =  cr_sec_nc_anu_n_grv (E, W, nu, m);
	break;
    }
}


return wynik/cm2*1e38 ;
//return wynik ;

}


