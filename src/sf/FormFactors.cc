#include "FormFactors.h"
#include "params.h"
extern NUWRO::params *p1;
double FA(const double q2) //axial form factor
{ 
	double MA=p1->qel_cc_axial_mass;
	double MA2=MA*MA;
	double a = 1.0/(1.0 - q2/MA2);
	return gA*a*a; 
}

double FP(const double q2) //pseudoscalar form factor
{
	return 2.0*M2*FA(q2)/(piMass2 - q2); 
}


//-------------------------------------Dipole-------------------------------------------------------
double DipoleGE(const double q2) //dipole electric form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	a = a*a;

	return 1.0/a;
}

double DipoleGM(const double q2) //dipole magnetic form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	a = a*a;

	return (mu_p - mu_n)/a;
}

double ProtonDipoleGE(const double q2) //dipole proton electric form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	a = a*a;

	return 1.0/a;
}

double ProtonDipoleGM(const double q2) //dipole proton magnetic form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	a = a*a;

	return mu_p/a;
}

double NeutronDipoleGE(const double q2) //dipole neutron electric form factor G_E^V
{
	return 0.0;
}

double NeutronDipoleGM(const double q2) //dipole neutron magnetic form factor G_E^V
{
	double a = 1.0 - q2/MV2;
	a = a*a;

	return mu_n/a;
}

//-------------------BBBA_05------------------------------------------------------------------------

double BBBA05_GE(const double q2) //electic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);

	const double GEp = ( 1.0-0.0578*tau )/( 1.0 + (11.1+(13.6+33.0*tau)*tau)*tau );
	
	const double GEn = ( 1.25+1.30*tau )*tau/( 1.0 +(-9.86+(305.0+(-758.0+802.0*tau)*tau)*tau)*tau );

	return (GEp - GEn);
}

double BBBA05_GM(const double q2) //magnetic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);

	const double GMp = ( 1.0+0.15*tau )/( 1.0 + (11.1+(19.6+7.54*tau)*tau)*tau );
	
	const double GMn = ( 1.0+1.81*tau )/( 1.0 + (14.1+(20.7+68.7*tau)*tau)*tau );

	return (mu_p*GMp - mu_n*GMn);
}

double ProtonBBBA05_GE(const double q2) //proton electic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);

	const double GEp = ( 1.0-0.0578*tau )/( 1.0 + (11.1+(13.6+33.0*tau)*tau)*tau );
	
	return GEp;
}

double ProtonBBBA05_GM(const double q2) //proton magnetic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);

	const double GMp = ( 1.0+0.15*tau )/( 1.0 + (11.1+(19.6+7.54*tau)*tau)*tau );

	return mu_p*GMp;
}

double NeutronBBBA05_GE(const double q2) //neutron electic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);

	const double GEn = ( 1.25+1.30*tau )*tau/( 1.0 +(-9.86+(305.0+(-758.0+802.0*tau)*tau)*tau)*tau );

	return GEn;
}

double NeutronBBBA05_GM(const double q2) //neutron magnetic form factor
{
	//arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV

	const double tau = -q2/(4*M2);
	
	const double GMn = ( 1.0+1.81*tau )/( 1.0 + (14.1+(20.7+68.7*tau)*tau)*tau );

	return mu_n*GMn;
}


//--------------------BBA_03------------------------------------------------------------------------

double BBA03_GE(const double q2) //electic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV

	const double Q2(-q2);
	const double tau( Q2/(4*M2) );

	double GEp = 1.0 + ( 3.253+( 1.422+( 0.08582+( 0.3318+(-0.09371+0.01076*Q2)*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GEp = 1.0/GEp;
	
	double GEn =  -0.942*mu_n*tau/(1 + 4.61*tau)/(1 + Q2*1.0e+6/MV2)/(1 + Q2*1.0e+6/MV2);

	return (GEp - GEn);
}

double BBA03_GM(const double q2) //magnetic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV

	const double Q2(-q2);

	double GMp= 1.0 + ( 3.104+( 1.428+( 0.1112+( -0.006981+(0.0003705-0.7063e-5*Q2)*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GMp = mu_p/GMp;
	
	double GMn = 1.0 + ( 3.043+( 0.8548+( 0.6806+( -0.1287+0.008912*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GMn = mu_n/GMn;

	return (GMp - GMn);
}

double ProtonBBA03_GE(const double q2) //proton electic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV2

	const double Q2(-q2);
	const double tau( Q2/(4*M2) );

	double GEp = 1.0 + ( 3.253+( 1.422+( 0.08582+( 0.3318+(-0.09371+0.01076*Q2)*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GEp = 1.0/GEp;
	
	return GEp;
}

double ProtonBBA03_GM(const double q2) //proton magnetic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV

	const double Q2(-q2);

	double GMp= 1.0 + ( 3.104+( 1.428+( 0.1112+( -0.006981+(0.0003705-0.7063e-5*Q2)*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GMp = mu_p/GMp;
	
	return GMp;
}

double NeutronBBA03_GE(const double q2) //neutron electic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV

	const double Q2(-q2);
	const double tau( Q2/(4*M2) );

	const double GEn =  -0.942*mu_n*tau/(1 + 4.61*tau)/(1 + Q2*1.0e+6/MV2)/(1 + Q2*1.0e+6/MV2);

	return GEn;
}

double NeutronBBA03_GM(const double q2) //neutron magnetic form factor
{
	//hep-ex/0308005 BBA-2003 for Q2<6 GeV

	const double Q2(-q2);

	double GMn = 1.0 + ( 3.043+( 0.8548+( 0.6806+( -0.1287+0.008912*Q2 )*Q2 )*Q2 )*Q2 )*Q2;
	GMn = mu_n/GMn;

	return GMn;
}


//--------------------JLab---------------------------------------------------------------------------

double JLabGE(const double q2) //electic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	
	//PHYSICAL REVIEW C, VOLUME 65, 051001(R)
	double GEp = 1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2;
	GEp = (1.0 - 0.13*(Q2 - 0.04))/GEp;
	
	//PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
	double GEn =  -1.25*mu_n*tau/(1 + 18.3*tau)/(1 + Q2*1.0e+6/MV2)/(1 + Q2*1.0e+6/MV2);

	return (GEp - GEn);	
}

double JLabGM(const double q2) //magnetic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	//PHYSICAL REVIEW C, VOLUME 65, 051001(R)
	double GMp = 1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2;
	GMp =  mu_p/GMp; 

	//PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
	double GMn = 1.0 - (1.74 + 7.63*Q2)*Q + (9.29 + 4.63*Q2)*Q2;
	GMn =  mu_n/GMn;

	return (GMp - GMn);
}

double ProtonJLabGE(const double q2) //proton electic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	
	//PHYSICAL REVIEW C, VOLUME 65, 051001(R)
	double GEp = 1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2;
	GEp = (1.0 - 0.13*(Q2 - 0.04))/GEp;

	return GEp;	
}

double ProtonJLabGM(const double q2) //proton magnetic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	//PHYSICAL REVIEW C, VOLUME 65, 051001(R)
	double GMp = 1.0 + (0.116 + (0.241 + 0.345*Q2)*Q2)*Q + (2.874 + 1.006*Q2)*Q2;
	GMp =  mu_p/GMp; 

	return GMp;
}

double NeutronJLabGE(const double q2) //neutron electic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	//PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
	double GEn =  -1.25*mu_n*tau/(1 + 18.3*tau)/(1 + Q2*1.0e+6/MV2)/(1 + Q2*1.0e+6/MV2);

	return GEn;	
}

double NeutronJLabGM(const double q2) //neutron magnetic form factor
{
	const double Q2(-q2);
	const double Q (sqrt(Q2)); 
	const double tau( Q2/(4*M2) );
	
	//PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
	double GMn = 1.0 - (1.74 + 7.63*Q2)*Q + (9.29 + 4.63*Q2)*Q2;
	GMn =  mu_n/GMn;

	return GMn;
}
