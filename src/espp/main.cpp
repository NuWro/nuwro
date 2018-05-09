#include "crossection.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                 CHANNELS:                          ///////////////////////////////
///////     1: e + p -> e + p + \pi^\0                                     ///////////////////////////////
///////     2: e + n -> e + p + \pi^-                                      ///////////////////////////////
///////     3: e + p -> e + n + \pi^+                                      ///////////////////////////////
///////     4: e + n -> e + n + \pi^0                                      ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;
//usually you pass this as C3V,C4V,C5V
	double vect[3];

const D4V<double> nuclmom(M,0,0,0);
D4V<double> lprime(0,0,0,0);
double En=0;

double dsigma_domega_dEprime_dcosphi( double cospi)
{
	double result=0;
	//cerr<<En<<" "<<lprime(0)<<" "<<lprime(1)<<" "<<lprime(2)<<" "<<lprime(3)<<" "<<cospi<<endl;
	result+=dsigma_dq0_dOmegal_dOmegaCMS(En,lprime, cospi, 0, nuclmom, vect);
	result+=dsigma_dq0_dOmegal_dOmegaCMS(En,lprime, cospi, 0.5*Pi, nuclmom, vect);
	result+=dsigma_dq0_dOmegal_dOmegaCMS(En,lprime, cospi, Pi, nuclmom, vect);
	result+=dsigma_dq0_dOmegal_dOmegaCMS(En,lprime, cospi, 1.5*Pi, nuclmom, vect);
	result*=0.5*Pi;
	return result;
}

void dosomecoolplots()
{
	double min =0.2*GeV;
	double step=0.001*GeV;
	double E=0.730*GeV;
	En=E;
	double cose=cos(37.1*degree);
	std::stringstream ss1;
	ss1<<"A_DELFF"<<ffset<<"_PV"<<PV<<"_PIFF"<<FP<<"_.txt"<<std::flush;
	//we compare the two procedures for nucleon at rest
	
	std::ofstream H(ss1.str().c_str());
	for(int i = 0; i<530 ; i++)
	{
		double W=min+i*step;
		double Eprime=E-W;
		lprime(0)=Eprime;
		lprime(1)=Eprime*sqrt(1.0-cose*cose);
		lprime(2)=0;
		lprime(3)=Eprime*cose;
		double a=0;
		double b=0;
		std::cerr<<"\n\n q^0 = "<<W/GeV<<"\n";
		H<<W/GeV;
		chan=1;
		a=dsigma_domega_dEprime_(E,Eprime,cose,vect)*1e33*GeV;//
		b=calg20(dsigma_domega_dEprime_dcosphi,-1,1,1)*1e33*GeV;//
		H<<" "<<a<<" "<<b<<" "<<(a-b)/a;
		chan=2;
		a=dsigma_domega_dEprime_(E,Eprime,cose,vect)*1e33*GeV;//
		b=calg20(dsigma_domega_dEprime_dcosphi,-1,1,1)*1e33*GeV;//
		H<<" "<<a<<" "<<b<<" "<<(a-b)/a;
		chan=3;
		a=dsigma_domega_dEprime_(E,Eprime,cose,vect)*1e33*GeV;//
		b=calg20(dsigma_domega_dEprime_dcosphi,-1,1,1)*1e33*GeV;
		H<<" "<<a<<" "<<b<<" "<<(a-b)/a;
		chan=4;
		a=dsigma_domega_dEprime_(E,Eprime,cose,vect)*1e33*GeV;//
		b=calg20(dsigma_domega_dEprime_dcosphi,-1,1,1)*1e33*GeV;
		H<<" "<<a<<" "<<b<<" "<<(a-b)/a;
		H<<"\n";
	}
	H.close();	
}

int main()
{
	// czynnik postaci pionu
	FP=0;
	//Sprzezenie pseudowektorowe pion-nukleon?
	PV=1;
	//which ff set 0-2 our fits,3-Olga
	ffset=2;
	dosomecoolplots();
//	std::cout<<DM(1)/2<<std::endl;
	return 0;
}
