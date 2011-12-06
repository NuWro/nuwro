#ifndef C_MOM_DISTRIB_O16_Ben_H
#define C_MOM_DISTRIB_O16_Ben_H

#include "GConstants.h"
#include "CLorentzDistrib.h"
#include "dirs.h"


//O. Benhar
//Momentum distribution normalized to 1 at interval [0:800 MeV/c]
// notation: k = momentum in fm^-1
//			 p = momentum in MeV 
class CMomDistrib_O16_Ben: public CMomDistrib
{
public:
	CMomDistrib_O16_Ben()
		:m_pRes(202),
		 m_pStep(4.0),
		 m_coeff( 1.0/(4.0*pi/fermi3)  ),
		 m_lorCenterV(147.0),
		 m_lorHalfWidth(90.0),
		 m_lorCoeff(0.042)
	{
		m_Distrib = new double[m_pRes];
		m_lorentz = new CLorentzDistrib(m_lorCenterV, m_lorHalfWidth);

		std::ifstream inputFile;
		open_data_file(inputFile, "sf/O16_Ben.dat" );
		
		if (inputFile.fail())
		{
			std::cerr<<"Indispensable file 'O16_Ben.dat' not found"<<std::endl;
			return;
		}

		double p(0), val(0.0);

		for (int iCnt(0); iCnt<m_pRes; ++iCnt)
		{
			inputFile>>p>>val;
			m_Distrib[iCnt] = val*(4*pi);
		}

	};

	double MF(const double p) const;
	double Corr(const double p) const;
	double Tot(const double p) const;
    double TotMean() const {return 167.83 ;}
	double generate() const;


private:
	const int m_pRes;
	const double m_pStep;
	const double m_coeff;

	mutable double* m_Distrib;
	
	const double m_lorCenterV;
	const double m_lorHalfWidth;
	const double m_lorCoeff;


	const CLorentzDistrib *m_lorentz;

	double evalDist(const double p) const;
	
};


double CMomDistrib_O16_Ben::evalDist(const double p) const
{
	//p = m_pStep*(np-0.5) + pR
	if ( p >= 798.0 )
		return m_coeff*0.085*exp(-0.37*p*p*fermi2);

	const int np( (int)floor((p+0.5*m_pStep)/m_pStep) ); 
	
	const double pR( p-m_pStep*(np-0.5) );
	
	const double c0( m_Distrib[np]   );
	const double c1( m_Distrib[np+1] );

	double f( (c1 - c0)*pR/m_pStep + c0 );

	return m_coeff*f;
}

double CMomDistrib_O16_Ben::MF(const double p) const
{
	//av momentum 167.83 MeV/c
	const double k2(p*p*fermi2);

	if ( p <= 360.0 )
		return evalDist(p) - m_coeff*( 0.326*exp(-1.0*k2) + 0.0263*exp(-0.22*k2) );

	return  m_coeff*( 2.74*exp(-3.3*k2)*(1.0 + 6.66*k2) );
	
}

double CMomDistrib_O16_Ben::Corr(const double p) const
{

	//Corr part is 24.65% of normalization for protons and neutrons
	const double k2(p*p*fermi2);

	if ( p <= 360.0 )
		return m_coeff*( 0.326*exp(-1.0*k2) + 0.0263*exp(-0.22*k2) );
			
	return  evalDist(p) - m_coeff*( 2.74*exp(-3.3*k2)*(1.0 + 6.66*k2) ); 

}

double CMomDistrib_O16_Ben::Tot(const double p) const
{
	return evalDist(p);
}


double CMomDistrib_O16_Ben::generate() const
{
	while (true)
	{
		const double p( m_lorentz->generate() );
		
		if ( p<0.0 || p>800.0 )
			continue;

		if ( frandom()*m_lorCoeff*m_lorentz->eval(p) < p*p*evalDist(p) )
			return p;
	}
}

#endif
