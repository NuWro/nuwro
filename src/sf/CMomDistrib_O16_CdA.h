#ifndef C_MOM_DISTRIB_O16_CdA_H
#define C_MOM_DISTRIB_O16_CdA_H

#include "GConstants.h"
#include "CLorentzDistrib.h"

//Ciofi degli Atti
//Momentum distribution normalized to 1 at interval [0:5 fm^{-1}]
// notation: k = momentum in fm^-1
//			 p = momentum in MeV 
class CMomDistrib_O16_CdA: public CMomDistrib
{
public:
	CMomDistrib_O16_CdA()
		:m_coeff( 1.0021/(4*pi/fermi3) ),
		 m_lorCenterV(146.0),
		 m_lorHalfWidth(80.0),
		 m_lorCoeff(0.040)
	{
		m_lorentz = new CLorentzDistrib(m_lorCenterV, m_lorHalfWidth);
	};

	double MF(const double p) const;
	double Corr(const double p) const;
	double Tot(const double p) const;
    double TotMean()const {return  162.121;}
	double generate() const;

private:
	const double m_coeff;
	
	const double m_lorCenterV;
	const double m_lorHalfWidth;
	const double m_lorCoeff;

	const CLorentzDistrib *m_lorentz;
	
};


double CMomDistrib_O16_CdA::MF(const double p) const
{
	const double k2(p*p*fermi2);

	return m_coeff*( 2.74*exp(-3.33*k2)*(1.0 + 6.66*k2) );
	
}

double CMomDistrib_O16_CdA::Corr(const double p) const
{
	//Corr part is 19.92% of normalization for protons and neutrons
	const double k2(p*p*fermi2);

	return m_coeff*( 0.326*exp(-1.40*k2) + 0.0263*exp(-0.22*k2) );

	
}

double CMomDistrib_O16_CdA::Tot(const double p) const
{
	const double k2(p*p*fermi2);

	return m_coeff*( 2.74*exp(-3.33*k2)*(1.0 + 6.66*k2) + 0.326*exp(-1.40*k2) + 0.0263*exp(-0.22*k2) );

}

double CMomDistrib_O16_CdA::generate() const
{
	while (true)
	{
		const double p( m_lorentz->generate() );

		const double k( p*fermi);
		
		if ( k<0.0 || k>5.0 )
			continue;
		
		if ( frandom()*m_lorCoeff*m_lorentz->eval(p) < p*p*Tot(p) )
			return p;
	}
}

#endif
