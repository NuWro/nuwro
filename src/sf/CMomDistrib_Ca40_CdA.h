#ifndef C_MOM_DISTRIB_CA40_CdA_H
#define C_MOM_DISTRIB_CA40_CdA_H

#include "GConstants.h"
#include "CLorentzDistrib.h"

//Ciofi degli Atti
//Momentum distribution normalized to 1 at interval [0:5 fm^{-1}]
// notation: k = momentum in fm^-1
//			 p = momentum in MeV 
class CMomDistrib_Ca40_CdA: public CMomDistrib
{
public:
	CMomDistrib_Ca40_CdA()
		:m_coeff( 20.0/20.0101/(4*pi/fermi3) ),
		 m_lorCenterV( 168.0 ),
		 m_lorHalfWidth( 80.0 ),
		 m_lorCoeff( 0.0385 )
	{
		m_lorentz = new CLorentzDistrib(m_lorCenterV, m_lorHalfWidth);
	};

	double Tot(const double p) const;
	double MF(const double p) const;
	double Corr(const double p) const;
    double TotMean()const {return 177.59 ;}
	double generate() const;

private:
	const double m_coeff;
	
	const double m_lorCenterV;
	const double m_lorHalfWidth;
	const double m_lorCoeff;

	const CLorentzDistrib *m_lorentz;
	
};



double CMomDistrib_Ca40_CdA::Tot(const double p) const
{

	const double k2(p*p*fermi2);

	return m_coeff*( 3.24*exp(-3.72*k2)*(1.0 + 11.1*k2*k2) + 0.419*exp(-1.77*k2) + 0.0282*exp(-0.22*k2) );
}


double CMomDistrib_Ca40_CdA::MF(const double p) const
{
	const double k2(p*p*fermi2);

	return m_coeff*( 3.24*exp(-3.72*k2)*(1.0 + 11.1*k2*k2) );	
}


double CMomDistrib_Ca40_CdA::Corr(const double p) const
{
	//Corr part is 19.8% of normalization for protons and neutrons
	const double k2(p*p*fermi2);

	return m_coeff*( 0.419*exp(-1.77*k2) + 0.0282*exp(-0.22*k2) );	
}



double CMomDistrib_Ca40_CdA::generate() const
{
	while (true)
	{
		const double p( m_lorentz->generate() );

		const double k( p*fermi );
		
		if ( k<0.0 || k>5.0 )
			continue;
		
		if ( frandom()*m_lorCoeff*m_lorentz->eval(p) < p*p*Tot(p) )
			return p;
	}

	return 0.0;
}

#endif
