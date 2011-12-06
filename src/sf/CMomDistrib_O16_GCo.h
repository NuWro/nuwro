#ifndef C_MOM_DISTRIB_O16_GCo_H
#define C_MOM_DISTRIB_O16_GCo_H

#include "GConstants.h"
#include "CLorentzDistrib.h"

//G. Co'
//Momentum distribution normalized to 1 at interval [0:5 fm^{-1}]
// notation: k = momentum in fm^-1
//			 p = momentum in MeV 
class CMomDistrib_O16_GCo: public CMomDistrib
{
public:
	CMomDistrib_O16_GCo()
		:m_kRes(122),
		 m_kStep(0.03),
		 m_coeff(1.020*16.0/8.0/pow(2.0*pi,3.0)*fermi3 ),
		 m_lorCenterV(154.0),
		 m_lorHalfWidth(86.0),
		 m_lorCoeff(0.042)
	{
		m_Distrib = new double[m_kRes];
		m_lorentz = new CLorentzDistrib(m_lorCenterV, m_lorHalfWidth);


		std::ifstream inputFile( "sf/O16_GCo.dat" );
		
		if (inputFile.fail())
		{
			std::cerr<<"Indispensable file 'O16_GCo.dat' not found"<<std::endl;
			return;
		}

		double k(0), val(0.0);

		for (int iCnt=0;iCnt<m_kRes;++iCnt)
		{
			inputFile>>k>>val;
			m_Distrib[iCnt] = val;
		}

	};

	double MF(const double p) const;
	double Corr(const double p) const;
	double Tot(const double p) const;
    double TotMean()const {return  174.4;}

	double generate() const;

private:
	const int m_kRes;
	const double m_kStep;
	const double m_coeff;

	mutable double* m_Distrib;

	const double m_lorCenterV;
	const double m_lorHalfWidth;
	const double m_lorCoeff;


	const CLorentzDistrib *m_lorentz;

  	double evalDist(const double k) const;

	
};


double CMomDistrib_O16_GCo::evalDist(const double k) const
{
	const double k2(k*k);

	if ( k >= 3.5850 )
		return 1.021*m_coeff*( 0.1678*exp(-0.241*k2) );

	const int nk( (int)floor((k+0.5*m_kStep)/m_kStep) ); 
	
	const double kR( k-m_kStep*(nk-0.5) );
	
	const double c0( m_Distrib[nk]   );
	const double c1( m_Distrib[nk+1] );

	double f( (c1 - c0)*kR/m_kStep + c0 );

	return m_coeff*f;
}

double CMomDistrib_O16_GCo::MF(const double p) const
{
	//average: 174.4 MeV
	//const double oxygenAlpha(1.028*fermi2) );

	const double k(p*fermi);
	const double k2(k*k);

	if (k >= 2.0250)
		return 0.0;

	return evalDist(k) - m_coeff*( 2.128*exp(-1.40*k2) + 0.1427*exp(-0.226*k2) );
	
}

double CMomDistrib_O16_GCo::Corr(const double p) const
{
	//Corr part is 12.0% of normalization for protons and neutrons
	const double k(p*fermi);
	const double k2(k*k);

	if (k <= 2.0250)
		return m_coeff*( 2.128*exp(-1.40*k2) + 0.1427*exp(-0.226*k2) );

	return evalDist(k);

	
}


double CMomDistrib_O16_GCo::Tot(const double p) const
{
	return evalDist(p*fermi);
}

double CMomDistrib_O16_GCo::generate() const
{
	while (true)
	{
		const double p( m_lorentz->generate() );

		const double k( p*fermi );
		
		if ( k<0.0 || k>5.0 )
			continue;
		
		if ( frandom()*m_lorCoeff*m_lorentz->eval(p) < p*p*evalDist(k) )
			return p;
	}
}

#endif
