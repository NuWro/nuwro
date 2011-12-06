#ifndef C_MOM_DISTRIB_CA40_GCo_H
#define C_MOM_DISTRIB_CA40_GCo_H

#include "GConstants.h"
#include "CLorentzDistrib.h"

//G. Co'
//Momentum distribution normalized to 1 at interval [0:5 fm^{-1}]
// notation: k = momentum in fm^-1
//			 p = momentum in MeV 
class CMomDistrib_Ca40_GCo: public CMomDistrib
{
public:
	CMomDistrib_Ca40_GCo(const IsospinOfSF i_isospin)
		:m_kRes(122),
		 m_kStep(0.03),
		 m_isospin(i_isospin),
		 m_coeff(1.037*40.0/20.0/pow(2.0*pi,3.0)*fermi3 ),
		 m_lorCenterV( (m_isospin == proton) ? 179.0 : 177.0 ),
		 m_lorHalfWidth( 90.0 ),
		 m_lorCoeff( 0.0423 )
	{
		m_Distrib = new double[m_kRes];
		m_lorentz = new CLorentzDistrib(m_lorCenterV, m_lorHalfWidth);
         
        const char * fname=(m_isospin == proton) ? "sf/Ca40P_GCo.dat" : "sf/Ca40N_GCo.dat" ;
		std::ifstream inputFile;
		open_data_file(inputFile,fname);
		
		if (inputFile.fail())
		{
			std::cerr<<"Indispensable file '"<<fname<<"' not found"<<std::endl;
			return;
		}

		double k(0), val(0.0);

		for (int iCnt=0;iCnt<m_kRes;++iCnt)
		{
			inputFile>>k>>val;
			m_Distrib[iCnt] = val;
		}

	};

	double Tot(const double p) const;
	double MF(const double p) const;
	double Corr(const double p) const;
	double TotMean()const {return   (m_isospin == proton) ? 189.07: 187.01;}


	double generate() const;

private:
	const int m_kRes;
	const double m_kStep;
	const IsospinOfSF m_isospin;
	const double m_coeff;

	mutable double* m_Distrib;

	const double m_lorCenterV;
	const double m_lorHalfWidth;
	const double m_lorCoeff;

	const CLorentzDistrib *m_lorentz;

};



double CMomDistrib_Ca40_GCo::Tot(const double p) const
{
	const double k(p*fermi);
	const double k2(k*k);

	if ( k >= 3.5850 )
		return m_coeff*( (m_isospin==proton)*0.2282*exp(-0.2580*k2) 
		              + (m_isospin==neutron)*0.2451*exp(-0.2648*k2) );

	const int nk( static_cast<int>((k+0.5*m_kStep)/m_kStep) ); 
	
	const double kR( k-m_kStep*(nk-0.5) );
	
	const double c0( m_Distrib[nk]   );
	const double c1( m_Distrib[nk+1] );

	double f( (c1 - c0)*kR/m_kStep + c0 );

	return m_coeff*f;
}


double CMomDistrib_Ca40_GCo::MF(const double p) const
{
	//const double calciumAlpha(0.8382*fermi2);//GCo 40Ca protons  (momentum: 189.1 MeV)
	//const double calciumAlpha(0.8560*fermi2);//GCo 40Ca neutrons (momentum: 187.1 MeV)

	const double k(p*fermi);
	const double k2(k*k);

	if (k >= 2.0250)
		return 0.0;

	return Tot(p) - m_coeff*( (m_isospin==proton)*( 4.2150*exp(-1.77*k2) + 0.1940*exp(-0.226*k2) ) 
								+ (m_isospin==neutron)*( 4.270*exp(-1.77*k2) + 0.1855*exp(-0.2142*k2) ) );
	
}

double CMomDistrib_Ca40_GCo::Corr(const double p) const
{
	//Corr part is 16.2% of normalization for protons and neutrons
	const double k(p*fermi);
	const double k2(k*k);

	if (k <= 2.0250)
		return m_coeff*( (m_isospin==proton)*( 4.2150*exp(-1.77*k2) + 0.1940*exp(-0.226*k2) ) 
						+ (m_isospin==neutron)*(4.270*exp(-1.77*k2) + 0.1855*exp(-0.2142*k2) ) );

	return Tot(p);

	
}



double CMomDistrib_Ca40_GCo::generate() const
{
	while (true)
	{
		const double p( m_lorentz->generate() );

		const double k( p*fermi );
		
		if ( k<0.0 || k>5.0 )
			continue;
		
		if ( frandom()*m_lorCoeff*m_lorentz->eval(p) < p*p*Tot(p) )
		{
			return p;
		}
	}

	return 0.0;
}


#endif
