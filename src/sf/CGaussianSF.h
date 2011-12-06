#ifndef C_GAUSSIAN_SF_H
#define C_GAUSSIAN_SF_H

#include "CSpectralFunc.h"

#include "CMomDistrib.h"
#include "GConstants.h"
#include "generatormt.h"

//Gaussian approximation of the spectral function



class CGaussianSF: public CSpectralFunc
{
public:
	CGaussianSF(const MomDistribs i_momDistrib,
			const IsospinOfSF i_isospin,
			TargetData atom,
//			int ZZ, int NN, double E2, double Pauli,
			double pmin, double pmax, double emin, double emax,
			const double *levP,const double * levN);


	double evalSF(const double p, const double removE) const;
	double norm() const;
	double generateE(const double p) const;


private:
	const IsospinOfSF m_isospin;
    int atomZ,atomN;
	const double m_targetBeta;
	const double m_targetE2;

	mutable double m_coeff;
	
	const double *levels;
    
};

inline double Gauss(double m,double s,const double x)
	{   static const double m_coeff=sqrt(8.0/M_PI);
		double t=(x-m)/s ;
		return m_coeff*exp( -8.0*t*t )/s;
	};

inline double GaussInt(double m,double s,double x)
	{   const double k=2*sqrt(2);
		return erf((x-m)/s*k)/2;
	};


#endif
