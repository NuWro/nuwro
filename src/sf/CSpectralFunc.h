#ifndef C_SPECTRAL_FUNC_H
#define C_SPECTRAL_FUNC_H

#include "GConstants.h"
#include "generatormt.h"
#include "CMomDistrib.h"

class CSpectralFunc
{
public:
 
	CSpectralFunc():m_momDistrib(0){}
	virtual ~CSpectralFunc(){delete m_momDistrib;};
	
	virtual double evalSF(const double p, const double removE) const = 0;
	virtual double generateE(const double p) const = 0;
	inline CMomDistrib* MomDist() const;
	
	virtual inline double get_pBlock() const;
	virtual inline double get_pMin() const;
	virtual inline double get_pMax() const;
	virtual inline double get_eMin() const;
	virtual inline double get_eMax() const;
	virtual inline double get_targetAlpha() const;
    void inline set_targetAlpha(double m_targetBeta);
	
	CMomDistrib * m_momDistrib;
	double m_eStart;
	double m_pStart;
	double m_eStop;
	double m_pStop;
	double m_pBlock;
	double m_targetAlpha;
};

CMomDistrib* CSpectralFunc::MomDist() const
{ 
	return m_momDistrib;
}


double CSpectralFunc::get_pBlock() const
{
	return m_pBlock;
}


double CSpectralFunc::get_pMin() const
{
	return m_pStart;
}


double CSpectralFunc::get_pMax() const
{
	return m_pStop;
}


double CSpectralFunc::get_eMin() const
{
	return m_eStart;
}


double CSpectralFunc::get_eMax() const
{
	return m_eStop;
}


double CSpectralFunc::get_targetAlpha() const
{
	return m_targetAlpha;
}

void CSpectralFunc::set_targetAlpha(double m_targetBeta)
{		       
	m_targetAlpha = f34(m_targetBeta,MomDist()->TotMean());
}



CSpectralFunc* createSF(const TargetElement i_target,
					    const MomDistribs i_momDistrib,
					    const IsospinOfSF i_isospin);

#endif
