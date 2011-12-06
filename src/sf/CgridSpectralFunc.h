#ifndef C_GRID_SPECTRAL_FUNC_H
#define C_GRID_SPECTRAL_FUNC_H

#include "CSpectralFunc.h"

#include <fstream>

#include "CMomDistrib.h"
#include "gridfun2d.h"

//Spectral function normalized to 1
//value at boundries is zero
using namespace std;

class CgridSpectralFunc: public CSpectralFunc
{
public:
	CgridSpectralFunc(const char* gridfile,CMomDistrib* md, double normal,double pauli,double Beta):
    m_normal(normal)

    {
	m_pBlock=pauli;
	m_targetAlpha=alpha;
    m_momDistrib=md;
    
	f.load(gridfile);
	m_eStart=f.Emin();
	m_eStop=f.Emax();
	m_pStart=f.Pmin();
	m_pStop=f.Pmax();
    set_targetAlpha(Beta);
    }

    

    double evalSF(const double p, const double e) const
    {
	return m_normal*f.value(p,e);
    }
    
    double generateE(const double p) const
    {
	return f.generateE(p);
    }

private:
	
	const double m_normal;
	gridfun2d f;

};

#endif
