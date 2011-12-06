#ifndef C_MOM_DISTRIB_H
#define C_MOM_DISTRIB_H

#include <fstream>
#include <iostream>
#include "GConstants.h"
#include "generatormt.h"

class CMomDistrib
{
public:
	virtual ~CMomDistrib() {};	
	virtual double MF(const double p) const = 0;
	virtual double Corr(const double p) const = 0;
	virtual double Tot(const double p) const = 0;
	virtual double generate() const = 0;
	virtual double  TotMean() const=0;
		
};

CMomDistrib* createDistrib(const MomDistribs i_md, const IsospinOfSF i_isospin);

#endif
