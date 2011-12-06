#ifndef _gridFun2d_h_
#define _gridFun2d_h_

class gridfun2d
{
public:
    gridfun2d():table(0),erow(0),prow(0){} 
    ~gridfun2d(){delete []table; delete [] erow; delete []prow;} 
	
	bool load(const char*);

	double value(const double p, const double removE) const;
	double generateE(const double p) const;
	double generateP() const;
	
    double Pmin(){return pMin;}
    double Pmax(){return pMax;}
    double Emin(){return eMin;}
    double Emax(){return eMax;}

private:
	int eRes;
	double eMin;
	double eMax;

	int pRes;
	double pMin;
	double pMax;

	double *table;
	double *erow;
	double *prow;
	
    double val(int i, int j) const
	{if (0<=i && i<pRes && 0<=j && j<eRes)
	   return table[i*eRes+j];
	 return 0;
	}   
};

#endif
