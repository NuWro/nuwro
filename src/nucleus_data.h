#ifndef _densities_h_
#define _densities_h_
#include <cstdlib>
#include "generatormt.h"
#include "jednostki.h"
#include "elements.h"

double Meff(double kf);                  ///< efective Mass from mean Fermi momentum

class nucleus_data
{
	private:
	
	int _p;
	int _n;
	double _r;
	double (*dens_fun)(double[],double);
	double *dens_params;
	double _max_rr_dens;
	double _kF;
	double _Mf;
	unsigned char _model_tag;            

    public:
   
    nucleus_data(int p0, 
                 int n0, 
		 double (*fun)(double[],double),
		 double * par=NULL,
                 unsigned char model_tag = 0
		        )
    : _p(p0), _n(n0), _r(0), dens_fun(fun), dens_params(par),
      _max_rr_dens(0), _kF(0), _Mf(0), _model_tag(model_tag)
	{}
	
	const  char* name();
	int p()		{return _p;}
	int n()		{return _n;}
	double A()	{return _p+_n;}
	double r();
	double random_r();
	double dens(double r);
	double kF();
	double Mf();
	unsigned char model_tag() const { return _model_tag; } 

	private:
	double kf_helper(double r)
	{   
		double den=dens(r);
		return cbrt(3*Pi*Pi*den/2)*r*r*den;
	}
	double dens_helper(double r)
	{
		return r*r*dens(r);
	}
/*	double mf_helper(double r)
	{
		double den=dens(r);
		return Meff(cbrt(3*Pi*Pi*den/2))*r*r*den;
	}
*/
};

extern nucleus_data dens_data[] ;

nucleus_data* best_data(int p, int n);
// model-aware selector
nucleus_data* best_data(int p, int n, unsigned char model_tag);

double density(double r,int p, int n);

inline double FermiMomentum(double density)
{
	return cbrt(1.5*Pi*Pi*density);
}

#endif
