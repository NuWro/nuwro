#ifndef _seek_h_
#define _seek_h_

class Seek
{   
	int i;
	double a;
  public:
	Seek(const double *E,double Ek)
	{ 
		seek(E,Ek);
	}
	/// find i such that E[i]<=Ek<E[i+1]
	Seek & seek(const double E[],double Ek)
	{	 
		i=0;
		while(Ek> E[i+1]) i++;
		a=(Ek-E[i])/(E[i+1]-E[i]);
		return *this;
	}
	double val(const double t[])
	{
		return t[i]*a+t[i+1]*(1-a);
	}
	double lval(const double t[])
	{
		return t[i];
	}
		
	double val2(const double t1[],const double t2[])
	{
		return (t1[i]+t2[i])*a+(t1[i+1]+t2[i+1])*(1-a);
	}
	double lval2(const double t1[],const double t2[])
	{
		return t1[i]+t2[i];
	}
	double val3(const double t1[],const double t2[],const double t3[])
	{
		return (t1[i]+t2[i]+t3[i])*a+(t1[i+1]+t2[i+1]+t3[i+1])*(1-a);
	}	
	double lval3(const double t1[],const double t2[],const double t3[])
	{
		return t1[i]+t2[i]+t3[i];
	}	
};

#endif
