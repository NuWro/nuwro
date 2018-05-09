#ifndef DM_h_
#define DM_h_
#include <cmath>
#include <complex>


//4x4 Complex matrix class

typedef std::complex<double> comp;

class DM
{
	comp 
	a0,a1,a2,a3,
	b0,b1,b2,b3,
	c0,c1,c2,c3,
	d0,d1,d2,d3;	
	
public:

    DM(const double a, const double b, const double c, const double d, const double e=0):
		a0(e+a,0),a1(0,0),  a2(-d,0), a3(-b,c),
		b0(0,0),  b1(e+a,0),b2(-b,-c),b3(d,0),
		c0(d,0),  c1(b,-c), c2(e-a,0),c3(0,0),
		d0(b,c),  d1(-d,0), d2(0,0),  d3(e-a,0)
    {}
    
	DM(comp A0, comp A1, comp A2, comp A3, 
	   comp B0, comp B1, comp B2, comp B3, 
	   comp C0, comp C1, comp C2, comp C3, 
	   comp D0, comp D1, comp D2, comp D3):
	   a0(A0), a1(A1), a2(A2), a3(A3), 
	   b0(B0), b1(B1), b2(B2), b3(B3), 
	   c0(C0), c1(C1), c2(C2), c3(C3), 
	   d0(D0), d1(D1), d2(D2), d3(D3)
	 {}  

	DM(comp x):	
	   a0(x), a1(0), a2(0), a3(0), 
	   b0(0), b1(x), b2(0), b3(0), 
	   c0(0), c1(0), c2(x), c3(0), 
	   d0(0), d1(0), d2(0), d3(x)
	{}

	DM(double x):	
	   a0(x), a1(0), a2(0), a3(0), 
	   b0(0), b1(x), b2(0), b3(0), 
	   c0(0), c1(0), c2(x), c3(0), 
	   d0(0), d1(0), d2(0), d3(x)
	{}
	 
	inline friend DM operator+ (const DM  &m1,const DM &m2)
	{
	   return DM(
			m1.a0+m2.a0,m1.a1+m2.a1,m1.a2+m2.a2,m1.a3+m2.a3,
			m1.b0+m2.b0,m1.b1+m2.b1,m1.b2+m2.b2,m1.b3+m2.b3,
			m1.c0+m2.c0,m1.c1+m2.c1,m1.c2+m2.c2,m1.c3+m2.c3,
			m1.d0+m2.d0,m1.d1+m2.d1,m1.d2+m2.d2,m1.d3+m2.d3);
	}         

	inline friend DM operator- (const DM  &m1,const DM &m2)
	{
	   return DM(
			m1.a0-m2.a0,m1.a1-m2.a1,m1.a2-m2.a2,m1.a3-m2.a3,
			m1.b0-m2.b0,m1.b1-m2.b1,m1.b2-m2.b2,m1.b3-m2.b3,
			m1.c0-m2.c0,m1.c1-m2.c1,m1.c2-m2.c2,m1.c3-m2.c3,
			m1.d0-m2.d0,m1.d1-m2.d1,m1.d2-m2.d2,m1.d3-m2.d3);
	}         


	inline DM& operator+= (const DM& m)
	{
		return *this=*this+m;
	}


	inline DM& operator-= (const DM& m)
	{
		return *this=*this-m;
	}


inline friend DM operator* (const DM  &m1,const DM &m2)
{
    DM res=0;
    const comp *a=&m1.a0;
    const comp *b=&m2.a0;
    comp *c=&res.a0;
    for(int i=0;i<4;i++)
    for(int k=0;k<4;k++)
    if(a[4*i+k]!=0.0)
    for(int j=0;j<4;j++)
		c[4*i+j]+=a[4*i+k]*b[4*k+j];
    return res;
}
	
	//multiplication A=A*B
	inline DM& operator*= (const DM&  m)
	{
		return *this=*this*m;

	}

	//multiplication A=B*A
	inline DM& operator&= (const DM&  m)
	{
		return *this=m* *this;
	}


	inline friend DM operator* (const DM  &m1,const double &d)
	{
	   return DM(
			m1.a0*d,m1.a1*d,m1.a2*d,m1.a3*d,
			m1.b0*d,m1.b1*d,m1.b2*d,m1.b3*d,
			m1.c0*d,m1.c1*d,m1.c2*d,m1.c3*d,
			m1.d0*d,m1.d1*d,m1.d2*d,m1.d3*d);
	}

	inline friend DM operator* (const DM  &m1,const comp &d)
	{
	   return DM(
			m1.a0*d,m1.a1*d,m1.a2*d,m1.a3*d,
			m1.b0*d,m1.b1*d,m1.b2*d,m1.b3*d,
			m1.c0*d,m1.c1*d,m1.c2*d,m1.c3*d,
			m1.d0*d,m1.d1*d,m1.d2*d,m1.d3*d);
	}
	
	template <class A>
	inline DM& operator*= (A d)
	{
		return *this=*this*d;
	}


	inline friend DM operator* (double d,const DM &m1)
	{
		return m1*d;
	}

	inline friend DM operator* (comp d,const DM &m1)
	{
		return m1*d;
	}
 
	inline friend DM operator/ (const DM  &m1,double d)
	{
	    return m1*(1.0/d);
	}

	inline friend DM operator/ (const DM  &m1,comp d)
	{
	    return m1*(1.0/d);
	}

	inline DM& operator/= (double d)
	{
		return *this=*this*(1.0/d);
	}

	inline DM& operator/= (comp d)
	{
		return *this=*this*(1.0/d);
	}

	void hermit();
	void transp();
	void conj();


	inline std::complex<double> operator() (unsigned short row, unsigned short col) const
	{
		return (&a0)[4*(row%4)+col%4];
	}

	inline std::complex<double>& operator() (unsigned short row, unsigned short col)
	{
		
			return (&a0)[4*(row%4)+col%4];
	}

  	//real part of trace
	inline double Trace()
	{
		return a0.real()+b1.real()+c2.real()+d3.real();

	}

	//complex trace
	inline std::complex<double> CTrace()
	{
		return a0+b1+c2+d3;

	}
	
	inline std::complex<double> BTrace(const DM & m2)
	{
		return a0*m2.a0 + a1*m2.b0
		     + b0*m2.a1 + b1*m2.b1
		     + c0*m2.a2 + c1*m2.b2
		     + d0*m2.a3 + d1*m2.b3;

	}
	
	friend std::ostream& operator<<(std::ostream& os, const DM& mat)
	{
		os<<mat.a0 <<"   "<<mat.a1 <<"   "<<mat.a2 <<"   "<<mat.a3 <<"   "<<std::endl;
		os<<mat.b0 <<"   "<<mat.b1 <<"   "<<mat.b2 <<"   "<<mat.b3 <<"   "<<std::endl;
		os<<mat.c0 <<"   "<<mat.c1 <<"   "<<mat.c2 <<"   "<<mat.c3 <<"   "<<std::endl;
		os<<mat.d0 <<"   "<<mat.d1 <<"   "<<mat.d2 <<"   "<<mat.d3 <<"   "<<std::endl;
		return os;
	}
};

DM Gamma(unsigned short i);

#include "D4V.h"
#include "D4T.h"

#endif
