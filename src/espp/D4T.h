//generic 4x4 Lorentz tensor

#ifndef _D4T_H_
#define _D4T_H_
#include "D4V.h"

const int dim=4;
template< class T >
class D4T
{
		D4V<T> t0_;
		D4V<T> t1_;
		D4V<T> t2_;
		D4V<T> t3_;
	public:
		//constructor as a tensor product of two vectors
		D4T(const D4V<T> &t,const D4V<T> &u):   t0_(t(0)*u(0), t(0)*u(1), t(0)*u(2), t(0)*u(3)),
							t1_(t(1)*u(0), t(1)*u(1), t(1)*u(2), t(1)*u(3)),
							t2_(t(2)*u(0), t(2)*u(1), t(2)*u(2), t(2)*u(3)),
							t3_(t(3)*u(0), t(3)*u(1), t(3)*u(2), t(3)*u(3)){}
		//element-by element constructor
		D4T(const T &t00,const T &t01, const T &t02, const T &t03,
		    const T &t10,const T &t11, const T &t12, const T &t13,
		    const T &t20,const T &t21, const T &t22, const T &t23,
		    const T &t30,const T &t31, const T &t32, const T &t33):
				t0_( t00, t01, t02, t03),
				t1_( t10, t11, t12, t13),
				t2_( t20, t21, t22, t23),
				t3_( t30, t31, t32, t33) {}
	


	inline T operator() (unsigned short i, unsigned short j) const
	{
		//return (&t0_)[i](j);
		switch(i)
		{
			case 0:	return t0_(j);
			case 1:	return t1_(j);
			case 2:	return t2_(j);
			case 3:	return t3_(j);
			default:
			{
				double failsafe=0;
				return t1_(0)*failsafe;
			}
		}
	}

	inline T& operator() (unsigned short i, unsigned short j)
	{
		//return (&t0_)[i](j);
		switch(i)
		{
			case 0:	return t0_(j);
			case 1:	return t1_(j);
			case 2:	return t2_(j);
			case 3:	return t3_(j);
			default:return t1_(0);
		}
	}


	inline D4V<T> operator() (unsigned short i)
	{
		//return (&t0_)[i];
		switch(i)
		{
			case 0:	return t0_;
			case 1:	return t1_;
			case 2:	return t2_;
			case 3:	return t3_;
			default:
			{
				double failsafe=0;
				return t1_*failsafe;
			}
		}
	}
	
//transposition
	void transp()
	{
		T tmp=(*this)(0,0);
		for (int a=1; a<4; a++)
			for (int b=0; b<a; b++)
			{
				tmp=(*this)(a,b);
				(*this)(a,b)=(*this)(b,a);
				(*this)(b,a)=tmp;
			}
		
	}

inline void hermit()
	{
		T tmp=(*this)(0,0);
		for (int a=1; a<4; a++)
			for (int b=0; b<a; b++)
			{
				tmp=(*this)(a,b);
				tmp.hermit();
				(*this)(a,b)=(*this)(b,a);
				(*this)(a,b).hermit();
				(*this)(b,a)=tmp;
			}
		(*this)(0,0).hermit();
		(*this)(1,1).hermit();
		(*this)(2,2).hermit();
		(*this)(3,3).hermit();
		/*for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b).hermit();
			}*/
	
}


	//funny addition procedures: imagine you have a tensor ~g^{\mu\nu} and want to add a piece ~p^\mu p^\nu without
	//creating a temporary tensor. you .Add(p^\mu*something,p^\nu)	Saves processor time
	template <class A,class B>
	inline void Add(const D4V<A> &t,const D4V<B> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	template <class A>
	inline void Add(const D4V<A> &t,const D4V<double> &u)
	{
		for (unsigned short b=0; b<dim; b++)
		if(u(b)!=0)
		for (unsigned short a=0; a<dim; a++)
		{
					(*this)(a,b)+=t(a)*u(b);
		}
	}
    
	
	template <class A,class B>
	inline void Replace(const D4V<A> &t,const D4V<B> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)=t(a)*u(b);
		}
	}



	inline friend D4T<T> operator+ (const D4T<T> &t1,const D4T<T> &t2)
	{
		return D4T<T>(

					t1.t0_ +t2.t0_,
					t1.t1_ +t2.t1_,
					t1.t2_ +t2.t2_,
					t1.t3_ +t2.t3_);
			
	}


	inline D4T<T>& operator+= (const D4T<T> &t1)
	{
		t0_ += t1.t0_;
		t1_ += t1.t1_;
		t2_ += t1.t2_;
		t3_ += t1.t3_;
		return *this;

	}

	inline friend D4T<T> operator- (const D4T<T> &t1,const D4T<T> &t2)
	{
		return D4T<T>(
					t1.t0_ -t2.t0_,
					t1.t1_ -t2.t1_,
					t1.t2_ -t2.t2_,
					t1.t3_ -t2.t3_);

	}

	inline D4T<T>& operator-= (const D4T<T> &t1)
	{
		t0_ -= t1.t0_;
		t1_ -= t1.t1_;
		t2_ -= t1.t2_;
		t3_ -= t1.t3_;
		return *this;

	}

	//multiplication T1^\mu_\alpha T2^{\alpha \nu}
	inline friend D4T<T> operator* (const D4T<T> &T1,const D4T<T> &T2)
	{
		D4T<T> result(T1.t0_*0,T1.t0_*0);
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		//T2.transp();
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)=(T1(a)*tmptransp(b));
				
			}


		/*
		result.t0_.t0_=T1.t0_.t0_*T2.t0_.t0_ - T1.t0_.t1_*T2.t1_.t0_ - T1.t0_.t2_*T2.t2_.t0_ - T1.t0_.t3_*T2.t3_.t0_;
		result.t0_.t1_=T1.t0_.t0_*T2.t0_.t1_ - T1.t0_.t1_*T2.t1_.t1_ - T1.t0_.t2_*T2.t2_.t1_ - T1.t0_.t3_*T2.t3_.t1_;
		result.t0_.t2_=T1.t0_.t0_*T2.t0_.t2_ - T1.t0_.t1_*T2.t1_.t2_ - T1.t0_.t2_*T2.t2_.t2_ - T1.t0_.t3_*T2.t3_.t2_;
		result.t0_.t3_=T1.t0_.t0_*T2.t0_.t3_ - T1.t0_.t1_*T2.t1_.t3_ - T1.t0_.t2_*T2.t2_.t3_ - T1.t0_.t3_*T2.t3_.t3_;

		result.t1_.t0_=T1.t1_.t0_*T2.t0_.t0_ - T1.t1_.t1_*T2.t1_.t0_ - T1.t1_.t2_*T2.t2_.t0_ - T1.t1_.t3_*T2.t3_.t0_;
		result.t1_.t1_=T1.t1_.t0_*T2.t0_.t1_ - T1.t1_.t1_*T2.t1_.t1_ - T1.t1_.t2_*T2.t2_.t1_ - T1.t1_.t3_*T2.t3_.t1_;
		result.t1_.t2_=T1.t1_.t0_*T2.t0_.t2_ - T1.t1_.t1_*T2.t1_.t2_ - T1.t1_.t2_*T2.t2_.t2_ - T1.t1_.t3_*T2.t3_.t2_;
		result.t1_.t3_=T1.t1_.t0_*T2.t0_.t3_ - T1.t1_.t1_*T2.t1_.t3_ - T1.t1_.t2_*T2.t2_.t3_ - T1.t1_.t3_*T2.t3_.t3_;

		result.t2_.t0_=T1.t2_.t0_*T2.t0_.t0_ - T1.t2_.t1_*T2.t1_.t0_ - T1.t2_.t2_*T2.t2_.t0_ - T1.t2_.t3_*T2.t3_.t0_;
		result.t2_.t1_=T1.t2_.t0_*T2.t0_.t1_ - T1.t2_.t1_*T2.t1_.t1_ - T1.t2_.t2_*T2.t2_.t1_ - T1.t2_.t3_*T2.t3_.t1_;
		result.t2_.t2_=T1.t2_.t0_*T2.t0_.t2_ - T1.t2_.t1_*T2.t1_.t2_ - T1.t2_.t2_*T2.t2_.t2_ - T1.t2_.t3_*T2.t3_.t2_;
		result.t2_.t3_=T1.t2_.t0_*T2.t0_.t3_ - T1.t2_.t1_*T2.t1_.t3_ - T1.t2_.t2_*T2.t2_.t3_ - T1.t2_.t3_*T2.t3_.t3_;

		result.t3_.t0_=T1.t3_.t0_*T2.t0_.t0_ - T1.t3_.t1_*T2.t1_.t0_ - T1.t3_.t2_*T2.t2_.t0_ - T1.t3_.t3_*T2.t3_.t0_;
		result.t3_.t1_=T1.t3_.t0_*T2.t0_.t1_ - T1.t3_.t1_*T2.t1_.t1_ - T1.t3_.t2_*T2.t2_.t1_ - T1.t3_.t3_*T2.t3_.t1_;
		result.t3_.t2_=T1.t3_.t0_*T2.t0_.t2_ - T1.t3_.t1_*T2.t1_.t2_ - T1.t3_.t2_*T2.t2_.t2_ - T1.t3_.t3_*T2.t3_.t2_;
		result.t3_.t3_=T1.t3_.t0_*T2.t0_.t3_ - T1.t3_.t1_*T2.t1_.t3_ - T1.t3_.t2_*T2.t2_.t3_ - T1.t3_.t3_*T2.t3_.t3_;*/
		
		return result;
		
	}
	
	inline friend D4T operator* (const D4T<T> &T1,const double &d)
	{
		D4T<T> result=T1;
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)*=d;
			}
		return result;
		
	}
	
	inline friend D4T operator* (const D4T<T> &T1,const DM &d)
	{
		D4T<T> result=T1;
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)*=d;
			}
		return result;
		
	}
	
	//multiplication V_\alpha T^{\alpha \nu}
	
	friend D4V<T> operator* (const D4V<double> &T1,const D4T<T> &T2)
	{
		D4V<T> result(T2(0,0),T2(0,0),T2(0,0),T2(0,0));
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1*tmptransp(a));
		return result;
		
	}
///Warning!!! this does not recognize right types. I do not know, why. Yet
	friend D4V<T> operator* (D4T<T> &T1,const D4V<double> &T2)
	{
		DM zero=0;
		D4V<DM> result(zero,zero,zero,zero);
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1(a)*T2);
		return result;
		
	}

	friend D4V<DM> operator* (const D4V<DM> &T1,const D4T<T> &T2)
	{
		DM zero=0;
		D4V<DM> result(zero,zero,zero,zero);
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1*tmptransp(a));
		return result;
		
	}
///Warning!!!
	friend D4V<DM> operator* (D4T<T> &T1,const D4V<DM> &T2)
	{
		DM zero=0;
		D4V<DM> result(zero,zero,zero,zero);
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1(a)*T2);
		return result;
		
	}

	//multiplication T1^\mu_\alpha T2^{\alpha \nu}
	inline D4T& operator*= (const D4T<T> &T2)
	{
		D4T<T> result(t0_.t0_*T2.t0_.t0_ -  t0_.t1_*T2.t1_.t0_ -  t0_.t2_*T2.t2_.t0_ -  t0_.t3_*T2.t3_.t0_,
		t0_.t0_*T2.t0_.t1_ -  t0_.t1_*T2.t1_.t1_ -  t0_.t2_*T2.t2_.t1_ -  t0_.t3_*T2.t3_.t1_,
		t0_.t0_*T2.t0_.t2_ -  t0_.t1_*T2.t1_.t2_ -  t0_.t2_*T2.t2_.t2_ -  t0_.t3_*T2.t3_.t2_,
		t0_.t0_*T2.t0_.t3_ -  t0_.t1_*T2.t1_.t3_ -  t0_.t2_*T2.t2_.t3_ -  t0_.t3_*T2.t3_.t3_,
		t1_.t0_*T2.t0_.t0_ -  t1_.t1_*T2.t1_.t0_ -  t1_.t2_*T2.t2_.t0_ -  t1_.t3_*T2.t3_.t0_,
		t1_.t0_*T2.t0_.t1_ -  t1_.t1_*T2.t1_.t1_ -  t1_.t2_*T2.t2_.t1_ -  t1_.t3_*T2.t3_.t1_,
		t1_.t0_*T2.t0_.t2_ -  t1_.t1_*T2.t1_.t2_ -  t1_.t2_*T2.t2_.t2_ -  t1_.t3_*T2.t3_.t2_,
		t1_.t0_*T2.t0_.t3_ -  t1_.t1_*T2.t1_.t3_ -  t1_.t2_*T2.t2_.t3_ -  t1_.t3_*T2.t3_.t3_,
		t2_.t0_*T2.t0_.t0_ -  t2_.t1_*T2.t1_.t0_ -  t2_.t2_*T2.t2_.t0_ -  t2_.t3_*T2.t3_.t0_,
		t2_.t0_*T2.t0_.t1_ -  t2_.t1_*T2.t1_.t1_ -  t2_.t2_*T2.t2_.t1_ -  t2_.t3_*T2.t3_.t1_,
		t2_.t0_*T2.t0_.t2_ -  t2_.t1_*T2.t1_.t2_ -  t2_.t2_*T2.t2_.t2_ -  t2_.t3_*T2.t3_.t2_,
		t2_.t0_*T2.t0_.t3_ -  t2_.t1_*T2.t1_.t3_ -  t2_.t2_*T2.t2_.t3_ -  t2_.t3_*T2.t3_.t3_,
		t3_.t0_*T2.t0_.t0_ -  t3_.t1_*T2.t1_.t0_ -  t3_.t2_*T2.t2_.t0_ -  t3_.t3_*T2.t3_.t0_,
		t3_.t0_*T2.t0_.t1_ -  t3_.t1_*T2.t1_.t1_ -  t3_.t2_*T2.t2_.t1_ -  t3_.t3_*T2.t3_.t1_,
		t3_.t0_*T2.t0_.t2_ -  t3_.t1_*T2.t1_.t2_ -  t3_.t2_*T2.t2_.t2_ -  t3_.t3_*T2.t3_.t2_,
		t3_.t0_*T2.t0_.t3_ -  t3_.t1_*T2.t1_.t3_ -  t3_.t2_*T2.t2_.t3_ -  t3_.t3_*T2.t3_.t3_);
		*this=result;
		return *this;
	}
	
	inline D4T& operator&= (const D4T<T> &T1)
	{
		D4T<T> resu(T1.t0_.t0_*t0_.t0_ - T1.t0_.t1_*t1_.t0_ - T1.t0_.t2_*t2_.t0_ - T1.t0_.t3_*t3_.t0_,
					T1.t0_.t0_*t0_.t1_ - T1.t0_.t1_*t1_.t1_ - T1.t0_.t2_*t2_.t1_ - T1.t0_.t3_*t3_.t1_,
					T1.t0_.t0_*t0_.t2_ - T1.t0_.t1_*t1_.t2_ - T1.t0_.t2_*t2_.t2_ - T1.t0_.t3_*t3_.t2_,
					T1.t0_.t0_*t0_.t3_ - T1.t0_.t1_*t1_.t3_ - T1.t0_.t2_*t2_.t3_ - T1.t0_.t3_*t3_.t3_,
					T1.t1_.t0_*t0_.t0_ - T1.t1_.t1_*t1_.t0_ - T1.t1_.t2_*t2_.t0_ - T1.t1_.t3_*t3_.t0_,
					T1.t1_.t0_*t0_.t1_ - T1.t1_.t1_*t1_.t1_ - T1.t1_.t2_*t2_.t1_ - T1.t1_.t3_*t3_.t1_,
					T1.t1_.t0_*t0_.t2_ - T1.t1_.t1_*t1_.t2_ - T1.t1_.t2_*t2_.t2_ - T1.t1_.t3_*t3_.t2_,
					T1.t1_.t0_*t0_.t3_ - T1.t1_.t1_*t1_.t3_ - T1.t1_.t2_*t2_.t3_ - T1.t1_.t3_*t3_.t3_,
					T1.t2_.t0_*t0_.t0_ - T1.t2_.t1_*t1_.t0_ - T1.t2_.t2_*t2_.t0_ - T1.t2_.t3_*t3_.t0_,
					T1.t2_.t0_*t0_.t1_ - T1.t2_.t1_*t1_.t1_ - T1.t2_.t2_*t2_.t1_ - T1.t2_.t3_*t3_.t1_,
					T1.t2_.t0_*t0_.t2_ - T1.t2_.t1_*t1_.t2_ - T1.t2_.t2_*t2_.t2_ - T1.t2_.t3_*t3_.t2_,
					T1.t2_.t0_*t0_.t3_ - T1.t2_.t1_*t1_.t3_ - T1.t2_.t2_*t2_.t3_ - T1.t2_.t3_*t3_.t3_,
					T1.t3_.t0_*t0_.t0_ - T1.t3_.t1_*t1_.t0_ - T1.t3_.t2_*t2_.t0_ - T1.t3_.t3_*t3_.t0_,
					T1.t3_.t0_*t0_.t1_ - T1.t3_.t1_*t1_.t1_ - T1.t3_.t2_*t2_.t1_ - T1.t3_.t3_*t3_.t1_,
					T1.t3_.t0_*t0_.t2_ - T1.t3_.t1_*t1_.t2_ - T1.t3_.t2_*t2_.t2_ - T1.t3_.t3_*t3_.t2_,
					T1.t3_.t0_*t0_.t3_ - T1.t3_.t1_*t1_.t3_ - T1.t3_.t2_*t2_.t3_ - T1.t3_.t3_*t3_.t3_);
		*this=resu;
		return *this;
		
	}

	inline D4T& operator*= (const std::complex<double> &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}

	inline D4T& operator*= (const double &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}


	inline D4T& operator*= (const DM &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}

	inline D4T& operator&= (const DM &d)
	{
		for (unsigned short a=0; a<4; a++)
			for (unsigned short b=0; b<4; b++)
			{
				(*this)(a,b)=d*(*this)(a,b);
			}
		return *this;
	}

	//returns the trace as a matrix
	inline D4T< std::complex<double> > CT(const D4T<DM> &D )
	{
		D4V< std::complex<double> > dummy(std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0));
		D4T< std::complex<double> > result(dummy,dummy);

		for (unsigned short a=0; a<dim; a++)
			for(unsigned short b=0; b<dim;b++)
			result(a,b)=D(a,b).CTrace();
		return result;
	}

	/*inline D4T<std::complex<double> > CT1(const D4T<DM> &D )
	{
		D4T< std::complex<double> > result((t(0)*u(0)).CTrace(), (t(0)*u(1)).CTrace(), (t(0)*u(2)).CTrace(), t(0)*u(3),
											t(1)*u(0), t(1)*u(1), t(1)*u(2), t(1)*u(3),
											t(2)*u(0), t(2)*u(1), t(2)*u(2), t(2)*u(3),
											t(3)*u(0), t(3)*u(1), t(3)*u(2), t(3)*u(3))
		return result;
	}*/


	//antisymmetrization by -\epsilon^{\mu\nu\alpha\beta}T_{\alpha\beta}
	inline D4T Ants(const D4T<T> &t)
	{
		D4T<T> result=t;
		double zero=0;
		double minus=-1.0;
		result.t0_.t0_*=zero;
		result.t1_.t1_*=zero;
		result.t2_.t2_*=zero;
		result.t3_.t3_*=zero;

		result.t0_.t1_=t.t2_.t3_-t.t3_.t2_;
		result.t1_.t0_=(result.t0_.t1_*minus);

		result.t0_.t2_=t.t3_.t1_-t.t1_.t3_;
		result.t2_.t0_=(result.t0_.t2_*minus);

		result.t0_.t3_=t.t1_.t2_-t.t2_.t1_;
		result.t3_.t0_=(result.t0_.t3_*minus);

		result.t1_.t2_=t.t3_.t0_-t.t0_.t3_;
		result.t2_.t1_=(result.t1_.t2_*minus);

		result.t2_.t3_=t.t1_.t0_-t.t0_.t1_;
		result.t3_.t2_=(result.t2_.t3_*minus);

		result.t3_.t1_=t.t2_.t0_-t.t0_.t2_;
		result.t1_.t3_=(result.t3_.t1_*minus);
		//result*=minus;
		return result;
	}
	//adds antisymmetric part -\epsilon^{\mu\nu\alpha\beta}k_\alpha k'_\beta
	inline void AddAnts(const D4V<double> &u,const D4V<double> &v)
	{
		double aaa=(u(2)*v(3)-u(3)*v(2));
		
		(*this)(0,1)+=aaa;
		(*this)(1,0)-=aaa;
		
		aaa=(u(3)*v(1)-u(1)*v(3));
		(*this)(0,2)+=aaa;
		(*this)(2,0)-=aaa;

		aaa=(u(1)*v(2)-u(2)*v(1));
		(*this)(0,3)+=aaa;
		(*this)(3,0)-=aaa;

		aaa=(u(3)*v(0)-u(0)*v(3));
		(*this)(1,2)+=aaa;
		(*this)(2,1)-=aaa;

		aaa=(u(1)*v(0)-u(0)*v(1));
		(*this)(2,3)+=aaa;
		(*this)(3,2)-=aaa;

		aaa=(u(2)*v(0)-u(0)*v(2));
		(*this)(3,1)+=aaa;
		(*this)(1,3)-=aaa;
	}
	
	//T1^{\mu\nu}T2_{\mu\nu} for the same type
	inline T contraction(const D4T<T> &t1,const D4T<T> &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}	
	
	//type of the second tensor is <double>
	inline T contrd (const D4T<T> &t1,const D4T<double> &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}
	//type of the second tensor is complex <double>
	
	inline T contrc (const D4T<T> &t1,const D4T<std::complex<double> > &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}
};


#endif
