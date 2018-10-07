//generic template of a Lorentz 4-vector

#ifndef _D4V_H_
#define _D4V_H_

template< class T>
class D4V
{
	//private:
		T t0_;
		T t1_;
		T t2_;
		T t3_;

	public:
	D4V(T t0,T t1, T t2, T t3): t0_(t0), t1_(t1), t2_(t2), t3_(t3)
	{}
	
	inline T operator() (unsigned short j) const
	{
		switch(j)
		{
			case 0: 	return t0_;
		    case 1:		return t1_;
		    case 2:		return t2_;
		    case 3:		return t3_;
			default:
			{
				double failsafe=0;
				return t1_*failsafe;
			} 
		}	
	}

	inline T& operator() (unsigned short j)
	{
		switch(j)
		{
			case 0: 	return t0_;
		    case 1:		return t1_;
		    case 2:		return t2_;
		    case 3:		return t3_;
			default:	return t1_;
			 
		} 
	}


	//again, bunch of overloaded operators
	inline friend D4V operator+ (const D4V<T> &v1,const D4V<T> &v2)
	{
		return D4V<T>(v1.t0_ + v2.t0_, v1.t1_ + v2.t1_, v1.t2_ + v2.t2_, v1.t3_ + v2.t3_);
	}

	
	inline D4V& operator+= (const D4V<T> &v1)
	{
		t0_+=v1.t0_;
		t1_+=v1.t1_;
		t2_+=v1.t2_;
		t3_+=v1.t3_;
		return *this;
	}

	inline friend D4V operator- (const D4V<T> &v1,const D4V<T> &v2)
	{
		return D4V<T>(v1.t0_ - v2.t0_, v1.t1_ - v2.t1_, v1.t2_ - v2.t2_, v1.t3_ - v2.t3_);

	}

	inline D4V& operator-= (const D4V<T> &v1)
	{
		t0_-=v1.t0_;
		t1_-=v1.t1_;
		t2_-=v1.t2_;
		t3_-=v1.t3_;
		return *this;
	}

	inline friend D4V<DM> operator* (const D4V<T> &v1, const  DM &d)
	{
		return D4V<DM> (v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
	}
	
	inline friend D4V<DM> operator* (const  DM &d, const D4V<T> &v1)
	{
		return D4V<DM>(d*v1.t0_, d*v1.t1_, d*v1.t2_, d*v1.t3_);
	}

	
	inline friend D4V<T> operator* (const D4V<T> &v1, const std::complex<double> &d)
	{
		return D4V<T> (v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
	}

	inline friend D4V<T> operator/ (const D4V<T> &v1, const std::complex<double> &d)
	{
		return D4V<T> (v1.t0_/d, v1.t1_/d, v1.t2_/d, v1.t3_/d);
	}
		
	inline friend D4V<T> operator* (const D4V<T> &v1, const double &d)
	{
		return D4V<T> (v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
	}

	inline friend D4V<T> operator/ (const D4V<T> &v1, const double &d)
	{
		return D4V<T> (v1.t0_/d, v1.t1_/d, v1.t2_/d, v1.t3_/d);
	}

	//vector-vector covariant multiplication. notice metric= (+,-,-,-)
	inline friend T operator* (const D4V<T> &v1,const D4V<double> &v2)
	{
		return v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3);
	}

	inline friend T operator* (D4V<T> v1,D4V< std::complex<double> > v2)
	{
		return v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3);
	}

	inline friend DM operator* (D4V<T> v1,D4V< DM > v2)
	{
		DM result=v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3);
		return result;

	}

	inline D4V& operator*= (const double &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}
	
	inline D4V& operator*= (const std::complex<double> &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}

	inline D4V& operator*= (const DM &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}
//left-handed multiplication by a matrix
	inline D4V& operator&= (const DM &d)
	{
		t0_=d*t0_;
		t1_=d*t1_;
		t2_=d*t2_;
		t3_=d*t3_;
		
		return *this;

	}


	inline void hermit()
	{
		(*this)(0).hermit();
		(*this)(1).hermit();
		(*this)(2).hermit();
		(*this)(3).hermit();
	}
	
	
};


#endif
