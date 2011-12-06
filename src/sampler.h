#ifndef _sampler_h_
#define _sampler_h_

#include <iostream>
/*
 * Adaptive random number sampler
 * Objective: reduce Variance (weight) while preserving Expected value(weight)
*/
template <int N> 
class sampler
{ 
 public:
    string name;
 private:
    bool active;
    double _cnt[N];
    double _sum[N];
    double _max[N];
    double _val[N];
    double _acc[N];           //
    double _Sum;
    double _Cnt;
    int curbin;
    double X;
 public:	
    sampler(string name0,bool a0=1):name(name0),active(a0)
    {
		for(int i=0;i<N;i++)
		{
		   _val[i]=0;
		   _cnt[i]=0;
		   _sum[i]=0;
		   _max[i]=0;
		  }
		 _Sum=N;
		 _Cnt=0;
		 X=1;
		 curbin=-1;
	}
    
    double random()
    {   if(active)
        {    
			double sum=0;
			int i=0;
			while(i<N)
			   sum+=_val[i++];
			_Sum=sum;
			i=0;
			double x=10000.0/(_Cnt+1);
			sum*=frandom()*(1+x);
			double s=0;
			while(sum>(s+=_val[i]+_Sum*x/N))
			  i++;
			curbin=i;
			X=_Sum*(1+x)/(_val[i]+_Sum*x/N)/N;
			return (curbin+frandom())/N;
	    }
	    else
	      { double res=frandom();
			curbin=res*N;
			return res;
		   }	
	    return frandom();
	}
 
    void report(double & weight,double w=-1)  /// correct the final weight for the 
                                  /// nonuniform sampling
    {   _Cnt++;
  	    _cnt[curbin]++;
  	    if(w==-1) w=weight;
		_sum[curbin]+=w;
		if(_max[curbin]<w)
		   {_max[curbin]=w;
	       }
	       
		if(active && curbin>=0)
		{ //  weight*=(_Sum/_val[curbin])/N; 
		     weight*=X; 
		    _val[curbin]=_sum[curbin]/_cnt[curbin];
//		    _val[curbin]=_max[curbin];
			curbin=-1;
		}
			
	}
    
	void dump()
	{  
		cout<<" Sampler "<<name<<": ";
		for(int i=0;i<N;i++)
		{
			cout<<" ";
			cout<<_val[i]<<" ";
		}
		for(int i=0;i<N;i++)
		{
			cout<<" ";
			cout<<_max[i]<<" ";
		}
		cout<<endl;
	} 
	~sampler()
	{
		dump();
	}
};

////////////////////////////////////////////////////////////////////////
//                     Implementation     
////////////////////////////////////////////////////////////////////////

#endif
