#ifndef _chooser_h_
#define _chooser_h_

#include <iostream>


template <int N> 
class chooser
{ 
 private:
 
    double Maxx;                ///< global max 
    double maxx[N];             ///< max wx   in i-th bin
    int    n[N];                 ///< count    in i-th bin
    double sb[N];                ///< sum b    in i-th bin
    double sw[N];                ///< sum w    in i-th bin
    double sx[N];                ///< sum x    in i-th bin
    double sxw[N];               ///< sum xw   in i-th bin
    double sx2w[N];              ///< sum x^2w in i-th bin

    double W[N];                ///< weight of  i-th bin
    double Wacc[N];             ///< accumalated weights
    double Desired[N];          ///< #events for i-th bin
    double Ready[N];            ///< #events for i-th bin
    void do_distrib();    
 public:
	int size(){return N;}       ///<number of chanels
    void reset(bool active[N]);
    void add(int i,double x, double bias=1);   ///< add x to i-th bin
    int choose();			                ///< choose bin
    void set_weights_to_avg();        
    bool accept(int i,double x,double bias=1); /// return 1 with probability x/max(i)
                                            /// increase maxwx and Maxwx if needed                                  
    double weight(int i);                   /// 
    double count(int i);
    double ratio(int i);                    /// probability of choosing i-th bin
    double avg (int i);
    double var (int i);
    double sigma (int i);
    double efficiency(int i);
    void report();
    void short_report(ostream &f);	///< write calculated total cross sections for each channel to file 
	void calculate_counts(int ilosc);
	int desired(int i);
	int ready(int i);
};

////////////////////////////////////////////////////////////////////////
//                     Implementation     
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::do_distrib()
	{
	  double prev=0;
	  for(int i=0;i<N;i++)
		 Wacc[i]=prev=prev+W[i];
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::reset(bool active[N])
    {double prev=0;
	 for(int i=0;i<N;i++)
        {
         W[i]=active[i];
         Wacc[i]=prev=prev+W[i];
		 sx[i]=0;
		 sb[i]=0;
		 sw[i]=0;
		 sxw[i]=0;
		 sx2w[i]=0;
		 n[i]=0;
		 maxx[i]=0;
		 Desired[i]=0;
		 Ready[i]=0;
	     } 
	  if(Wacc[N-1]==0)
	    {cerr<<"No active dynamics - chooser invalid"<<endl;
	     exit(17);
	    }
	   else
	    std::cout<<"chooser created"<<std::endl;
    }
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::add(int i,double x,double bias)
    {double w=1;
     x/=bias;
	 if(!(x==x)) // refuse to add 'nan'
	    return;
	 double xw=x*w;
     n[i]+=1;
     sb[i]+=1/bias;
     sw[i]+=w;
     sx[i]+=x;
     sxw[i]+=xw;
     sx2w[i]+=x*xw;
     if(x>maxx[i])
        {maxx[i]=x;
		 if(x>Maxx)
		   Maxx=x;
		}   
    }
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::choose()//wybiera dynamike
    { 
      double x=frandom()*Wacc[N-1];
      int i=0;
	  while(x>=Wacc[i]) 
	      {++i;
	       if(i==N) return N-1;
	      }
      return i;      
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::set_weights_to_avg()
    {
    cout<<"Updating"<<endl;
     for(int i=0;i<N;i++)
       {
		W[i]=(n[i] ? avg(i) : 0);
		//maxx[i]=0;
	   }
	  do_distrib(); 
    }

////////////////////////////////////////////////////////////////////////
/*
template <int N> 
inline void chooser<N>::set_weights_to_max()
    {
     cout<<"Updating"<<endl;
     for(int i=0;i<N;i++)
			W[i]=max(i);
	 do_distrib();
    }
*/
////////////////////////////////////////////////////////////////////////
/// return 1 with probability x/max(i)
/// increase maxwx and Maxwx if needed
template <int N> 
inline bool chooser<N>::accept(int i,double x,double bias)
	{  double w=1;
		double prevmax=maxx[i];
		double prevready=Ready[i];
	   	add(i,x,bias);
	   	if(prevmax==maxx[i])
		   if(x/bias>frandom()*maxx[i])
		     {Ready[i]+=1;
		       return 1;
		     }
		    else 
		       return 0;  
		else
		  {
           Ready[i]*=prevmax/maxx[i];
  	       cout.precision(6);
/*	       cout << "Dyn[" << i << "]"
//	       <<"   New max =" << setw(11)<< max(i) 
//		   <<"   Prev max=" << setw(11)<< prevmax 
		   <<"   " <<setw(11)<<prevready-Ready[i] << " events deleted." 
		   <<" Efficiency "<<setw(11)<<sxw[i]/n[i]/maxx[i]*100 <<" %."<< endl;
		   */
           Ready[i]*=prevmax/maxx[i];
		   Ready[i]+=1;
		   return 1;
	      }
	}                                      
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::weight(int i)
    {
       return W[i];
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::count(int i)
    {
       return n[i];
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::ratio(int i)/// probability of choosing i-th bin
    {
       return W[i]/Wacc[N-1];
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::avg (int i)
    {
      if(n[i]==0)
	return 0;
      return sxw[i]/sw[i]/sb[i]*n[i];
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::var (int i)
    {
      if(n[i]==0)
	return 0;
      return (sx2w[i]/sw[i]- sxw[i]*sxw[i]/sw[i]/sw[i])/sb[i]*n[i]/sb[i]*n[i];
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::sigma (int i)
    {
      if(n[i]<=1)
	return 0;
      return sqrt(var(i)/(n[i]-1));
    }

////////////////////////////////////////////////////////////////////////
/// write calculated total cross sections for each channel to file 
template <int N> 
inline void chooser<N>::short_report(ostream &f)
    {
	  f <<"dyn  n   ratio     sigma[cm2] " <<  endl;
	  for (int k = 0; k < N; k++)
		{
		  f << k <<" "<< setw(5)  << desired(k) 
				<<" "<<  setw(10) << ratio(k)
				<<" "<<  setw(10) << avg(k) << endl;
		}
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::report()
    {
 cout<<"--------------------------------------------------------------";
 cout<<"-------------------------------------------"<<endl;
	  cout <<"dyn|"
	       <<"  weight      |"
	       <<"  ratio       |"
	       <<"  efficiency    |"
	       <<"  mean_value     |" 
	       <<"  deviation      |" 
	       <<"  sigma          |" 
	       <<  endl;
 cout<<"--------------------------------------------------------------";
 cout<<"-------------------------------------------"<<endl;
      for (int j = 0; j < N; j++)
        cout <<setw(2)<< j<<setprecision(6) <<" |"
             <<setw(12)<< weight (j)  <<            "  |"
             <<setw(12)<< ratio (j)  <<             "  |"
             <<setw(12)<< sx[j]/n[j]/maxx[j]*100 << " %  | "
             <<setw(12)<< avg (j)    <<             " cm2| "
             <<setw(12)<< sqrt(var(j))<<            " cm2| "
             <<setw(12)<< sqrt(var(j)/n[j])<<       " cm2| "
             << endl;
 cout<<"--------------------------------------------------------------";
 cout<<"-------------------------------------------"<<endl;
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::calculate_counts(int ilosc)
	{
      double total=1.0;
      for (int k = 0; k < N-1; k++)
      {  Ready[k]=0;
         double frac=ratio(k);
         if(frac==0)
	        Desired[k]=0;
         else
	        Desired[k] = int(frac/total*ilosc+0.5 );
         ilosc-=Desired[k];
         total-=frac;
      }
      Desired[N-1]=ilosc;
	 }
	 
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::desired(int i)
	{
      return Desired[i];
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::ready(int i)
	{
      return Ready[i];
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::efficiency(int i)
	{
      return sxw[i]/n[i]/maxx[i];
	}
 
////////////////////////////////////////////////////////////////////////

#endif
