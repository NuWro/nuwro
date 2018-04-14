#ifndef _chooser_h_
#define _chooser_h_

#include <iostream>
#include "params.h"

using namespace std;

template <int N> 
class chooser
{ 
    class bin
    {
	public: 
        double maxxw; ///< max xw   in i-th bin
        int n;        ///< count    in i-th bin
        double sw;    ///< sum w    in i-th bin
        double sxw;   ///< sum xw   in i-th bin
        double sx2w;  ///< sum x^2w in i-th bin

        bool active;    ///< the bin is active
        double W;       ///< weight of  i-th bin
        double Wacc;    ///< accumalated weights
        double Desired; ///< #events for i-th bin
        double Ready;   ///< #events for i-th bin
        int dyn;        ///< code for the dynamics
        double *Maxxw;  /// Total Max

        bin(double *Maxxwp=NULL,int code=0,bool active0=true):Maxxw(Maxxwp)
        {
            dyn = code;
            active = active0;
            W = active;
            sw = 0;
            sxw = 0;
            sx2w = 0;
            n = 0;
            maxxw = 0;
            Desired = 0;
            Ready = 0;
        }

        void add(double x,double bias=1)
        {
            if (!(x == x)) // refuse to add 'nan'
                return;
            double w = 1 / bias;
            double xw = x * w;
            n += 1;
            sw += w;
            sxw += xw;
            sx2w += x * xw;
            if (xw > maxxw)
            {
                maxxw = xw;
                if (xw > *Maxxw)
                    *Maxxw = xw;
            }
        }
	double avg (){ return n ? sxw/sw : 0;}

        double var ()  
	{
	    if(n==0)
		return 0;
	    return sx2w/sw - sxw*sxw/sw/sw;
	}
        double sigma ()
	{
	    if(n<=1)
		return 0;
	    return sqrt(var()/(n-1));
	}
	double efficiency()
	{ if(n>0 && maxxw>0)
		return sxw/n/maxxw;
	  else
		return 0;
	}
        
    };
 private:
    double maxxw;
    bin proc[N];
    void do_distrib();    
 public:
    int size(){return N;}       ///<number of chanels
    void reset(params &p);
    void add(int i,double x, double bias=1);   ///< add x to i-th bin
    int choose();			       ///< choose bin
    void set_weights_to_avg();        
    bool accept(int i,double x,double bias=1); /// return 1 with probability x/max(i)
                                            /// increase maxwx and Maxwx if needed                                  
    bool active(int i){return proc[i].active;}
    int dyn(int i);                         ///      
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
    double total(){return proc[N-1].Wacc;}

};

////////////////////////////////////////////////////////////////////////
//                     Implementation     
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::do_distrib()
	{
	  double prev=0;
	  for(int i=0;i<N;i++)
		 proc[i].Wacc=prev=prev+proc[i].W;
	}
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::reset(params &p)
    {
        double prev=0;
        proc[0]=bin(&maxxw,0,p.dyn_qel_cc);
        proc[1]=bin(&maxxw,1,p.dyn_qel_nc);
        proc[2]=bin(&maxxw,2,p.dyn_res_cc);
        proc[3]=bin(&maxxw,3,p.dyn_res_nc);
        proc[4]=bin(&maxxw,4,p.dyn_dis_cc);
        proc[5]=bin(&maxxw,5,p.dyn_dis_nc);
        proc[6]=bin(&maxxw,6,p.dyn_coh_cc);
        proc[7]=bin(&maxxw,7,p.dyn_coh_nc);
        proc[8]=bin(&maxxw,8,p.dyn_mec_cc);
        proc[9]=bin(&maxxw,9,p.dyn_mec_nc);
        do_distrib();
	    if(proc[N-1].Wacc==0)
	    {
            cerr<<"No active dynamics - chooser invalid"<<endl;
	        exit(19);
	    }
	   else
	    std::cout<<"chooser created"<<std::endl;
    }
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline void chooser<N>::add(int i,double x,double bias)
    {
        proc[i].add(x,bias);
    }
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::choose()//wybiera dynamike
    {
        double x = frandom() * proc[N-1].Wacc;
        int i = 0;
        while (x >= proc[i].Wacc)
        {
            ++i;
            if (i == N)
                return N - 1;
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
		proc[i].W=(proc[i].n ? avg(i) : 0);
	}
	do_distrib(); 
    }


////////////////////////////////////////////////////////////////////////
/// return 1 with probability x/max(i)
/// increase maxwx and Maxwx if needed
template <int N> 
inline bool chooser<N>::accept(int i,double x,double bias)
	{ // double w=1;
		double prevmax=proc[i].maxxw;
		double prevready=proc[i].Ready;
	   	add(i,x,bias);
	   	if(prevmax==proc[i].maxxw)
		   if(x/bias>frandom()*proc[i].maxxw)
		     {proc[i].Ready+=1;
		       return 1;
		     }
		    else 
		       return 0;  
		else
		  {
           proc[i].Ready*=prevmax/proc[i].maxxw;
  	       cout.precision(6);
/*	       cout << "Dyn[" << i << "]"
//	       <<"   New max =" << setw(11)<< max(i) 
//		   <<"   Prev max=" << setw(11)<< prevmax 
		   <<"   " <<setw(11)<<prevready-Ready[i] << " events deleted." 
		   <<" Efficiency "<<setw(11)<<sxw[i]/n[i]/maxxw[i]*100 <<" %."<< endl;
		   */
           proc[i].Ready*=prevmax/proc[i].maxxw;
           proc[i].Ready += 1;
           return 1;
	      }
	}                                      
    
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::dyn(int i)
    {
       return proc[i].dyn;
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::weight(int i)
    {
       return proc[i].W;
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::count(int i)
    {
       return proc[i].n;
    }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::ratio(int i)/// probability of choosing i-th bin
    {
       return proc[i].W/proc[N-1].Wacc;
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::avg (int i)
    {
        return proc[i].avg();
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::var (int i)
    {
      
      return proc[i].var();
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::sigma (int i)
    {
      return proc[i].sigma();
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
        cout <<setw(2)<<  proc[j].dyn<<setprecision(6) <<" |"
             <<setw(12)<< weight (j)  <<            "  |"
             <<setw(12)<< ratio (j)  <<             "  |"
             <<setw(12)<< efficiency(j)*100 << " %  | "
             <<setw(12)<< avg (j)    <<             " cm2| "
             <<setw(12)<< sqrt(var(j))<<            " cm2| "
             <<setw(12)<< sigma(j)<<       " cm2| "
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
      {  
          proc[k].Ready=0;
          double frac=ratio(k);
          if(frac==0)
               proc[k].Desired=0;
          else
               proc[k].Desired = int(frac/total*ilosc+0.5);
         ilosc-=proc[k].Desired;
         total-=frac;
      }
      proc[N-1].Desired=ilosc;
      }

////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::desired(int i)
    {
      return proc[i].Desired;
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline int chooser<N>::ready(int i)
    {
      return proc[i].Ready;
    }
////////////////////////////////////////////////////////////////////////
template <int N> 
inline double chooser<N>::efficiency(int i)
    { 
        return proc[i].efficiency();
    }
 
////////////////////////////////////////////////////////////////////////

#endif
