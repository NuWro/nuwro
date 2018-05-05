#ifndef _chooser_h_
#define _chooser_h_

#include <iostream>
#include <vector>
#include "params.h"

using namespace std;

class Dyn
{
public: 
    double Mxw;     ///< max x*w   
    int n;          ///< count    
    double Sw;      ///< sum w    
    double Sxw;     ///< sum x*w   
    double Sxxw;    ///< sum x*x*w 

    bool active;    ///< is the dynamics enabled?
    double W;       ///< weight ~ probabability of choosing this channel
    double Wacc;    ///< accumalated weights 
    double Desired; ///< number of events to be generated
    double Ready;   ///< number of accepted events
    int dyn;        ///< code for the dynamics
    string label;   ///< label for the dynamics

    Dyn(int dyn0, string label0="",bool active0=true)
    {
        dyn = dyn0;
        label=label0;
        active = active0;
        W = active;
        Sw = Sxw = Sxxw = Mxw = Desired = Ready = n = 0;
    }

    void add(double x,double bias=1)
    {
        if (!(x == x)) // refuse to add 'nan'
            return;
        double w = 1 / bias;
        double xw = x * w;
        if (xw > Mxw)
            Mxw = xw;
        n ++;
        Sw += w;
        Sxw += xw;
        Sxxw += x * xw;
    }
    
	double avg (){ return Sw==0 ? 0 : Sxw/Sw;}

    double var (){ return Sw==0 ? 0 : Sxxw/Sw - Sxw*Sxw/Sw/Sw;}
    
    double sigma (){ return n<=1 ? 0 : sqrt(var()/(n-1));}
	
    double efficiency()	{ return n==0 || Mxw==0 ? 0 : Sxw/n/Mxw;}

};


class chooser
{ 
 private:
    int N;             // number of channels
    vector<Dyn> proc;  // dynamical channels
    void do_distrib(); // initialize    
 public:
    chooser():N(0){}
    void reset(params &p);                     ///< initialize active dynamics
    int size(){return N;}                      ///< number of chanels
    void add(int i,double x, double bias=1);   ///< add x to i-th channel
    int choose();			                   ///< choose bin number (according to bin weight)
    Dyn& chooseDyn(){return proc[choose()];}   ///< choose bin (according to bin weight)
    void set_weights_to_avg();                 ///< make weights proportional to cross sections
    bool accept(int i,double x,double bias=1); ///< add x and return true with probability xw/Mxw 
    bool active(int i){return proc[i].active;} ///< is the channel active?
    int dyn(int i) {return proc[i].dyn;}       ///< code for the dynamics
    double weight(int i) {return proc[i].W;}   ///< weight of i-th channel
    double total() {return proc[N-1].Wacc;}    ///< Sum of channel weigths 
    double count(int i) {return proc[i].n;}    ///< Number af additions to the bin
    double ratio(int i) {return proc[i].W/proc[N-1].Wacc;}///< bin probability 
    double avg (int i) {return proc[i].avg();} ///< weighted average
    double var (int i) {return proc[i].var();} ///< variance 
    double sigma (int i){return proc[i].sigma();} ///< sigma
    double efficiency(int i) {return proc[i].efficiency();} /// efficincy of i-th channel
    void report();                  ///< report active channels characteristics
    void short_report(ostream &f);	///< write calculated total cross sections for each channel to file 
    void calculate_counts(int ilosc); ///< calculate how many events to generate for each channel
    int desired(int i) {return proc[i].Desired;} ///< how many events to generate for each channel
    int ready(int i) {return proc[i].Ready;}; ///< how many events got accepted (and writen to file)
    
};

////////////////////////////////////////////////////////////////////////
//                     Implementation     
////////////////////////////////////////////////////////////////////////

inline void chooser::do_distrib()
{
    double prev=0;
    for(int i=0;i<N;i++)
        proc[i].Wacc=prev=prev+proc[i].W;
}

////////////////////////////////////////////////////////////////////////
inline void chooser::reset(params &p)
{
    proc.clear();
    if(p.dyn_qel_cc) proc.push_back(Dyn(0,"QELcc"));
    if(p.dyn_qel_nc) proc.push_back(Dyn(1,"QELnc"));
    if(p.dyn_res_cc) proc.push_back(Dyn(2,"REScc"));
    if(p.dyn_res_nc) proc.push_back(Dyn(3,"RESnc"));
    if(p.dyn_dis_cc) proc.push_back(Dyn(4,"DIScc"));
    if(p.dyn_dis_nc) proc.push_back(Dyn(5,"DISnc"));
    if(p.dyn_coh_cc) proc.push_back(Dyn(6,"COHcc"));
    if(p.dyn_coh_nc) proc.push_back(Dyn(7,"COHnc"));
    if(p.dyn_mec_cc) proc.push_back(Dyn(8,"MECcc"));
    if(p.dyn_mec_nc) proc.push_back(Dyn(9,"MECnc"));
    if(p.dyn_e_el  ) proc.push_back(Dyn(20,"EEL  "));
    N=proc.size();
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
inline void chooser::add(int i, double x, double bias)
{
    proc[i].add(x,bias);
}

////////////////////////////////////////////////////////////////////////
inline int chooser::choose() 
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
inline void chooser::set_weights_to_avg()
{
	cout<<"Updating"<<endl;
	for(int i=0;i<N;i++)
		proc[i].W=avg(i);
	do_distrib(); 
}

////////////////////////////////////////////////////////////////////////
inline bool chooser::accept(int i, double x, double bias)
{
    /// return 1 with probability x/max(i)
    /// increase Mwx if needed
    double prevmax = proc[i].Mxw;
    double prevready = proc[i].Ready;
    add(i, x, bias);
    if (prevmax == proc[i].Mxw)
    {
        if (x / bias > frandom() * proc[i].Mxw)
        {
            proc[i].Ready++;
            return true;
        }
        else 
            return false;
    }
    else
    {

/*	       
        cout.precision(6);
        cout << "Dyn[" << dyn(i) << "]"
	    <<"   New max =" << setw(11)<< proc[i].maxxw
        <<"   Prev max=" << setw(11)<< prevmax 
        <<"   " <<setw(11)<<prevready-ready(i) << " events deleted." 
        <<" Efficiency "<<setw(11)<<efficiency(i)*100 <<" %."<< endl;
*/

        proc[i].Ready*=prevmax/proc[i].Mxw;// "remove" some events
        proc[i].Ready++;
        return true;
    }
}

////////////////////////////////////////////////////////////////////////
inline void chooser::short_report(ostream &f)
{
/// write calculated total cross sections for each channel to file
    f << "dyn  n        ratio          sigma[cm2] " << endl;
    for (int k = 0; k < N; k++)
    {
        f << dyn(k) << " " << setw(5) << desired(k)
          << " " << setw(15) << ratio(k)
          << " " << setw(15) << avg(k) << endl;
    }
}

////////////////////////////////////////////////////////////////////////
inline void chooser::report()
{
 cout<<"--------------------------------------------------------------";
 cout<<"---------------------------------------------------"<<endl;
	  cout <<"dyn| label |"
	       <<"  weight      |"
	       <<"  ratio       |"
	       <<"  efficiency    |"
	       <<"  mean_value     |" 
	       <<"  deviation      |" 
	       <<"  sigma          |" 
	       <<  endl;
 cout<<"--------------------------------------------------------------";
 cout<<"---------------------------------------------------"<<endl;
      for (int j = 0; j < N; j++)
        cout <<setw(2)<<  proc[j].dyn<<setprecision(6) << " | "
             <<proc[j].label  << " |"
             <<setw(12)<< weight (j) << "  |"
             <<setw(12)<< ratio (j) << "  |"
             <<setw(12)<< efficiency(j)*100 << " %  | "
             <<setw(12)<< avg (j) << " cm2| "
             <<setw(12)<< sqrt(var(j))<< " cm2| "
             <<setw(12)<< sigma(j)<< " cm2| "
             << endl;
 cout<<"--------------------------------------------------------------";
 cout<<"---------------------------------------------------"<<endl;
}

////////////////////////////////////////////////////////////////////////
inline void chooser::calculate_counts(int Total)
{
    double total=1.0;
    for (int k = 0; k < N-1; k++)
    {  
        proc[k].Ready=0;
        double frac=ratio(k);
        if(frac==0)
            proc[k].Desired=0;
        else
            proc[k].Desired = int(frac/total*Total+0.5);
        Total-=proc[k].Desired;
        total-=frac;
    }
    proc[N-1].Desired=Total;
}

#endif
