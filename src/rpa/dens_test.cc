#include "../jednostki.h"
#include "../nucleus_data.h"
#include "../calg5.h"
#include "../util2.h"
#include <sstream>
#include <iostream>
#include <algorithm>
	
using namespace std;

int main()
{
	cout<<"<table border=1>"<<endl;
	//if(false)
	for(int j=0;dens_data[j].p()!=0;j++)
		for(int i=0;dens_data[i+1].p()!=0;i++)
		{
			nucleus_data &a(dens_data[i]);
			nucleus_data &b(dens_data[i+1]);
			if(a.p()*1000+a.n()	>b.p()*1000+b.n())
			{
				nucleus_data c=a;
				a=b;b=c;
			}
		 
		}
	int lp=0,lup=1;
	for(nucleus_data *i=dens_data;i->p();++i)
	{   
		int p=i->p();
		int n=i->n();
		if(lp)
			if(i[-1].p()!=i->p() ||i[-1].n()!=i->n())
				lup++;
		double r=i->r();
		stringstream s;

		s<<"dens-"<<p<<','<<n<<'-'<<el[p].symbol<<p+n<<'-'<<i->name()<<".dat\0"<<flush;
		const char* name =s.str().c_str();
		auto dens=[=](double x)
		{ 
			return i->dens(x);
		};
		plot(dens,0.,r,name,100,fermi,1/fermi3);

		auto rrdens=[=](double x)
		{ 
			return i->dens(x)*x*x*4*M_PI;
		};
		double masa=calg5a(rrdens,0.,i->r());	

		stringstream s2;
		s2<<"rrdens-"<<p<<','<<n<<'-'<<el[p].symbol<<p+n<<'-'<<i->name()<<".dat\0"<<flush;
		const char* name2 =s2.str().c_str();
		plot(rrdens,0.,r,name2,100,fermi,1/fermi);
		
		
		auto kfrrdens=[=](double x)
		{
			double ro=i->dens(x);
			return FermiMomentum(ro)*ro*x*x*4*M_PI;
		};
		double kf=calg5a(kfrrdens,0.,r);

		auto prrdens=[=](double x)
		{
			double ro=i->dens(x);
			return FermiMomentum(ro*2*p/(p+n))*ro*x*x*4*M_PI;
		};
		double kfp=calg5a(prrdens,0.,r);

		auto nrrdens=[=](double x)
		{
			double ro=i->dens(x);
			return FermiMomentum(ro*2*n/(p+n))*ro*x*x*4*M_PI;
		};			
		double kfn=calg5a(nrrdens,0.,r);

		auto rrrrdens=[=](double x)
		{
			double ro=i->dens(x);
			return x*x*ro*x*x*4*M_PI;
		};
		double rms=calg5a(rrrrdens,0.,r);
		cout<<"<tr><td>"<< ++lp<<"<td>"<< lup
			<<"<td>"<<p<<','<<n
			<<"<td> "<<el[p].symbol<<p+n<<'\t'
			<<"<td>"<<i->name()
			<<"<td> rms="<<sqrt(rms/masa)/fermi	<<'\t'
			<<"<td> kf= "<<i->kF()
//			<<"<td>  "<<kf/masa-i->kF()<<'\t'
//			<<kfp/masa<<'\t'
//			<<kfn/masa<<'\t'
			<<"<td> Mf= "<<i->Mf()
			<<"<td> norm/A="<<masa/(p+n)<<'\t'
			<<"<td> A/norm="<<(p+n)/masa<<'\t'
			<<"<td> r="<<i->r()/fermi<<'\t'
			<<"</tr>"<<endl;
	}
	cout<<"</table"<<endl;
//	for(int i=0;i<1000;i++)
		cout<<best_data(6,6)->kF()<<endl;
		cout<<best_data(6,6)->Mf()<<endl;

		cout<<best_data(18,22)->kF()<<endl;
		cout<<best_data(18,22)->Mf()<<endl;
	
}
