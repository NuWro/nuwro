#include <iostream>

#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "hist.h"

using namespace std;

struct ev1
{
	double weight,q2,costheta;
};


vector<ev1> v;

void makehist(string parname,string filename,double ev1::*val,double xunit,double yunit,int nbins=100)
{
	double a=0,b=0,x=0;
	int n=0;
	for(ev1 e:v)
	{
		if(n++==0)
		    a=b=e.*val;
		else
		{
			x=e.*val;
			if(x<a) a=x;
			if(x>b) b=x;
		}
	}		
	hist h(parname,a,b,nbins);
	
	for(auto e:v)
	{
		
		if(e.costheta<0.999 and e.costheta>-2)
			h.insert_value(e.*val,e.weight);
		else
			h.insert_value(e.*val,0);
	}

	h.plot(filename,xunit,yunit);
	cout<<"file \""<< filename<< "\" created"<<endl;

}

int main (int argc, char* argv[])
{

	event *e=new event;

	for(int j=1;j<argc;j++)
	{
		TFile * f1 = new TFile(argv[j]);
		TTree * t1 = (TTree*)f1->Get("treeout");
		t1->SetBranchAddress("e", &e);
		int n = t1->GetEntries();
		for (int i=0;i<n;i++)
		{
			if((i+1)%10000 == 0) 
				cout<< argv[j]<<": "<<(i+1)*100/n<<"%        \r"<<flush;
			t1->GetEntry(i);

			ev1 e1;
			e1.weight=e->weight*cm2;
			e1.q2=e->q2();
			e1.costheta=e->costheta();

			v.push_back(e1);			
		}
		f1->Close();
		delete f1;
		cout<<endl;
	}
	
	makehist("q2","whist.q2.txt",&ev1::q2,GeV2,cm2/GeV2);
	makehist("cos","whist.cos1.txt",&ev1::costheta,1,cm2/4/M_PI,100);
	makehist("cos","whist.cos2.txt",&ev1::costheta,1,cm2/4/M_PI,500);
    return 0;
}
