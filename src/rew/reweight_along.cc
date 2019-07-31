#include "../event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Reweighters.h"
#include "../nucleus.h"
#include <iostream>

void SetupSPP(params&); // defined in rewRES.cc

static string name(RewParam &a, double k1, RewParam&b, double k2)
{
	if(k2==0)
		return a.name+(k1>0?"Up":"Down");
	else
		return a.name+(k1>0?"Up_":"Down_")+b.name+(k2>0?"Up":"Down");
}

int main(int argc, char *argv[])
{
	if(argc<4)
	{
		cerr<<"[INFO] Usage: " 
		"reweight_along <nuwro_output.root> <results_filename.root> rewpar1 rewpar2 ...\n";
		exit(-1);
	}

    /// setup input tree 
	TFile * f1 = new TFile(argv[1]);
	if (f1->IsZombie()) 
	{
 	    std::cout << "Error opening file" << std::endl;
   		exit(-1);
	}

	TTree * t1 = dynamic_cast<TTree *>(f1->Get("treeout"));
	if (!t1) 
	{
	    std::cerr << "[ERROR]: Couldn't find TTree (\"treeout\") in file " << argv[1] << "." << std::endl;
   		exit(-1);
  	}	

	event *e=new event;
	t1->SetBranchAddress("e", &e);

	vector<RewParam*> args;
	int nargs=0;
	for(int i=3;i<argc;i++)
	{
		RewParam& p=rew(argv[i]);
		if(p.name=="")
		{	

			cerr <<"[Error] parameter \""<<argv[i]<<"\" can not be reweighted. Try one of:\n";
			rew.list(cerr);
			exit(1);
		}
		else 
		{  
           args.push_back(&p);
           REW(p.engine).active=true;
           nargs++;
		}
	}
    /// create output file
	TFile * f2 = new TFile(argv[2],"recreate");
	TTree * t2 = new TTree("weights","Tree of weights");

	vector<string> names;
	vector<double> weights(500);
	int dir=0;
    // Setup branches 
	for(int i=0;i<nargs;i++)
	{
		for(int k1=-1;k1<=1;k1+=2)
		{
			///setup branch for args[i].twk=k1*2;
			names.push_back(name(*args[i],k1,rew.End,0));
			t2->Branch(names[dir].c_str(),&(weights[dir]),(names[dir]+"/D").c_str());
			dir++;

			for(int j=i+1;j<nargs;j++)
			{
				for(int k2=-1;k2<=1;k2+=2)
				{
					///setup branch for args[i].twk=k1*2 and args[j].twk=k2*2;
					names.push_back(name(*args[i],k1,*args[j],k2));
					t2->Branch(names[dir].c_str(),&(weights[dir]),string(names[dir]+"/D").c_str());
		  		   	dir++;
				}

			}
		}
	}

	dir=0;
    // Calculate and save weights for each event 
	int n = t1->GetEntries();
	//if(false)
	for (int ie=0;ie<n;ie++)
	{
		t1->GetEntry(ie);
		REW.init(e->par);
		nucleus t(e->par);
		if(ie==0)
			SetupSPP(e->par);
		double nominal=REW.weight(*e,e->par,t);
		dir=0;
		for(int i=0;i<nargs;i++)
		{
			for(int k1=-1;k1<=1;k1+=2)
			{
				//calc weight for args[i].twk=k1*2;
				args[i]->setTwk(k1*2);
				//cout<<args[i]->name<<"="<<args[i]->val<<endl;
				/// change params if needed
				weights[dir]=REW.weight(*e,e->par,t)/nominal;
				dir++;
				for(int j=i+1;j<nargs;j++)
				{
					for(int k2=-1;k2<=1;k2+=2)
					{
			  		   //calc weight for args[i].twk=k1*2 and args[j].twk=k2*2;
						args[j]->setTwk(k2*2);
						/// change params if needed
						weights[dir]=REW.weight(*e,e->par,t)/nominal;
						dir++;
					}
					args[j]->setTwk(0);
				}
				args[i]->setTwk(0);
			}
		}

		t2->Fill();
	}
	f1->Close();
	delete f1;
	f2->Write();
	f2->Close();
	delete e;
	delete f2;
	cout<<"Output file: \""<<argv[2]<<"\""<<endl;	
	return 0;
}