#include "../event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Reweighters.h"
#include "../nucleus.h"
#include <iostream>

void SetupSPP(params&); // defined in rewRES.cc




int main(int argc, char *argv[])
{
	if(argc<4)
	{
		cerr<<"[INFO] Usage: " 
		"reweight_to <nuwro_output.root> -o <results_filename.root> -p par1 val1 -p par2 val2 ...\n";
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
	vector<double> vals;
	int nargs=0;
	char* outname=NULL;
	for(int i=2;i<argc;i++)
	{

		if(string(argv[i])=="-o")
		{
			outname=argv[++i];
			continue;
		}
		else
		if(string(argv[i])=="-p")
		{
			 
			RewParam& p=rew(argv[++i]);
			if(p.name=="")
			{	

				cerr <<"[Error] parameter \""<<argv[i]<<"\" can not be reweighted. Try one of:\n";
				rew.list(cerr);
				exit(1);
			}
			else 
			{  
			   args.push_back(&p);
			   stringstream s(argv[++i]);
			   double x;
			   if(i<argc && (s>>x))
			   {
				   vals.push_back(x);
				   REW(p.engine).active=true;
				   nargs++;
			   }
			   else
			   {
				   cerr <<"[Error] parameter \""<<argv[i-1]<<"\" must have numeric value (not \""<<argv[i]<<"\").\n";
				   exit(1);
			   }
			}
		}
		else
		{
		   cerr <<"[Error] -p or -o expexted instead of\""<<argv[i]<<"\".\n";
			cerr<<"[INFO] Usage: " 
			"reweight_to <nuwro_output.root> -o <results_filename.root> -p par1 val1 -p par2 val2 ...\n";
			exit(1);
			
			
		}
			
	}
	
	if(nargs==0 or outname==NULL)
	{
		   cerr <<"[Error] -p or -o .\n";
			cerr<<"[INFO] Usage: " 
			"reweight_to <nuwro_output.root> -o <results_filename.root> -p par1 val1 -p par2 val2 ...\n";
			exit(1);
		
	}	
	
    /// create output file
	TFile * f2 = new TFile(outname,"recreate");
	TTree * t2 = new TTree("weights","Tree of weights");

	//event *e=new event;
	TFile * f3 = new TFile("out.root","recreate");
	TTree * t3 = new TTree("treeout","Tree of events");
	t3->Branch("e","event",&e);


	double weight;

    // Setup branches 
	
	t2->Branch("weight",&weight,"weight/D");

    // Calculate and save weights for each event 
	int n = t1->GetEntries();
	//~ for(int j=0;j<nargs;j++)
		//~ args[j]->set(vals[j]);

	//~ rew.list_vals(cout);
	
	for (int ie=0;ie<n;ie++)
	{
		t1->GetEntry(ie);
		REW.init(e->par);
		nucleus t(e->par);
		if(ie==0)
			SetupSPP(e->par);
			
		double nominal=REW.weight(*e,e->par,t);

		for(int j=0;j<nargs;j++)
			args[j]->set(vals[j]);
		ff_configure(e->par);
        //cout<<rew.qel_cc_axial_mass.val<<endl;
        
        
		weight=REW.weight(*e,e->par,t)/nominal;
		cout<<weight<< " ";//endl;
		t2->Fill();
		e->weight*=weight;
		t3->Fill();
	}
	
	f1->Close();
	delete f1;
	f2->Write();
	f2->Close();
	f3->Write();
	f3->Close();
	delete e;
	delete f2;
	delete f3;
	cout<<"Output file: \""<<outname<<"\""<<endl;	
	return 0;
}
