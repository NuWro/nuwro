#include <iostream>

#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;
using namespace NUWRO;
int main (int argc, char* argv[])
{

	event *e=new event;
	TFile * f2 = new TFile("out.root","recreate");
	TTree * t2 = new TTree("treeout","Tree of events");
	t2->Branch("e","event",&e);

	for(int j=1;j<argc;j++)
	{
		TFile * f1 = new TFile(argv[j]);
		TTree * t1 = (TTree*)f1->Get("treeout");
		t1->SetBranchAddress("e", &e);
		int n = t1->GetEntries();
		for (int i=0;i<n;i++)
		{
			if(i%50 == 0) printf("%s: Event: %d\n",argv[j],i);
			t1->GetEntry(i);
			t2->Fill();
		}
		f1->Close();
		delete f1;
	}
	f2->Write();
	f2->Close();
	delete e;
	delete f2;
	cout<<"Plik wynikowy: \"out.root\""<<endl;
    return 0;
}
