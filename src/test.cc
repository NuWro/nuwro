#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include "jednostki.h"
#include "data.h"
#include "nucleus.h"
#include "nucleusmaker.h"

using namespace std;

int main()
{	

	TFile *tf1 = new TFile("test.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();

	tt1->SetBranchAddress("e",&e);

	double E[2][50] = {{0}};

	for (int i = 0; i < 100000; i++)
	{
		tt1->GetEntry(i);
		
		double energy = 0;
		
		for (int k = 0; k < e->post.size(); k++)
			if (e->post[k].nucleon())
				energy += e->post[k].Ek();
			else
				energy += e->post[k].E();
		
		int bin = energy / 10;
					
		int ktory = 0;
		
		if (e->out.size() == 3)
			ktory = 1;
				
		if (bin >= 0 and bin < 50)		
			E[ktory][bin]++;
		
		if (e->in[0].E() < energy)
		{
			cout << "WARNING!!!!\n";
			cout << "Energy = " << energy << endl;
			cout << "fof(2212) = " << e->fof(2212) << endl;
			cout << "fof(2112) = " << e->fof(2112) << endl;
			cout << "nof(2212) = " << e->nof(2212) << endl;
			cout << "nof(2112) = " << e->nof(2112) << endl;
			cout << "Energy out[0] = " << e->out[0].E() << endl;
			cout << "Energy K out[1] = " << e->out[1].Ek() << endl;
			cout << "number_of_interactions() = " << e->number_of_interactions() << endl;
			cout << "number of pions = " << e->fof(211) + e->fof(-211) + e->fof(111) << endl;
			cout << "Tk (his_fermi) = "; for (int i = 1; i < e->post.size(); i++) cout << e->post[i].Ek() << " (" << e->post[i].his_fermi << ") "; cout << endl;
			cout << "--------------------------------------\n";
		}
		cout << 100*i/100000 << "%\r" << flush;
	}
		
	delete e;
	delete tt1;
	delete tf1;
	
	for (int i = 0; i < 50; i++)
		cout << i*10 << " " << E[0][i] << " " << E[1][i] << endl << (i+1)*10 << " " << E[0][i] << " " << E[1][i] << endl;
		
	return 1;
}
