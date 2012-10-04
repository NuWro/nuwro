#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include "jednostki.h"
#include "data.h"

using namespace std;

int main()
{	
	double E[3][20] = {{0}};
	
	double Enorm[3] = {0};
	
	TFile *tf1 = new TFile("eventsout.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < 100000; i++)
	{
		tt1->GetEntry(i);
		
		double En = 0;
		
		for (int k = 0; k < e1->post.size(); k++)
		{
			if(e1->post[k].nucleon())
				En += e1->post[k].Ek();
			else
				En += e1->post[k].E();
		}
				
		int a = (En-650)/20;
		
		if (a >= 20) a = 19;
		else if (a < 0) a = 0;
				
		int ile = e1->fof(2212) + e1->fof(2112) - 1;
		
		if (ile >= 0 && ile <= 2)
		{
			E[ile][a]++;
			Enorm[ile]++;
		}
		
		cout << 100*i/100000 << "%\r" << flush;
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream file("energy_1.txt");
	
	for (int i = 0; i < 20; i++)
		file << i*20 + 650 << " " << E[0][i]/Enorm[0] << " " << E[1][i]/Enorm[1] << " " << E[2][i]/Enorm[2] << endl;
		
	file.close();
	return 1;
}
