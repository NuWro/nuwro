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
	int n[10] = {0};
	int norma = 0;
	
	double dE0[100] = {0};
	double dE1[100] = {0};
	double dE2[100] = {0};
	double dE3[100] = {0};
	double dE4[100] = {0};
		
	TFile *tf1 = new TFile("kaskada.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < 100000; i++)
	{
		tt1->GetEntry(i);
		
		double E = 0;
		
		if (e1->fof(211) + e1->fof(-211) + e1->fof(111) == 0)
		{
			
			int nofn = e1->fof(2212) + e1->fof(2112);
			
			n[nofn]++;
			norma++;
			for (int k = 0; k < e1->post.size(); k++)
			{
				if(e1->post[k].nucleon())
					E += e1->post[k].Ek();
				else
					E += e1->post[k].E();
			}
			
			double temp = e1->out[0].E() - E;
					
			int a = temp/2;
			
			if (a > 100) a = 100;
			else if (a < 0) a = 0;
			
			int ile = e1->fof(2112) + e1->fof(2212);
			
			switch(ile)
			{
				case 0:	dE0[a]++; break;
				case 1:	dE1[a]++; break;
				case 2:	dE2[a]++; break;
				case 3:	dE3[a]++; break;
				default: dE4[a]++; break;
			}		
			
			cout << 100*i/100000 << "%\r" << flush;
		}
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream file("energy_test.txt");
	
	for (int i = 0; i < 100; i++)
	{
		file << i*2 << " " << dE0[i] << " " << dE1[i] << " " << dE2[i] << " " << dE3[i] << " " << dE4[i] << endl;
		file << (i+1)*2 << " " << dE0[i] << " " << dE1[i] << " " << dE2[i] << " " << dE3[i] << " " << dE4[i] << endl;
	}
	
	file.close();
	
	params p;
	p.nucleus_p = 6;
	p.nucleus_n = 6;
	
	nucleus* nucl= make_nucleus(p);
	
	ofstream plik("ef.txt");
			
	double const M = 0.5*(PDG::mass_proton+PDG::mass_neutron);

	for (int i = 0; i < 100; i++)
	{
		double kmom = nucl -> localkf_(2212,i*0.1*fm);
		plik << i*0.1 << " " << sqrt(kmom*kmom + M*M) - M << endl;
	}
		
	plik.close();
		
	for (int i = 0; i < 10; i++)
		cout << i << ": " << 100.0*n[i]/norma << endl;

	
	return 1;
}
