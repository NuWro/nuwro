#include <iostream>
#include <string>

#include "fsi.h"
#include "calculations.h"
#include "event1.h"
#include <TFile.h>
#include <TTree.h>

using namespace std;

const int nofb = 10;
const int nofa = 40;
const int nofp = 15;
const int events = 100000;

void multi (std::string filename, double vivi[][nofa][nofp], double norm[][nofa], double counter[][5])
{
	for (int l = 0; l < nofb; l++)
	{		
		for (int m = 0; m < nofa; m++)
		{	
			zero(vivi[l][m], nofp);
		}
	}
	
	for (int l = 0; l < nofb; l++) zero(counter[l], 5);

	for (int z = 0; z < nofb; z++) zero(norm[z], nofa);
		
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);

		int przed = e1->n();
		int po = e1->f();
				
		double table[nofp];
		
		int help;
		
		if (przed > 10) help = 11;
		else help = przed - 1;
		
		if (help != 11)
		{
			for (int j = 0; j < przed; j++)
			{
				switch (e1->out[j].pdg)
				{
					case 2212: counter[help][0]++; break;
					case 2112: counter[help][1]++; break;
					case  211: counter[help][2]++; break;
					case -211: counter[help][3]++; break;
					case  111: counter[help][4]++; break;
					default: break;
				}
			}
		}		
		
		table[0]  = e1->number_of_nucleon_elastic();
		table[1]  = e1->number_of_nucleon_spp();
		table[2]  = e1->number_of_nucleon_dpp();
		table[3]  = e1->number_of_pion_elastic();
		table[4]  = e1->number_of_pion_ce();
		table[5]  = e1->number_of_pion_spp();
		table[6]  = e1->number_of_pion_dpp();
		table[7]  = e1->number_of_pion_abs();
		table[8]  = e1->fof(2212); //protons
		table[9]  = e1->fof(2112); //neutrons
		table[10] = e1->fof(211);
		table[11] = e1->fof(-211);
		table[12] = e1->fof(111);
		table[13] = 0;
		table[14] = 0;
		
		int pomoc = table[8] + table[10] + table[11];
		
		for (int k = 0; k < po; k++)
		{
			if (e1->post[k].pdg == 2212 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 370 && e1->post[k].momentum() < 700) table[13]++;
			else if (e1->post[k].pdg == -211 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 350 && e1->post[k].momentum() < 800) table[14]++;
		}
		
		for (int k = 0; k < nofb; k++)
		{
			if (k+1 == przed)
			{
				for (int l = 0; l < nofa; l++)
				{
					if (l+1 == pomoc)  //<-change it
 					{
						norm[k][l]++;
						
						for (int x = 0; x < nofp; x++)
						{
							vivi[k][l][x] += table[x];
						}
					}
				}
			}
		}

		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
		
	cout<<filename<<": done"<<endl<<endl;
	
	for (int i = 0; i < nofb; i++)
	{
		
		for (int k = 0; k < nofa; k++)
		{
			if (norm[i][k] != 0)
			{
				for (int l = 0; l < nofp; l++)
				{
					vivi[i][k][l] /= norm[i][k];
				}
			}
		}
	}
	
	double norm2[nofb];
	
	zero(norm2, nofb);
	
	for (int k = 0; k < nofb; k++)
	{
		for (int l = 0; l < nofa; l++)
		{
			norm2[k] += norm[k][l];
		}
	}
	
	for (int i = 0; i < nofb; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			if (norm2[i] != 0) counter[i][j] /= norm2[i];
		}
	}
}

void maketex (string file, double vivi[][nofa][nofp], double norm[][nofa], double counter[][5], string name, int fz)
{
	ofstream plik(file.c_str());
		
	plik<<setprecision (2);
	plik<<"\\documentclass[titlepage]{article}"<<endl<<endl<<"\\usepackage[margin = 0in, tmargin=0.5in, landscape]{geometry}"<<endl<<"\\pagestyle{empty}"<<endl;
	plik<<"\\title{"<<name<<"}"<<endl;
	plik<<"\\author{"<<fzname[fz]<<"}"<<endl;	
	plik<<endl<<"\\begin{document}"<<endl<<"\\maketitle"<<endl;
	
	int ile = 0;
	
	if (strcmp(name.c_str(), "Nomad") == 0) ile = nofb;
	else if (strcmp(name.c_str(), "Pion (3GeV) - Barium") == 0 or strcmp(name.c_str(), "Proton (1GeV) - Barium") == 0) ile = 1;
	
	for (int j = 0; j < ile; j++)
	{
		plik<<"\\begin{table}[!ht]"<<endl;
		plik<<"\\begin{center}"<<endl;
		plik<<"\\begin{tabular}{|c||c||c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}"<<endl<<"\\hline"<<endl;
		plik<<"\\#particles & \\%events & \\#NEL & \\#NSPP & \\#NDPP & \\#PEL & \\#PCEX & \\#PSPP & \\#PDPP & \\#PABS & \\#proton & \\#neutrons & \\#$\\pi^{+}$ & \\#$\\pi^{-}$ & \\#$\\pi^{0}$ & \\#Bp & \\#B$\\pi^{-}$ \\\\"<<endl;
		plik<<"\\hline \\hline"<<endl;
		
		for (int i = 0; i < nofa; i++)
		{
			plik<<i+1<<" & "<<100.0*norm[j][i]/events<<"\\% ";
			
			for (int k = 0; k < nofp; k++)
			{
				plik<<"& "<<vivi[j][i][k]<<" ";
			}
			
			plik<<"\\\\ \\hline"<<endl;
		}
		
		plik<<"\\end{tabular}"<<endl<<"\\end{center}"<<endl<<"\\caption{Number of particles before cascade: "<<j+1<<" ("<<counter[j][0]<<"*protons, "<<counter[j][1]<<"*neutrons, "<<counter[j][2]<<"*$\\pi^{+}$, "<<counter[j][3]<<"*$\\pi^{-}$, "<<counter[j][4]<<"*$\\pi^{0}$)}";
		plik<<endl<<"\\end{table}"<<endl<<endl<<"\\newpage"<<endl<<endl;
	}
	
	plik<<"\\end{document}"<<endl;

	plik.close();
}
	

int viviNomad (int fz)
{
	cout<<endl<<endl<<"Making vivisection for Nomad ("<<fzname[fz]<<"): "<<endl<<endl;
	
	string file = "root_files/Nomad_100k_";		
	file += fzwork[fz] + string("*.root");
	file = find_last(file);
	
	if (noFile(file))
	{
		logfile<<"Vivisection - Nomad ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"Vivisection - Nomad ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}

	get_date();
	
	double vivi[nofb][nofa][nofp];  //number of particles before cascade - number of particles after cascade - number of processes
	double norm[nofb][nofa];
	double counter[nofb][5];

	multi(file, vivi, norm, counter);
	
	string help = "vivisection/Nomad/";
	run(string("mkdir -p ") + help);
	help += string("tex/");
	run(string("mkdir -p ") + help);
	help += string("nomad_") + fzwork[fz] + string(".tex");
	
	maketex(help, vivi, norm, counter, string("Nomad"), fz);
		
	string pomoc = "nomad_";
	pomoc += fzwork[fz];
	
	string dvifile = pomoc + string(".dvi");
	string psfile  = pomoc + string(".ps");
	string pdffile = pomoc + string(".pdf");

	run(string("latex ") + help);
	run(string("dvips ") + dvifile);
	run(string("ps2pdf ") + psfile);
	run(string("mv ") + pdffile + string(" vivisection/Nomad/")); 
	run(string("mv nomad_* vivisection/Nomad/tex/")); 
	
	return 1;
}

int viviMultiplicity (int fz)
{
	cout<<endl<<endl<<"Making vivisection for Multiplicity ("<<fzname[fz]<<"): "<<endl<<endl;
	
	string file = "root_files/multiplicity_pip_3GeV_barium_";		
	file += fzwork[fz] + string("*.root");
	file = find_last(file);

	string file2 = "root_files/multiplicity_proton_1GeV_barium_";		
	file2 += fzwork[fz] + string("*.root");
	file2 = find_last(file2);
	
	if (noFile(file))
	{
		logfile<<"Vivisection - Multiplicity ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"Vivisection - Multiplicity ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}

	get_date();
	
	double vivi[nofb][nofa][nofp];  //number of particles before cascade - number of particles after cascade - number of processes
	double norm[nofb][nofa];
	double counter[nofb][5];
	
	double vivi2[nofb][nofa][nofp];  //number of particles before cascade - number of particles after cascade - number of processes
	double norm2[nofb][nofa];
	double counter2[nofb][5];
	
	multi(file, vivi, norm, counter);
	multi(file2, vivi2, norm2, counter2);
	
	string help = "vivisection/Multiplicity/";
	run(string("mkdir -p ") + help);
	help += string("tex/");
	run(string("mkdir -p ") + help);
	
	string help2 = help + string("multiplicity_proton_") + fzwork[fz] + string(".tex");
	
	help += string("multiplicity_pion_") + fzwork[fz] + string(".tex");
	
	maketex(help, vivi, norm, counter, string("Pion (3GeV) - Barium"), fz);
	maketex(help2, vivi2, norm2, counter2, string("Proton (1GeV) - Barium"), fz);
		
	string pomoc = "multiplicity_pion_";
	string pomoc2 = "multiplicity_proton_";
	
	pomoc += fzwork[fz];
	pomoc2 += fzwork[fz];
	
	string dvifile = pomoc + string(".dvi");
	string psfile  = pomoc + string(".ps");
	string pdffile = pomoc + string(".pdf");
	string dvifile2 = pomoc2 + string(".dvi");
	string psfile2 = pomoc2 + string(".ps");
	string pdffile2 = pomoc2 + string(".pdf");

	run(string("latex ") + help);
	run(string("dvips ") + dvifile);
	run(string("ps2pdf ") + psfile);
	
	run(string("latex ") + help2);
	run(string("dvips ") + dvifile2);
	run(string("ps2pdf ") + psfile2);
		
	run(string("mv ") + pdffile + string(" vivisection/Multiplicity/")); 
	run(string("mv ") + pdffile2 + string(" vivisection/Multiplicity/")); 

	run(string("mv multiplicity_* vivisection/Multiplicity/tex/"));
	
	return 1;
}
