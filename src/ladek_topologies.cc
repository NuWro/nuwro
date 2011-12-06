#include <iostream>
#include <stdio.h>
#include <vector>
#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "pdg.h"
#include "nucleus.h"
#include <fstream>
#include "args.h"
#include "dirs.h"

int main(int argc, char **argv)
{
	int Npi0, Npip, Npim, Npi0_fsi, Npip_fsi, Npim_fsi;
	int before_fsi, after_fsi, x, y;
	int pions[12][12];
	
	for (int i = 0; i<12; i++)
	{
		for (int j = 0; j<12; j++)
		{
			pions[i][j] = 0;
		}	
	}
	
	/*
	0 - no pions
	1 - pi0
	2 - pi+
	3 - pi-
	4 - 2pi0
	5 - 2pi+
	6 - 2pi-
	7 - pi0pi+
	8 - pi0pi-	 
	9 - pi+pi-
	10 - > 2 pions
	 
	example: pions[2][3] means that before fsi was pi+ and after fsi was pi-
	*/
	set_dirs(argv[0]);    
    args a;
    a.read (argc, argv);	
    event *e = new event;
    TFile *f = new TFile (a.output);
    TTree *t = (TTree *) f->Get ("treeout");
    t->SetBranchAddress ("e", &e);
    int n = t->GetEntries ();
   
    for (int i = 0; i < n; i++)
    {   delete e;
        e=new event;
    	t->GetEntry (i);
    	Npi0 = 0; Npip = 0; Npim = 0; Npi0_fsi = 0; Npip_fsi = 0; Npim_fsi = 0;
        
    	for (int l = 0; l<e->out.size(); l++)
	    {	
			switch(e->out[l].pdg)
			{
				case 211: Npip++; break;
				case -211: Npim++; break;
				case 111: Npi0++; break;
			}
		}
	
		for (int l = 0; l<e->post.size(); l++)
		{
			switch(e->post[l].pdg)
			{
				case 211: Npip_fsi++; break;
				case -211: Npim_fsi++; break;
				case 111: Npi0_fsi++; break;
			}
		}
		
		before_fsi = 100*Npim + 10*Npip + Npi0;
		after_fsi  = 100*Npim_fsi + 10*Npip_fsi + Npi0_fsi;
		
		switch (before_fsi)
		{
			case 0: x = 0; break;
			case 1: x = 1; break;
			case 10: x = 2; break;
			case 100: x = 3; break;
			case 2: x = 4; break;
			case 20: x = 5; break;
			case 200: x = 6; break;
			case 11: x = 7; break;
			case 101: x = 8; break;
			case 110: x = 9; break;
			default: x = 10; break;
		}
		
		switch (after_fsi)
		{
			case 0: y = 0; break;
			case 1: y = 1; break;
			case 10: y = 2; break;
			case 100: y = 3; break;
			case 2: y = 4; break;
			case 20: y = 5; break;
			case 200: y = 6; break;
			case 11: y = 7; break;
			case 101: y = 8; break;
			case 110: y = 9; break;
			default: y = 10; break;
		}
		
		pions[x][y]++;
		
	cout<<i<<" event ready.\r";
	}
	
	for (int i = 0; i<11; i++)
	{
		for (int j = 0; j<11; j++)
		{
			pions[i][11] = pions[i][11] + pions[i][j];
			pions[11][i] = pions[11][i] + pions[j][i];
		}
	}    
    
    for (int i = 0; i<11; i++)
    {
    	pions[11][11] = pions[11][11] + pions[11][i];
	}
    
    ofstream plik;
    plik.open("ladek.tex");
	plik<<"\\documentclass{article}"<<endl<<endl<<"\\usepackage[a4paper,landscape]{geometry}"<<endl<<"\\addtolength{\\leftskip}{-1.0in}"<<endl<<endl<<"\\begin{document}"<<endl<<"\\begin{center}"<<endl<<"\\begin{tabular}{||c||c|c|c|c|c|c|c|c|c|c|c||c||}"<<endl<<"\\hline"<<endl;
    plik<<" & \\multicolumn{11}{|c||}{Primary Hadronic System} & \\\\"<<endl<<"\\hline"<<endl;
    plik<<"Final State & $0\\pi$ & $\\pi^0$ & $\\pi^+$ & $\\pi^-$ & $2\\pi^0$ & $2\\pi^+$ & $2\\pi^-$ & $\\pi^0\\pi^+$ & $\\pi^0\\pi^-$ & $\\pi^+\\pi^-$ & $\\geq 3 \\pi$ & Total \\\\"<<endl<<"\\hline\\hline"<<endl;
    plik<<"$0\\pi$ & "<<pions[0][0]<<" & "<<pions[1][0]<<" & "<<pions[2][0]<<" & "<<pions[3][0]<<" & "<<pions[4][0]<<" & "<<pions[5][0]<<" & "<<pions[6][0]<<" & "<<pions[7][0]<<" & "<<pions[8][0]<<" & "<<pions[9][0]<<" & "<<pions[10][0]<<" & "<<pions[11][0]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^0$ & "<<pions[0][1]<<" & "<<pions[1][1]<<" & "<<pions[2][1]<<" & "<<pions[3][1]<<" & "<<pions[4][1]<<" & "<<pions[5][1]<<" & "<<pions[6][1]<<" & "<<pions[7][1]<<" & "<<pions[8][1]<<" & "<<pions[9][1]<<" & "<<pions[10][1]<<" & "<<pions[11][1]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^+$ & "<<pions[0][2]<<" & "<<pions[1][2]<<" & "<<pions[2][2]<<" & "<<pions[3][2]<<" & "<<pions[4][2]<<" & "<<pions[5][2]<<" & "<<pions[6][2]<<" & "<<pions[7][2]<<" & "<<pions[8][2]<<" & "<<pions[9][2]<<" & "<<pions[10][2]<<" & "<<pions[11][2]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^-$ & "<<pions[0][3]<<" & "<<pions[1][3]<<" & "<<pions[2][3]<<" & "<<pions[3][3]<<" & "<<pions[4][3]<<" & "<<pions[5][3]<<" & "<<pions[6][3]<<" & "<<pions[7][3]<<" & "<<pions[8][3]<<" & "<<pions[9][3]<<" & "<<pions[10][3]<<" & "<<pions[11][3]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$2\\pi^0$ & "<<pions[0][4]<<" & "<<pions[1][4]<<" & "<<pions[2][4]<<" & "<<pions[3][4]<<" & "<<pions[4][4]<<" & "<<pions[5][4]<<" & "<<pions[6][4]<<" & "<<pions[7][4]<<" & "<<pions[8][4]<<" & "<<pions[9][4]<<" & "<<pions[10][4]<<" & "<<pions[11][4]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$2\\pi^+$ & "<<pions[0][5]<<" & "<<pions[1][5]<<" & "<<pions[2][5]<<" & "<<pions[3][5]<<" & "<<pions[4][5]<<" & "<<pions[5][5]<<" & "<<pions[6][5]<<" & "<<pions[7][5]<<" & "<<pions[8][5]<<" & "<<pions[9][5]<<" & "<<pions[10][5]<<" & "<<pions[11][5]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$2\\pi^-$ & "<<pions[0][6]<<" & "<<pions[1][6]<<" & "<<pions[2][6]<<" & "<<pions[3][6]<<" & "<<pions[4][6]<<" & "<<pions[5][6]<<" & "<<pions[6][6]<<" & "<<pions[7][6]<<" & "<<pions[8][6]<<" & "<<pions[9][6]<<" & "<<pions[10][6]<<" & "<<pions[11][6]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^0\\pi^+$ & "<<pions[0][7]<<" & "<<pions[1][7]<<" & "<<pions[2][7]<<" & "<<pions[3][7]<<" & "<<pions[4][7]<<" & "<<pions[5][7]<<" & "<<pions[6][7]<<" & "<<pions[7][7]<<" & "<<pions[8][7]<<" & "<<pions[9][7]<<" & "<<pions[10][7]<<" & "<<pions[11][7]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^0\\pi^-$ & "<<pions[0][8]<<" & "<<pions[1][8]<<" & "<<pions[2][8]<<" & "<<pions[3][8]<<" & "<<pions[4][8]<<" & "<<pions[5][8]<<" & "<<pions[6][8]<<" & "<<pions[7][8]<<" & "<<pions[8][8]<<" & "<<pions[9][8]<<" & "<<pions[10][8]<<" & "<<pions[11][8]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\pi^+\\pi^-$ & "<<pions[0][9]<<" & "<<pions[1][9]<<" & "<<pions[2][9]<<" & "<<pions[3][9]<<" & "<<pions[4][9]<<" & "<<pions[5][9]<<" & "<<pions[6][9]<<" & "<<pions[7][9]<<" & "<<pions[8][9]<<" & "<<pions[9][9]<<" & "<<pions[10][9]<<" & "<<pions[11][9]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"$\\geq 3 \\pi$ & "<<pions[0][10]<<" & "<<pions[1][10]<<" & "<<pions[2][10]<<" & "<<pions[3][10]<<" & "<<pions[4][10]<<" & "<<pions[5][10]<<" & "<<pions[6][10]<<" & "<<pions[7][10]<<" & "<<pions[8][10]<<" & "<<pions[9][10]<<" & "<<pions[10][10]<<" & "<<pions[11][10]<<" \\\\ "<<endl<<"\\hline\\hline"<<endl;
    plik<<"Total & "<<pions[0][11]<<" & "<<pions[1][11]<<" & "<<pions[2][11]<<" & "<<pions[3][11]<<" & "<<pions[4][11]<<" & "<<pions[5][11]<<" & "<<pions[6][11]<<" & "<<pions[7][11]<<" & "<<pions[8][11]<<" & "<<pions[9][11]<<" & "<<pions[10][11]<<" & "<<pions[11][11]<<" \\\\ "<<endl<<"\\hline"<<endl;
    plik<<"\\end{tabular}"<<endl<<"\\end{center}"<<endl<<"\\end{document}";
    
        
    plik.close();
    delete e;
    delete f;
//    delete t;
   
    return 0; 
}
