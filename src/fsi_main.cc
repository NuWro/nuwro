#include "fsi.h"
#include "plots.h"
#include "mplots.h"
#include "calculations.h"
#include "vivisection.h"
#include "simulations.h"
#include "event1.h"
#include "dirs.h"

#include <TFile.h>
#include <TTree.h>

void make_simulations (int expr, int fz, int xs)
{
	switch (expr){
		case 0: K2K(fz, xs); break;
		case 1: MB(fz, xs, 0);break; //1); break;
		case 2: MB(fz, xs, 0); break;
		case 3: PNS(fz, xs); break;
		case 4: PrT(fz, xs, (char *)"high_energy"); break;
		case 5: PrT(fz, xs, (char *)"low_energy"); break;
		case 6: PiT(fz, xs, (char *)"high_energy"); break;
		case 7: PiT(fz, xs, (char *)"low_energy"); break;
		case 8: Nomad(fz, xs); break;
		case 9: AtmNu(fz, xs); break;
		case 10: Multiplicity(fz, xs); break;
		case 11: MBCC(fz, xs); break;
		case 12: simMBCCtotal(fz, xs); break;
		case 13: simSBCCtotal(fz, xs); break;
		case 14: simNOMADCCtotal (fz, xs); break;
		case 15: simMINOSCCtotal (fz, xs); break;
		case 16: simMBCCratio(fz, xs); break;
		case 17: break;//simnuintTr(fz, xs); break;
		case 18: simMBback(fz, xs, 1); break;
		default: break;
	}
}

void make_calculations (int expr, int fz, int xs)
{
	switch (expr){
		case 0: calcK2K(fz, xs); break;
		case 1: {calcMB(fz, xs, 0);} break;// calcMB(fz, xs, 1);} break;
		case 2: calcSB(fz, xs); break;
		case 3: calcPNS(fz, xs); break;
		case 4: calcPrThe(fz, xs); break;
		case 5: calcPrTle(fz, xs); break;
		case 6: calcPiThe(fz, xs); break;
		case 7: calcPiTle(fz, xs); break;
		case 8: calcNomad(fz, xs); break;
		case 9: calcAtmNu(fz, xs); break;
		case 11: calcMBCC(fz, xs); break;
		case 12: calcMBCCtotal(fz, xs); break;
		case 13: calcSBCCtotal(); break;
		case 14: calcNOMADCCtotal(); break;
		case 15: calcMINOSCCtotal(); break;
		case 16: calcMBCCrat(fz, xs); break;
		case 17: calcnuintTr(fz, xs); break;
		case 18: {calcMBback(fz, xs, 0); calcMBback(fz, xs, 1);} break;
		default: break;
	}
}

void make_plots (int expr, int fz, int xs)
{
	switch (expr){
		case 0: plotK2K(fz, xs); break;
		case 1: plotMB(fz, xs); break;
		case 2: plotSB(fz, xs); break;
		case 3: plotPNS(fz, xs); break;
		case 4: plotPrThe(fz, xs); break;
		case 5: plotPrTle(fz, xs); break;
		case 6: plotPiThe(fz, xs); break;
		case 7: plotPiTle(fz, xs); break;
		case 8: plotNomad(fz, xs); break;
		case 9: plotAtmNu(fz, xs); break;
		case 11: plotMBCC(fz, xs); break;
		case 12: plotMBCCtotal(fz, xs); break;
		case 13: break;
		case 14: break;
		case 15: break;
		case 16: plotMBCCrat(fz, xs); break;
		case 17: plotnuintPr(fz, xs); break;
		default: break;
	}
}

void make_mplots (int expr)
{
	switch (expr){
		case 0: mplotK2K(); break;
    	case 1: mplotMB(); break;
		case 2: mplotSB(); break;
		case 3: mplotPNS(); break;
		case 4: mplotPrThe(); break;
		case 5: mplotPrTle(); break;
		case 6: mplotPiThe(); break;
		case 7: mplotPiTle(); break;
		case 8: mplotNomad(); break;
		case 9: mplotAtmNu(); break;
		case 11: mplotMBCC(); break;
		case 12: break;
		case 13: break;
		case 14: break;
		case 15: break;
		case 16: break;
		case 17: break;
		default: break;
	}
}

void make_vivisection (int expr, int fz)
{
	switch (expr){
		case 0: break;
		case 1: break;
		case 2: break;
		case 3: break;
		case 4: break;
		case 5: break;
		case 6: break;
		case 7: break;
		case 8: viviNomad(fz); break;
		case 9: break;
		case 10: viviMultiplicity(fz); break;
		case 11: break;
		case 12: break;
		case 13: break;
		case 14: break;
		case 15: break;
		case 16: break;
		case 17: break;
		default: break;
	}
}

void poster()
{
	const int events = 1000000;
	const int bins = 20;
	const double min = 0;
	const double max = 900;
	double mom[bins];

	double counter[6][bins]; //0 - pi+ before fsi, 1 - 0pi, 2 - pi+, 3 - pi-, 4 - pi0, 5 - more pi
	
	for	(int i = 0; i < 6; i++) zero(counter[i], bins);
	
	for (int i = 0; i < bins; i++) mom[i] = min + (i+1.0)*(max - min)/bins;
		
	TFile *tf1 = new TFile("eventsout.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
		
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion    = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		int pionfsi = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		 
		if (pion == 100)
		{
			double val;
			double help;
			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg==211)
				{
					val = e1->out[k].momentum();
				}
			}
			
			put(val, mom, counter[0], help, bins);
			
			if (pionfsi == 0) put(val, mom, counter[1], help, bins);
			else if (pionfsi == 100) put(val, mom, counter[2], help, bins);
			else if (pionfsi == 10) put(val, mom, counter[3], help, bins);
			else if (pionfsi == 1) put(val, mom, counter[4], help, bins);
			else put(val, mom, counter[5], help, bins);			
		}
	}
		
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int i = 0; i < bins; i++) mom[i] = min + (i+0.5)*(max - min)/bins;
	
	ofstream file("pitop.txt");
	
	for (int i = 0; i < bins; i++)
	{
		file<<mom[i];
		for (int k = 0; k < 6; k++) file<<" "<<counter[k][i];
		file<<endl;
	}
	file.close();
}

int main(int argc, char** argv)
{
    set_dirs(argv[0]);
	cout.precision(2);
	
	//poster();

	int space = 50 - allinone.length();
	
	for (int i = 0; i < nof_expr; i++) expr_on[i] = 0;
	for (int i = 0; i < nof_fz;  i++) fz_on[i]  = 0;
	for (int i = 0; i < nof_opt; i++) options_on[i] = 0;
	for (int i = 0; i < nof_xsec; i++) xsec_on[i] = 0;

	char type;
	
	do{	
		if(system("clear"));
		cout<<"Select comparisions you want to do (type '0' to select all of them). Press 'enter' to go to the next options."<<endl<<endl;
			
		for (int i = 0; i < nof_expr; i++)
		{
			int spaces = 50 - expr[i].length();
			
			cout<<((char)(97+i))<<" - "<<expr[i]; for (int k = 0; k < spaces; k++) cout<<" "; if (expr_on[i]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		cout<<endl<<((char)(97+nof_expr))<<" - "<<allinone; for (int k = 0; k < space; k++) cout<<" "; if (allinone_on) cout<<"*"; cout<<endl; 
		
		type = getch();
		
		if (type == '0') for (int i = 0; i < nof_expr; i++) expr_on[i] = 1;
		else if (type == ((char)(97+nof_expr))) allinone_on = !allinone_on;
		else for (int i = 0; i < nof_expr; i++) if (type == ((char)(97+i))) expr_on[i] = !expr_on[i];
	
	}while(type != 10);
	
	do{	
		if(system("clear"));
		cout<<"Select xsec model you want to use (type '0' to select all of them). Press 'enter' to go to the next options."<<endl<<endl;
			
		for (int i = 0; i < nof_xsec; i++)
		{
			int spaces = 50 - xsec[i].length();
			
			cout<<((char)(97+i))<<" - "<<xsec[i]; for (int k = 0; k < spaces; k++) cout<<" "; if (xsec_on[i]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		type = getch();
		
		if (type == '0') for (int i = 0; i < nof_xsec; i++) xsec_on[i] = 1;
		else for (int i = 0; i < nof_xsec; i++) if (type == ((char)(97+i))) xsec_on[i] = !xsec_on[i];
	
	}while(type != 10);
	
	do{	
		if(system("clear"));
		cout<<"Select formation zone models you want to use (type '0' to select all of them). Press 'enter' to go to the next options."<<endl<<endl;
			
		for (int i = 0; i < nof_fz; i++)
		{
			int spaces = 50 - fzname[i].length();
			
			cout<<((char)(97+i))<<" - "<<fzname[i];	for (int k = 0; k < spaces; k++) cout<<" ";	if (fz_on[i]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		type = getch();
		
		if (type == '0') for (int i = 0; i < nof_fz; i++) fz_on[i] = 1;
		else for (int i = 0; i < nof_fz; i++) if (type == ((char)(97+i))) fz_on[i] = !fz_on[i];
	
	}while(type != 10);
		
	do{	
		if(system("clear"));
		cout<<"Select operations you want to do (type '0' to select all of them). Press 'enter' to start."<<endl<<endl;
			
		for (int i = 0; i < nof_opt; i++)
		{
			int spaces = 50 - options[i].length();
			
			cout<<((char)(97+i))<<" - "<<options[i]; for (int k = 0; k < spaces; k++) cout<<" "; if (options_on[i]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		cout<<endl<<endl; for (int k = 0; k < 50; k++) cout<<"-";
		cout<<endl<<"Note:"<<endl<<" - you need gnuplot installed to make plots"<<endl<<" - to generate pdf-table you need: latex, dvips and ps2pdf installed"<<endl;
		for (int k = 0; k < 50; k++) cout<<"-";
		cout<<endl<<endl;
		
		type = getch();
		
		if (type == '0') for (int i = 0; i < nof_opt; i++) options_on[i] = 1;
		else for (int i = 0; i < nof_opt; i++) if (type == ((char)(97+i))) options_on[i] = !options_on[i];
	
	}while(type != 10);
	
	logfile.open("compare_log");
	run("mkdir -p tmp/");
	run("mkdir -p root_files/");
	run("mkdir -p results/");
	run("mkdir -p plots/");
	run("mkdir -p multiplots/");
	run("mkdir -p vivisection/");
		
	if (options_on[0])
	{
		for (int i = 0; i < nof_expr; i ++)
		{
			if (expr_on[i])
			{
				for (int k = 0; k < nof_fz; k++)
				{
					if (fz_on[k])
					{
						for (int j = 0; j < nof_xsec; j++)
						{
							if (xsec_on[j])	make_simulations(i, k, j);
						}
					}						
				}
			}
		}
		
		if(system("clear"));
		cout<<"Simulations done!"<<endl<<endl;
	}
	
	if (options_on[1])
	{
		for (int i = 0; i < nof_expr; i ++)
		{
			if (expr_on[i])
			{
				for (int k = 0; k < nof_fz; k++)
				{
					if (fz_on[k])
					{
						for (int j = 0; j < nof_xsec; j++)
						{
							if (xsec_on[j]) make_calculations(i, k, j);
						}
					}
				}
			}
		}
		
		cout<<"Calculations done!"<<endl<<endl;
	}			
	
	if (options_on[2])
	{
		for (int i = 0; i < nof_expr; i ++)
		{
			if (expr_on[i])
			{
				for (int k = 0; k < nof_fz; k++)
				{
					if (fz_on[k])
					{
						for (int j = 0; j < nof_xsec; j++)
						{
							if (xsec_on[j]) make_plots(i, k, j);
						}
					}
				}
			}
		}
		
		cout<<"Plots done!"<<endl<<endl;
	}
	
	if (options_on[3])
	{
		for (int i = 0; i < nof_expr; i ++)
		{
			if (expr_on[i])
			{
				make_mplots(i);
			}
		}
		
		cout<<"Multiplots done!"<<endl<<endl;
	}
	
	if (options_on[4])
	{
		for (int i = 0; i < nof_expr; i ++)
		{
			if (expr_on[i])
			{
				for (int k = 0; k < nof_fz; k++)
				{
					if (fz_on[k]) make_vivisection(i, k);
				}
			}
		}
		
		cout<<"Multiplicity done!"<<endl<<endl;
	}
	
	if (allinone_on) plotFZ();
	
	//simsfg();
	//calcsfg();
	
	//calcbg();
	
	//roman_sim();
	//roman_calc();
	//roman_3b();	
	//roman_3d();
	//roman_7b();
	//roman_extra();
	
	//oset_sim();
	//oset_calc();
	
	//t2k_anal_sim();
	//t2k_anal_calc();
	//t2k_anal_plot();
	
	//hayato_sim1();
	//hayato_sim2();
	//hayato_sim3();
	//hayato_sim4();	
	
	//hayato_calc();
	
	//dens_test_sim();
	//dens_test_calc();
	
	//ptsim();
	//ptcalc();
	
	//angle_test();
	//ang_calc();
	
	viviNomad_new(0, 1);
	//viviNomad_new(8, 1);
	
	run("rm -r tmp/");
	logfile.close();
	return 0;	
}
