#include "niwg.h"
#include "niwg_tech.h"
#include "simulations.h"
#include "vivisection.h"
#include "calculations.h"
#include "event1.h"
#include "dirs.h"
#include <TFile.h>
#include <TTree.h>


double mom_bin[40] = {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000};
double ang_bin[20] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};


void tech_sim()
{
	for (int n = 0; n < 4; n++)
	{
		string sim[3];
		
		sim[0] = get_bin_dir()+"nuwro -o 'NIWG/tech/root_files/" + nuname[n] + "_hydrogen.root' " + par + hydrogen + nuflux[n] + "-p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
		sim[1] = get_bin_dir()+"nuwro -o 'NIWG/tech/root_files/" + nuname[n] + "_oxygen_fg.root' " + par + oxygen_fg + nuflux[n];
		sim[2] = get_bin_dir()+"nuwro -o 'NIWG/tech/root_files/" + nuname[n] + "_oxygen_sf.root' " + par + oxygen_sf + nuflux[n];
		
	//	for (int i = 0; i < 3; i++) run(sim[i]);	
	}
	
	string com = get_bin_dir()+"nuwro -o 'NIWG/tech/root_files/" + nuname[0] + "_hydrogen_stat.root' " + par_stat + hydrogen + nuflux[0] + "-p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
	run(com);
	com = get_bin_dir()+"nuwro -o 'NIWG/tech/root_files/" + nuname[0] + "_oxygen_stat.root' " + par_stat + oxygen_fg + nuflux[0];
	run(com);
	
}

void make_vivi(string filename, double *result, double *extra, bool nucl)
{
	zero(result, 27);
	zero(extra, 32);
	
	cout<<endl; 
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();
		
	tt1->SetBranchAddress("e",&e);

	for (int i = 0; i < 100000; i++)
	{
		tt1->GetEntry(i);
		
		int electron = 0;
		int positron = 0;
		int muon = 0;
		int antimuon = 0;
		int proton = 0;
		int pip = 0;
		int pim = 0;
		int pi0 = 0;
		int neutrino = 0;

		if (e->flag.nc) neutrino++;
		else
		{
			if (e->in[0].pdg == 12) electron++;
			else if (e->in[0].pdg == -12) positron++;
			else if (e->in[0].pdg ==  14 and e->out[0].momentum() > 150) muon++;
			else if (e->in[0].pdg == -14) antimuon++;
		}
		
		int ile = e->n();
		
		if (nucl) ile = e->f();
				
		for (int k = 0; k < ile; k++)
		{
			int pdg = e->out[k].pdg;
			double momentum = e->out[k].momentum();
			
			if (nucl)
			{
				pdg = e->post[k].pdg;
				momentum = e->post[k].momentum();
			}
			
			if (pdg == pdg_proton and momentum > 1090) proton++;
			else if (pdg == pdg_piP and momentum > 180) pip++;
			else if (pdg == -pdg_piP and momentum > 180) pim++;
			else if (pdg == pdg_pi) pi0++;
		}
		
		int pic = pip + pim;
		
		if (muon == 1)
		{
			if (proton + pic + pi0 == 0) result[0]++;
			else if (proton == 1 and pic + pi0 == 0) result[1]++;
			else if (pi0 == 1 and proton + pic == 0) result[2]++;
			else if (pic == 1 and proton + pi0 == 0) result[3]++;
			else if (proton == 1 and pi0 == 1 and pic == 0) result[4]++;
			else if (proton == 1 and pic == 1 and pi0 == 0) result[5]++;
			else if (pi0 == 1)
			{
				result[6]++;
				if (pic == 1 and proton == 0) extra[0]++;
				else if (pic == 1 and proton == 1) extra[1]++;
				else if (pic == 2 and proton == 0) extra[2]++;
				else if (pic == 2 and proton == 1) extra[3]++;
				else extra[4]++;
			}
			else
			{
				result[7]++;
				if (pi0 == 2 and proton + pic == 0) extra[5]++;
				else if (pi0 == 2 and pic == 0 and proton == 1) extra[6]++;
				else if (pic == 2 and proton + pi0 == 0) extra[7]++;
				else if (pic == 2 and pi0 == 0 and proton == 1) extra[8]++;
				else if (pi0 == 2 and pic == 1 and proton == 0) extra[9]++;
				else if (pi0 == 2 and pic == 1 and proton == 1) extra[10]++;
				else if (pi0 == 3 and pic + proton == 0) extra[11]++;
				else if (pi0 == 3 and pic == 0 and proton == 1) extra[12]++;
				else if (pic == 3 and pi0 + proton == 0) extra[13]++;
				else if (pic == 3 and pi0 == 0 and proton == 1) extra[14]++;
				else extra[15]++;
			}			
		}
		else if (antimuon == 1)
		{
			if (proton + pic + pi0 == 0) result[8]++;
			else result[9]++;
		}
		else if (electron == 1)
		{
			if (proton + pic + pi0 == 0) result[10]++;
			else if (proton == 1 and pic + pi0 == 0) result[11]++;
			else if (pi0 == 1 and proton + pic == 0) result[12]++;
			else if (pic == 1 and proton + pi0 == 0) result[13]++;
			else if (proton == 1 and pi0 == 1 and pic == 0) result[14]++;
			else if (proton == 1 and pic == 1 and pi0 == 0) result[15]++;
			else if (pi0 == 1) result[16]++;
			else result[17]++;
		}
		else if (neutrino == 1)
		{
			if (proton + pic + pi0 == 0) result[18]++;
			else if (proton == 1 and pic + pi0 == 0) result[19]++;
			else if (pi0 == 1 and proton + pic == 0) result[20]++;
			else if (pic == 1 and proton + pi0 == 0) result[21]++;
			else if (proton == 1 and pi0 == 1 and pic == 0) result[22]++;
			else if (proton == 1 and pic == 1 and pi0 == 0) result[23]++;
			else if (pi0 == 1)
			{
				result[24]++;
				if (pic == 1 and proton == 0) extra[16]++;
				else if (pic == 1 and proton == 1) extra[17]++;
				else if (pic == 2 and proton == 0) extra[18]++;
				else if (pic == 2 and proton == 1) extra[19]++;
				else extra[20]++;
			}
			else
			{
				result[25]++;
				if (pi0 == 2 and proton + pic == 0) extra[21]++;
				else if (pi0 == 2 and pic == 0 and proton == 1) extra[22]++;
				else if (pic == 2 and proton + pi0 == 0) extra[23]++;
				else if (pic == 2 and pi0 == 0 and proton == 1) extra[24]++;
				else if (pi0 == 2 and pic == 1 and proton == 0) extra[25]++;
				else if (pi0 == 2 and pic == 1 and proton == 1) extra[26]++;
				else if (pi0 == 3 and pic + proton == 0) extra[27]++;
				else if (pi0 == 3 and pic == 0 and proton == 1) extra[28]++;
				else if (pic == 3 and pi0 + proton == 0) extra[29]++;
				else if (pic == 3 and pi0 == 0 and proton == 1) extra[30]++;
				else extra[31]++;
			}
		}
		else result[26]++;
	
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}	
	
	delete e;
	delete tt1;
	delete tf1;

	cout<<endl;
	
}

void distribution(string filename, double mom[2][15][41], double ang[2][15][20], bool nucl)
{
	for (int i = 0; i < 15; i++)
	{
		zero(mom[0][i], 41);
		zero(ang[0][i], 20);
		zero(mom[1][i], 41);
		zero(ang[1][i], 20);
	}
		
	cout<<endl; 
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();
		
	tt1->SetBranchAddress("e",&e);

	for (int i = 0; i < 1000000; i++)
	{
		tt1->GetEntry(i);
		
		int electron = 0;
		int positron = 0;
		int muon = 0;
		int antimuon = 0;
		int proton = 0;
		int pip = 0;
		int pim = 0;
		int pi0 = 0;
		int neutrino = 0;

		if (e->flag.nc) neutrino++;
		else
		{
			if (e->in[0].pdg == 12) electron++;
			else if (e->in[0].pdg == -12) positron++;
			else if (e->in[0].pdg ==  14 and e->out[0].momentum() > 150) muon++;
			else if (e->in[0].pdg == -14) antimuon++;
		}
		
		int ile = e->n();
		
		if (nucl) ile = e->f();
				
		for (int k = 0; k < ile; k++)
		{
			int pdg = e->out[k].pdg;
			double momentum = e->out[k].momentum();
			
			if (nucl)
			{
				pdg = e->post[k].pdg;
				momentum = e->post[k].momentum();
			}
			
			if (pdg == pdg_proton and momentum > 1090) proton++;
			else if (pdg == pdg_piP and momentum > 180) pip++;
			else if (pdg == -pdg_piP and momentum > 180) pim++;
			else if (pdg == pdg_pi) pi0++;
		}
		
		int pic = pip + pim;
		int dyn = 1;
		
		if (neutrino == 0) dyn = 0;
		
		if (pi0 == 1 and pic == 1 and proton == 0)
		{
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				double rest;
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if (pdg == pdg_pi)
				{
					put(momentum, mom_bin, mom[dyn][0], mom[dyn][0][40], 40);
					put(angle, ang_bin, ang[dyn][0], rest, 20);
				}
				else if ((pdg == pdg_piP or pdg == -pdg_piP) and momentum > 180)
				{
					put(momentum, mom_bin, mom[dyn][1], mom[dyn][1][40], 40);
					put(angle, ang_bin, ang[dyn][1], rest, 20);
				}
			}
		}
		else if (pi0 == 1 and pic == 1 and proton == 1)
		{
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				double rest;
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if (pdg == pdg_pi)
				{
					put(momentum, mom_bin, mom[dyn][2], mom[dyn][2][40], 40);
					put(angle, ang_bin, ang[dyn][2], rest, 20);
				}
				else if ((pdg == pdg_piP or pdg == -pdg_piP) and momentum > 180)
				{
					put(momentum, mom_bin, mom[dyn][3], mom[dyn][3][40], 40);
					put(angle, ang_bin, ang[dyn][3], rest, 20);
				}
				else if (pdg == pdg_proton and momentum > 1090)
				{
					put(momentum, mom_bin, mom[dyn][4], mom[dyn][4][40], 40);
					put(angle, ang_bin, ang[dyn][4], rest, 20);
				}
			}
		}
		else if (pi0 == 2 and pic + proton == 0)
		{
			double mom1 = 0;
			double mom2 = 0;
			double ang1 = 0;
			double ang2 = 0;
			double rest = 0;
			int next = 0;
			
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if (pdg == pdg_pi)
				{
					if (next == 0)
					{
						mom1 = momentum;
						ang1 = angle;
						next++;
					}
					else
					{
						mom2 = momentum;
						ang2 = angle;
					}
				}
				
				//if (mom1 >= mom2)
				//{
					put(mom1, mom_bin, mom[dyn][5], mom[dyn][5][40], 40);
					put(mom2, mom_bin, mom[dyn][6], mom[dyn][6][40], 40);
					put(ang1, ang_bin, ang[dyn][5], rest, 20);
					put(ang2, ang_bin, ang[dyn][6], rest, 20);
				//}
				//else
				//{
				//	put(mom1, mom_bin, mom[dyn][6], mom[dyn][6][20], 20);
				//	put(mom2, mom_bin, mom[dyn][5], mom[dyn][5][20], 20);
				//	put(ang1, ang_bin, ang[dyn][6], rest, 20);
				//	put(ang2, ang_bin, ang[dyn][5], rest, 20);
				//}
			}
		}
		else if (pi0 == 2 and pic == 0 and proton == 1)
		{
			double mom1 = 0;
			double mom2 = 0;
			double ang1 = 0;
			double ang2 = 0;
			double rest = 0;
			int next = 0;
			
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if (pdg == pdg_pi)
				{
					if (next == 0)
					{
						mom1 = momentum;
						ang1 = angle;
						next++;
					}
					else
					{
						mom2 = momentum;
						ang2 = angle;
					}
				}
				else if (pdg == pdg_proton and momentum > 1090)
				{
					put(momentum, mom_bin, mom[dyn][9], mom[dyn][9][40], 40);
					put(angle, ang_bin, ang[dyn][9], rest, 20);
				}
				
				//if (mom1 >= mom2)
				//{
					put(mom1, mom_bin, mom[dyn][7], mom[dyn][7][40], 40);
					put(mom2, mom_bin, mom[dyn][8], mom[dyn][8][40], 40);
					put(ang1, ang_bin, ang[dyn][7], rest, 20);
					put(ang2, ang_bin, ang[dyn][8], rest, 20);
				//}
				//else
				//{
				//	put(mom1, mom_bin, mom[dyn][8], mom[dyn][8][20], 20);
				//	put(mom2, mom_bin, mom[dyn][7], mom[dyn][7][20], 20);
				//	put(ang1, ang_bin, ang[dyn][8], rest, 20);
				//	put(ang2, ang_bin, ang[dyn][7], rest, 20);
				//}
			}
		}
		else if (pic == 2 and pi0 + proton == 0)
		{
			double mom1 = 0;
			double mom2 = 0;
			double ang1 = 0;
			double ang2 = 0;
			double rest = 0;
			int next = 0;
			
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if ((pdg == pdg_piP or pdg == -pdg_piP) and momentum > 180)
				{
					if (next == 0)
					{
						mom1 = momentum;
						ang1 = angle;
						next++;
					}
					else
					{
						mom2 = momentum;
						ang2 = angle;
					}
				}
				
				//if (mom1 >= mom2)
				//{
					put(mom1, mom_bin, mom[dyn][10], mom[dyn][10][40], 40);
					put(mom2, mom_bin, mom[dyn][11], mom[dyn][11][40], 40);
					put(ang1, ang_bin, ang[dyn][10], rest, 20);
					put(ang2, ang_bin, ang[dyn][11], rest, 20);
				//}
				//else
				//{
				//	put(mom1, mom_bin, mom[dyn][11], mom[dyn][11][20], 20);
				//	put(mom2, mom_bin, mom[dyn][10], mom[dyn][10][20], 20);
				//	put(ang1, ang_bin, ang[dyn][11], rest, 20);
				//	put(ang2, ang_bin, ang[dyn][10], rest, 20);
				//}
			}
		}
		else if (pic == 2 and pi0 == 0 and proton == 1)
		{
			double mom1 = 0;
			double mom2 = 0;
			double ang1 = 0;
			double ang2 = 0;
			double rest = 0;
			int next = 0;
			
			for (int k = 0; k < ile; k++)
			{
				int pdg = e->out[k].pdg;
				double momentum = e->out[k].momentum();
				double angle = e->out[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
		
				if (nucl)
				{
					pdg = e->post[k].pdg;
					momentum = e->post[k].momentum();
					angle = e->post[k].p()*e->in[0].p()/momentum/e->in[0].momentum();
				}
				
				if ((pdg == pdg_piP or pdg == -pdg_piP) and momentum > 180)
				{
					if (next == 0)
					{
						mom1 = momentum;
						ang1 = angle;
						next++;
					}
					else
					{
						mom2 = momentum;
						ang2 = angle;
					}
				}
				else if (pdg == pdg_proton and momentum > 1090)
				{
					put(momentum, mom_bin, mom[dyn][14], mom[dyn][14][40], 40);
					put(angle, ang_bin, ang[dyn][14], rest, 20);
				}
				
				//if (mom1 >= mom2)
				//{
					put(mom1, mom_bin, mom[dyn][12], mom[dyn][12][40], 40);
					put(mom2, mom_bin, mom[dyn][13], mom[dyn][13][40], 40);
					put(ang1, ang_bin, ang[dyn][12], rest, 20);
					put(ang2, ang_bin, ang[dyn][13], rest, 20);
				//}
				//else
				//{
				//	put(mom1, mom_bin, mom[dyn][13], mom[dyn][13][20], 20);
				//	put(mom2, mom_bin, mom[dyn][12], mom[dyn][12][20], 20);
				//	put(ang1, ang_bin, ang[dyn][13], rest, 20);
				//	put(ang2, ang_bin, ang[dyn][12], rest, 20);
				//}
			}
		}
		
		cout<<filename<<": "<<100*i/1000000<<"%\r"<<flush;
	}	
	
	delete e;
	delete tt1;
	delete tf1;

	cout<<endl;
		
}

void nu_mix(double in[4][3][27], double ex[4][3][32], double *H, double *Ofg, double * Osf, double *fg, double *sf, double *fgex, double *sfex)
{
	for (int i = 0; i < 27; i++)
	{		
		for (int n = 0; n < 4; n++)
		{
			double fermi;
			double spectral;
			
			fermi = in[n][0][i]*H[n]*2.0 + in[n][1][i]*Ofg[n]*16.0;
			spectral = in[n][0][i]*H[n]*2.0 + in[n][2][i]*Osf[n]*16.0;
			
			fg[i] += fermi*nupart[n];
			sf[i] += spectral*nupart[n];
		}
	}
	
	for (int i = 0; i < 32; i++)
	{
		for (int n = 0; n < 4; n++)
		{
			double fex = ex[n][0][i]*H[n]*2.0 + ex[n][1][i]*Ofg[n]*16.0;
			double sex = ex[n][0][i]*H[n]*2.0 + ex[n][2][i]*Osf[n]*16.0;
				
			fgex[i] += fex*nupart[n];
			sfex[i] += sex*nupart[n];
		}
	}
	
}

void nu_norm(double *table, double norm, int bins)
{
	double sum = 0;
	
	for (int i = 0; i < bins; i++) sum += table[i];
	for (int i = 0; i < bins; i++) table[i] *= norm/sum;
	
}

void make_table(double *fg, double *sf, double tot_fg, double tot_sf)
{
	ofstream file("tech_table.tex");
	
	file<<fixed;
	file<<setprecision (1);
	file<<"\\documentclass[titlepage]{article}"<<endl<<endl<<"\\usepackage[margin = 0in, tmargin=0.5in]{geometry}"<<endl<<"\\pagestyle{empty}"<<endl;
	file<<"\\usepackage{array}"<<endl<<"\\newcolumntype{j}{m{50pt}}"<<"\\newcolumntype{R}{>{\\hfill\\arraybackslash}m{2cm}<{}}"<<endl;
	
	file<<endl<<"\\begin{document}"<<endl;

	file<<"\\begin{table}[!ht]"<<endl;
	file<<"\\begin{center}"<<endl;
	file<<"\\begin{tabular}{c|RRRR}"<<endl;
	file<<"FS Paticles & NEUT & GENIE & NuWro (FG) & NuWro(SF) \\\\"<<endl;
	file<<"\\hline \\hline"<<endl;

	for (int l = 0; l < 26; l++)
	{
		file<<nt[l];
		
		if (neut[l] < 0.1) file << " $<$ 0.1\\% & ";
		else file<<neut[l]<<" \\% & ";
		
		if (genie[l] < 0.1) file << " $<$ 0.1\\% & ";
		else file<<genie[l]<<" \\% & ";
		
		if (fg[l] < 0.1) file << " $<$ 0.1\\% & ";
		else file<<fg[l]<<" \\% & ";
		
		if (sf[l] < 0.1) file << " $<$ 0.1\\% \\\\ ";
		else file<<sf[l]<<" \\% \\\\ ";
		
		if (l == 0 or l == 7 or l == 9 or l == 17 or l == 25) file << " \\hline ";
		
		file<<endl;
	}

	file<<"Total & "<<tot_neut<<" & "<<tot_genie<<" & "<<tot_fg<<" & "<<tot_sf<<" \\\\ \\hline \\hline"<<endl; 
	
	file<<"\\end{tabular}"<<endl<<"\\end{center}"<<endl<<"\\end{table}"<<endl<<"\\end{document}"<<endl;
	
	file.close();
	
	run(string("latex tech_table.tex"));
	run(string("dvips tech_table.dvi"));
	run(string("ps2pdf tech_table.ps"));
	run(string("mv tech_table.pdf NIWG/tech/results/")); 
	run(string("mv tech_table* NIWG/tech/results/tex/"));
	
}

void make_table_ex(double *fg, double *sf, double *all_fg, double *all_sf)
{
	
	ofstream file("tech_table_ex.tex");
	
	file<<fixed;
	file<<setprecision (1);
	file<<"\\documentclass[titlepage]{article}"<<endl<<endl<<"\\usepackage[margin = 0in, tmargin=0.5in]{geometry}"<<endl<<"\\pagestyle{empty}"<<endl;
	file<<"\\usepackage{array}"<<endl<<"\\newcolumntype{j}{m{50pt}}"<<"\\newcolumntype{R}{>{\\hfill\\arraybackslash}m{2cm}<{}}"<<endl;
	
	file<<endl<<"\\begin{document}"<<endl;

	file<<"\\begin{table}[!ht]"<<endl;
	file<<"\\begin{center}"<<endl;
	file<<"\\begin{tabular}{l|RR}"<<endl;
	file<<"FS Paticles & NuWro (FG) & NuWro(SF) \\\\"<<endl;
	file<<"\\hline \\hline"<<endl;

	file<<nt[6]<<all_fg[6]<<" \\% & "<<all_sf[6]<<" \\% \\\\ \\hline\\hline"<<endl;
	
	for (int i = 0; i < 32; i++)
	{
		file<<nt_ex[i]<<" & ";
		
		if (fg[i] < 0.1) file << " $<$ 0.1\\% & ";
		else file<<fg[i]<<" \\% & ";
		
		if (sf[i] < 0.1) file << " $<$ 0.1\\% \\\\ ";
		else file<<sf[i]<<" \\% \\\\ ";
						
		if (i == 4 or i == 15 or i == 20) file<<"\\hline\\hline"<<endl;
	
		if (i == 4) file<<nt[7]<<all_fg[7]<<" \\% & "<<all_sf[7]<<" \\% \\\\ \\hline\\hline"<<endl;
	
		if (i == 15) file<<nt[24]<<all_fg[24]<<" \\% & "<<all_sf[24]<<" \\% \\\\ \\hline\\hline"<<endl;
	
		if (i == 20) file<<nt[25]<<all_fg[25]<<" \\% & "<<all_sf[25]<<" \\% \\\\ \\hline\\hline"<<endl;
	}
	
	file<<"\\hline\\hline"<<endl;		
	file<<"\\end{tabular}"<<endl<<"\\end{center}"<<endl<<"\\end{table}"<<endl<<"\\end{document}"<<endl;

	file.close();
	
	run(string("latex tech_table_ex.tex"));
	run(string("dvips tech_table_ex.dvi"));
	run(string("ps2pdf tech_table_ex.ps"));
	run(string("mv tech_table_ex.pdf NIWG/tech/results/")); 
	run(string("mv tech_table_ex* NIWG/tech/results/tex/"));
	
}

void norm_ex(double *tab, double norm, int from, int to)
{
	double sum = 0;
	
	for (int i = from; i <= to; i++)
	{
		sum += tab[i];
	}
	
	for (int i = from; i <= to; i++)
	{
		tab[i] *= norm/sum;
	}
}

void make_mom_plot(string filename, string p1, string p2, string p3, int k, bool nc)
{
	string in  = "NIWG/tech/results/dis/txt/" + filename;
	string out = "NIWG/tech/results/dis/plots/" + filename + ".eps";
	
	bool trzy = false;
	if (strcmp(p3.c_str(), "null")) trzy = true;
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	double ar = height[k];
	if (nc) ar = height[k+6];
	
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<out<<"'"<<endl;
	gfile<<"set xrange[0:2100]"<<endl;
	gfile<<"set yrange[0:"<<ar<<"]"<<endl;	
	gfile<<"set title 'momentum distribution for "<<p1<<p2;
	if (trzy) gfile<<p3;
	gfile<<" events";
	if (nc) gfile<<" (NC)";
	else gfile<<" (CC)";
	gfile<<"' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel '{momentum [MeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set arrow 1 from 2000,0 to 2000,"<<ar/2.0<<" nohead"<<endl;
	gfile<<"set arrow 2 from 2000,"<<ar*9.0/20.0<<" to 2050,"<<ar*9.0/20.0<<endl;

	gfile<<"plot '"<<in<<"' using 1:2 with lines lt 1 lc rgb 'black' lw 5 title '";
	
	if (strcmp(p1.c_str(), p2.c_str()) == 0) gfile<<"first "<<p1;
	else gfile<<p1;
	
	gfile<<"', '"<<in<<"' using 1:3 with lines lt 1 lc rgb 'blue' lw 5 title '";
	
	if (strcmp(p1.c_str(), p2.c_str()) == 0) gfile<<"second "<<p2;
	else gfile<<p2;
	
	gfile<<"'";
	
	if (trzy) gfile<<", '"<<in<<"' using 1:4 with lines lt 1 lc rgb 'red' lw 5 title '"<<p3<<"'";
	
	if (!trzy) gfile<<", '"<<in<<"' using 4:5 with lines lt 1 lc rgb 'black' lw 5 notitle, '"<<in<<"' using 4:6 with lines lt 1 lc rgb 'blue' lw 5 notitle";
	else gfile<<", '"<<in<<"' using 5:6 with lines lt 1 lc rgb 'black' lw 5 notitle, '"<<in<<"' using 5:7 with lines lt 1 lc rgb 'blue' lw 5 notitle, '"<<in<<"' using 5:8 with lines lt 1 lc rgb 'red' lw 5 notitle";
	
	gfile.close();	
	run("gnuplot tmp/gnuplot.gnu");
	
}

void make_ang_plot(string filename, string p1, string p2, string p3, int k, bool nc)
{
	string in  = "NIWG/tech/results/dis/txt/" + filename;
	string out = "NIWG/tech/results/dis/plots/" + filename + ".eps";
	
	bool trzy = false;
	if (strcmp(p3.c_str(), "null")) trzy = true;
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	double ar = height[k];
	if (nc) ar = height[k+6];
	
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<out<<"'"<<endl;
	gfile<<"set xrange[-1:1]"<<endl;
	//gfile<<"set yrange[0:"<<ar<<"]"<<endl;	
	gfile<<"set title 'angle distribution for "<<p1<<p2;
	if (trzy) gfile<<p3;
	gfile<<" events";
	if (nc) gfile<<" (NC)";
	else gfile<<" (CC)";
	gfile<<"' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'cos{/Symbol Q}' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;

	gfile<<"plot '"<<in<<"' using 1:2 with lines lt 1 lc rgb 'black' lw 5 title '";
	
	if (strcmp(p1.c_str(), p2.c_str()) == 0) gfile<<"first "<<p1;
	else gfile<<p1;
	
	gfile<<"', '"<<in<<"' using 1:3 with lines lt 1 lc rgb 'blue' lw 5 title '";
	
	if (strcmp(p1.c_str(), p2.c_str()) == 0) gfile<<"second "<<p2;
	else gfile<<p2;
	
	gfile<<"'";
	
	if (trzy) gfile<<", '"<<in<<"' using 1:4 with lines lt 1 lc rgb 'red' lw 5 title '"<<p3<<"'";
		
	gfile.close();	
	run("gnuplot tmp/gnuplot.gnu");
	
}
	
void tech_calc()
{
	double nof[4][2];
	double vivi[4][3][27]; //neutrino, target, vivi
	double extra[4][3][32];
	double crossH[4];
	double crossOfg[4];
	double crossOsf[4];
	
	double momH[2][15][41];
	double angH[2][15][20];
	double momO[2][15][41];
	double angO[2][15][20];
	double mom[2][15][41];
	double ang[2][15][20];
	
	for (int n = 0; n < 4; n++)
	{
		string files[3];
		
		files[0] = "NIWG/tech/root_files/" + nuname[n] + "_hydrogen.root";
		files[1] = "NIWG/tech/root_files/" + nuname[n] + "_oxygen_fg.root";
		files[2] = "NIWG/tech/root_files/" + nuname[n] + "_oxygen_sf.root";
		
		crossH[n]   = crosssection(files[0] + ".txt");
		crossOfg[n] = crosssection(files[1] + ".txt");
		crossOsf[n] = crosssection(files[2] + ".txt");
		
		nof[n][0] = (nofH*crossH[n] + 16.0*nofO*crossOfg[n])*nufield[n];
		nof[n][1] = (nofH*crossH[n] + 16.0*nofO*crossOsf[n])*nufield[n];
		
		make_vivi(files[0], vivi[n][0], extra[n][0], false);
		make_vivi(files[1], vivi[n][1], extra[n][1], true);
		make_vivi(files[2], vivi[n][2], extra[n][2], true);
		
		if (n == 0)
		{
			string fstat[2];
			
			fstat[0] = "NIWG/tech/root_files/" + nuname[0] + "_hydrogen_stat.root";
			fstat[1] = "NIWG/tech/root_files/" + nuname[0] + "_oxygen_stat.root";

			distribution(fstat[0], momH, angH, 0);
			distribution(fstat[1], momO, angO, 1);
	
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 15; j++)
				{
					for (int k = 0; k < 41; k++)
					{
						mom[i][j][k] = crossH[0]*2.0*momH[i][j][k] + crossOfg[0]*16.0*momO[i][j][k];
						if (k < 20) ang[i][j][k] = crossH[0]*2.0*angH[i][j][k] + crossOfg[0]*16.0*angO[i][j][k];
					}
				}
			}
		}
						
	}

	ofstream noffile("NIWG/tech/results/number_of_interactions.txt");
	noffile<<"Neutrino / #interactions (fg/sf)/cross section (H/O(fg)/O(sf) / flux"<<endl<<endl;
	
	double nof_all[2] = {0, 0};
	
	for (int n = 0; n < 4; n++)
	{
		noffile<<nuname[n]<<": "<<nof[n][0]<<" "<<nof[n][1]<<" "<<crossH[n]<<" "<<crossOfg[n]<<" "<<crossOsf[n]<<" "<<nufield[n]<<endl<<endl;
		nof_all[0] += nof[n][0];
		nof_all[1] += nof[n][1];
	}
	
	noffile<<"--------------------"<<endl<<endl;
	noffile<<"all: "<<nof_all[0]<<" "<<nof_all[1]<<endl<<endl;
	
	noffile.close();
	
	double res_fg[27]; zero(res_fg, 27);
	double res_sf[27]; zero(res_sf, 27);

	double res_fg_ex[32]; zero(res_fg_ex, 32);
	double res_sf_ex[32]; zero(res_sf_ex, 32);
	
	nu_mix(vivi, extra, crossH, crossOfg, crossOsf, res_fg, res_sf, res_fg_ex, res_sf_ex);
	
	nu_norm(res_fg, 100.0, 27);
	nu_norm(res_sf, 100.0, 27);
	
	norm_ex(res_fg_ex, res_fg[6], 0, 4);
	norm_ex(res_fg_ex, res_fg[7], 5, 15);
	norm_ex(res_fg_ex, res_fg[24], 16, 20);
	norm_ex(res_fg_ex, res_fg[25], 21, 31);

	norm_ex(res_sf_ex, res_sf[6], 0, 4);
	norm_ex(res_sf_ex, res_sf[7], 5, 15);
	norm_ex(res_sf_ex, res_sf[24], 16, 20);
	norm_ex(res_sf_ex, res_sf[25], 21, 31);
				
	make_table(res_fg, res_sf, nof_all[0], nof_all[1]);
	make_table_ex(res_fg_ex, res_sf_ex, res_fg, res_sf);
	
	int help[15] = {0, 0, 1, 1, 1, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8};
	int help2[15] = {16, 16, 17, 17, 17, 21, 21, 22, 22, 22, 23, 23, 24, 24, 24};
	
	for (int i = 0; i < 15; i++)
	{
		nu_norm(mom[0][i], res_fg_ex[help[i]]*nof[0][0]/100.0, 41);
		nu_norm(ang[0][i], res_fg_ex[help[i]]*nof[0][0]/100.0, 20);
		
		nu_norm(mom[1][i], res_fg_ex[help2[i]]*nof[0][0]/100.0, 41);
		nu_norm(ang[1][i], res_fg_ex[help2[i]]*nof[0][0]/100.0, 20);
	}
	
	ofstream momfiles[12];
	ofstream angfiles[12];
	
	for (int i = 0; i < 6; i++)
	{
		string f = "NIWG/tech/results/dis/txt/mu_mom_" + channels[i];
		momfiles[i].open(f.c_str());
		
		f = "NIWG/tech/results/dis/txt/nu_mom_" + channels[i];
		momfiles[i+6].open(f.c_str());
		
		f = "NIWG/tech/results/dis/txt/mu_ang_" + channels[i];
		angfiles[i].open(f.c_str());
		
		f = "NIWG/tech/results/dis/txt/nu_ang_" + channels[i];
		angfiles[i+6].open(f.c_str());
	}
	
	int help3[6] = {0, 2, 5, 7, 10, 12};
	
	for (int i = 0; i < 6; i += 2)
	{
		for (int k = 0; k < 40; k++)
		{
			momfiles[i] << mom_bin[k] - 25 << " " << mom[0][help3[i]][k] << " " << mom[0][help3[i]+1][k];
			
			if (k == 38) momfiles[i] << " " << mom_bin[39] << " " << mom[0][help3[i]][40] << " " << mom[0][help3[i]+1][40];
			if (k == 39) momfiles[i] << " " << mom_bin[39] + 100 << " " << mom[0][help3[i]][40] << " " << mom[0][help3[i]+1][40];
			
			momfiles[i] << endl;
			
			if (k < 20) angfiles[i] << ang_bin[k] - 0.05 << " " << ang[0][help3[i]][k] << " " <<ang[0][help3[i]+1][k] << endl;
			
			momfiles[i+6] << mom_bin[k] - 25 << " " << mom[1][help3[i]][k] << " " << mom[1][help3[i]+1][k];
			
			if (k == 38) momfiles[i+6] << " " << mom_bin[39] << " " << mom[1][help3[i]][40] << " " << mom[1][help3[i]+1][40];
			if (k == 39) momfiles[i+6] << " " << mom_bin[39] + 100 << " " << mom[1][help3[i]][40] << " " << mom[1][help3[i]+1][40];
			
			momfiles[i+6] << endl;
			
			if (k < 20) angfiles[i+6] << ang_bin[k] - 0.05 << " " << ang[1][help3[i]][k] << " " <<ang[1][help3[i]+1][k] << endl;
		}
	}
	
	for (int i = 1; i < 6; i += 2)
	{
		for (int k = 0; k < 40; k++)
		{
			momfiles[i] << mom_bin[k] - 25 << " " << mom[0][help3[i]][k] << " " << mom[0][help3[i]+1][k] << " " << mom[0][help3[i]+2][k];
			
			if (k == 38) momfiles[i] << " " << mom_bin[39] << " " << mom[0][help3[i]][40] << " " << mom[0][help3[i]+1][40] << " " << mom[0][help3[i]+2][40];
			if (k == 39) momfiles[i] << " " << mom_bin[39] + 100 << " " << mom[0][help3[i]][40] << " " << mom[0][help3[i]+1][40] << " " << mom[0][help3[i]+2][40];
						
			momfiles[i] << endl;
			
			if (k < 20) angfiles[i] << ang_bin[k] - 0.05 << " " << ang[0][help3[i]][k] << " " << ang[0][help3[i]+1][k] << " " << ang[0][help3[i]+2][k] << endl;
			
			momfiles[i+6] << mom_bin[k] - 25 << " " << mom[1][help3[i]][k] << " " << mom[1][help3[i]+1][k] << " " << mom[1][help3[i]+2][k];
			
			if (k == 38) momfiles[i+6] << " " << mom_bin[39] << " " << mom[1][help3[i]][40] << " " << mom[1][help3[i]+1][40] << " " << mom[1][help3[i]+2][40];
			if (k == 39) momfiles[i+6] << " " << mom_bin[39] + 100 << " " << mom[1][help3[i]][40] << " " << mom[1][help3[i]+1][40] << " " << mom[1][help3[i]+2][40];
			
			momfiles[i+6] << endl;
			
			if (k < 20) angfiles[i+6] << ang_bin[k] - 0.05 << " " << ang[1][help3[i]][k] << " " <<ang[1][help3[i]+1][k] << " " << ang[1][help3[i]+2][k] << endl;
		}
	}
	
	for (int i = 0; i < 6; i++)
	{
		string f1 = "mu_mom_" + channels[i];
		string f2 = "nu_mom_" + channels[i];
		
		make_mom_plot(f1, part[top[i][0]], part[top[i][1]], part[top[i][2]], i, false);
		make_mom_plot(f2, part[top[i][0]], part[top[i][1]], part[top[i][2]], i, true);

		f1 = "mu_ang_" + channels[i];
		f2 = "nu_ang_" + channels[i];

		make_ang_plot(f1, part[top[i][0]], part[top[i][1]], part[top[i][2]], i, false);
		make_ang_plot(f2, part[top[i][0]], part[top[i][1]], part[top[i][2]], i, true);
	}
	
}
