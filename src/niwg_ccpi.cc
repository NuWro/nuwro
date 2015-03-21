#include "fsi.h"
#include "niwg.h"
#include "niwg_ccpi.h"
#include "calculations.h"
#include "event1.h"
#include <TFile.h>
#include <TTree.h>

int ccpi_sim_help(string *tarpar, string en)
{
	string energia = "-p 'beam_type = 0' -p 'beam_energy = " + en + string("' ");
	
	for (int l = 0; l < 4; l++)
	{		
		string command = get_bin_dir()+"nuwro -o NIWG/ccpi/root_files/e" + en + am0 + vm0 + target[l] + string(".root ") + ccpi + tarpar[l] + pam0 + pvm0 + energia;
		run(command);
			
		command = get_bin_dir()+"nuwro -o NIWG/ccpi/root_files/e" + en + am1 + vm0 + target[l] + string(".root ") + ccpi + tarpar[l] + pam1 + pvm0 + energia;
		run(command);

		command = get_bin_dir()+"nuwro -o NIWG/ccpi/root_files/e" + en + am2 + vm0 + target[l] + string(".root ") + ccpi + tarpar[l] + pam2+ pvm0 + energia;
		run(command);

		command = get_bin_dir()+"nuwro -o NIWG/ccpi/root_files/e" + en + am0 + vm1 + target[l] + string(".root ") + ccpi + tarpar[l] + pam0 + pvm1 + energia;
		run(command);
			
		command = get_bin_dir()+"nuwro -o NIWG/ccpi/root_files/e" + en + am0 + vm2 + target[l] + string(".root ") + ccpi + tarpar[l] + pam0 + pvm2 + energia;
		run(command);
	}	

	return 1;
}

int ccpi_sim()
{
	string tarpar[4];	
	tarpar[0] = "-p 'nucleus_p = 0' -p 'nucleus_n = 1' -p 'nucleus_target = 0' -p 'kaskada_on = 1' -p 'pauli_blocking = 0' ";
	tarpar[1] = "-p 'nucleus_p = 8' -p 'nucleus_n = 8' -p 'nucleus_target = 1' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' ";
	tarpar[2] = "-p 'nucleus_p = 8' -p 'nucleus_n = 8' -p 'nucleus_target = 4' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' ";
	tarpar[3] = "-p 'nucleus_p = 26' -p 'nucleus_n = 30' -p 'nucleus_target = 1' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' ";
	
	for (int energy = 200; energy < 1000; energy += 50)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccpi_sim_help(tarpar, en);
	}

	for (int energy = 1000; energy <= 10000; energy += 200)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
	
		ccpi_sim_help(tarpar, en);
	}
	
	return 1;
}

void pifact(string file, double fac[][2])
{	
	int events = 10000;
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion[3][2];
		
		pion[0][0] = e1->nof(211);
		pion[1][0] = e1->nof(-211);
		pion[2][0] = e1->nof(111);
		
		pion[0][1] = e1->fof(211);
		pion[1][1] = e1->fof(-211);
		pion[2][1] = e1->fof(111);
		
		int check = e1->nof(211) + e1->nof(-211) + e1->nof(111);
		
		if (check == 1)
		{
			for (int l = 0; l < 3; l++) if (pion[l][0] == 1) fac[l][0]++;
		}

		check = e1->fof(211) + e1->fof(-211) + e1->fof(111);
		
		if (check == 1)
		{
			for (int l = 0; l < 3; l++) if (pion[l][1] == 1) fac[l][1]++;
		}		

		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	for (int l = 0; l < 3; l++)
	{
		fac[l][0] /= events;
		fac[l][1] /= events;
	}
	
	cout<<file<<": done"<<endl<<endl;
	
}

int ccpi_calc_help (ofstream files[][3][2], string en)
{
	double cross[5][3][2]; //0,1,2 - pi+/pi-/pi0; 0 - before fsi, 1 - after fsi
	
	for (int l = 0; l < 4; l++)
	{
		string filename[5];
		
		if (l == 1) l++; //skip fg oxygen
			
		filename[0] = "NIWG/ccpi/root_files/e" + en + am0 + vm0 + target[l] + string(".root");
		filename[1] = "NIWG/ccpi/root_files/e" + en + am1 + vm0 + target[l] + string(".root");
		filename[2] = "NIWG/ccpi/root_files/e" + en + am2 + vm0 + target[l] + string(".root");
		filename[3] = "NIWG/ccpi/root_files/e" + en + am0 + vm1 + target[l] + string(".root");
		filename[4] = "NIWG/ccpi/root_files/e" + en + am0 + vm2 + target[l] + string(".root");
		
		for (int k = 0; k < 5; k++)
		{
			double factor[3][2]; for (int i = 0; i < 3; i++) zero(factor[i], 2);
			pifact(filename[k], factor);
			
			for (int l = 0; l < 3; l++)
			{
				cross[k][l][0] = crosssection(filename[k] + ".txt")*factor[l][0];
				cross[k][l][1] = crosssection(filename[k] + ".txt")*factor[l][1];
			}
		}
		
		for (int q = 0; q < 3; q++)
		{
			for (int c = 0; c < 2; c++)
			{
						
				double result = cross[0][q][c];
						
				double errup1 = cross[0][q][c] - cross[1][q][c];
				double errup2 = cross[0][q][c] - cross[3][q][c];
				double errup  = sqrt(errup1*errup1 + errup2*errup2);
							
				double errdown1 = cross[0][q][c] - cross[2][q][c];
				double errdown2 = cross[0][q][c] - cross[4][q][c];
				double errdown  = sqrt(errdown1*errdown1 + errdown2*errdown2);
						
				double err;
				
				if (errup > errdown) err = errup;
				else err = errdown;
				
				double relerr = 0;
				
				if (result != 0) relerr = err/result;
						
				files[l][q][c] << en << " " << result << " " << err << " " << relerr << " " << errup << " " << errdown << endl;
			}
		}
	}
	
	return 1;
}

int ccpi_calc()
{	
	ofstream plik[4][3][2]; //target/pi/fsi
	
	for (int i = 0; i < 4; i++)
	{
		for (int k = 0; k < 3; k++)
		{			
			string help = ccpidir + kofpi[k] + target[i] + am0 + ".txt";
			plik[i][k][0].open(help.c_str());
			
			help = ccpidir + kofpi[k] + target[i] + am0 + "_FSI.txt";
			plik[i][k][1].open(help.c_str());
		}
	}			
		
	for (int i = 0; i < 4; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			plik[i][k][0] << "#energy / cross section / error / relerror / errorup / errordown" << endl;
			plik[i][k][1] << "#energy / cross section / error / relerror / errorup / errordown" << endl;
		}
	}
				
	for (int energy = 200; energy < 1000; energy += 50)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccpi_calc_help(plik, en);
	}
	
	for (int energy = 1000; energy <= 4600; energy += 200) 
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccpi_calc_help(plik, en);
	}
	
	for (int i = 0; i < 4; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			plik[i][k][0].close();
			plik[i][k][1].close();
		}
	}
			
	return 1;
}

int ccpi_plot()
{
	string plik[4][3][2]; //target/pi/fsi
	
	for (int i = 0; i < 4; i++)
	{
		for (int k = 0; k < 3; k++)
		{			
			plik[i][k][0] = ccpidir + kofpi[k] + target[i] + am0 + ".txt";
			plik[i][k][1] = ccpidir + kofpi[k] + target[i] + am0 + "_FSI.txt";
		}
	}
	
	string pi[3] = {"{/Symbol p}^{+}", "{/Symbol p}^{-}", "{/Symbol p}^{0}"};
	string pin[3] = {"pip", "pim", "pi0"};
		
	for (int l = 0; l < 4; l++)
	{
		if (l==1) l++;
		
		for (int q = 0; q < 3; q++)
		{
			ofstream gfile("tmp/gnuplot.gnu");
			
			gfile<<setprecision (2);
			gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
			gfile<<"set output 'NIWG/ccpi/plots/ccpi_nuwro_"<<target[l]+am0<<"_"<<pin[q]<<".eps'"<<endl;
			gfile<<"set xrange[0:5000]"<<endl;
			//gfile<<"set yrange[0:2e-38]"<<endl;
			gfile<<"set title '{/Symbol n}_{/Symbol m}-"<<target[l]<<" CC "<<pi[q]<<"' font 'Arial, 20'"<<endl;
			gfile<<"set xlabel '{/Symbol n}_{/Symbol m} energy [MeV]' font 'Arial, 16'"<<endl;
			gfile<<"set ylabel 'cross section [mb]' font 'Arial, 16'"<<endl;
			gfile<<"set key spacing 1.5"<<endl;
			gfile<<"set bar 0"<<endl;

			gfile<<"plot '"<<plik[l][q][0]<<"' using 1:2:3 with yerrorbars notitle lc rgb 'black' lw 2, '"<<plik[l][q][0]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'before FSI', '"<<plik[l][q][0]<<"' using 1:2 notitle with points 7, '"<<plik[l][q][1]<<"' using 1:2:3 with yerrorbars notitle lc rgb 'blue' lw 2, '"<<plik[l][q][1]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'blue' title 'after FSI', '"<<plik[l][q][1]<<"' using 1:2 notitle with points 7"<<endl;
			
			gfile.close();
			run("gnuplot tmp/gnuplot.gnu");
		}
	}
	
	return 1;
}


