#include "fsi.h"
#include "niwg.h"
#include "niwg_ccqe.h"
#include "calculations.h"

#include <sstream>
using namespace std;

int ccqe_sim_help(string *tarpar, string en)
{
	string energia = "-p 'beam_type = 0' -p 'beam_energy = " + en + string("' ");
		
	for (int l = 0; l < 4; l++)
	{		
		string command = get_bin_dir()+"nuwro -o NIWG/ccqe/root_files/e" + en + am0 + vm0 + target[l] + string(".root ") + qelcc + tarpar[l] + pam0 + pvm0 + energia;
		run(command);
			
		command = get_bin_dir()+"nuwro -o NIWG/ccqe/root_files/e" + en + am1 + vm0 + target[l] + string(".root ") + qelcc + tarpar[l] + pam1 + pvm0 + energia;
		run(command);

		command = get_bin_dir()+"nuwro -o NIWG/ccqe/root_files/e" + en + am2 + vm0 + target[l] + string(".root ") + qelcc + tarpar[l] + pam2 + pvm0 + energia;
		run(command);

		command = get_bin_dir()+"nuwro -o NIWG/ccqe/root_files/e" + en + am0 + vm1 + target[l] + string(".root ") + qelcc + tarpar[l] + pam0 + pvm1 + energia;
		run(command);
			
		command = get_bin_dir()+"nuwro -o NIWG/ccqe/root_files/e" + en + am0 + vm2 + target[l] + string(".root ") + qelcc + tarpar[l] + pam0 + pvm1 + energia;
		run(command);
	}	

	return 1;
}

int ccqe_sim()
{
	string tarpar[4];	
	tarpar[0] = "-p 'nucleus_p = 0' -p 'nucleus_n = 1' -p 'nucleus_target = 0' -p 'kaskada_on = 0' -p 'pauli_blocking = 0' ";
	tarpar[1] = "-p 'nucleus_p = 8' -p 'nucleus_n = 8' -p 'nucleus_target = 1' -p 'kaskada_on = 0' -p 'pauli_blocking = 1' ";
	tarpar[2] = "-p 'nucleus_p = 8' -p 'nucleus_n = 8' -p 'nucleus_target = 4' -p 'kaskada_on = 0' -p 'pauli_blocking = 1' ";
	tarpar[3] = "-p 'nucleus_p = 26' -p 'nucleus_n = 30' -p 'nucleus_target = 1' -p 'kaskada_on = 0' -p 'pauli_blocking = 1' ";
	
	for (int energy = 200; energy < 1000; energy += 50)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccqe_sim_help(tarpar, en);
	}

	for (int energy = 1000; energy <= 10000; energy += 200)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
	
		ccqe_sim_help(tarpar, en);
	}
	
	return 1;
}

int ccqe_calc_help (ofstream *files, string en)
{
	double cross[5];
	
	for (int l = 0; l < 4; l++) 
	{
		string filename[5];
			
		filename[0] = "NIWG/ccqe/root_files/e" + en + am0 + vm0 + target[l] + string(".root.txt");
		filename[1] = "NIWG/ccqe/root_files/e" + en + am1 + vm0 + target[l] + string(".root.txt");
		filename[2] = "NIWG/ccqe/root_files/e" + en + am2 + vm0 + target[l] + string(".root.txt");
		filename[3] = "NIWG/ccqe/root_files/e" + en + am0 + vm1 + target[l] + string(".root.txt");
		filename[4] = "NIWG/ccqe/root_files/e" + en + am0 + vm2 + target[l] + string(".root.txt");
		
		for (int k = 0; k < 5; k++)
		{
			cross[k] = crosssection(filename[k]);
		}
				
		double result = cross[0];
				
		double errup1 = cross[0] - cross[1];
		double errup2 = cross[0] - cross[3];
		double errup  = sqrt(errup1*errup1 + errup2*errup2);
					
		double errdown1 = cross[0] - cross[2];
		double errdown2 = cross[0] - cross[4];
		double errdown  = sqrt(errdown1*errdown1 + errdown2*errdown2);
				
		double err;
		
		if (errup > errdown) err = errup;
		else err = errdown;

		double relerr = err/result;
				
		files[l] << en << " " << result << " " << err << " " << relerr << " " << errup << " " << errdown << endl;
	}
	
	return 1;
}

int ccqe_sfg_ratio()
{
	string plik[3];
	
	plik[0] = "NIWG/ccqe/results/ccqe_nuwro_oxygen_FG" + am0 + ".txt";
	plik[1] = "NIWG/ccqe/results/ccqe_nuwro_oxygen_SF" + am0 + ".txt";
	plik[2] = "NIWG/ccqe/results/ccqe_nuwro_SFG_ratio" + am0 + ".txt";
	
	ofstream res;
	res.open(plik[2].c_str());
	
	res<<"#energy vs FG/SF"<<endl;

	ifstream Input1 (plik[0].c_str());
	ifstream Input2 (plik[1].c_str());
	
	double c1, c2;
	
	if (Input1 and Input2)
	{
		string help;
		
		do
		{
			getline (Input1, help);
			getline (Input2, help);

			Input1>>help;
			Input2>>help;

			Input1>>c1;
			Input2>>c2;
			
			if (Input1) res<<help<<" "<<c1/c2<<endl;	
					
		}while (Input1);
	}
	
	res.close();
		
	return 1;
}

int ccqe_calc()
{	
	ofstream plik[4];
	
	string help = "NIWG/ccqe/results/ccqe_nuwro_neutron" + am0 + ".txt";
	plik[0].open(help.c_str());
	help = "NIWG/ccqe/results/ccqe_nuwro_oxygen_FG" + am0 + ".txt";
	plik[1].open(help.c_str());
	help = "NIWG/ccqe/results/ccqe_nuwro_oxygen_SF" + am0 + ".txt";
	plik[2].open(help.c_str());
	help = "NIWG/ccqe/results/ccqe_nuwro_iron" + am0 + ".txt";
	plik[3].open(help.c_str());
	
	for (int i = 0; i < 3; i++) plik[i] << "#energy / cross section / error / relerror / errorup / errordown" << endl;
				
	for (int energy = 200; energy < 1000; energy += 50)
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccqe_calc_help(plik, en);
	}
	
	for (int energy = 1000; energy <= 10000; energy += 200) 
	{
		stringstream temp;
		string en;

		temp << energy;
		temp >> en;
		
		ccqe_calc_help(plik, en);
	}
	
	plik[0].close();
	plik[1].close();
	plik[2].close();
	plik[3].close();
	
	ccqe_sfg_ratio();
				
	return 1;
}

int ccqe_sfg_plot()
{
	string plik = "NIWG/ccqe/results/ccqe_nuwro_SFG_ratio" + am0 + ".txt";
	
	ofstream gfile("tmp/gnuplot.gnu");
		
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;		
	gfile<<"set output 'NIWG/ccqe/plots/ccqe_nuwro_sfg_ratio_"<<am0<<".eps'"<<endl;
	gfile<<"set xrange[0:10000]"<<endl;
	//gfile<<"set yrange[0:2e-38]"<<endl;
	gfile<<"set title '{/Symbol n}_{/Symbol m}-oxygen CCQE (axial mass = "<<am<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel '{/Symbol n}_{/Symbol m} energy [MeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'ratio' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;

	gfile<<"plot '"<<plik<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'Fermi Gas/Spectral function"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");	
	
	return 1;
}

int ccqe_plot()
{
	string plik[4];
	
	plik[0] = "NIWG/ccqe/results/ccqe_nuwro_neutron" + am0 + ".txt";
	plik[1] = "NIWG/ccqe/results/ccqe_nuwro_oxygen_FG" + am0 + ".txt";
	plik[2] = "NIWG/ccqe/results/ccqe_nuwro_oxygen_SF" + am0 + ".txt";
	plik[3] = "NIWG/ccqe/results/ccqe_nuwro_iron" + am0 + ".txt";
		
	for (int l = 0; l < 4; l++)
	{
		ofstream gfile("tmp/gnuplot.gnu");
		
		gfile<<setprecision (2);
		gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;		
		gfile<<"set output 'NIWG/ccqe/plots/ccqe_nuwro_"<<target[l]+am0<<".eps'"<<endl;
		gfile<<"set xrange[0:10000]"<<endl;
		//gfile<<"set yrange[0:2e-38]"<<endl;
		if (l == 1) gfile<<"set title '{/Symbol n}_{/Symbol m}-oxygen (FG) CCQE (axial mass = "<<am<<")' font 'Arial, 20'"<<endl;
		else if (l == 2) gfile<<"set title '{/Symbol n}_{/Symbol m}-oxygen (SF) CCQE (axial mass = "<<am<<")' font 'Arial, 20'"<<endl;
		else gfile<<"set title '{/Symbol n}_{/Symbol m}-"<<target[l]<<" CCQE (axial mass = "<<am<<")' font 'Arial, 20'"<<endl;
		gfile<<"set xlabel '{/Symbol n}_{/Symbol m} energy [MeV]' font 'Arial, 16'"<<endl;
		gfile<<"set ylabel 'cross section [mb]' font 'Arial, 16'"<<endl;
		gfile<<"set key spacing 1.5"<<endl;
		gfile<<"set bar 0"<<endl;

		gfile<<"plot '"<<plik[l]<<"' using 1:2:3 with yerrorbars notitle lc rgb 'black' lw 2, '"<<plik[l]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' notitle, '"<<plik[l]<<"' using 1:2 notitle with points 7"<<endl;
		
		gfile.close();
		run("gnuplot tmp/gnuplot.gnu");	
	}
	
	ccqe_sfg_plot();
	
	return 1;
}


