#include <fstream>
#include <iostream>

#include "fsi.h"
#include "data.h"
#include "plots.h"
#include "calculations.h"

using namespace std;

void make_data(const double *x, const double *y, const double *xerr, const double *yerr, const int bins, bool rest)
{
	ofstream datafile("tmp/data.txt");
	
	for (int i = 0; i < bins; i++)
	{
		datafile<<x[i]<<" "<<y[i]<<" "<<xerr[i]<<" "<<yerr[i]<<endl;
	}
	
	if (rest) datafile<<x[bins]<<" "<<y[bins]<<" "<<xerr[bins]<<" "<<yerr[bins];

	datafile.close();
}	

void plotlog(string text)
{
	logfile<<day<<" "<<fullm<<" "<<year<<" "<<hour<<":"<<minute<<endl<<text<<endl<<endl;
}

void findratio(string file, double &ratio, double &ratiofsi)
{
	ifstream input (file.c_str());
	string help;
	
	if (input)
	{		
		for (int k = 0; k < 4; k++)	input>>help;
		input>>ratio;
		for (int k = 0; k < 4; k++)	input>>help;
		input>>ratiofsi;
	}
}

void findnof(string file, double *back)
{
	ifstream input (file.c_str());
	string help;
	
	if (input)
	{		
		for (int k = 0; k < 7; k++)	input>>help;
		input>>back[0];
		input>>back[1];
		input>>back[2];
		input>>back[3];
	}
}

int plotK2K(int fz, int xs)
{
	get_date();
	string filename;
	const int bins = 8;
	
	filename = string("results/K2K/k2k_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename = find_last(filename);
	
	if (noFile(filename))
	{
		cout<<"K2K ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"K2K ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	double ratio;
	double ratiofsi;
	
	findratio(filename, ratio, ratiofsi);
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/K2K/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("k2k_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:900]"<<endl;
	gfile<<"set yrange[0:1300]"<<endl;
	gfile<<"set title '{/=24 K2K} NC{/Symbol p}^{0} production on H_{2}0 ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [MeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set arrow 1 from 800,0 to 800,500 nohead"<<endl;
	gfile<<"set arrow 2 from 800,400 to 850,400"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 350, 1000 font 'Arial Bold, 16'"<<endl;
	gfile<<"set label 2 'NC{/Symbol p}^{0}/CC ratio: ' at 550, 625"<<endl;
	gfile<<"set label 3 'Data:' at 550, 550"<<endl;
	gfile<<"set label 4 'Before FSI:' at 550, 500"<<endl;
	gfile<<"set label 5 'After FSI:' at 550, 450"<<endl;
	gfile<<"set label 6 '"<<K2Kratio<<"{/Symbol \261}"<<K2Kratioerr<<"{/Symbol \261}"<<K2Kratioerr2<<"' at 650, 550"<<endl;
	gfile<<"set label 7 '"<<ratio<<"' at 650, 500"<<endl;
	gfile<<"set label 8 '"<<ratiofsi<<"' at 650, 450"<<endl;

	make_data(K2Kx, K2Ky, K2Kxerr, K2Kyerr, bins, 1);
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'K2K data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:2 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename<<"' using 1:3 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename<<"' using 10:11 with lines lt 2 lw 5 lc rgb 'black' notitle, '";
	gfile<<filename<<"' using 10:12 with lines lt 1 lw 5 lc rgb 'black' notitle, '";
	gfile<<filename<<"' using 1:4 smooth csplines lt 2 lw 5 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename<<"' using 1:5 smooth csplines lt 1 lw 5 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename<<"' using 1:6 smooth csplines lt 2 lw 5 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename<<"' using 1:7 smooth csplines lt 1 lw 5 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename<<"' using 1:8 smooth csplines lt 2 lw 5 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename<<"' using 1:9 smooth csplines lt 1 lw 5 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'";
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

int plotMBangle(int fz, int xs)
{
	get_date();
	string filename[4]; //0 - neutrino, 1 - antineutrino
	const int bins  = 18;
	const int binsa = 10;
	
	filename[0] = string("results/MB/mb_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[1] = string("results/MB/mb_anti_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[2] = string("results/MB/mb_normalized_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[3] = string("results/MB/mb_normalized_anti_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");

	for (int i = 0; i < 4; i++) filename[i] = find_last(filename[i]);
	
	bool allfiles = true;
	for (int i = 0; i < 4; i++) if (noFile(filename[i])) allfiles = false;
	
	if (!allfiles)
	{
		cout<<"MB angle ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"MB angle ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/MB/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("mb_angle_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[-1:1]"<<endl;
	gfile<<"set yrange[0:1e-39]"<<endl;
	gfile<<"set xlabel 'cos{/Symbol O}_{/Symbol p_{0}}'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dcos{/Symbol O}_{/Symbol p_{0}} [10^{-39} cm^{2}/nucleon]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1,1.1"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl; 
	gfile<<"set key left"<<endl;
	gfile<<"set ytics ('0' 0, '0.2' 0.2e-39, '0.4' 0.4e-39, '0.6' 0.6e-39, '0.8' 0.8e-39)"<<endl;

	make_data(MBAnglex, MBAngley, MBAnglexerr, MBAngleyerr, bins, 0);
	run("cp tmp/data.txt tmp/data2.txt");
	make_data(MBantiAnglex, MBantiAngley, MBantiAnglexerr, MBantiAngleyerr, binsa, 0);
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set title 'neutrino mode'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[0]<<"' using 10:11 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[0]<<"' using 10:12 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)'"<<endl;

	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0"<<endl;
	gfile<<"set title 'neutrino mode - normalized'"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[2]<<"' using 10:11 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[2]<<"' using 10:12 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)'"<<endl;

	gfile<<"set label 1 '{/=16 MiniBooNE} NC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<") - "<<fzname[fz]<<" ' at -3.5, 0.65e-39 font 'Arial, 20'"<<endl;
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0.5"<<endl;
	gfile<<"set title 'anti-neutrino mode'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dcos{/Symbol O}_{/Symbol p_{0}} [10^{-39} cm^{2}/nucleon]'"<<endl;	
	gfile<<"set yrange[0:0.5e-39]"<<endl;
	gfile<<"set ytics ('0' 0, '0.1' 0.1e-39, '0.2' 0.2e-39, '0.3' 0.3e-39, '0.4' 0.4e-39)"<<endl;	

	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[1]<<"' using 10:11 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[1]<<"' using 10:12 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)'"<<endl;
 
    gfile<<"unset label 1"<<endl;

	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0"<<endl;
	gfile<<"set title 'anti-neutrino mode - normalized'"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[3]<<"' using 10:11 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[3]<<"' using 10:12 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)'"<<endl;
				
	gfile<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	for (int i = 0; i < 2; i++) help += filename[i] + (" ");
	plotlog(help);
	
	return 1;
}

int plotMB(int fz, int xs)
{
	get_date();
	string filename[4]; //0 - neutrino, 1 - antineutrino, 2,3 - normalized
	const int bins  = 11;
	const int binsa = 10;
	
	filename[0] = string("results/MB/mb_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[1] = string("results/MB/mb_anti_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[2] = string("results/MB/mb_normalized_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename[3] = string("results/MB/mb_normalized_anti_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");

	for (int i = 0; i < 4; i++) filename[i] = find_last(filename[i]);
	
	bool allfiles = true;
	for (int i = 0; i < 4; i++) if (noFile(filename[i])) allfiles = false;
	
	if (!allfiles)
	{
		cout<<"MB ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"MB ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/MB/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("mb_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:1000]"<<endl;
	gfile<<"set yrange[0:2e-42]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [MeV]'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-42} cm^{2}/MeV/nucleon]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1.5,1.2"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl; 
	gfile<<"set nokey"<<endl;
	gfile<<"set ytics ('0' 0, '0.5' 0.5e-42, '1.0' 1e-42, '1.5' 1.5e-42, '2.0' 2e-42)"<<endl;

	make_data(MBMomx, MBMomy, MBMomxerr, MBMomyerr, bins, 0);
	run("cp tmp/data.txt tmp/data2.txt");
	make_data(MBantiMomx, MBantiMomy, MBantiMomxerr, MBantiMomyerr, binsa, 0);
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set title 'neutrino'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[0]<<"' using 1:2 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[0]<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename[0]<<"' using 1:4 smooth csplines lt 2 lw 3 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename[0]<<"' using 1:5 smooth csplines lt 1 lw 3 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[0]<<"' using 1:6 smooth csplines lt 2 lw 3 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename[0]<<"' using 1:7 smooth csplines lt 1 lw 3 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[0]<<"' using 1:8 smooth csplines lt 2 lw 3 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename[0]<<"' using 1:9 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'"<<endl;

	gfile<<"set label 1 '{/=24 MiniBooNE} NC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<")' at -800, 3e-42 font 'Arial, 20'"<<endl;
	
	double cs;
	double csfsi;
	double csa;
	double csfsia;
	
	findratio(filename[0], cs, csfsi);
	findratio(filename[1], csa, csfsia);
	
	gfile<<"set label 2 'Neutrino cross section [10^{-40} cm^{2}/nucleon]: ' at 1200, 2e-42"<<endl;
	gfile<<"set label 3 'Data: ' at 1200, 1.75e-42"<<endl;
	gfile<<"set label 4 'Before FSI: ' at 1200, 1.5e-42"<<endl;
	gfile<<"set label 5 'After FSI:  ' at 1200, 1.25e-42"<<endl;
	gfile<<"set label 6 'Antineutrino cross section [10^{-40} cm^{2}/nucleon]: ' at 1200, 0.75e-42"<<endl;
	gfile<<"set label 7 'Data: ' at 1200, 0.5e-42"<<endl;
	gfile<<"set label 8 'Before FSI: ' at 1200, 0.25e-42"<<endl;
	gfile<<"set label 9 'After FSI:  ' at 1200, 0e-42"<<endl;
	gfile<<"set label 10 '"<<MBcross<<"{/Symbol \261}"<<MBcrosserr<<"{/Symbol \261}"<<MBcrosserr2<<"' at 1600, 1.75e-42"<<endl;
	gfile<<"set label 11 '"<<cs*1e40<<"' at 1600, 1.5e-42"<<endl;
	gfile<<"set label 12 '"<<csfsi*1e40<<"' at 1600, 1.25e-42"<<endl;
	gfile<<"set label 13 '"<<MBcrossanti<<"{/Symbol \261}"<<MBcrossantierr<<"{/Symbol \261}"<<MBcrossantierr2<<"' at 1600, 0.5e-42"<<endl;
	gfile<<"set label 14 '"<<csa*1e40<<"' at 1600, 0.25e-42"<<endl;
	gfile<<"set label 15 '"<<csfsia*1e40<<"' at 1600, 0.0e-42"<<endl;
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0.5"<<endl;
	gfile<<"set title 'neutrino - normalized'"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[2]<<"' using 1:2 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[2]<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename[2]<<"' using 1:4 smooth csplines lt 2 lw 3 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename[2]<<"' using 1:5 smooth csplines lt 1 lw 3 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[2]<<"' using 1:6 smooth csplines lt 2 lw 3 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename[2]<<"' using 1:7 smooth csplines lt 1 lw 3 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[2]<<"' using 1:8 smooth csplines lt 2 lw 3 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename[2]<<"' using 1:9 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'"<<endl;
	
	for (int i = 1; i <= 16; i++) gfile<<"unset label "<<i<<endl;
			
	gfile<<"set yrange[0:0.8e-42]"<<endl;
	gfile<<"set ytics ('0' 0, '2.0' 2e-43, '4.0' 4e-43, '6.0' 6e-43, '8.0' 8e-43)"<<endl;
	gfile<<"set title 'antineutrino'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-43} cm^{2}/MeV/nucleon]'"<<endl;

	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[1]<<"' using 1:2 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[1]<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename[1]<<"' using 1:4 smooth csplines lt 2 lw 3 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename[1]<<"' using 1:5 smooth csplines lt 1 lw 3 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[1]<<"' using 1:6 smooth csplines lt 2 lw 3 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename[1]<<"' using 1:7 smooth csplines lt 1 lw 3 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[1]<<"' using 1:8 smooth csplines lt 2 lw 3 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename[1]<<"' using 1:9 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'"<<endl;

	gfile<<"set size 1, 0.5"<<endl;
	gfile<<"set origin 0.5, 0"<<endl;
	gfile<<"set key out vert"<<endl;
	gfile<<"set label '"<<fzname[fz]<<"' at 1400,1e-43 font 'Arial Bold, 16'"<<endl;
	gfile<<"set title 'antineutrino - normalized"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename[3]<<"' using 1:2 smooth csplines lt 2 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[3]<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename[3]<<"' using 1:4 smooth csplines lt 2 lw 3 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename[3]<<"' using 1:5 smooth csplines lt 1 lw 3 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[3]<<"' using 1:6 smooth csplines lt 2 lw 3 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename[3]<<"' using 1:7 smooth csplines lt 1 lw 3 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename[3]<<"' using 1:8 smooth csplines lt 2 lw 3 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename[3]<<"' using 1:9 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'"<<endl;
	
	gfile<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	for (int i = 0; i < 4; i++) help += filename[i] + (" ");
	plotlog(help);
	
	plotMBangle(fz, xs);
	
	return 1;
}

int plotMBCC(int fz, int xs)
{
	get_date();
	
	string filename = "results/MBCC/mb_cc_" + fzwork[fz] + sep + xsec[xs] + string("*.txt");
	filename = find_last(filename);

	string filename2 = "results/MBCC/mb_cc_normalized_" + fzwork[fz] + sep + xsec[xs] + string("*.txt");
	filename2 = find_last(filename2);

	const int bins  = 11;
		
	if (noFile(filename) or noFile(filename2))
	{
		cout<<"MB CC ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"MB CC ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/MBCC/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("mb_cc_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:1.4]"<<endl;
	gfile<<"set yrange[0:4.5e-38]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [GeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-38}cm^{2}/GeV/CH_{2}]' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set ytics ('0' 0, '0.5' 0.5e-38, '1.0' 1e-38, '1.5' 1.5e-38, '2.0' 2e-38, '2.5' 2.5e-38, '3.0' 3e-38, '3.5' 3.5e-38)"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 0.8, 2.0e-38 font 'Arial Bold, 16'"<<endl;
	gfile<<"set size 2,1.2"<<endl;
	gfile<<"set multiplot"<<endl; 
	gfile<<"set label 2 '{/=24 MB} CC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<")' at 1, 5e-38 font 'Arial, 20'"<<endl;
	
	make_data(MBCCmom, MBCCxsecnuance, MBCCmomerr, MBCCxsecerr, bins, 0);
	run(string("cp tmp/data.txt tmp/nuance.txt"));
	
	double nu_norm[bins];
	for (int i = 0; i < bins; i++) nu_norm[i] = MBCCxsecnuance[i];
	double fact = factor(nu_norm, MBCCxsec, bins);
	for (int i = 0; i < bins; i++) nu_norm[i] *= fact;
	
	make_data(MBCCmom, nu_norm, MBCCmomerr, MBCCxsecerr, bins, 0);
	run(string("cp tmp/data.txt tmp/nu_norm.txt"));
	
	make_data(MBCCmom, MBCCxsec, MBCCmomerr, MBCCxsecerr, bins, 0);
	
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set size 1,1"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:2 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'red' title 'NuWro (after FSI)', ";
	gfile<<"'tmp/nuance.txt' using 1:2 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'Nuance'"<<endl;
	
	gfile<<"set ylabel 'arbitrary units' font 'Arial, 16'"<<endl;	
	gfile<<"unset label 2"<<endl;
	gfile<<"set origin 1,0"<<endl;
	gfile<<"set size 1,1"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename2<<"' using 1:2 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename2<<"' using 1:3 smooth csplines lt 1 lw 3 lc rgb 'red' title 'NuWro (after FSI)', ";
	gfile<<"'tmp/nu_norm.txt' using 1:2 smooth csplines lt 1 lw 3 lc rgb 'blue' title 'Nuance'"<<endl;
	
	gfile<<"unset multiplot"<<endl;
	
	gfile.close();

	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

int plotMBCCtotal(int fz, int xs)
{
	get_date();
	string filename;
	const int bins  = 14;
	
	filename = string("results/MBCC/mb_cc_total_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");

	filename = find_last(filename);
		
	if (noFile(filename))
	{
		cout<<"MB CC ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"MB CC ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/MBCC/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("mb_cc_total_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0.5:2]"<<endl;
	gfile<<"set yrange[0:30e-39]"<<endl;
	gfile<<"set xlabel '{/Symbol n}_{/Symbol m} Energy [GeV]'"<<endl;
	gfile<<"set ylabel '{/Symbol s} [10^{-39} cm^{2}/CH_{2}]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set title 'MB CC{/Symbol p}^{0} on CH_{2}'"<<endl;
	gfile<<"set key left top"<<endl;
	
	
	double energy[bins+1];
	for (int l = 0; l < bins; l++) energy[l+1] = MBCCnuen[l] + MBCCnuenerr[l];
	energy[0] = MBCCnuen[0] - MBCCnuenerr[0];
	
	make_data(energy, MBCCtotalnuance, MBCCnuenerr, MBCCtotalerr, bins+1, 0);
	run(string("cp tmp/data.txt tmp/data2.txt"));
	make_data(MBCCnuen, MBCCtotal, MBCCnuenerr, MBCCtotalerr, bins, 0);
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:2 smooth frequency with fstep lt 1 lw 3 lc rgb 'blue' title 'NuWro', ";
	gfile<<"'tmp/data2.txt' using 1:2 smooth frequency with fstep lt 1 lw 3 lc rgb 'red' title 'Nuance'"<<endl;

	gfile.close();

	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}


int plotSB(int fz, int xs)
{
	get_date();
	string filename;
	const int bins = 9;
	
	filename = string("results/SB/sb_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename = find_last(filename);
	
	if (noFile(filename))
	{
		cout<<"SB ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"SB ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	double ratio;
	double ratiofsi;
	
	findratio(filename, ratio, ratiofsi);
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/SB/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("sb_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:800]"<<endl;
	gfile<<"set yrange[0:0.3]"<<endl;
	gfile<<"set title '{/=24 SciBooNE} NC{/Symbol p}^{0} production on C_{8}H_{8} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [MeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set arrow 1 from 720,0 to 720,0.15 nohead"<<endl;
	gfile<<"set arrow 2 from 720,0.125 to 760,0.125"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 350, 0.25 font 'Arial Bold, 16'"<<endl;
	gfile<<"set label 2 'NC{/Symbol p}^{0}/CC ratio: ' at 450, 0.175"<<endl;
	gfile<<"set label 3 'Data:' at 450, 0.16"<<endl;
	gfile<<"set label 4 'Before FSI:' at 450, 0.15"<<endl;
	gfile<<"set label 5 'After FSI:' at 450, 0.14"<<endl;
	gfile<<"set label 6 '"<<SBratio<<"0{/Symbol \261}"<<SBratioerr<<"{/Symbol \261}"<<SBratioerrdown<<"' at 550, 0.16"<<endl;
	gfile<<"set label 7 '"<<ratio<<"' at 550, 0.15"<<endl;
	gfile<<"set label 8 '"<<ratiofsi<<"' at 550, 0.14"<<endl;

	make_data(SBx, SBy, SBxerr, SByerr, bins, 1);
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'SB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:2 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename<<"' using 1:3 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (after FSI)', '";
	gfile<<filename<<"' using 10:11 with lines lt 2 lw 5 lc rgb 'black' notitle, '";
	gfile<<filename<<"' using 10:12 with lines lt 1 lw 5 lc rgb 'black' notitle, '";
	gfile<<filename<<"' using 1:4 smooth csplines lt 2 lw 5 lc rgb 'red' title '{/Symbol p}^{0} -> no {/Symbol p}', '";
	gfile<<filename<<"' using 1:5 smooth csplines lt 1 lw 5 lc rgb 'red' title 'no {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename<<"' using 1:6 smooth csplines lt 2 lw 5 lc rgb 'green' title '{/Symbol p}^{0} -> charged {/Symbol p}', '";
	gfile<<filename<<"' using 1:7 smooth csplines lt 1 lw 5 lc rgb 'green' title 'charged {/Symbol p} -> {/Symbol p}^{0}', '";
	gfile<<filename<<"' using 1:8 smooth csplines lt 2 lw 5 lc rgb 'blue' title '{/Symbol p}^{0} -> more {/Symbol p}', '";
	gfile<<filename<<"' using 1:9 smooth csplines lt 1 lw 5 lc rgb 'blue' title 'more {/Symbol p} -> {/Symbol p}^{0}'";
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

int plotPNS(int fz, int xs)
{
	get_date();
	string filenameC, filenameI;
	const int bins = 10;
	
	filenameC = string("results/PNS/pns_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameC = find_last(filenameC);

	filenameI = string("results/PNS/pns_iron_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameI = find_last(filenameI);
	
	if (noFile(filenameC) or noFile(filenameI))
	{
		cout<<"Pion-Nucleus ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"Pion-Nucleus ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/PNS/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("pns_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:550]"<<endl;
	gfile<<"set yrange[0:500]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{+} kinetic energy [MeV/c^2]'"<<endl;
	gfile<<"set ylabel 'cross section [mb]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1.15, 0.5"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl;
	gfile<<"set key out vert bot right"<<endl;
	//gfile<<"set label 1 '"<<fzname[fz]<<"' at 600, 600 font 'Arial Bold, 16'"<<endl;
	//gfile<<"set label 2 '{/Symbol p}^{+}N scattering ("<<day<<" "<<fullm<<" "<<year<<")' at 500, 700 font 'Arial, 20'"<<endl;	
	
	double err6[6] = {0, 0, 0, 0, 0, 0};
	make_data(CarbonEnergyAshery, CarbonReacAshery, err6, err6, 6, 0);
	run("cp tmp/data.txt tmp/car.txt");
	make_data(CarbonEnergyAshery, CarbonAbsAshery, err6, CarbonAbsAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/caa.txt");
	make_data(CarbonEnergyAshery, CarbonSCXAshery, err6, CarbonSCXAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/cas.txt");
	make_data(CarbonEnergyAshery, CarbonInelAshery, err6, CarbonInelAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/cai.txt");
	double err1[1] = {0};
	make_data(CarbonEnergyNavon, CarbonAbsNavon, err1, CarbonAbsNavonErr, 1, 0);
	run("cp tmp/data.txt tmp/cna.txt");
	double err4[4] = {0, 0, 0, 0};
	make_data(CarbonEnergyJones, CarbonAbsJones, err4, CarbonAbsJonesErr, 4, 0);
	run("cp tmp/data.txt tmp/cja.txt");
	make_data(CarbonEnergyJones, CarbonSCXJones, err4, err4, 4, 0);
	run("cp tmp/data.txt tmp/cjs.txt");
	double err2[2] = {0, 0};
	make_data(CarbonEnergyJones2, CarbonReacJones, err2, err2, 2, 0);
	run("cp tmp/data.txt tmp/cjr.txt");
	make_data(CarbonEnergyJones2, CarbonInelJones, err2, err2, 2, 0);
	run("cp tmp/data.txt tmp/cji.txt");
	
	gfile<<"set size 0.65, 0.5"<<endl;
	gfile<<"set origin 0, 0"<<endl;
	gfile<<"set title 'Carbon'"<<endl;
	gfile<<"plot 'tmp/car.txt' using 1:2:3:4 notitle with xyerrorbars pt 7 lc rgb 'black', ";
	gfile<<"'tmp/caa.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/caa.txt' using 1:2 title 'Ashery' lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/cas.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', ";
	gfile<<"'tmp/cai.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'green', ";
	gfile<<"'tmp/cna.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'red', ";
	gfile<<"'tmp/cna.txt' using 1:2 title 'Navon' lt 1 pt 6 lc rgb 'red', ";
	gfile<<"'tmp/cjr.txt' using 1:2:3:4 notitle with xyerrorbars pt 5 lc rgb 'black', ";
	gfile<<"'tmp/cja.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 5 lc rgb 'red', ";
	gfile<<"'tmp/cja.txt' using 1:2 title 'Jones' lt 1 pt 5 lc rgb 'red', ";
	gfile<<"'tmp/cjs.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 5 lc rgb 'blue', ";
	gfile<<"'tmp/cji.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 5 lc rgb 'green', '";
	gfile<<filenameC<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 3 title 'Reaction', '";
	gfile<<filenameC<<"' using 1:3 smooth csplines lt 1 lc rgb 'red' lw 3 title 'Absorption', '";
	gfile<<filenameC<<"' using 1:4 smooth csplines lt 1 lc rgb 'blue' lw 3 title 'CEX', '";
	gfile<<filenameC<<"' using 1:5 smooth csplines lt 1 lc rgb 'green' lw 3 title 'Inelastic'"<<endl;

	make_data(IronEnergyAshery, IronReacAshery, err6, err6, 6, 0);
	run("cp tmp/data.txt tmp/iar.txt");
	make_data(IronEnergyAshery, IronAbsAshery, err6, IronAbsAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/iaa.txt");
	make_data(IronEnergyAshery, IronSCXAshery, err6, IronSCXAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/ias.txt");
	make_data(IronEnergyAshery, IronInelAshery, err6, IronInelAsheryErr, 6, 0);
	run("cp tmp/data.txt tmp/iai.txt");
	make_data(IronEnergyNavon, IronAbsNavon, err1, IronAbsNavonErr, 1, 0);
	run("cp tmp/data.txt tmp/ina.txt");
	
	gfile<<"unset label 1"<<endl;
	gfile<<"unset label 2"<<endl;
	gfile<<"set nokey"<<endl;
	gfile<<"set size 0.5,0.5"<<endl;
	gfile<<"set origin 0.65,0"<<endl;
	gfile<<"set xrange[0:400]"<<endl;	
	gfile<<"set title 'Iron'"<<endl;
	gfile<<"set yrange[0:1400]"<<endl;
	gfile<<"plot 'tmp/iar.txt' using 1:2:3:4 notitle with xyerrorbars pt 7 lc rgb 'black', ";
	gfile<<"'tmp/iaa.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/iaa.txt' using 1:2 title 'Ashery' lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/ias.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', ";
	gfile<<"'tmp/iai.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'green', ";
	gfile<<"'tmp/ina.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'red', ";
	gfile<<"'tmp/ina.txt' using 1:2 title 'Navon' lt 1 pt 6 lc rgb 'red', '";
	gfile<<filenameI<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 3 title 'Reaction', '";
	gfile<<filenameI<<"' using 1:3 smooth csplines lt 1 lc rgb 'red' lw 3 title 'Absorption', '";
	gfile<<filenameI<<"' using 1:4 smooth csplines lt 1 lc rgb 'blue' lw 3 title 'CEX', '";
	gfile<<filenameI<<"' using 1:5 smooth csplines lt 1 lc rgb 'green' lw 3 title 'Inelastic'"<<endl;
	
	gfile<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filenameC + string(" and ") + filenameI;
	plotlog(help);
	
	return 1;
}

int plotPrThe(int fz, int xs)
{
	get_date();
	string filenameC, filenameI;
	const int bins = 20;
		
	filenameC = string("results/PrTrans/prtrans_he_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameC = find_last(filenameC);
	filenameI = string("results/PrTrans/prtrans_he_iron_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameI = find_last(filenameI);
	
	if (noFile(filenameC) or noFile(filenameI))
	{
		cout<<"ProtonTransparency (high-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"ProtonTransparency (high-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/PrTrans/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("prtrans_he_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:7]"<<endl;
	gfile<<"set yrange[0.2:1]"<<endl;
	gfile<<"set title '{/=24 Proton Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'Q^{2} [GeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 2, 0.9 font 'Arial Bold, 16'"<<endl;

	double err5[5] = {0, 0, 0, 0, 0};
	double err4[4] = {0, 0, 0, 0};
	
	double CarbonBatesTerr[4]; for (int i = 0; i < 4; i++) CarbonBatesTerr[i] = 0.032*CarbonBatesT[i];
	double IronBatesTerr[4]; for (int i = 0; i < 4; i++) IronBatesTerr[i] = 0.032*IronBatesT[i];
	
	make_data(CarbonNE18Q2, CarbonNE18T, err5, CarbonNE18Terr, 5, 0);
	run("cp tmp/data.txt tmp/cn.txt");
	make_data(CarbonBatesQ2, CarbonBatesT, err4, CarbonBatesTerr, 4, 0);
	run("cp tmp/data.txt tmp/cb.txt");
	make_data(IronNE18Q2, IronNE18T, err5, IronNE18Terr, 5, 0);
	run("cp tmp/data.txt tmp/in.txt");
	make_data(IronBatesQ2, IronBatesT, err4, IronBatesTerr, 4, 0);
	run("cp tmp/data.txt tmp/ib.txt");
	
	gfile<<"plot 'tmp/cn.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/cn.txt' using 1:2 title 'NE18 data' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/cb.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 8 lc rgb 'black', '";
	gfile<<"tmp/cb.txt' using 1:2 title 'Bates data' lt 1 pt 8 lc rgb 'black', '";
	gfile<<"tmp/in.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', '";
	gfile<<"tmp/ib.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 8 lc rgb 'red', '";
	gfile<<filenameC<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'Carbon', '";
	gfile<<filenameI<<"' using 1:2 smooth csplines lt 1 lc rgb 'red' lw 5 title 'Iron'"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filenameC + string(" and ") + filenameI;
	plotlog(help);
	
	return 1;
}

int plotPrTle(int fz, int xs)
{
	get_date();
	string filenameL, filenameC, filenameA;
	const int bins = 10;
		
	filenameL = string("results/PrTrans/prtrans_le_lithium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameL = find_last(filenameL);
	filenameC = string("results/PrTrans/prtrans_le_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameC = find_last(filenameC);
	filenameA = string("results/PrTrans/prtrans_le_aluminium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameA = find_last(filenameA);
	
	if (noFile(filenameL) or noFile(filenameC) or noFile(filenameA))
	{
		cout<<"ProtonTransparency (low-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"ProtonTransparency (low-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/PrTrans/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("prtrans_le_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[125:300]"<<endl;
	gfile<<"set yrange[0.2:1]"<<endl;
	gfile<<"set title '{/=24 Proton Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'Kinetic energy [MeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 200, 0.35 font 'Arial Bold, 16'"<<endl;

	double err7[7] = {0, 0, 0, 0, 0, 0, 0};
	
	make_data(PrTrLeEn, PrTrLeLi, err7, PrTrLeLierr, 7, 0);
	run("cp tmp/data.txt tmp/li.txt");
	make_data(PrTrLeEn, PrTrLeC, err7, PrTrLeCerr, 7, 0);
	run("cp tmp/data.txt tmp/c.txt");
	make_data(PrTrLeEn, PrTrLeAl, err7, PrTrLeAlerr, 7, 0);
	run("cp tmp/data.txt tmp/al.txt");
	
	gfile<<"plot 'tmp/li.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/c.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', '";
	gfile<<"tmp/al.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', '";
	gfile<<filenameL<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'Lithium', '";
	gfile<<filenameC<<"' using 1:2 smooth csplines lt 1 lc rgb 'red' lw 5 title 'Carbon', '";
	gfile<<filenameA<<"' using 1:2 smooth csplines lt 1 lc rgb 'blue' lw 5 title 'Aluminium'"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filenameL + string(" and ") + filenameC + string(" and" ) + filenameA;
	plotlog(help);
	
	return 1;
}

int plotPiTle(int fz, int xs)
{
	get_date();
	string filenameLp, filenameCp, filenameAp, filenameL0, filenameC0, filenameA0;
	const int bins = 10;
	
	filenameLp = string("results/PiTrans/pitrans_le_pip_lithium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameLp = find_last(filenameLp);
	filenameCp = string("results/PiTrans/pitrans_le_pip_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameCp = find_last(filenameCp);
	filenameAp = string("results/PiTrans/pitrans_le_pip_aluminium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameAp = find_last(filenameAp);
	filenameL0 = string("results/PiTrans/pitrans_le_pi0_lithium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameL0 = find_last(filenameL0);
	filenameC0 = string("results/PiTrans/pitrans_le_pi0_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameC0 = find_last(filenameC0);
	filenameA0 = string("results/PiTrans/pitrans_le_pi0_aluminium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameA0 = find_last(filenameA0);
		
	if (noFile(filenameLp) or noFile(filenameCp) or noFile(filenameAp) or noFile(filenameL0) or noFile(filenameC0) or noFile(filenameA0))
	{
		cout<<"PionTransparency (low-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"PionTransparency (low-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/PiTrans/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("pitrans_le_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:500]"<<endl;
	gfile<<"set yrange[0:1]"<<endl;
	gfile<<"set label 2 '{/=24 Pion Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20' at 300, 3.2"<<endl;
	gfile<<"set xlabel 'Kinetic energy [MeV]' font 'Arial, 12'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 12'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1,1.2"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 400, 3 font 'Arial Bold, 16'"<<endl;

	double err9[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	make_data(PiTrLeEnp, PiTrLeLip, err9, PiTrLeLiperr, 9, 0);
	run("cp tmp/data.txt tmp/lip.txt");
	make_data(PiTrLeEnp, PiTrLeCp, err9, PiTrLeCperr, 9, 0);
	run("cp tmp/data.txt tmp/cp.txt");
	make_data(PiTrLeEnp, PiTrLeAlp, err9, PiTrLeAlperr, 9, 0);
	run("cp tmp/data.txt tmp/alp.txt");

	make_data(PiTrLeEn0, PiTrLeLi0, err9, PiTrLeLi0err, 9, 0);
	run("cp tmp/data.txt tmp/li0.txt");
	make_data(PiTrLeEn0, PiTrLeC0, err9, PiTrLeC0err, 9, 0);
	run("cp tmp/data.txt tmp/c0.txt");
	make_data(PiTrLeEn0, PiTrLeAl0, err9, PiTrLeAl0err, 9, 0);
	run("cp tmp/data.txt tmp/al0.txt");
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set nokey"<<endl;
	gfile<<"set title 'Lithium' font 'Arial, 12'"<<endl;
	gfile<<"plot 'tmp/lip.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/li0.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'black', '";
	gfile<<"tmp/lip.txt' using 1:2 title 'Charged Pions' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/li0.txt' using 1:2 title 'Neutral Pions' lt 1 pt 6 lc rgb 'black', '";	
	gfile<<filenameLp<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'NuWro (charged)', '";
	gfile<<filenameL0<<"' using 1:2 smooth csplines lt 2 lc rgb 'black' lw 5 title 'NuWro (neutral)'"<<endl;
	
	gfile<<"unset label 1"<<endl;
	gfile<<"unset label 2"<<endl;
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set title 'Carbon' font 'Arial, 12'"<<endl;
	gfile<<"plot 'tmp/cp.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/c0.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'black', '";
	gfile<<"tmp/cp.txt' using 1:2 title 'Charged Pions' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/c0.txt' using 1:2 title 'Neutral Pions' lt 1 pt 6 lc rgb 'black', '";	
	gfile<<filenameCp<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'NuWro (charged)', '";
	gfile<<filenameC0<<"' using 1:2 smooth csplines lt 2 lc rgb 'black' lw 5 title 'NuWro (neutral)'"<<endl;
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0.5"<<endl;
	gfile<<"set key at 450, -0.8"<<endl;

	gfile<<"set title 'Aluminium' font 'Arial, 12'"<<endl;
	gfile<<"plot 'tmp/alp.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al0.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'black', '";
	gfile<<"tmp/alp.txt' using 1:2 title 'Charged Pions' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al0.txt' using 1:2 title 'Neutral Pions' lt 1 pt 6 lc rgb 'black', '";	
	gfile<<filenameAp<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'NuWro (charged)', '";
	gfile<<filenameA0<<"' using 1:2 smooth csplines lt 2 lc rgb 'black' lw 5 title 'NuWro (neutral)'"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filenameLp + string(" and ") + filenameCp + string(" and ") + filenameAp + string(" and ") + filenameL0 + string(" and ") + filenameC0 + string(" and ") + filenameA0;
	plotlog(help);
	
	return 1;
}

void plotFZ ()
{
	calcFZ();
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output 'multiplots/fz_comp.eps'"<<endl;
	gfile<<"set xrange[0:5000]"<<endl;
	gfile<<"set yrange[0:3]"<<endl;
	gfile<<"set xlabel 'Momentum [MeV]'"<<endl;
	gfile<<"set ylabel 'Formation length [fm]'"<<endl;
	gfile<<"set key spacing 1"<<endl;
	gfile<<"set multiplot layout 2,2 title 'Formation Models'"<<endl;

	int nof = 0;
	int ile = 0;
	
	for (int k = 1; k < nof_fz; k++) if (fz_on[k]) nof++;
	
	if (nof == 0)
	{
		for (int k = 1; k < nof_fz; k++) fz_on[k] = true;
		nof = nof_fz - 1;
	}
	
	gfile<<"set title 'Nucleon Elastic'"<<endl;
	gfile<<"plot ";
	
	for (int i = 1; i < nof_fz; i++)
	{
		if (fz_on[i])
		{
			gfile<<"'tmp/nucleon_qel.txt' using 1:"<<i+2<<" title '"<<fzname[i]<<"' smooth csplines lw 5 ";//ls "<<i;
			ile++;
			if (ile < nof) gfile<<", ";
		}
	}
	
	ile = 0;
	
	gfile<<endl<<"set title 'Nucleon Inelastic'"<<endl;
	gfile<<"plot ";
	
	for (int i = 1; i < nof_fz; i++)
	{
		if (fz_on[i])
		{
			gfile<<"'tmp/nucleon_inel.txt' using 1:"<<i+2<<" title '"<<fzname[i]<<"' smooth csplines lw 5";// ls "<<i;
			ile++;
			if (ile < nof) gfile<<", ";
		}
	}
	
	ile = 0;

	gfile<<endl<<"set title 'Pion Elastic'"<<endl;
	gfile<<"plot ";
	
	for (int i = 1; i < nof_fz; i++)
	{
		if (fz_on[i])
		{
			gfile<<"'tmp/pion_qel.txt' using 1:"<<i+2<<" title '"<<fzname[i]<<"' smooth csplines lw 5 ";//ls "<<i;
			ile++;
			if (ile < nof) gfile<<", ";
		}
	}
	
	ile = 0;

	gfile<<endl<<"set title 'Pion Inelastic'"<<endl;
	gfile<<"plot ";
	
	for (int i = 1; i < nof_fz; i++) 
	{
		if (fz_on[i])
		{
			gfile<<"'tmp/pion_inel.txt' using 1:"<<i+2<<" title '"<<fzname[i]<<"' smooth csplines lw 5 ";//ls "<<i;
			ile++;
			if (ile < nof) gfile<<", ";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	
	cout<<endl<<"multiplots/fz_comp.eps created."<<endl<<endl;
}

int plotNomad(int fz, int xs)
{
	get_date();
	string prfilename;
	string pifilename;
	
	const int prbins = 8;
	const int pibins = 6;
	
	double evepr[4];
	double evepi[4];
		
	prfilename = string("results/Nomad/nomad_protons_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	prfilename = find_last(prfilename);

	pifilename = string("results/Nomad/nomad_pions_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	pifilename = find_last(pifilename);
	
	if (noFile(prfilename) or noFile(pifilename))
	{
		cout<<"Nomad ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"Nomad ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	findnof(prfilename, evepr);
	findnof(pifilename, evepi);
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/Nomad/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("Nomad_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	//gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1, 1.05"<<endl;
	gfile<<"set multiplot layout 2,2 title '{/=20 Average number of backward moving particles}'"<<endl;
	gfile<<"set xlabel 'Q^{2} [GeV]'"<<endl;
	gfile<<"set xrange [0:100]"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 90,-0.1 font 'Arial Bold, 16'"<<endl;
	gfile<<"set label 2 'Number of events with x backwards protons: ' at 5, -0.2"<<endl;
	gfile<<"set label 3 'NOMAD: ' at 20, -0.25"<<endl;
	gfile<<"set label 4 'NuWro: ' at 60, -0.25"<<endl;
	gfile<<"set label 5 '0: ' at 0, -0.3"<<endl;
	gfile<<"set label 6 '1: ' at 0, -0.35"<<endl;
	gfile<<"set label 7 '2: ' at 0, -0.4"<<endl;
	gfile<<"set label 8 '3: ' at 0, -0.45"<<endl;
	gfile<<"set label 9 '"<<NOMADmultiProton[0]<<"' at 20, -0.3"<<endl;
	gfile<<"set label 10 '"<<NOMADmultiProton[1]<<"' at 20, -0.35"<<endl;
	gfile<<"set label 11 '"<<NOMADmultiProton[2]<<"' at 20, -0.4"<<endl;
	gfile<<"set label 12 '"<<NOMADmultiProton[3]<<"' at 20, -0.45"<<endl;
	gfile<<"set label 13 '"<<evepr[0]<<"' at 60, -0.3"<<endl;
	gfile<<"set label 14 '"<<evepr[1]<<"' at 60, -0.35"<<endl;
	gfile<<"set label 15 '"<<evepr[2]<<"' at 60, -0.4"<<endl;
	gfile<<"set label 16 '"<<evepr[3]<<"' at 60, -0.45"<<endl;
	
	gfile<<"set label 17 'Number of events with x backwards pions: ' at 145, -0.2"<<endl;
	gfile<<"set label 18 'NOMAD: ' at 160, -0.25"<<endl;
	gfile<<"set label 19 'NuWro: ' at 200, -0.25"<<endl;
	gfile<<"set label 20 '0: ' at 140, -0.3"<<endl;
	gfile<<"set label 21 '1: ' at 140, -0.35"<<endl;
	gfile<<"set label 22 '2: ' at 140, -0.4"<<endl;
	gfile<<"set label 23 '3: ' at 140, -0.45"<<endl;
	gfile<<"set label 24 '"<<NOMADmultiPion[0]<<"' at 160, -0.3"<<endl;
	gfile<<"set label 25 '"<<NOMADmultiPion[1]<<"' at 160, -0.35"<<endl;
	gfile<<"set label 26 '"<<NOMADmultiPion[2]<<"' at 160, -0.4"<<endl;
	gfile<<"set label 27 '"<<NOMADmultiPion[3]<<"' at 160, -0.45"<<endl;
	gfile<<"set label 28 '"<<evepi[0]<<"' at 200, -0.3"<<endl;
	gfile<<"set label 29 '"<<evepi[1]<<"' at 200, -0.35"<<endl;
	gfile<<"set label 30 '"<<evepi[2]<<"' at 200, -0.4"<<endl;
	gfile<<"set label 31 '"<<evepi[3]<<"' at 200, -0.45"<<endl;
		
	gfile<<"set title 'Protons' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel '<#Bp>'"<<endl;
	gfile<<"set yrange [0:0.3]"<<endl;	

	make_data(NOMADprotonsQ2, NOMADprotonsBack, NOMADprotonsQ2err, NOMADprotonsBackerr, prbins, 0);
	run ("cp tmp/data.txt tmp/data2.txt");
	make_data(NOMADpionsQ2, NOMADpionsBack, NOMADpionsQ2err, NOMADpionsBackerr, pibins, 0);
	
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'NOMAD' pt 7 lc rgb 'black', '";
	gfile<<prfilename<<"' using 1:2 title 'NuWro' smooth csplines lt 1 lc rgb 'blue' lw 5"<<endl;
	
	gfile<<"unset label 1"<<endl;
	
	gfile<<"set title 'Pions' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel '<#B{/Symbol p}^{-}>'"<<endl;
	gfile<<"set yrange [0:0.04]"<<endl;		
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'NOMAD' pt 7 lc rgb 'black', '";
	gfile<<pifilename<<"' using 1:2 title 'NuWro' smooth csplines lt 1 lc rgb 'blue' lw 5"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + prfilename + string(" and ") + pifilename;
	plotlog(help);
	
	return 1;
}

int plotPiThe(int fz, int xs)
{
	get_date();
	string filenameC, filenameA, filenameM;
	const int bins = 20;
		
	filenameC = string("results/PiTrans/pitrans_he_carbon_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameC = find_last(filenameC);
	filenameA = string("results/PiTrans/pitrans_he_aluminium_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameA = find_last(filenameA);
	filenameM = string("results/PiTrans/pitrans_he_copper_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filenameM = find_last(filenameM);
	
	if (noFile(filenameC) or noFile(filenameA) or noFile(filenameM))
	{
		cout<<"PionTransparency (high-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"PionTransparency (high-energy) ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/PiTrans/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("pitrans_he_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[2.5:5]"<<endl;
	gfile<<"set yrange[0.2:1]"<<endl;
	gfile<<"set title '{/=24 Pion Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'p_{/Symbol p} [GeV/c]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 2, 0.9 font 'Arial Bold, 16'"<<endl;

	double err5[5] = {0, 0, 0, 0, 0};
	
	make_data(PTppi, PTC, err5, PTCerr, 5, 0);
	run("cp tmp/data.txt tmp/c.txt");
	make_data(PTppi, PTAl, err5, PTAlerr, 5, 0);
	run("cp tmp/data.txt tmp/al.txt");
	make_data(PTppi, PTCu, err5, PTCuerr, 5, 0);
	run("cp tmp/data.txt tmp/cu.txt");
	
	gfile<<"plot 'tmp/c.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', '";
	gfile<<"tmp/cu.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', '";
	gfile<<filenameC<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'Carbon', '";
	gfile<<filenameA<<"' using 1:2 smooth csplines lt 1 lc rgb 'red' lw 5 title 'Aluminium', '";
	gfile<<filenameM<<"' using 1:2 smooth csplines lt 1 lc rgb 'blue' lw 5 title 'Copper'"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filenameC + string(" and ") + filenameA + string(" and ") + filenameM;
	plotlog(help);
	
	return 1;
}

int plotAtmNu(int fz, int xs)
{
	get_date();
	string filename;
	const int bins = 8;
	string names[bins] = {"0pi", "pi+", "pi0", "pi-pi+", "2pi0", "pi+(>0*pi0)", "2pi+(>=0*pi0)", "pi-"};
	
	filename = string("results/AtmNu/atmnu_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename = find_last(filename);
	
	if (noFile(filename))
	{
		cout<<"Atmospheric Neutrinos ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"Atmospheric Neutrinos ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/Atmnu/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("atmnu_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set style data histogram"<<endl;
	gfile<<"set style histogram cluster gap 1"<<endl;
	gfile<<"set style fill solid noborder"<<endl;
	gfile<<"set xtic rotate by -45 scale 0"<<endl;
	gfile<<"set yrange[0:1]"<<endl;
	gfile<<"set title '{/=24 Atmospheric Neutrinos} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set ylabel 'Number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set label 1 '"<<fzname[fz]<<"' at 1, 0.8 font 'Arial Bold, 16'"<<endl;

	gfile<<"set xtics ('' -1";
	for (int i = 0; i < bins; i++) gfile<<", '"<<names[i]<<"' "<<i;
	gfile<<")"<<endl;

	ofstream gdata("tmp/data.txt");
	for (int i = 0; i < bins; i++)
	{
		gdata<<AtmosphericDeuterium[i]<<" "<<AtmosphericNeon[i]<<endl;
	}
	
	gfile<<"plot 'tmp/data.txt' using 1 title 'Deuterium' lc rgb 'red' fill solid 0.50, '";
	gfile<<filename<<"' using 2 title 'NuWro (before FSI)' lc rgb 'green' fill solid 0.50, '";
	gfile<<"tmp/data.txt' using 2 title 'Neon' lc rgb 'red', '";
	gfile<<filename<<"' using 3 title 'NuWro (after FSI)' lc rgb 'green'";
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

int plotMBCCrat(int fz, int xs)
{
	get_date();
	string filename;
	const int bins = 13;
	
	filename = string("results/MBratio/mb_ratio_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename = find_last(filename);
	
	if (noFile(filename))
	{
		cout<<"MB ccpi+/ccqe ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"MB ccpi+/ccqe ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
		
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/MBratio/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("mb_ratio_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set size 1, 2"<<endl;
	gfile<<"set multiplot"<<endl;
	gfile<<"set xrange[0:2.5]"<<endl;
	gfile<<"set xlabel 'E_{/Symbol n} [GeV]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set nokey"<<endl;

	gfile<<"set origin 0, 0"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:1.2]"<<endl;
	gfile<<"set title 'MiniBooNE CC{/Symbol p}^{+}-like/CCQE-like on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CC{/Symbol p}^{+}-like}/{/Symbol s}_{CCQE-like}'"<<endl;

	make_data(mbraten, mbrat, mbratenerr, mbraterr, bins, 0);
	run("cp tmp/data.txt tmp/data1.txt");
	
	gfile<<"plot 'tmp/data1.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data1.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:4 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;
	
	gfile<<"set origin 0.5, 0"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:10e-38]"<<endl;
	gfile<<"set title 'MiniBooNE CCQE-like on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CCQE-like} [10^{-38}cm^{2}]'"<<endl;
	gfile<<"set ytics ('0' 0, '1' 1e-38, '2' 2e-38, '3' 3e-38, '4' 4e-38, '5' 5e-38, '6' 6e-38, '7' 7e-38, '8' 8e-38, '9' 9e-38)"<<endl;

	make_data(mbccqeen, mbccqe, mbccqeenerr, mbccqeerr, 14, 0);
	run("cp tmp/data.txt tmp/data2.txt");
	
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:2 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;
	
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:15e-38]"<<endl;
	gfile<<"set title 'MiniBooNE CC{/Symbol p}^{+}-like on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CC{/Symbol p}^{+}-like}'"<<endl;
	gfile<<"set ytics ('0' 0, '2' 2e-38, '4' 4e-38, '6' 6e-38, '8' 8e-38, '10' 10e-38, '12' 12e-38, '14' 14e-38)"<<endl;

	make_data(mbccpien, mbccpi, mbccpienerr, mbccpierr, 27, 0);
	run("cp tmp/data.txt tmp/data3.txt");
	
	gfile<<"plot 'tmp/data3.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data3.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:3 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;
	
	gfile<<"set origin 0.5, 0.5"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:3e-38]"<<endl;
	gfile<<"set title 'MiniBooNE CC{/Symbol p}^{0} on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CC{/Symbol p}^{0}}'"<<endl;
	gfile<<"set ytics ('0' 0, '0.5' 0.5e-38, '1' 1e-38, '1.5' 1.5e-38, '2' 2e-38, '2.5' 2.5e-38, '3' 3e-38)"<<endl;

	make_data(mbccpi0en, mbccpi0, mbccpi0enerr, mbccpi0err, 14, 0);
	run("cp tmp/data.txt tmp/data4.txt");
	
	gfile<<"plot 'tmp/data4.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data4.txt' using 1:2 title 'MB data' pt 7 lc rgb 'black', '";
	gfile<<filename<<"' using 1:5 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;
	
	gfile<<"set origin 0, 1"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:1]"<<endl;
	gfile<<"set title 'MiniBooNE CCQE-like/CCtotal on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CCQE-like}/{/Symbol s}_{CCtotal}'"<<endl;
	gfile<<"unset ytics"<<endl;
	gfile<<"set ytics"<<endl;
	
	gfile<<"plot '"<<filename<<"' using 1:6 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;
	
	gfile<<"set origin 0.5, 1"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0:0.5]"<<endl;
	gfile<<"set title 'MiniBooNE CC{/Symbol p}^{+}-like/CCtotal on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CC{/Symbol p}^{+}-like}/{/Symbol s}_{CCtotal}'"<<endl;
	
	gfile<<"plot '"<<filename<<"' using 1:7 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;

	gfile<<"set origin 0, 1.5"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[1:1.2]"<<endl;
	gfile<<"set title 'MiniBooNE CCQE-like/CCQE-true on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CCQE-like}/{/Symbol s}_{CCQE-true}'"<<endl;
	
	gfile<<"plot '"<<filename<<"' using 1:8 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;	
	
	gfile<<"set origin 0.5, 1.5"<<endl;
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set yrange[0.5:1]"<<endl;
	gfile<<"set title 'MiniBooNE CC{/Symbol p}^{+}-like/CC{/Symbol p}^{+}-true on CH_{2}'"<<endl;
	gfile<<"set ylabel '{/Symbol s}_{CC{/Symbol p}^{+}-like}/{/Symbol s}_{CC{/Symbol p}^{+}-true}'"<<endl;
	
	gfile<<"plot '"<<filename<<"' using 1:9 smooth csplines lt 2 lw 5 lc rgb 'black' title 'NuWro'"<<endl;	
		
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

int plotnuintPr(int fz, int xs)
{
	get_date();
	string filename;
	const int bins = 13;
	
	filename = string("results/NUINT/proton_") + fzwork[fz] + sep + xsec[xs] + ("*.txt");
	filename = find_last(filename);
	
	if (noFile(filename))
	{
		cout<<"NUINT pr("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		logfile<<"NUINT pr ("<<fzname[fz]<<"): there are no txt files, make calculation first!"<<endl<<endl;
		return 1;
	}
		
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "plots/NUINT/";
	run(string("mkdir -p ") + plotfile);
	plotfile += string("proton_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".eps");
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:2000]"<<endl;
	gfile<<"set yrange[0:0.7]"<<endl;
	gfile<<"set title 'p^{12}C' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'proton momentum [MeV/c]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'Probability' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1.5, 1"<<endl;
	gfile<<"set key out"<<endl;
	//gfile<<"set label 1 '"<<fzname[fz]<<"' at 1.5, 0.2 font 'Arial Bold, 16'"<<endl;

	make_data(mbraten, mbrat, mbratenerr, mbraterr, bins, 0);
	
	gfile<<"plot '"<<filename<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'there was interaction', '";
	gfile<<filename<<"' using 1:3 smooth csplines lt 1 lw 5 lc rgb 'blue' title 'there was one {/Symbol p}', '";
	gfile<<filename<<"' using 1:4 smooth csplines lt 2 lw 5 lc rgb 'blue' title 'there was one {/Symbol p}^{+}', '";
	gfile<<filename<<"' using 1:5 smooth csplines lt 3 lw 5 lc rgb 'blue' title 'there was one {/Symbol p}^{-}', '";
	gfile<<filename<<"' using 1:6 smooth csplines lt 4 lw 5 lc rgb 'blue' title 'there was one {/Symbol p}^{0}', '";
	gfile<<filename<<"' using 1:7 smooth csplines lt 1 lw 5 lc rgb 'green' title 'there was more than one {/Symbol p}'";
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ") + filename;
	plotlog(help);
	
	return 1;
}

void t2k_anal_plot()
{
	run("mkdir -p plots/t2k_anal/");
	
	string f1[7] = {"results/t2k_anal/CCQE_", "results/t2k_anal/CCPiP_", "results/t2k_anal/CCPiM_", "results/t2k_anal/CCPi0_", "results/t2k_anal/NCPiP_", "results/t2k_anal/NCPiM_", "results/t2k_anal/NCPi0_"};
	string f2[2] = {"before_fsi_", "after_fsi_"};
	string f3[3] = {"neutrino_energy_", "lepton_momentum_", "lepton_angle_"};
	string f4[2] = {"carbon.txt", "oxygen.txt"};

	string of1[7] = {"plots/t2k_anal/CCQE_", "plots/t2k_anal/CCPiP_", "plots/t2k_anal/CCPiM_", "plots/t2k_anal/CCPi0_", "plots/t2k_anal/NCPiP_", "plots/t2k_anal/NCPiM_", "plots/t2k_anal/NCPi0_"};
	string of2[2] = {"before_fsi_", "after_fsi_"};
	string of3[3] = {"neutrino_energy_", "lepton_momentum_", "lepton_angle_"};
	string of4[2] = {"carbon.eps", "oxygen.eps"};
	string of4b[2] = {"carbon_ratio.eps", "oxygen_ratio.eps"};
	
	string tit1[7] = {"CCQE ", "CC{/Symbol p}^{+} ", "CC{/Symbol p}^{-} ", "CC{/Symbol p}^{0} ", "NC{/Symbol p}^{+} ", "NC{/Symbol p}^{-} ", "NC{/Symbol p}^{0} "};
	string tit2[2] = {"carbon ", "oxygen "};
	string tit3[2] = {"(before FSI)", "(after FSI)"};
	
	string lab[3] = {"incident neutrino energy [MeV]", "outgoing lepton momentum [MeV]", "cos{/Symbol Q}"};
	
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 2; l++)
				{
					string fname = f1[i] + f2[j] + f3[k] + f4[l];
					string pname = of1[i] + of2[j] + of3[k] + of4[l];
					string rname = of1[i] + of2[j] + of3[k] + of4b[l];
					
					ofstream gfile("tmp/gnuplot.gnu");

					gfile << "set terminal postscript eps enhanced color 'Arial' 16 " << endl;
					gfile << "set output '" << pname << "'" << endl;
					gfile << "set title '" << tit1[i] + tit2[l] + tit3[j] << "' " << endl;
					gfile << "set xlabel '" << lab[k] << "' " << endl;
					if (k != 2) gfile << "set xrange [200:1400]" << endl;
					gfile << "set ylabel 'cross section [cm^{2}]'" << endl;
					
					gfile << "plot '" << fname << "' u 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'Fermi Gas', '";
					gfile << fname << "' u 1:3 smooth csplines lt 1 lw 5 lc rgb 'blue' title 'Spectral Function'" << endl;
					
					gfile << "set output '" << rname << "'" << endl;
					gfile << "set ylabel ''" << endl;
					gfile << "set yrange [1:1.5]" << endl;					
					gfile << "plot '" << fname << "' u 1:4 smooth csplines lt 1 lw 5 lc rgb 'black' title 'FG/SF ratio'";
					
					gfile.close();
					run("gnuplot tmp/gnuplot.gnu");
				}
			}
		}
	}					
}
