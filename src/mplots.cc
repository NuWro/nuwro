#include "fsi.h"
#include "data.h"
#include "plots.h"
#include "mplots.h"
#include <string>

std::string const colours[11] = {"black", "red", "blue", "green", "purple", "yellow", "#006400", "brown", "coral", "#00FFFF", "#800000"};

int mplotK2K()
{
	get_date();
	const int bins = 8;
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}	
	
	string filename[nof_fz];
	
	for (int k = 0; k < nof_fz; k++)
	{
		filename[k] = string("results/K2K/k2k_") + fzwork[k] + ("*.txt");
		filename[k] = find_last(filename[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and noFile(filename[k]))
		{
			cout<<"multiK2K ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"multiK2K ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_k2k.eps";
	
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

	make_data(K2Kx, K2Ky, K2Kxerr, K2Kyerr, bins, 1);
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'K2K data' with points 7, '";
	gfile<<filename[first]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[first]<<"' using 10:11 with lines lt 1 lw 5 lc rgb 'black' notitle, '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[k]<<"' using 1:3 title '"<<fzname[k]<<"' smooth csplines lw 5 lt "<<k+1<<", '";
			gfile<<filename[k]<<"' using 10:12 with lines notitle lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
				
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) help += filename[k];
	}
	plotlog(help);
	
	return 1;
}

int mplotMBCC()
{
	get_date();
	const int bins = 11;
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}	
	
	string filename[nof_fz];
	string filename2[nof_fz];
	
	for (int k = 0; k < nof_fz; k++)
	{
		filename[k] = string("results/MBCC/mb_cc_") + fzwork[k] + ("*.txt");
		filename[k] = find_last(filename[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		filename2[k] = string("results/MBCC/mb_cc_normalized_") + fzwork[k] + ("*.txt");
		filename2[k] = find_last(filename2[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and (noFile(filename[k]) or noFile(filename2[k])))
		{
			cout<<"multiMBCC ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"multiMBCC ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_mbcc.eps";
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:1.2]"<<endl;
	gfile<<"set yrange[0:4.5e-38]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [GeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-38}cm^{2}/GeV/CH_{2}]' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set ytics ('0' 0, '0.5' 0.5e-38, '1.0' 1e-38, '1.5' 1.5e-38, '2.0' 2e-38, '2.5' 2.5e-38, '3.0' 3e-38, '3.5' 3.5e-38)"<<endl;
	gfile<<"set size 2,1.2"<<endl;
	gfile<<"set multiplot"<<endl; 
	gfile<<"set label 2 '{/=24 MB} CC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<")' at 1, 5e-38 font 'Arial, 20'"<<endl;

	make_data(MBCCmom, MBCCxsec, MBCCmomerr, MBCCxsecerr, bins, 0);
	
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set size 1,1"<<endl;
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'K2K data' with points 7, '";
	gfile<<filename[first]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[k]<<"' using 1:3 title '"<<fzname[k]<<"' smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
		
	gfile<<endl<<"set ylabel 'arbitrary units' font 'Arial, 16'"<<endl;	
	gfile<<"unset label 2"<<endl;
	gfile<<"set origin 1,0"<<endl;
	gfile<<"set size 1,1"<<endl;
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'K2K data' with points 7, '";
	gfile<<filename2[first]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	
	ile=0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename2[k]<<"' using 1:3 title '"<<fzname[k]<<"' smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile<<endl<<"unset multiplot"<<endl;
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) help += filename[k] + string(" ") + filename2[k];
	}
	plotlog(help);
	
	return 1;
}

int mplotMBangle()
{
	get_date();
	
	int nof = 0;
	int ile = 0;
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}

	string filename[2][nof_fz]; //0 - neutrino, 1 - antineutrino, 2,3 - normalized

	const int bins  = 18;
	const int binsa = 10;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			filename[0][k] = string("results/MB/mb_normalized_") + fzwork[k] + ("*.txt");
			filename[1][k] = string("results/MB/mb_normalized_anti_") + fzwork[k] + ("*.txt");

			for (int i = 0; i < 2; i++) filename[i][k] = find_last(filename[i][k]);
		}
	}	
	
	bool allfiles = true;
	
	for (int k = 0; k <nof_fz; k++)
	{
		
		for (int i = 0; i < 2; i++) if (noFile(filename[i][k])) allfiles = false;
		if (fz_on[k] and !allfiles)
		{
			cout<<"multiMB angle ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"multiMB angle ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
		
		allfiles = true;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_mb_angle.eps";
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[-1:0]"<<endl;
	gfile<<"set yrange[0:0.3e-39]"<<endl;
	gfile<<"set xlabel 'cos{/Symbol O}_{/Symbol p_{0}}'"<<endl;
	//gfile<<"set ylabel 'd{/Symbol s}/dcos{/Symbol O}_{/Symbol p_{0}} [10^{-39} cm^{2}/nucleon]'"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;	
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1,0.6"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl; 
	gfile<<"set key left"<<endl;
	gfile<<"set ytics ('0' 0, '0.05' 0.05e-39, '0.1' 0.1e-39, '0.15' 0.15e-39, '0.2' 0.2e-39)"<<endl;
	
	make_data(MBAnglex, MBAngley, MBAnglexerr, MBAngleyerr, bins, 0);
	run("cp tmp/data.txt tmp/data2.txt");
	make_data(MBantiAnglex, MBantiAngley, MBantiAnglexerr, MBantiAngleyerr, binsa, 0);
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0"<<endl;
	gfile<<"set title 'neutrino mode'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 with xyerrorbars lw 2 lc rgb 'black' notitle, '";
	gfile<<"tmp/data2.txt' using 1:2 with points pt 7 lc rgb 'black' title 'MB data', '";
	gfile<<filename[0][first]<<"' using 10:11 smooth csplines title 'NuWro (before FSI)' lt 1 lw 3 lc rgb 'black' , '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[0][k]<<"' using 10:12 ";
			gfile<<"smooth csplines lw 5 lt "<<k+1<<" ";
			gfile<<"title '"<<fzname[k]<<"' ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
			
	gfile<<endl<<"set label 1 '{/=16 MiniBooNE} NC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<")' at -3, 0.65e-39 font 'Arial, 20'"<<endl;
	
	gfile<<"set yrange[0:0.15e-39]"<<endl;
	gfile<<"set ytics ('0' 0, '0.02' 0.02e-39, '0.04' 0.04e-39, '0.06' 0.06e-39, '0.08' 0.08e-39)"<<endl;
			
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0"<<endl;
	gfile<<"set title 'anti-neutrino mode'"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'MB data' with points pt 7 lc rgb 'black', '";
	gfile<<filename[1][first]<<"' using 10:11 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";

	ile = 0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[1][k]<<"' using 10:12 ";
			gfile<<"smooth csplines lw 5 lt "<<k+1<<" ";
			gfile<<"title '"<<fzname[k]<<"' ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile<<endl<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{ 
				for (int i = 0; i < 2; i++) help += filename[i][k] + (" ");
		}
	}
	
	plotlog(help);
	
	return 1;
}


int mplotMB()
{
	/*get_date();
	
	int nof = 0;
	int ile = 0;
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}

	string filename[4][nof_fz]; //0 - neutrino, 1 - antineutrino, 2,3 - normalized

	const int bins  = 11;
	const int binsa = 10;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			filename[0][k] = string("results/MB/mb_") + fzwork[k] + ("*.txt");
			filename[1][k] = string("results/MB/mb_anti_") + fzwork[k] + ("*.txt");
			filename[2][k] = string("results/MB/mb_normalized_") + fzwork[k] + ("*.txt");
			filename[3][k] = string("results/MB/mb_normalized_anti_") + fzwork[k] + ("*.txt");

			for (int i = 0; i < 4; i++) filename[i][k] = find_last(filename[i][k]);
		}
	}
	
	
	bool allfiles = true;
	
	for (int k = 0; k <nof_fz; k++)
	{
		
		for (int i = 0; i < 4; i++) if (noFile(filename[i][k])) allfiles = false;
		if (fz_on[k] and !allfiles)
		{
			cout<<"multiMB ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"multiMB ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
		
		allfiles = true;
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_mb.eps";
	
	gfile<<setprecision (3);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:1000]"<<endl;
	gfile<<"set yrange[0:2e-42]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{0} momentum [MeV]'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-42} cm^{2}/MeV/nucleon]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1,1.2"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl; 
	//gfile<<"set nokey"<<endl;
	gfile<<"set ytics ('0' 0, '0.5' 0.5e-42, '1.0' 1e-42, '1.5' 1.5e-42, '2.0' 2e-42)"<<endl;

	make_data(MBMomx, MBMomy, MBMomxerr, MBMomyerr, bins, 0);
	run("cp tmp/data.txt tmp/data2.txt");
	make_data(MBantiMomx, MBantiMomy, MBantiMomxerr, MBantiMomyerr, binsa, 0);
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set title 'neutrino'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'MB data' with points 7, '";
	gfile<<filename[0][first]<<"' using 1:2 smooth csplines lt 1 lw 3 lc rgb 'black' title 'NuWro (before FSI)', '";
	
	int leg;
	
	if (nof < 5) leg = 100;
	else leg = nof/2;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[0][k]<<"' using 1:3 ";
			if (k < leg) gfile<<"title '"<<fzname[k]<<"' ";
			else gfile<<"notitle ";
			gfile<<"smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
			
	gfile<<endl<<"set label 1 '{/=24 MiniBooNE} NC{/Symbol p}^{0} production on CH_{2} ("<<day<<" "<<fullm<<" "<<year<<")' at -1100, 3e-42 font 'Arial, 20'"<<endl;
			
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0.5"<<endl;
	gfile<<"set title 'neutrino - normalized'"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 notitle with points 7, '";

	ile = 0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[2][k]<<"' using 1:3 ";
			if (k >= leg) gfile<<"title '"<<fzname[k]<<"' ";
			else gfile<<"notitle ";
			gfile<<"smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}

	gfile<<endl<<"unset label 1"<<endl;
	
	gfile<<"set yrange[0:0.8e-42]"<<endl;
	gfile<<"set ytics ('0' 0, '2.0' 2e-43, '4.0' 4e-43, '6.0' 6e-43, '8.0' 8e-43)"<<endl;
	gfile<<"set title 'antineutrino'"<<endl;
	gfile<<"set ylabel 'd{/Symbol s}/dp [10^{-43} cm^{2}/MeV/nucleon]'"<<endl;

	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 notitle with points 7, '";
	
	ile = 0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[1][k]<<"' using 1:3 notitle smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile<<endl<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0.5, 0"<<endl;
	//gfile<<"set key out vert"<<endl;
	gfile<<"set title 'antineutrino - normalized"<<endl;
	gfile<<"set ylabel 'arbitrary units'"<<endl;
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 notitle with points 7, '";

	ile = 0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[3][k]<<"' using 1:3 notitle smooth csplines lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile<<endl<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{ 
				for (int i = 0; i < 4; i++) help += filename[i][k] + (" ");
		}
	}
	
	plotlog(help);*/
	
	mplotMBangle();
	
	return 1;
}

int mplotSB()
{
	get_date();

	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
	
	string filename[nof_fz];
	const int bins = 9;
	
	for (int k = 0; k < nof_fz; k++)
	{
		filename[k] = string("results/SB/sb_") + fzwork[k] + ("*.txt");
		filename[k] = find_last(filename[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and noFile(filename[k]))
		{
			cout<<"SB ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"SB ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
		
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_sb.eps";
	
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

	make_data(SBx, SBy, SBxerr, SByerr, bins, 1);
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'SB data' with points 7, '";
	gfile<<filename[first]<<"' using 1:2 smooth csplines lt 1 lw 5 lc rgb 'black' title 'NuWro (before FSI)', '";
	gfile<<filename[first]<<"' using 10:11 with lines lt 1 lw 5 lc rgb 'black' notitle, '";
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[k]<<"' using 1:3 title '"<<fzname[k]<<"' smooth csplines lw 5 lt "<<k+1<<", '";
			gfile<<filename[k]<<"' using 10:12 with lines notitle lw 5 lt "<<k+1<<" ";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) help += filename[k];
	}

	plotlog(help);
	
	return 1;
}

int mplotPNS()
{
	get_date();
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
		
	string filenameC[nof_fz], filenameI[nof_fz];
	
	const int bins = 10;
	
	for (int k = 0; k < nof_fz; k++)
	{	
		filenameC[k] = string("results/PNS/pns_carbon_") + fzwork[k] + ("*.txt");
		filenameC[k] = find_last(filenameC[k]);

		filenameI[k] = string("results/PNS/pns_iron_") + fzwork[k] + ("*.txt");
		filenameI[k] = find_last(filenameI[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and (noFile(filenameC[k]) or noFile(filenameI[k])))
		{
			cout<<"Pion-Nucleus ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"Pion-Nucleus ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_pns.eps";
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:550]"<<endl;
	gfile<<"set yrange[0:500]"<<endl;
	gfile<<"set xlabel '{/Symbol p}^{+} kinetic energy [MeV]'"<<endl;
	gfile<<"set ylabel 'cross section [mb]'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 0.8, 1.4"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl;
	gfile<<"set key out vert bot right"<<endl;
	gfile<<"set label 2 '{/Symbol p}^{+}N scattering ("<<day<<" "<<fullm<<" "<<year<<")' at 200, 700 font 'Arial, 20'"<<endl;	
	
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
	
	gfile<<"set size 0.65, 0.6"<<endl;
	gfile<<"set origin 0, 0.6"<<endl;
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
	gfile<<filenameC[first]<<"' using 1:2 smooth csplines lt "<<first+1<<" lc rgb 'black' lw 3 title 'Reaction', '";
	gfile<<filenameC[first]<<"' using 1:3 smooth csplines lt "<<first+1<<" lc rgb 'red' lw 3 title 'Absorption', '";
	gfile<<filenameC[first]<<"' using 1:4 smooth csplines lt "<<first+1<<" lc rgb 'blue' lw 3 title 'CEX', '";
	gfile<<filenameC[first]<<"' using 1:5 smooth csplines lt "<<first+1<<" lc rgb 'green' lw 3 title 'Inelastic'";
	
	//ile++;
	
	if (nof > 1) gfile<<", '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and k != first)
		{
			gfile<<filenameC[k]<<"' using 1:2 smooth csplines lt "<<k+1<<" lc rgb 'black' lw 3 notitle, '";
			gfile<<filenameC[k]<<"' using 1:3 smooth csplines lt "<<k+1<<" lc rgb 'red' lw 3 notitle, '";
			gfile<<filenameC[k]<<"' using 1:4 smooth csplines lt "<<k+1<<" lc rgb 'blue' lw 3 notitle, '";
			gfile<<filenameC[k]<<"' using 1:5 smooth csplines lt "<<k+1<<" lc rgb 'green' lw 3 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
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
	
	gfile<<endl<<"set size 0.78,0.6"<<endl;
	gfile<<"unset label 2"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set xrange[0:400]"<<endl;	
	gfile<<"set title 'Iron'"<<endl;
	gfile<<"set yrange[0:1400]"<<endl;
	gfile<<"plot 'tmp/iar.txt' using 1:2:3:4 notitle with xyerrorbars pt 7 lc rgb 'black', ";
	gfile<<"'tmp/iaa.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/iaa.txt' using 1:2 notitle lt 1 pt 7 lc rgb 'red', ";
	gfile<<"'tmp/ias.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', ";
	gfile<<"'tmp/iai.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'green', ";
	gfile<<"'tmp/ina.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'red', ";
	gfile<<"'tmp/ina.txt' using 1:2 notitle lt 1 pt 6 lc rgb 'red', '";
	
	ile = 0;
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{		
			gfile<<filenameI[k]<<"' using 1:2 smooth csplines lt "<<k+1<<" lc rgb 'black' lw 3 title '"<<fzname[k]<<"', '";
			gfile<<filenameI[k]<<"' using 1:3 smooth csplines lt "<<k+1<<" lc rgb 'red' lw 3 notitle, '";
			gfile<<filenameI[k]<<"' using 1:4 smooth csplines lt "<<k+1<<" lc rgb 'blue' lw 3 notitle, '";
			gfile<<filenameI[k]<<"' using 1:5 smooth csplines lt "<<k+1<<" lc rgb 'green' lw 3 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile<<endl<<"unset multiplot"<<endl;
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) help += filenameC[k] + string(" and ") + filenameI[k];
	}
	plotlog(help);
	
	return 1;
}

int mplotPrThe()
{
	get_date();
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
	
	string filenameC[nof_fz], filenameI[nof_fz];
	
	const int bins = 20;
	
	for (int k = 0; k < nof_fz; k++)
	{		
		filenameC[k] = string("results/PrTrans/prtrans_he_carbon_") + fzwork[k] + ("*.txt");
		filenameC[k] = find_last(filenameC[k]);
		filenameI[k] = string("results/PrTrans/prtrans_he_iron_") + fzwork[k] + ("*.txt");
		filenameI[k] = find_last(filenameI[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and (noFile(filenameC[k]) or noFile(filenameI[k])))
		{
			cout<<"ProtonTransparency (high-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"ProtonTransparency (high-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_prtrans_he.eps";
	
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
	gfile<<"set key out vert bot right"<<endl;

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
	gfile<<filenameC[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'black' lw 5 title 'Carbon', '";
	gfile<<filenameI[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'red' lw 5 title 'Iron', '";
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameC[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'black' lw 5 title '"<<fzname[k]<<"', '";
			gfile<<filenameI[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'red' lw 5 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
			
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += filenameC[k] + string(" and ") + filenameI[k];
	plotlog(help);
	
	return 1;
}

int mplotPrTle()
{
	get_date();

	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
	
	string filenameL[nof_fz], filenameC[nof_fz], filenameA[nof_fz];
	const int bins = 10;
				
	for (int k = 0; k < nof_fz; k++)
	{		
		filenameL[k] = string("results/PrTrans/prtrans_le_lithium_") + fzwork[k] + ("*.txt");
		filenameL[k] = find_last(filenameL[k]);
		filenameC[k] = string("results/PrTrans/prtrans_le_carbon_") + fzwork[k] + ("*.txt");
		filenameC[k] = find_last(filenameC[k]);
		filenameA[k] = string("results/PrTrans/prtrans_le_aluminium_") + fzwork[k] + ("*.txt");
		filenameA[k] = find_last(filenameA[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and (noFile(filenameL[k]) or noFile(filenameC[k]) or noFile(filenameA[k])))
		{
			cout<<"ProtonTransparency (low-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"ProtonTransparency (low-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
			
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_prtrans_le.eps";
	
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
	gfile<<"set key out vert bot right"<<endl;

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
	gfile<<filenameL[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'black' lw 5 title 'Lithium', '";
	gfile<<filenameC[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'red' lw 5 title 'Carbon', '";
	gfile<<filenameA[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'blue' lw 5 title 'Aluminium', '";

		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameL[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'black' lw 5 title '"<<fzname[k]<<"', '";
			gfile<<filenameC[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'red' lw 5 notitle, '";
			gfile<<filenameA[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'blue' lw 5 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += filenameL[k] + string(" and ") + filenameC[k] + string(" and ") + filenameA[k];
	plotlog(help);
	
	return 1;
}

int mplotPiThe()
{
	get_date();
	string filenameC[nof_fz], filenameA[nof_fz], filenameM[nof_fz];
	const int bins = 20;
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
					
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]){		
		filenameC[k] = string("results/PiTrans/pitrans_he_carbon_") + fzwork[k] + ("*.txt");
		filenameC[k] = find_last(filenameC[k]);
		filenameA[k] = string("results/PiTrans/pitrans_he_aluminium_") + fzwork[k] + ("*.txt");
		filenameA[k] = find_last(filenameA[k]);
		filenameM[k] = string("results/PiTrans/pitrans_he_copper_") + fzwork[k] + ("*.txt");
		filenameM[k] = find_last(filenameM[k]);
		}
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and (noFile(filenameC[k]) or noFile(filenameA[k]) or noFile(filenameM[k])))
		{
			cout<<"PionTransparency (high-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"PionTransparency (high-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
		
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_pitrans_he.eps";
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[1:5]"<<endl;
	gfile<<"set yrange[0.2:1]"<<endl;
	gfile<<"set title '{/=24 Pion Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set xlabel 'Q^{2} [GeV]' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set key out vert bot right"<<endl;

	double err5[5] = {0, 0, 0, 0, 0};
	
	make_data(PTQ2, PTC, PTQ2err, PTCerr, 5, 0);
	run("cp tmp/data.txt tmp/c.txt");
	make_data(PTQ2, PTAl, PTQ2err, PTAlerr, 5, 0);
	run("cp tmp/data.txt tmp/al.txt");
	make_data(PTQ2, PTCu, PTQ2err, PTCuerr, 5, 0);
	run("cp tmp/data.txt tmp/cu.txt");
	
	gfile<<"plot 'tmp/c.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'red', '";
	gfile<<"tmp/cu.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'blue', '";
	gfile<<filenameC[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'black' lw 5 title 'Carbon', '";
	gfile<<filenameA[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'red' lw 5 title 'Aluminium', '";
	gfile<<filenameM[first]<<"' using 1:2 smooth csplines lt "<<first + 1<<" lc rgb 'blue' lw 5 title 'Copper', '";

		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameC[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'black' lw 5 title '"<<fzname[k]<<"', '";
			gfile<<filenameA[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'red' lw 5 notitle, '";
			gfile<<filenameM[k]<<"' using 1:2 smooth csplines lt "<<k + 1<<" lc rgb 'blue' lw 5 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += filenameC[k] + string(" and ") + filenameA[k] + string(" and ") + filenameM[k];
	plotlog(help);
		
	return 1;
}

int mplotPiTle()
{
	get_date();
	string filenameLp[nof_fz], filenameCp[nof_fz], filenameAp[nof_fz], filenameL0[nof_fz], filenameC0[nof_fz], filenameA0[nof_fz];
	const int bins = 10;
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		filenameLp[k] = string("results/PiTrans/pitrans_le_pip_lithium_") + fzwork[k] + ("*.txt");
		filenameLp[k] = find_last(filenameLp[k]);
		filenameCp[k] = string("results/PiTrans/pitrans_le_pip_carbon_") + fzwork[k] + ("*.txt");
		filenameCp[k] = find_last(filenameCp[k]);
		filenameAp[k] = string("results/PiTrans/pitrans_le_pip_aluminium_") + fzwork[k] + ("*.txt");
		filenameAp[k] = find_last(filenameAp[k]);
		filenameL0[k] = string("results/PiTrans/pitrans_le_pi0_lithium_") + fzwork[k] + ("*.txt");
		filenameL0[k] = find_last(filenameL0[k]);
		filenameC0[k] = string("results/PiTrans/pitrans_le_pi0_carbon_") + fzwork[k] + ("*.txt");
		filenameC0[k] = find_last(filenameC0[k]);
		filenameA0[k] = string("results/PiTrans/pitrans_le_pi0_aluminium_") + fzwork[k] + ("*.txt");
		filenameA0[k] = find_last(filenameA0[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{		
		if (fz_on[k] and (noFile(filenameLp[k]) or noFile(filenameCp[k]) or noFile(filenameAp[k]) or noFile(filenameL0[k]) or noFile(filenameC0[k]) or noFile(filenameA0[k])))
		{
			cout<<"PionTransparency (low-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"PionTransparency (low-energy) ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_pitrans_le.eps";
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set xrange[0:500]"<<endl;
	gfile<<"set yrange[0:1]"<<endl;
	gfile<<"set label 2 '{/=24 Pion Transparency} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20' at 200, 5"<<endl;
	gfile<<"set xlabel 'Kinetic energy [MeV]' font 'Arial, 12'"<<endl;
	gfile<<"set ylabel 'Transparency' font 'Arial, 12'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1,1.7"<<endl;
	gfile<<"set origin 0,0"<<endl;
	gfile<<"set multiplot"<<endl;

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
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameLp[k]<<"' using 1:2 smooth csplines lt 1 lc rgb '"<<colours[k]<<"' lw 5 title 'NuWro (charged)', '";
			gfile<<filenameL0[k]<<"' using 1:2 smooth csplines lt 2 lc rgb '"<<colours[k]<<"' lw 5 title 'NuWro (neutral)'";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	ile = 0;
	
	gfile<<endl<<"unset label 2"<<endl;
	
	gfile<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 0.5"<<endl;
	gfile<<"set title 'Carbon' font 'Arial, 12'"<<endl;
	gfile<<"plot 'tmp/cp.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/c0.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'black', '";
	gfile<<"tmp/cp.txt' using 1:2 title 'Charged Pions' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/c0.txt' using 1:2 title 'Neutral Pions' lt 1 pt 6 lc rgb 'black', '";	
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameCp[k]<<"' using 1:2 smooth csplines lt 1 lc rgb '"<<colours[k]<<"' lw 5 title 'NuWro (charged)', '";
			gfile<<filenameC0[k]<<"' using 1:2 smooth csplines lt 2 lc rgb '"<<colours[k]<<"' lw 5 title 'NuWro (neutral)'";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	ile = 0;
		
	gfile<<endl<<"set size 0.5, 0.5"<<endl;
	gfile<<"set origin 0, 1"<<endl;
	gfile<<"set key horizontal"<<endl;
	gfile<<"set key at 1200, 0.5"<<endl;
	
	gfile<<"set title 'Aluminium' font 'Arial, 12'"<<endl;
	gfile<<"plot 'tmp/alp.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al0.txt' using 1:2:3:4 notitle with xyerrorbars lt 1 pt 6 lc rgb 'black', '";
	gfile<<"tmp/alp.txt' using 1:2 title 'Charged Pions' lt 1 pt 7 lc rgb 'black', '";
	gfile<<"tmp/al0.txt' using 1:2 title 'Neutral Pions' lt 1 pt 6 lc rgb 'black', '";	
	gfile<<filenameAp[first]<<"' using 1:2 smooth csplines lt 1 lc rgb 'black' lw 5 title 'NuWro (charged)', '";
	gfile<<filenameA0[first]<<"' using 1:2 smooth csplines lt 2 lc rgb 'black' lw 5 title 'NuWro (neutral)', '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filenameAp[k]<<"' using 1:2 smooth csplines lt 1 lc rgb '"<<colours[k]<<"' lw 5 title '"<<fzname[k]<<"', '";
			gfile<<filenameA0[k]<<"' using 1:2 smooth csplines lt 2 lc rgb '"<<colours[k]<<"' lw 5 notitle";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += filenameLp[k] + string(" and ") + filenameCp[k] + string(" and ") + filenameAp[k] + string(" and ") + filenameL0[k] + string(" and ") + filenameC0[k] + string(" and ") + filenameA0[k];
	plotlog(help);
	
	return 1;
}

int mplotNomad()
{
	get_date();
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}

	string prfilename[nof_fz];
	string pifilename[nof_fz];
	
	const int prbins = 8;
	const int pibins = 6;
	
	for (int k = 0; k < nof_fz; k++)
	{
		prfilename[k] = string("results/Nomad/nomad_protons_") + fzwork[k] + ("*.txt");
		prfilename[k] = find_last(prfilename[k]);

		pifilename[k] = string("results/Nomad/nomad_pions_") + fzwork[k] + ("*.txt");
		pifilename[k] = find_last(pifilename[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		
		if (fz_on[k] and (noFile(prfilename[k]) or noFile(pifilename[k])))
		{
			cout<<"Nomad ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"Nomad ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
		
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_nomad.eps";
	
	//gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 12 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;
	gfile<<"set size 1, 1"<<endl;
	gfile<<"set multiplot layout 2,2 title '{/=20 Average number of backward moving particles}'"<<endl;
	gfile<<"set xlabel 'Q^{2} [GeV]'"<<endl;
	gfile<<"set xrange [0:100]"<<endl;
	gfile<<"set key at 200, -0.15"<<endl;
		
	gfile<<"set title 'Protons' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel '<#Bp>'"<<endl;
	gfile<<"set yrange [0:0.3]"<<endl;	

	make_data(NOMADprotonsQ2, NOMADprotonsBack, NOMADprotonsQ2err, NOMADprotonsBackerr, prbins, 0);
	run ("cp tmp/data.txt tmp/data2.txt");
	make_data(NOMADpionsQ2, NOMADpionsBack, NOMADpionsQ2err, NOMADpionsBackerr, pibins, 0);
	
	gfile<<"plot 'tmp/data2.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data2.txt' using 1:2 title 'NOMAD' with points 7, '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<prfilename[k]<<"' using 1:2 smooth csplines lt 1 lc rgb '"<<colours[k]<<"' lw 5 title '"<<fzname[k]<<"'";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	ile = 0;
	
	gfile<<endl<<"unset label 1"<<endl;
	gfile<<endl<<"set nokey"<<endl;
	
	gfile<<"set title 'Pions' font 'Arial, 16'"<<endl;
	gfile<<"set ylabel '<#B{/Symbol p}^{-}>'"<<endl;
	gfile<<"set yrange [0:0.04]"<<endl;		
	
	gfile<<"plot 'tmp/data.txt' using 1:2:3:4 notitle with xyerrorbars lc rgb 'black' lw 2, '";
	gfile<<"tmp/data.txt' using 1:2 title 'NOMAD' with points 7, '";
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<pifilename[k]<<"' using 1:2 smooth csplines lt 1 lc rgb '"<<colours[k]<<"' lw 5 title '"<<fzname[k]<<"'";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}	
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += prfilename[k] + string(" and ") + pifilename[k];
	plotlog(help);
	
	return 1;
}

int mplotAtmNu()
{
	get_date();
	
	int nof = 0;
	int ile = 0;
	
	int first = 0;
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) nof++;
		
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k]) first = k;
		break;
	}
	
	string filename[nof_fz];
	
	const int bins = 8;
	string names[bins] = {"0pi", "pi+", "pi0", "pi-pi+", "2pi0", "pi+(>0*pi0)", "2pi+(>=0*pi0)", "pi-"};
	
	for (int k = 0; k < nof_fz; k++)
	{
		filename[k] = string("results/AtmNu/atmnu_") + fzwork[k] + ("*.txt");
		filename[k] = find_last(filename[k]);
	}
	
	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k] and noFile(filename[k]))
		{
			cout<<"Atmospheric Neutrinos ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			logfile<<"Atmospheric Neutrinos ("<<fzname[k]<<"): there are no txt files, make calculation first!"<<endl<<endl;
			return 1;
		}
	}
	
	ofstream gfile("tmp/gnuplot.gnu");
	
	string plotfile = "multiplots/multi_atmnu.eps";
	
	gfile<<setprecision (2);
	gfile<<"set terminal postscript eps enhanced color 'Arial' 16 "<<endl;
	gfile<<"set output '"<<plotfile<<"'"<<endl;
	gfile<<"set style data histogram"<<endl;
	gfile<<"set style histogram cluster gap 1"<<endl;
	gfile<<"set style fill solid noborder"<<endl;
	gfile<<"set xtic rotate by -45 scale 0"<<endl;
	gfile<<"set yrange[0:1]"<<endl;
	gfile<<"set title '{/=24 Atmospheric Neutrinos - neon} ("<<day<<" "<<fullm<<" "<<year<<")' font 'Arial, 20'"<<endl;
	gfile<<"set ylabel 'Number of events' font 'Arial, 16'"<<endl;
	gfile<<"set key spacing 1.5"<<endl;
	gfile<<"set bar 0"<<endl;

	gfile<<"set xtics ('' -1";
	for (int i = 0; i < bins; i++) gfile<<", '"<<names[i]<<"' "<<i;
	gfile<<")"<<endl;

	ofstream gdata("tmp/data.txt");
	for (int i = 0; i < bins; i++)
	{
		gdata<<AtmosphericDeuterium[i]<<" "<<AtmosphericNeon[i]<<endl;
	}
		
	gfile<<"plot 'tmp/data.txt' using 2 title 'Data' lc rgb 'black' fill solid 0.50, '";

	for (int k = 0; k < nof_fz; k++)
	{
		if (fz_on[k])
		{
			gfile<<filename[k]<<"' using 3 title '"<<fzname[k]<<"' lc rgb '"<<colours[k]<<"'";
			ile++;
			if (ile < nof) gfile<<", '";
		}
	}
	
	gfile.close();
	run("gnuplot tmp/gnuplot.gnu");
	cout<<endl<<plotfile<<" created."<<endl<<endl;
	
	string help = plotfile + string(" created from ");
	
	for (int k = 0; k < nof_fz; k++) if (fz_on[k]) help += filename[k] + string(" and ");
	plotlog(help);
	
	return 1;
}
