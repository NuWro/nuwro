#include "ganalysis.h"

histogram :: histogram (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int typnormy, double norma)
				: name (name_of_hist), data (data_file), title (plot_title), xlabel (x_name), ylabel (y_name), bins_x (nof_bins), norm_type (typnormy), normalization (norma), begin_x(begin)
{
	finalized = false;
	saved = false;
	ploted = false;
	
	width_x = (end - begin) / bins_x;

	for (int i = 0; i < nof_dyn; i++)
	{
		tnorm[i] = new int [bins_x];

		for (int j = 0; j < bins_x; j++)
			tnorm[i][j] = 0;
	}
}	
histogram :: ~histogram ()
{
	for (int i = 0; i < nof_dyn; i++)
		delete [] tnorm[i];
}

void histogram :: put (double x, int dyn)
{
	int a = (x - begin_x) / width_x;
	
	if (a >= 0 and a < bins_x)
		tnorm [dyn][a]++;
}

hist1D :: hist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int typnormy, double norma)
			: histogram (name_of_hist, data_file, plot_title, x_name, y_name, nof_bins, begin, end, typnormy, norma)			
{
	result = new double [bins_x];
	
	for (int i = 0; i < bins_x; i++)
		result[i] = 0;
	
	for (int i = 0; i < nof_dyn; i++)
	{
		part_result[i] = new double[bins_x];
		for (int j = 0; j < bins_x; j++)
			part_result[i][j] = 0;
	}			
}

hist1D :: ~hist1D ()
{
	finalize ();
	save ();
	plot ();
	
	for (int i = 0; i < nof_dyn; i++)
		delete [] part_result[i];
	
	delete [] result;
}

void hist1D :: finalize ()
{
	if (finalized)
		return;
	
	for (int i = 0; i < nof_dyn; i++)
		for (int j = 0; j < bins_x; j++)
			if (norm_type == 6 and tnorm[i][j] != 0)
				result[j] += part_result[i][j] / tnorm[i][j];
			else if (norm_type == 5)
				result[j] += part_result[i][j] / xsec[i];
			else if (norm[i] != 0)
				result[j] += part_result[i][j] / norm[i] / width_x;

	if (norm_type == 1)
	{
		double sum = 0;
		
		for (int i = 0; i < bins_x; i++)
			sum += result[i];
		
		for (int i = 0; i < bins_x; i++)
			result[i] *= normalization / sum;
	}
	
	finalized = true;
}

void hist1D :: save ()
{
	if (saved)
		return;
	
	string out = "analysis/" + name + ".txt";
	
	ofstream file(out.c_str());
	
	for (int i = 0; i < bins_x; i++)
		file << begin_x + width_x * i << " " << result[i] << endl << begin_x + width_x * (i + 1) << " " << result[i] << endl;
	
	file.close();
	
	saved = true;
}

void hist1D :: plot ()
{
	if (ploted)
		return;
	
	string in = "analysis/" + name + ".txt";
	string out = "analysis/" + name + ".eps";
	string gs = "analysis/" + name + ".gnu";
	
	ofstream gnu(gs.c_str());
	
	gnu << "set terminal postscript eps enhanced color 'ArialBold' 18" << endl;
	gnu << "set border lw 3" << endl;
	gnu << "set output '" << out << "'" << endl;
	gnu << "set title '" << title << "'" << endl;
	gnu << "set xlabel '" << xlabel << "'" << endl;
	gnu << "set ylabel '" << ylabel << "'" << endl;
	
	gnu << "plot '" << in << "' u 1:2 with lines lw 5 lt 1 title 'NuWro'";
		
	if (strcmp(data.c_str(), "") != 0)
		gnu << ", \\" << endl<< "'" << data << "' u 1:3 pt 7 lc rgb 'black' title 'data', \\" << endl << "'" << data << "' u 1:3:2:4 with xyerrorbars pt 7 lc rgb 'black' notitle";
			
	gnu.close();
	
	run_command ("gnuplot " + gs);
	
	ploted = true;
}

void hist1D :: put (double x, int dyn, double weight)
{
	int a = (x - begin_x) / width_x;
	
	if (a >= 0 and a < bins_x)
		part_result [dyn][a] += weight;
}

void hist1D :: put_over (double x, int dyn, double weight)
{
	int a = (x - begin_x) / width_x;
	
	if (a >= bins_x)
		a = bins_x - 1;
	
	if (a >= 0)
		part_result [dyn][a] += weight;
}

hist2D :: hist2D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, string z_name, int nof_bins_x, double begin_x, double end_x, int nof_bins_y, double begin_y, double end_y, int typnormy, double norma)
			: histogram (name_of_hist, data_file, plot_title, x_name, y_name, nof_bins_x, begin_x, end_x, typnormy, norma), zlabel (z_name), bins_y (nof_bins_y), begin_y(begin_y)
{
	max = 0;
	
	width_y = (end_y - begin_y) / bins_y;
	
	result = new double *[bins_x];
	
	for (int i = 0; i < bins_x; i++)
	{
		result [i] = new double [bins_y];
		
		for (int j = 0; j < bins_y; j++)
			result [i][j] = 0;
	}
	
	for (int i = 0; i < nof_dyn; i++)
	{
		part_result [i] = new double *[bins_x];
		
		for (int j = 0; j < bins_x; j++)
		{
			part_result [i][j] = new double [bins_y];
			
			for (int k = 0; k < bins_y; k++)
				part_result [i][j][k] = 0;
		}
	}
}

hist2D :: ~hist2D ()
{
	finalize ();
	save ();
	plot ();	

	for (int i = 0; i < bins_x; i++)
		delete [] result[i];
		
	delete [] result;
	
	for (int i = 0; i < nof_dyn; i++)
	{
		for (int j = 0; j < bins_x; j++)
			delete [] part_result [i][j];
			
		delete [] part_result [i];
	}
}

void hist2D :: finalize ()
{
	if (finalized)
		return;
	
	for (int i = 0; i < nof_dyn; i++)
		for (int j = 0; j < bins_x; j++)
			for (int k = 0; k < bins_y; k++)
				if (norm_type == 5 and xsec[i] != 0)
					result[j][k] += part_result[i][j][k] / xsec[i];
				else if (norm[i] != 0)
					result[j][k] += part_result[i][j][k] / norm[i] / width_x / width_y;

	if (norm_type == 1)
	{
		double sum = 0;
	
		for (int i = 0; i < bins_x; i++)
			for (int j = 0; j < bins_y; j++)
				sum += result[i][j];
				
		for (int i = 0; i < bins_x; i++)
			for (int j = 0; j < bins_y; j++)
				result[i][j] *= normalization / sum;
	}
	
	finalized = true;
}

void hist2D :: save ()
{
	if (saved)
		return;
	
	string out1 = "analysis/" + name + ".txt";
	string out2 = "analysis/" + name + "_box.txt";
	
	ofstream file1(out1.c_str());
	ofstream file2(out2.c_str());
		
	for (int i = 0; i < bins_x; i++)
	{	
		for (int j = 0; j < bins_y; j++)
		{
			if (max < result[i][j]) max = result[i][j];
			
			file1 << begin_x + width_x * i << " " << begin_y + width_y * j << " " << result[i][j] << endl;
			file1 << begin_x + width_x * i << " " << begin_y + width_y * (j + 1) << " " << result[i][j] << endl;
		}
		
		file1 << endl;
		
		for (int j = 0; j < bins_y; j++)
		{
			file1 << begin_x + width_x * (i + 1) << " " << begin_y + width_y * j << " " << result[i][j] << endl;
			file1 << begin_x + width_x * (i + 1) << " " << begin_y + width_y * (j + 1) << " " << result[i][j] << endl;
		}
		
		file1 << endl;
	}
	
	for (int i = 0; i < bins_x; i++)
	{
		for (int j = 0; j < bins_y; j++)
		{
			double x1 = begin_x + width_x * i;
			double x2 = begin_x + width_x * (i + 1);
			
			double y1 = begin_y + width_y * j;
			double y2 = begin_y + width_y * (j + 1);
			
			double factor = result[i][j] / max;
			
			double xshift = (x2 - x1) * (1 - factor) / 2.0;
			double yshift = (y2 - y1) * (1 - factor) / 2.0;
			
			x1 += xshift;
			x2 -= xshift;
			
			y1 += yshift;
			y2 -= yshift;
			
			file2 << x1 << " " << y1 << endl;
			file2 << x2 << " " << y1 << endl;
			file2 << x2 << " " << y2 << endl;
			file2 << x1 << " " << y2 << endl;
			file2 << x1 << " " << y1 << endl;
			file2 << endl;
		}
	}
		
	file1.close();
	file2.close();

	saved = true;
}

void hist2D :: plot ()
{
	if (ploted)
		return;
	
	string in1 = "analysis/" + name + ".txt";
	string in2 = "analysis/" + name + "_box.txt";
	
	string out1 = "analysis/" + name + ".eps";
	string out2 = "analysis/" + name + "_map.eps";
	string out3 = "analysis/" + name + "_box.eps";
	
	string gs = "analysis/" + name + ".gnu";
	
	ofstream gnu(gs.c_str());
	
	gnu << "set terminal postscript eps enhanced color 'ArialBold' 18" << endl;
	gnu << "set border lw 3" << endl;
	gnu << "set title '" << title << "'" << endl;
	gnu << "set xlabel '" << xlabel << "' offset -3, 0" << endl;
	gnu << "set ylabel '" << ylabel << "' offset -4, -1" << endl;

	gnu << "set output '" << out3 << "'" << endl;
	gnu << "plot '" << in2 << "' with lines lw 2 lt 1 lc rgb 'black' title 'NuWro (max = " << max << ")'";

	gnu << endl << endl;

	gnu << "set zlabel '" << zlabel << "' rotate left offset -0.5, 0" << endl;
	gnu << "set hidden3d" << endl;
	gnu << "set grid ztics" << endl;
	gnu << "set ticslevel 0" << endl;
	gnu << "set border 4095 lw 0.5" << endl;

	gnu << "set output '" << out1 << "'" << endl;
	gnu << "splot '" << in1 << "' with lines lw 2 lt 1 lc rgb 'black' title 'NuWro'";

	gnu << endl << endl;

	gnu << "set output '" << out2 << "'" << endl;
	gnu << "set pm3d map" << endl;
	gnu << "splot '" << in1 << "' title 'NuWro'";

	gnu.close();
	
	run_command ("gnuplot " + gs);	
	
	ploted = true;
}

void hist2D :: put (double x, double y, int dyn, double weight)
{
	int a = (x - begin_x) / width_x;
	int b = (y - begin_y) / width_y;
	
	if (a < bins_x and b < bins_y and a >= 0 and b >= 0)
		part_result [dyn][a][b] += weight;
}

mhist1D :: mhist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int nof_cases, int typnormy, double norma)
			: histogram (name_of_hist, data_file, plot_title, x_name, y_name, nof_bins, begin, end, typnormy, norma), cases (nof_cases)
{
	cnames = new string [cases];
	
	result = new double * [cases];
	
	for (int i = 0; i < cases; i++)
	{
		result [i] = new double [bins_x];
		
		for (int j = 0; j < bins_x; j++)
			result [i][j] = 0;
	}
	
	for (int i = 0; i < nof_dyn; i++)
	{
		part_result [i] = new double *[cases];
		
		for (int j = 0; j < cases; j++)
		{
			part_result [i][j] = new double [bins_x];
			
			for (int k = 0; k < bins_x; k++)
				part_result [i][j][k] = 0;
		}
	}
}

mhist1D :: ~mhist1D ()
{
	finalize ();
	save ();
	plot ();	

	for (int i = 0; i < cases; i++)
		delete [] result[i];
		
	delete [] result;
	
	for (int i = 0; i < nof_dyn; i++)
	{
		for (int j = 0; j < cases; j++)
			delete [] part_result [i][j];
			
		delete [] part_result [i];
	}

	delete [] cnames;
}

void mhist1D :: finalize ()
{
	if (finalized)
		return;
	
	for (int i = 0; i < nof_dyn; i++)
		for (int j = 0; j < cases; j++)
			for (int k = 0; k < bins_x; k++)
				if (norm_type == 6 and tnorm[i][k] != 0)
						result [j][k] += part_result [i][j][k] / tnorm[i][k];
				else if (norm_type == 5 and xsec[i] != 0)
						result [j][k] += part_result [i][j][k] / xsec[i];
				else if (norm[i] != 0)
						result [j][k] += part_result [i][j][k] / norm[i] / width_x;
						
	if (norm_type == 2)
	{
		double *sum = new double [cases];

		for (int i = 0; i < cases; i++)
		{
			sum [i] = 0;
			
			for (int j = 0; j < bins_x; j++)
				sum [i] += result [i][j];
		}
		
		for (int i = 0; i < cases; i++)
			for (int j = 0; j < bins_x; j++)
				result [i][j] *= normalization / sum [i];
				
		delete sum;
	}
	else if (norm_type == 3 or norm_type == 4)
	{
		double sum = 0;

		if (norm_type == 3)
			for (int i = 0; i < cases; i++)
				for (int j = 0; j < bins_x; j++)
					sum += result [i][j];
		else if (norm_type == 4)
			for (int j = 0; j < bins_x; j++)
				sum += result [0][j];
	
		for (int i = 0; i < cases; i++)
			for (int j = 0; j < bins_x; j++)
				result [i][j] *= normalization / sum;
	}
	
	finalized = true;
}

void mhist1D :: save ()
{
	if (saved)
		return;
	
	string out = "analysis/" + name + ".txt";
	
	ofstream file(out.c_str());
	
	for (int i = 0; i < bins_x; i++)
	{
		file << begin_x + width_x * i;
		
		for (int j = 0; j < cases; j++)
			file << " " << result [j][i];
			
		file << endl << begin_x + width_x * (i + 1);
		
		for (int j = 0; j < cases; j++)
			file << " " << result [j][i];
			
		file << endl;
	}
		
	file.close();
	
	saved = true;
}

void mhist1D :: plot ()
{
	if (ploted)
		return;
	
	string in = "analysis/" + name + ".txt";
	string out = "analysis/" + name + ".eps";
	string gs = "analysis/" + name + ".gnu";
	
	ofstream gnu(gs.c_str());
	
	gnu << "set terminal postscript eps enhanced color 'ArialBold' 18" << endl;
	gnu << "set border lw 3" << endl;
	gnu << "set output '" << out << "'" << endl;
	gnu << "set title '" << title << "'" << endl;
	gnu << "set xlabel '" << xlabel << "'" << endl;
	gnu << "set ylabel '" << ylabel << "'" << endl;
	
	gnu << "plot ";
	
	for (int i = 0; i < cases; i++)
	{
		gnu << "'" << in << "' u 1:" << i + 2 << " with lines lw 5 lt " << i + 1 << " title '" << cnames[i] << "'";
		
		if (i != cases - 1)
					gnu << ", \\" << endl;
	}
		
	if (strcmp(data.c_str(), "") != 0)
		gnu << ", \\" << endl<< "'" << data << "' u 1:3 pt 7 lc rgb 'black' title 'data', \\" << endl << "'" << data << "' u 1:3:2:4 with xyerrorbars pt 7 lc rgb 'black' notitle";
	
	gnu.close();
	
	run_command ("gnuplot " + gs);
	
	ploted = true;
}

void mhist1D :: put (double x, int dyn, double weight, int cas)
{
	int a = (x - begin_x) / width_x;
	
	if (a >= 0 and a < bins_x)
		part_result [dyn][cas][a] += weight;
}
void pattern :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
			
	N = new NuWro ();
	
	N -> set (P);
	
	if(P.dyn_dis_nc or P.dyn_res_nc  or P.dyn_dis_cc or P.dyn_res_cc)
		singlepion(P);
		
	run ();
	
	delete N;
}

void pattern :: run ()
{	
	for (int i = 0; i < events; i++)
	{
		event *e = new event();
		e->dyn = N->proces();
				
		N -> makeevent(e,P);
		N -> finishevent(e,P);
						
		if (e->weight >= 0)
		{
			norm[e->dyn]++;
			xsec[e->dyn] += e->weight;
			calculate (e);
		}
			
		delete e;
		
		cout << name << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << endl;
}

void mb_ncpi0 :: set_params ()
{
	P.read("data/beam/MBnumu.txt");
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void mb_ccpi0 :: set_params ()
{
	P.read("data/beam/MBccpi0.txt");
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
}

void mb_ccpip :: set_params ()
{
	P.read("data/beam/MBnumu.txt");
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
}

void mb_ncpi0 :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e->fof(pdg_piP) + 10 * e->fof(-pdg_piP) + e->fof(pdg_pi);

	if (pion == 1)
	{
		double mom = 0;
		double cos = 0;
		
		for (int i = 0; i < e->post.size(); i++)
			if (e->post[i].pdg == pdg_pi)
			{
				mom = e->post[i].momentum() / 1000.0;
				cos = e->post[i].p().z / mom / 1000.0;
				break;
			}
			
		h1 -> put (mom, e->dyn, e->weight);
		h2 -> put (cos, e->dyn, e->weight);
	}
}

void mb_ccpi0 :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e->fof(pdg_piP) + 10 * e->fof(-pdg_piP) + e->fof(pdg_pi);
	
	h1 -> histogram :: put (e -> in[0].E() / 1000.0, e->dyn);
	
	if (pion == 1)
	{
		double mom = 0;
		double cos = 0;
		
		for (int i = 0; i < e->post.size(); i++)
			if (e->post[i].pdg == pdg_pi)
			{
				mom = e->post[i].momentum() / 1000.0;
				cos = e->post[i].p().z / mom / 1000.0;
				break;
			}
			
		h3 -> put (mom, e->dyn, 14.08 * e->weight);
		h4 -> put (cos, e->dyn, 14.08 * e->weight);
			
		h1 -> put (e -> in[0].E() / 1000.0, e->dyn, 14.08 * e->weight);
		h2 -> put (-e -> q2() / 1000000.0, e->dyn, 14.08 * e->weight);
		h5 -> put (e -> out[0].Ek() / 1000.0, e->dyn, 14.08 * e->weight);
	}
}

void mb_ccpip :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e->fof(pdg_piP) + 10 * e->fof(-pdg_piP) + e->fof(pdg_pi);
	
	h1 -> histogram :: put (e -> in[0].E(), e->dyn);
	
	if (pion == 100)
	{
		h1 -> put (e -> in[0].E(), e->dyn, 14.08 * e->weight);
		h2 -> put (-e -> q2() / 1000000.0, e->dyn, 14.08 * e->weight);
		h4 -> put (e -> out[0].Ek(), e->dyn, 14.08 * e->weight);
		
		double ek = 0;
		
		for (int i = 0; i < e->post.size(); i++)
			if (e->post[i].pdg == pdg_piP)
			{
				ek = e->post[i].Ek();
				break;
			}
			
		h3 -> put (ek, e->dyn, 14.08 * e->weight);
	}
}

void A :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "2000";
	P.read("data/target/C.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void A :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e -> fof (pdg_piP) + 10 * e -> fof (-pdg_piP) + e -> fof (pdg_pi);
	
	double nu = e -> q0();
	double total = 0;
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_pi)
			total += e->post[i].E();
			
	
	if (e -> flag.nc)
		h4 -> put(total/nu, e -> dyn, e->weight);

	if (e -> flag.nc and e->fof(pdg_pi) > 0)
	{
		h4b -> put(total/nu, e -> dyn, e->weight);
		if (e -> flag.res)
			h4r -> put(total/nu, e -> dyn, e->weight);
		else if (e->flag.dis)
			h4d -> put(total/nu, e -> dyn, e->weight);
		else if (e->flag.coh)			
			h4c -> put(total/nu, e -> dyn, e->weight);
	}
	
	if (e -> flag.cc and pion == 0)
	{
		int Nproton = e -> fof (pdg_proton);
	
		double Emu = e -> out[0].E();
		double cos = e -> out[0].p().z / e -> out[0].momentum();
		
		double M = mass_proton;
		double m = mass_mu;
		
		double Erec = (2.0 * M * Emu - m * m) / (2.0 * (M - Emu + sqrt(Emu * Emu - m * m) * cos)); 
			
		if (Nproton < h1 -> bins_x)
			h1 -> part_result [e -> dyn][Nproton] += e -> weight;
		
		h2 -> put (Erec, e->dyn, e->weight);
	}
	
	if (e -> flag.nc and e -> fof (pdg_pi) > 0)
	{
		double energy = 0;
		double cos;
		
		for (int i = 0; i < e -> post.size(); i++)
			if (e -> post[i].pdg == pdg_pi and e -> post[i].E() > energy)
			{
				energy = e -> post[i].E();
				cos = e -> post[i].p().z / e -> post[i].momentum();
			}
		
		h3 -> put (energy, cos, e->dyn, e->weight * 12.0);
		h3b -> put (energy, cos, e->dyn, e->weight);
	}
	
}

void antiA :: set_params ()
{
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "2000";
	P.read("data/target/C.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void B :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "5000";
	P.read("data/target/Fe.txt");
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void B2 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "3000";
	P.read("data/target/Fe.txt");
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void antiB2 :: set_params ()
{
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "3000";
	P.read("data/target/Fe.txt");
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;
}

void B :: calculate (event *e)
{
	using namespace PDG;

	double nu = e -> q0 ();
	
	if (nu > 0 and e -> fof(pdg_pi) > 0)
	{
		
		for (int i = 0; i < e -> post.size(); i++)
			if (e -> post[i].pdg == pdg_pi)
			{
				double x = e -> post[i].E() / nu;
				
				h1 -> put (x, e->dyn, e->weight);
			}
	}
}

void B2 :: calculate (event *e)
{
	using namespace PDG;
	
	double total = 0;

	for (int i = 0; i < e -> post.size(); i++)
	{
		int pdg = e -> post[i].pdg;
		double energy = e -> post[i].E();
		
		if (pdg == pdg_pi)
			energy *= 1.3;
		else if (pdg == pdg_proton)
		{
			energy = e -> post[i].Ek();
			
			if (energy < 150) energy = 0;
		}
		else if (pdg == pdg_neutron)
		{
			energy = e-> post[i].Ek();
			
			if (energy < 300) energy = 0;
			else energy *= 0.5;
		}
		
		total += energy;
	}
				
	h1 -> put (total, e->dyn, e->weight);
}

void antiB :: set_params ()
{
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "5000";
	P.read("data/target/Fe.txt");
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;	
}

void C1 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "0 10000";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void C1 :: calculate (event *e)
{
	using namespace PDG;
	
	double energy = e -> in[0].E() / 1000.0;
	
	h1 -> histogram :: put (energy, e->dyn);
		
	int pion = 100 * e -> fof (pdg_piP) + 10 * e -> fof (-pdg_piP) + e -> fof (pdg_pi);
	
	int a = energy / h1 -> width_x;
		
	if (a < h1 -> bins_x)
	{
		h1 -> part_result [e->dyn][0][a] += e -> weight * 40.0;
				
		if (pion == 0)
			h1 -> part_result [e->dyn][1][a] += e -> weight * 40.0;
	}
}

void antiC1 :: set_params ()
{
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "0 10000";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void C2 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "2500";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void C2 :: calculate (event *e)
{
	using namespace PDG;
	
	int N1 = e -> fof (pdg_proton);
	
	if (N1 < h1 -> bins_x)
		h1 -> part_result [e -> dyn][0][N1] += e -> weight;
	
	int N2 = 0;
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_proton and e -> post[i].Ek() > 50)
			N2++;
			
	if (N2 < h1 -> bins_x)
		h1 -> part_result [e->dyn][1][N2] += e -> weight;	
}

void antiC2 :: set_params ()
{
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "2500";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void D :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "3000";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void D :: calculate (event *e)
{
	using namespace PDG;
	
	double totalE = 0;
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg != pdg_neutron)
		{
			if (e -> post[i].pdg == pdg_proton)
				totalE += e -> post[i].Ek()/1000.0;
			else
				totalE += e -> post[i].E()/1000.0;
		}
			
	h4 -> put(totalE, e->dyn, e->weight);
				
	int pion = 100 * e -> fof (pdg_piP) + 10 * e -> fof (-pdg_piP) + e -> fof (pdg_pi);
	
	if (pion == 0)
	{
		int N1 = e -> fof (pdg_proton);
		int N2 = 0;
		
		for (int i = 0; i < e -> post.size(); i++)
			if (e -> post[i].pdg == pdg_proton and e -> post[i].Ek() > 50)
				N2++;
				
		if (N1 < h1 -> bins_x and N1 > 0)
			h1 -> part_result [e->dyn][0][N1] += e -> weight;

		if (N2 < h1 -> bins_x)
			h1 -> part_result [e->dyn][1][N2] += e -> weight;
		
		if (N2 > 0)
		{
			double mom = 0;
			double cos;

			for (int i = 0; i < e -> post.size(); i++)
				if (e -> post[i].pdg == pdg_proton and e -> post[i].Ek() > 50 and e -> post[i].momentum() > mom)
				{
					mom = e -> post[i].momentum();
					cos = e -> post[i].p().z / mom;
				}
				
			if (mom != 0)
			{
				int a = mom / h2 -> width_x;
				int b = cos / h3 -> width_x + 10;
				
				if (N2 > 4) N2 = 4;
				
				if (a < h2 -> bins_x)
					h2 -> part_result [e->dyn][N2-1][a] += e -> weight;
					
				if (b < h3 -> bins_x)
					h3 -> part_result [e->dyn][N2-1][b] += e -> weight;
			}
		}	
		if (N1 > 0)
		{
			double mom = 0;
			double cos;

			for (int i = 0; i < e -> post.size(); i++)
				if (e -> post[i].pdg == pdg_proton and e -> post[i].momentum() > mom)
				{
					mom = e -> post[i].momentum();
					cos = e -> post[i].p().z / mom;
				}
				
			if (mom != 0)
			{
				int a = mom / h2b -> width_x;
				
				if (N1 > 4) N1 = 4;
				
				if (a < h2 -> bins_x)
					h2b -> part_result [e->dyn][N1-1][a] += e -> weight;
			}
		}	
	}
}

void E :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "1000";
	P.read("data/target/Ar.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void F :: set_params ()
{
	P.read("data/beam/nuintF.txt");
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 0;	
}

void F :: calculate (event *e)
{
	double energy = e -> in[0].E() / 1000.0;

	int a = energy / h1 -> width_x;
	
	if (e -> in[0].pdg == PDG::pdg_nu_mu)
	{
		h1 -> histogram :: put (energy, e -> dyn);
		
		if (a < h1 -> bins_x)
		{
			if (e -> flag.cc)
			{
				h1 -> part_result [e -> dyn][0][a] += e -> weight;
				
				if (!e -> flag.mec)			
					h1 -> part_result [e -> dyn][2][a] += e -> weight;
			}
			else
			{
				h1 -> part_result [e -> dyn][1][a] += e -> weight;
				
				if (!e -> flag.mec)			
					h1 -> part_result [e -> dyn][3][a] += e -> weight;				
			}
		}
	}
	else
	{
		h2 -> histogram :: put (energy, e -> dyn);
		
		if (a < h2 -> bins_x)
		{
			if (e -> flag.cc)
			{
				h2 -> part_result [e -> dyn][0][a] += e -> weight;
				
				if (!e -> flag.mec)			
					h2 -> part_result [e -> dyn][2][a] += e -> weight;
			}
			else
			{
				h2 -> part_result [e -> dyn][1][a] += e -> weight;
				
				if (!e -> flag.mec)			
					h2 -> part_result [e -> dyn][3][a] += e -> weight;				
			}
		}
	}
}

F :: ~F ()
{
	h1 -> finalize();
	h2 -> finalize();
				
	for (int j = 0; j < h1 -> cases; j++)
		for (int k = 0; k < h1 -> bins_x; k++)
			if (h2 -> result[j][k] != 0)
				h3 -> result [j][k] = h1 -> result[j][k] / h2 -> result[j][k];
	
	h3 -> finalized = true;
				
	delete h1;
	delete h2;
	delete h3;
}

void G1 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "600";
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G1 :: calculate (event *e)
{
	double ek = e -> out[0].Ek();
	double cos = e -> out[0].p().z / e -> out[0].momentum ();
	
	h2 -> put(ek, cos, e->dyn, e->weight);
	
	if (e -> flag.qel)
		h1 -> put(ek, cos, e->dyn, e->weight);
}

void G1_2 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G1_3 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "1200";
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G1_4 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "600";
	P.read("data/target/h2o.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G1_5 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/h2o.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G1_6 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "1200";
	P.read("data/target/h2o.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G2 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G2 :: calculate (event *e)
{
	double sum = 0;
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == PDG::pdg_proton)
			sum += e -> post[i].Ek();
			
	int a = sum / h1 -> width_x;
	
	if (a < h1 -> bins_x)
	{
		h1 -> part_result [e->dyn][0][a] += e -> weight;
		
		if (!e -> flag.mec)
			h1 -> part_result [e->dyn][1][a] += e -> weight;
	}
}

void G2_2 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/CH2.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G2_3 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/h2o.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void G2_4 :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_energy = "800";
	P.read("data/target/h2o.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0; //should be 1

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void H :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_energy = "5000";
	P.read("data/target/C.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void antiH :: set_params ()
{
	P.beam_type = 0;
	P.beam_particle = -PDG::pdg_nu_mu;
	P.beam_energy = "5000";
	P.read("data/target/C.txt");
	
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void H :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e -> fof(pdg_piP) + 10 * e -> fof(-pdg_piP) + e -> fof(pdg_pi);
	int nof_pi = e -> fof(pdg_piP) + e -> fof(-pdg_piP) + e -> fof(pdg_pi);
	
	double transfer = e -> q0 ()/ 1000.0;
	
	double lep = 0; //lead energy proton
	double len = 0;

	double sump = 0; //sum of protons energy
	double sumn = 0;
	
	double sumpi0 = 0;
	double sumpip = 0;
	double sumpim = 0;

	double lepi0 = 0;
	double lepip = 0;
	double lepim = 0;
	double lcpi0 = 0;
	double lcpip = 0;
	double lcpim = 0;
		
	for (int i = 0; i < e -> post.size(); i++)
	{
		int pdg = e -> post[i].pdg;
		double E = e -> post[i].E() / 1000.0;
		double Ek = e -> post[i].Ek() / 1000.0;
		double cos = e -> post[i].p().z / e -> post[i].momentum();
		
		if (pdg == pdg_proton)
		{
			sump += Ek;
			
			if (Ek > lep) lep = Ek;
		}
		else if (pdg == pdg_neutron)
		{
			sumn += Ek;
			
			if (Ek > len) len = Ek;
		}
		else if (pdg == pdg_pi)
		{
			sumpi0 += E;
			
			if (E > lepi0)
			{
				lepi0 = Ek;
				lcpi0 = cos;
			}
		}
		else if (pdg == pdg_piP)
		{
			sumpip += E;
			
			if (E > lepip)
			{
				lepip = Ek;
				lcpip = cos;
			}
		}
		else if (pdg == -pdg_piP)
		{
			sumpim += E;
			
			if (E > lepim)
			{
				lepim = Ek;
				lcpim = cos;
			}
		}
	}
	
	if (pion == 0)
	{
		h1 -> put (lep, e->dyn, e->weight);
		h2 -> put (sump, e->dyn, e->weight);
	}
	else if (pion == 100)
		h3 -> put (lepip, lcpip, e->dyn, 12.0 * e->weight);
	else if (pion == 10)
		h4 -> put (lepim, lcpim, e->dyn, 12.0 * e->weight);
	else if (pion == 1)
		h5 -> put (lepi0, lcpi0, e->dyn, 12.0 * e->weight);
	else if (nof_pi >= 2)
	{
		if (lepip != 0)
		{
			h6 -> put (lepip, lcpip, e->dyn, e->weight);
			h9 -> put (sumpip, e->dyn, e->weight);
		}

		if (lepim != 0)
		{
			h7 -> put (lepim, lcpim, e->dyn, e->weight);
			h10 -> put (sumpim, e->dyn, e->weight);
		}

		if (lepi0 != 0)
		{
			h8 -> put (lepi0, lcpi0, e->dyn, e->weight);
			h11 -> put (sumpi0, e->dyn, e->weight);
		}
	}

	if (e -> fof (pdg_pi) > 0) h13 -> put (transfer, sumpi0 / transfer, e->dyn, e->weight);
	if (e -> fof (pdg_piP) > 0)h19 -> put (transfer, sumpip / transfer, e->dyn, e->weight);
	if (e -> fof (-pdg_piP) > 0)h20 -> put (transfer, sumpim / transfer, e->dyn, e->weight);

	h14 -> put (transfer, e -> fof (pdg_pi), e->dyn, e->weight);
	h15 -> put (transfer, e -> fof (pdg_piP), e->dyn, e->weight);
	h16 -> put (transfer, e -> fof (-pdg_piP), e->dyn, e->weight);
	
	if (true)//e -> fof (pdg_neutron) > 0) 
	{
		h17 -> put (transfer, sumn / transfer, e->dyn, e->weight);
		if (nof_pi == 0) h17a -> put (transfer, sumn / transfer, e->dyn, e->weight);
		else if (nof_pi == 1) h17b -> put (transfer, sumn / transfer, e->dyn, e->weight);
		else if (nof_pi == 2) h17c -> put (transfer, sumn / transfer, e->dyn, e->weight);
		else if (nof_pi > 2) h17d -> put (transfer, sumn / transfer, e->dyn, e->weight);
	}
	if (true)//e -> fof (pdg_proton) > 0)
	{
		 h18 -> put (transfer, sump / transfer, e->dyn, e->weight);
	
		if (nof_pi == 0) h18a -> put (transfer, sump / transfer, e->dyn, e->weight);
		else if (nof_pi == 1) h18b -> put (transfer, sump / transfer, e->dyn, e->weight);
		else if (nof_pi == 2) h18c -> put (transfer, sump / transfer, e->dyn, e->weight);
		else if (nof_pi > 2) h18d -> put (transfer, sump / transfer, e->dyn, e->weight);
	}	
}

void test_mec_old_cc :: set_params ()
{
	P.read("data/beam/test_mec.txt");
	P.read("data/target/C.txt");
	
	P.mec_kind = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void test_mec_new_cc :: set_params ()
{
	P.read("data/beam/test_mec.txt");
	P.read("data/target/C.txt");
	
	P.mec_kind = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void test_mec_old_cc :: calculate (event *e)
{
	double energy = e -> in[0].E();

	
	if (e -> in[0].pdg == PDG::pdg_nu_mu)
	{
		h1 -> histogram :: put (energy, e -> dyn);
		h1 -> put (energy, e -> dyn, e->weight);
	}
	else
	{
		h2 -> histogram :: put (energy, e -> dyn);
		h2 -> put (energy, e -> dyn, e->weight);
	}
}

void test2_mec_old_cc :: set_params ()
{
	P.read("data/target/C.txt");
	
	P.beam_type = 0;
	P.beam_particle = 14;
	P.beam_energy = "1000";
	
	P.kaskada_on = 0;
	
	P.mec_kind = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void test2_mec_new_cc :: set_params ()
{
	P.read("data/target/C.txt");
	
	P.beam_type = 0;
	P.beam_particle = 14;
	P.beam_energy = "1000";
	
	P.kaskada_on = 0;
	
	P.mec_kind = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void test2_mec_old_cc :: calculate (event *e)
{
	using namespace PDG;
	
	int proton  = e -> nof (pdg_proton);
	int which = 2;
	
	if (proton == 2)
		which = 0;
	else if (proton == 0)
		which = 1;
		
	double mom = 0, cos;
	
	for (int i = 0; i < e -> out.size(); i++)
		if (e -> out[i].pdg == pdg_proton or e -> out[i].pdg == pdg_neutron)
			if (e -> out[i].momentum () > mom)
			{
				mom = e -> out[i].momentum ();
				cos = e -> out[i].p().z / mom;
			}
			
	h1 -> put (mom, e -> dyn, e -> weight, which);
	h2 -> put (cos, e -> dyn, e -> weight, which);
}

void test_mec_nc :: set_params ()
{
	P.read("data/target/C.txt");
	P.read("data/beam/nuintF.txt");
	
	P.kaskada_on = 0;
	
	P.mec_kind = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 1;	
}

void test_mec_nc :: calculate (event *e)
{
	double energy = e -> in[0].E() / 1000.0;

	int x = 2;
	
	if (e -> flag.qel)
		x = 1;
	
	if (e -> in[0].pdg == PDG::pdg_nu_mu)
	{
		h1 -> histogram :: put (energy, e -> dyn);
		h1 -> put (energy, e -> dyn, e->weight, 0);
		h1 -> put (energy, e -> dyn, e->weight, x);
	}
	else
	{
		h2 -> histogram :: put (energy, e -> dyn);
		h2 -> put (energy, e -> dyn, e->weight, 0);
		h2 -> put (energy, e -> dyn, e->weight, x);
	}
}

mb_nce :: mb_nce () : pattern ("MB NCE", 10000000)
{	
	run_command ("mkdir -p analysis/mb_nce/");

	nof_Ma = 91;
	nof_ds = 1;
			
	Ma = new int [nof_Ma];
	ds = new double [nof_ds];
	
	for (int i = 0; i < nof_Ma; i++)
		Ma [i] = 700 + i*10;
		
	for (int i  = 0; i < nof_ds; i++)
		ds [i] = 0;//-0.5 + i*0.1;
	
	chi_2d = new hist2D ("mb_nce/chi2_2d", "", "", "Axial Mass [MeV]", "{/Symbol D}s", "{/Symbol x}^{2}", nof_Ma, Ma[0], Ma[nof_Ma - 1], nof_ds, ds[0], ds[nof_ds - 1], 0);
	chi_2d -> finalized = true;

	chi_2d_mec = new hist2D ("mb_nce/chi2_2d_mec", "", "", "Axial Mass [MeV]", "{/Symbol D}s", "{/Symbol x}^{2}", nof_Ma, Ma[0], Ma[nof_Ma - 1], nof_ds, ds[0], ds[nof_ds - 1], 0);
	chi_2d_mec -> finalized = true;
	
	chi = new mhist1D * [nof_ds];
	
	for (int i = 0; i < nof_ds; i++)
	{
		stringstream temp1;
		string des;
	
		temp1 << ds[i];
		temp1 >> des;
			
		string name = "mb_nce/chi2_" + des;

		chi[i] = new mhist1D (name, "", "", "Axial Mass [MeV]", "{/Symbol x}^{2}", nof_Ma, Ma[0], Ma[nof_Ma - 1], 2, 0);
		chi[i] -> finalized = true;
		chi[i] -> cnames[0] = "without MEC";
		chi[i] -> cnames[1] = "with MEC";
	}
	
	for (int i = 0; i < 51; i++)
	{
		not_qel [i] = 0;
		mec     [i] = 0;
	}
}

mb_nce :: ~mb_nce ()
{	
	delete [] Ma;
	delete [] ds;
				
	for (int i = 0; i < nof_ds; i++)
		delete chi[i];
		
	delete [] chi;
	
	delete chi_2d;
	delete chi_2d_mec;
}

void mb_nce :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
	
	N = new NuWro ();

	N -> set(P);

	singlepion(P);

	cout << "QEL on Hydrogen" << endl;
	
	run ();

	cout << "RES on Carbon" << endl;

	P.read("data/target/C.txt");
	
	N -> refresh_target (P);
	
	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 1;
			
	N -> refresh_dyn (P);

	run ();

	cout << "DIS on Carbon" << endl;
	
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 1;
	
	N -> refresh_dyn (P);
	
	run ();
			
	cout << "MEC on Carbon" << endl;
		
	P.dyn_dis_nc = 0;
	P.dyn_mec_nc = 1;	

	N -> refresh_dyn (P);

	run ();
	
	//QEL on Carbon

	P.dyn_qel_nc = 1;
	P.dyn_mec_nc = 0;
	
	N -> refresh_dyn (P);	

	for (int i = 0; i < nof_Ma; i++)
		for (int j = 0; j < nof_ds; j++)
		{
			P.qel_nc_axial_mass = Ma[i];
			P.delta_s = ds[j];

			ff_configure (P);
			
			stringstream temp1, temp2;
			string mas, des;
	
			temp1 << ds[j];
			temp1 >> des;
			
			temp2 << Ma[i];
			temp2 >> mas;
			
			string name = "mb_nce/" + mas + "_" + des;
			string title = "Axial mass = " + mas + " MeV, {/Symbol D}s = " + des;

			qel = new mhist1D (name, "data/data/mb_ncel.txt", title, "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 2, 0);
			qel -> cnames[0] = "without MEC";
			qel -> cnames[1] = "with MEC";
			qel -> finalized = true;
			
			cout << "QEL on Carbon, Axial mass = " << mas << " MeV, Delta s = " << des << endl;
			
			run ();
			
			for (int z = 0; z < 51; z++)
			{
				qel -> result [0][z] += not_qel [z] + mb_nce_bg [z];
				qel -> result [1][z] += not_qel [z] + mec [z] + mb_nce_bg [z];
			}
			
			chi2 (i, j);
			
			delete qel;
		}
				
	delete N;
}

void mb_nce :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/H.txt");
	
	P.qel_nc_axial_mass = 1030;
	
	P.mec_kind = 0;
		
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void mb_nce :: calculate (event *e)
{		
	using namespace PDG;
	
	const double factor = 2.95e45 / 14.0 / events / 18.0;

	if (e -> flag.mec)
	{		
		double total_sn = 0;
		double total_sp1 = 0;
		double total_sp2 = 0;
		double total_m = 0;
		
		for (int i = 0; i < e -> post.size(); i++)
		{
			if (e -> post[i].pdg == pdg_proton)
			{
				if (e -> post[i].primary)
				{
					if (total_sp1 == 0)
						total_sp1 += e -> post[i].Ek();
					else
						total_sp2 += e -> post[i].Ek();
				}
				else
					total_m += e -> post[i].Ek();
			}
			else if (e -> post[i].pdg == pdg_neutron)
			{
				if (e -> post[i].primary)
					total_sn += e -> post[i].Ek();
				else
					total_m += e -> post[i].Ek();
			}
		}
		
		int scenario;		
		
		if (e -> in[1].pdg == pdg_neutron and e -> in[2].pdg == pdg_neutron)
			scenario = 0;
		else if (e -> in[1].pdg == pdg_proton and e -> in[2].pdg == pdg_proton)
		{
			if (total_sp1 != 0)
			{
				if (total_sp2 != 0)
					scenario = 4;
				else 
					scenario = 5;
			}
			else
				scenario = 6;
		}
		else
		{
			if (total_sp1 != 0)
				scenario = 1;
			else if (total_sn != 0)
				scenario = 2;
			else
				scenario = 3;
		}
		
		int bin = 0;
		int w;
		double waga = 1;
		
		switch (scenario)
		{
			case 0:
				bin = t2r (total_sn + total_m, 3, w);
				if (bin >= 0) waga = wagi [3][w];
				break;
			case 1:
				bin = t2r_mec (total_sp1, 1, total_m + total_sn, 3);
				break;
			case 2:
				bin = t2r_mec (total_m, 2, total_sn, 3);
				break;
			case 3:
				bin = t2r_mec (total_m / 2.0, 2, total_m / 2.0, 3);
				break;
			case 4:
				bin = t2r_mec (total_sp1, 1, total_sp2, 1);
				break;
			case 5:
				bin = t2r_mec (total_sp1, 1, total_m, 2);
				break;
			case 6:
				bin = t2r (total_m, 2, w);
				if (bin >= 0) waga = wagi [2][w];
				break;
			default: break;
		}
		
		if (bin > 50) bin = 50;
		
		if (bin >= 0)
			mec [bin] += e -> weight * factor * 12.0 * waga;
	}
	else
	{
		double total = 0;
			
		for (int i = 0; i < e -> post.size(); i++)
			if (e -> post[i].pdg == pdg_proton or e -> post[i].pdg == pdg_neutron)
				total += e -> post[i].Ek();
				
		if (P.nucleus_p == 1)
		{
			int w;
			int bin = t2r (total, 0, w);
			if (bin >= 0)
				not_qel [bin] += e -> weight * factor * 2.0 * wagi [0][w];
		}
		else if ( (e -> flag.res or e -> flag.dis) and (e -> fof (pdg_pi) + e -> fof (pdg_piP) + e -> fof (-pdg_piP) == 0))
		{
			int w;
			int bin = t2r (total, 4, w);
			if (bin >= 0)
				not_qel [bin] += e -> weight * factor * 12.0 * wagi [4][w];
		}
		else if (e -> flag.qel)
		{
			int scenario = 1;
			
			if (e -> in[1].pdg == pdg_neutron)
				scenario = 3;
			else if (e -> number_of_interactions ())
				scenario = 2;
			
			int w;
			int bin = t2r (total, scenario, w);
			if (bin >= 0)
			{
				qel -> result [0][bin] += e -> weight * factor * 12.0 * wagi [scenario][w];			
				qel -> result [1][bin] += e -> weight * factor * 12.0 * wagi [scenario][w];			
			}
		}	
	}
}

int mb_nce :: t2r_mec (double x1, int s1, double x2, int s2)
{
	int w1, w2;
	
	int bin1 = t2r (x1, s1, w1);
	int bin2 = t2r (x2, s2, w2);
	
	double waga1 = 0;
	double waga2 = 0;
	
	if (bin1 >= 0)
		waga1 = wagi [s1][w1];
	
	if (bin2 >= 0)
		waga2 = wagi [s2][w2];
		
	bool a1 = false;
	bool a2 = false;
	
	double los1 = frandom();
	double los2 = frandom();
	
	if (los1 < waga1)
		a1 = true;
	
	if (los2 < waga2)
		a2 = true;
		
	if (a1 and a2)
		return bin1 + bin2;
	else if (a1)
		return bin1;
	else if (a2)
		return bin2;
	else
		return -1;	
}

int mb_nce :: t2r (double x, int s, int &w)
{
	if (x <= 0)
		return -1;
	
	w = x / 18;
	
	if (w > 50)	w = 50;
	
	double los = frandom00 ();
		
	return find_bin (s, w, los);
}

int mb_nce :: find_bin (int s, int col, double x)
{
	int a = 0;
	int b = 50;
	int c = (b - a) / 2 + a;
			
	if (dystrR [s][0][col] == 1)
		return -1;
	
	if (x <= dystrR [s][0][col])
		return 0;
	
	do
	{
		if (x > dystrR [s][a][col] and x <= dystrR [s][c][col])
			b = c;
		else
			a = c;
			
		c = (b - a) / 2 + a;
	}
	while (b - a > 1);

	return b;
}		

void mb_nce :: chi2 (int a, int b)
{
	double c0 = calc_chi (qel -> result [0]);
	double c1 = calc_chi (qel -> result [1]);

	chi [b] -> result [0][a] = c0;
	chi [b] -> result [1][a] = c1;
				
	chi_2d -> result [a][b] = c0;
	chi_2d_mec -> result [a][b] = c1;
}

double mb_nce :: calc_chi (double *x)
{
	double res = 0;
	
	double dif [51];
		
	for (int i = 0; i < 51; i++)
		dif[i] = mb_nce_data[i] - x[i];
				
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			res += dif[i] * Mrev [i][j] * dif[j];
	
	return res;	
}

pattern * choose (int x)
{
	pattern *wsk;
	
	switch (x)
	{
		case  0: wsk = new mb_ncpi0; break;
		case  1: wsk = new mb_ccpi0; break;
		case  2: wsk = new mb_ccpip; break;
		case  3: wsk = new A; break;
		case  4: wsk = new antiA; break;
		case  5: wsk = new B; break;
		case  6: wsk = new antiB; break;
		case  7: wsk = new C1; break;
		case  8: wsk = new antiC1; break;
		case  9: wsk = new C2; break;
		case 10: wsk = new antiC2; break;
		case 11: wsk = new D; break;
		case 12: wsk = new E; break;
		case 13: wsk = new F; break;
		case 14: wsk = new G1; break;
		case 15: wsk = new G1_2; break;
		case 16: wsk = new G1_3; break;
		case 17: wsk = new G1_4; break;
		case 18: wsk = new G1_5; break;
		case 19: wsk = new G1_6; break;
		case 20: wsk = new G2; break;
		case 21: wsk = new G2_2; break;
		case 22: wsk = new G2_3; break;
		case 23: wsk = new G2_4; break;
		case 24: wsk = new H; break;
		case 25: wsk = new antiH; break;
		case 26: wsk = new B2; break;
		case 27: wsk = new antiB2; break;
		case 28: wsk = new test_mec_old_cc; break;
		case 29: wsk = new test_mec_new_cc; break;
		case 30: wsk = new test2_mec_old_cc; break;
		case 31: wsk = new test2_mec_new_cc; break;
		case 32: wsk = new test_mec_nc; break;
		case 33: wsk = new mb_nce; break;
		default: wsk = NULL;
	}
	
	return wsk;
}

int main (int argc, char **argv)
{
	set_dirs(argv[0]);
	init_genrand(time(NULL));
	
	run_command ("mkdir -p analysis/");
			
	for (int i = 0; i < nof_class; i++)
	{
		if (!active_class[i])
			continue;
			
		pattern *wsk = choose(i);
		
		if (wsk) 
		{
			wsk -> start();
			delete wsk;
		}
	}
			
	return 0;
}

