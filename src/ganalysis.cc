#include "ganalysis.h"

hist1D :: hist1D (string name_of_hist, string data_file, string x_name, string y_name, int nof_bins, double begin, double end, double norma)
			: histogram (name_of_hist, data_file, x_name, y_name, nof_bins, begin, end, norma)			
{
	result = new double[bins_x];
	
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
	for (int i = 0; i < nof_dyn; i++)
		delete [] part_result[i];
	
	delete [] result;
}

void hist1D :: finalize (double *norm)
{
	for (int i = 0; i < nof_dyn; i++)
		if (norm[i] != 0)
			for (int j = 0; j < bins_x; j++)
				result[j] += part_result[i][j] / norm[i] / width_x;

	if (normalization != 0)
	{
		double sum = 0;
		
		for (int i = 0; i < bins_x; i++)
			sum += result[i];
		
		for (int i = 0; i < bins_x; i++)
			result[i] *= normalization / sum;
	}
}

void hist1D :: save ()
{
	string out = "analysis/" + name + ".txt";
	
	ofstream file(out.c_str());
	
	for (int i = 0; i < bins_x; i++)
		file << width_x * i << " " << result[i] << endl << width_x * (i + 1) << " " << result[i] << endl;
	
	file.close();
}

void hist1D :: plot ()
{
	string in = "analysis/" + name + ".txt";
	string out = "analysis/" + name + ".eps";
	string gs = "analysis/" + name + ".gnu";
	
	ofstream gnu(gs.c_str());
	
	gnu << "set terminal postscript eps enhanced color 'Arial' 18" << endl;
	gnu << "set output '" << out << "'" << endl;
	gnu << "set xlabel '" << xlabel << "'" << endl;
	gnu << "set ylabel '" << ylabel << "'" << endl;
	
	gnu << "plot '" << in << "' u 1:2 with lines lw 5 lt 1 lc rgb 'black' title 'NuWro'";
		
	if (strcmp(data.c_str(), "") != 0)
		gnu << ", \\" << endl<< "'" << data << "' u 1:3 pt 7 lc rgb 'black' title 'data', \\" << endl << "'" << data << "' u 1:3:2:4 with xyerrorbars pt 7 lc rgb 'black' notitle";
		
	gnu.close();
	
	run_command ("gnuplot " + gs);
}

hist2D :: hist2D (string name_of_hist, string data_file, string x_name, string y_name, string z_name, int nof_bins_x, double begin_x, double end_x, int nof_bins_y, double begin_y, double end_y, double norma)
			: histogram (name_of_hist, data_file, x_name, y_name, nof_bins_x, begin_x, end_x, norma), zlabel (z_name), bins_y (nof_bins_y)
{
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

void hist2D :: finalize (double *norm)
{
	for (int i = 0; i < nof_dyn; i++)
		if (norm[i] != 0)
			for (int j = 0; j < bins_x; j++)
				for (int k = 0; k < bins_y; k++)
					result[j][k] += part_result[i][j][k] / norm[i] / width_x / width_y;

	if (normalization != 0)
	{
		double sum = 0;
	
		for (int i = 0; i < bins_x; i++)
			for (int j = 0; j < bins_y; j++)
				sum += result[i][j];
				
		for (int i = 0; i < bins_x; i++)
			for (int j = 0; j < bins_y; j++)
				result[i][j] *= normalization / sum;
	}
}

void hist2D :: save ()
{
	string out = "analysis/" + name + ".txt";
	
	ofstream file(out.c_str());
	
	for (int i = 0; i < bins_x; i++)
	{
		for (int j = 0; j < bins_y; j++)
		{
			file << width_x * i << " " << width_y * j << " " << result[i][j] << endl;
			file << width_x * i << " " << width_y * (j + 1) << " " << result[i][j] << endl;
		}
		
		file << endl;
		
		for (int j = 0; j < bins_y; j++)
		{
			file << width_x * (i + 1) << " " << width_y * j << " " << result[i][j] << endl;
			file << width_x * (i + 1) << " " << width_y * (j + 1) << " " << result[i][j] << endl;
		}
		
		file << endl;
	}
		
	file.close();
}

void hist2D :: plot ()
{
	string in = "analysis/" + name + ".txt";
	string out = "analysis/" + name + ".eps";
	string gs = "analysis/" + name + ".gnu";
	
	ofstream gnu(gs.c_str());
	
	gnu << "set terminal postscript eps enhanced color 'Arial' 18" << endl;
	gnu << "set output '" << out << "'" << endl;
	gnu << "set xlabel '" << xlabel << "' offset -3, 0" << endl;
	gnu << "set ylabel '" << ylabel << "' offset -4, -1" << endl;
	gnu << "set zlabel '" << zlabel << "' rotate left offset -0.5, 0" << endl;
	gnu << "set hidden3d" << endl;
	gnu << "set grid ztics" << endl;
	gnu << "set ticslevel 0" << endl;
	gnu << "set border 4095 lw 0.5" << endl;

	gnu << "splot '" << in << "' with lines lw 2 lt 1 lc rgb 'black' title 'NuWro'";
				
	gnu.close();
	
	run_command ("gnuplot " + gs);	
}

void pattern :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
			
	N = new NuWro (P);
		
	run ();
	
	delete N;
	
	finish();
}

void pattern :: run ()
{	
	if(P.dyn_dis_nc or P.dyn_res_nc  or P.dyn_dis_cc or P.dyn_res_cc)
		singlepion(P);
		
	for (int i = 0; i < events; i++)
	{
		event *e = new event();
		e->dyn = N->proces();
		
		N -> makeevent(e,P);
		N -> finishevent(e,P);
						
		if (e->weight >= 0)
			norm[e->dyn]++;
		
		calculate (e);
			
		delete e;
		
		cout << name << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << endl;
}

void mb_ncpi0 :: set_params ()
{
	P.read("data/beam/newMB.txt");
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

void mb_ncpi0 :: finish ()
{
	h1 -> finalize (norm);
	h1 -> save ();
	h1 -> plot ();
	h2 -> finalize (norm);
	h2 -> save ();
	h2 -> plot ();
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
				mom = e->post[i].momentum();
				cos = e->post[i].p().z / mom;
				break;
			}
			
		int a = mom / h1->width_x;
			
		if (a < h1->bins_x)
			h1->part_result[e->dyn][a] += e->weight;
			
		int b = mom / h2->width_x;
		int c  = cos / h2->width_y + 10;
		
		if (b < h2->bins_x and c < h2->bins_y)
			h2->part_result[e->dyn][b][c] += e->weight;
	}
}

pattern * choose (int x)
{
	pattern *wsk;
	
	switch (x)
	{
		case 0: wsk = new mb_ncpi0; break;
		case 1: wsk = new A; break;
		case 2: wsk = new antiA; break;
		default: wsk = NULL;
	}
	
	return wsk;
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
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 1;
}

void A :: finish ()
{
	h1 -> finalize (norm);
	h1 -> save ();
	h1 -> plot ();
	h2 -> finalize (norm);
	h2 -> save ();
	h2 -> plot ();
	h3 -> finalize (norm);
	h3 -> save ();
	h3 -> plot ();
}

void A :: calculate (event *e)
{
	using namespace PDG;
	
	int pion = 100 * e -> fof (pdg_piP) + 10 * e -> fof (-pdg_piP) + e -> fof (pdg_pi);
	
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
			
		int a = Erec / h2->width_x;
		
		if (a < h2 -> bins_x)
			h2 -> part_result [e -> dyn][a] += e -> weight;
	}
	else if (e -> fof (pdg_pi) > 0)
	{
		double energy = 0;
		double cos;
		
		for (int i = 0; i < e -> post.size(); i++)
			if (e -> post[i].pdg == pdg_pi and e -> post[i].E() > energy)
			{
				energy = e -> post[i].E();
				cos = e -> post[i].p().z / e -> post[i].momentum();
			}
			
		int a = energy / h3 -> width_x;
		int b = cos / h3 -> width_y + 10;
		
		if (a < h3 -> bins_x and b < h3 -> bins_y)
			h3 -> part_result [e -> dyn][a][b] += e -> weight;
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
	P.dyn_mec_cc = 1;

	P.dyn_qel_nc = 1;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 1;
	P.dyn_mec_nc = 1;
}

void run_command (string com)
{
	FILE *fp;
	fp = popen (com.c_str(), "w");
	pclose (fp);
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

