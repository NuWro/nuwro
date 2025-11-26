#include "ganalysis.h"
#include "beam.h"
#include "beam_uniform.h"

string intToStr(int n)
{
     string tmp;
     
     if(n < 0) {
          tmp = "-";
          n = -n;
     }
     
     if(n > 9)
          tmp += intToStr(n / 10);
     
     tmp += n % 10 + 48;
     
     return tmp;
}

void parabola (double &a, double &b, double &c, const double* x, const double& start, const double &step, const int bins)
{
	double w[8] = {0};
	
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < bins; j++)
		{
			if (x[j] != 0)
			{
				double temp = start + j*step;
				w[i] += pow(temp, i);
			}
		}
		
	for (int i = 0; i < bins; i++)
	{
		double temp = start + i*step;
		
		w[5] += x[i];
		w[6] += x[i]*temp;
		w[7] += x[i]*temp*temp;
	}
	
	double W  = w[2]*w[2]*w[2] + w[3]*w[3]*w[0] + w[4]*w[1]*w[1] - w[4]*w[2]*w[0] - w[3]*w[1]*w[2] - w[2]*w[3]*w[1];
	double Wc = w[2]*w[2]*w[7] + w[3]*w[3]*w[5] + w[4]*w[1]*w[6] - w[4]*w[2]*w[5] - w[3]*w[1]*w[7] - w[2]*w[3]*w[6];
	double Wb = w[2]*w[6]*w[2] + w[3]*w[7]*w[0] + w[4]*w[5]*w[1] - w[4]*w[6]*w[0] - w[3]*w[5]*w[2] - w[2]*w[7]*w[1];
	double Wa = w[5]*w[2]*w[2] + w[6]*w[3]*w[0] + w[7]*w[1]*w[1] - w[7]*w[2]*w[0] - w[6]*w[1]*w[2] - w[5]*w[3]*w[1];
	
	a = Wa/W;
	b = Wb/W;
	c = Wc/W;
}

void find_minimum (const hist2D& h, double &a, double &b, double &c)
{
	a = 0;
	b = 0;
	c = 1e10;
		
	for (int i = 0; i < h.bins_x; i++)
		for (int j = 0; j < h.bins_y; j++)
			if (c > h.result[i][j])
			{
				a = h.begin_x + i*h.width_x;
				b = h.begin_y + j*h.width_y;
				c = h.result[i][j];
			}
}
	
void find_sigma (const double& a, const double& b, const double &c, double& x1, double& x2, const double& minimum, const double& sigma)
{			
	double delta = b*b - 4.0*a*(c - 54.4045 - 1);//(c - minimum - sigma);

	if (delta >= 0)
	{
		delta = sqrt(delta);
	
		x1 = (-b - delta)/2.0/a;
		x2 = (-b + delta)/2.0/a;
	}
	else
		x1 = x2 = 0;
}

void make_smooth (const hist2D& in, hist2D& out, int ile, double start, double step, int bins, bool inverse)
{
	if (!inverse)
	{
		for (int i = 0; i < ile; i++)
		{
			double a,b,c;
			parabola (a, b, c, in.result[i], start, step, bins);
				
			for (int j = 0; j < bins; j++)
			{
				double x = start + j*step;
			
				out.result[i][j] = a*x*x+b*x+c;
			}
		}
	}
	else
	{	
		for (int i = 0; i < ile; i++)
		{
			double *chi_temp = new double[bins];
		
			for (int j = 0; j < bins; j++)
				chi_temp[j] = in.result[j][i];

			double a,b,c;
			parabola (a, b, c, chi_temp, start, step, bins);
								
			for (int j = 0; j < bins; j++)
			{
				double x = start + j * step;
			
				out.result[j][i] = a*x*x+b*x+c;			
			}
		
			delete chi_temp;
		}
	}
}

void expand (const hist2D& in, hist2D& out, bool inverse)
{
	if (!inverse)
	{
		for (int i = 0; i < in.bins_x; i++)
		{
			double a,b,c;
			parabola (a, b, c, in.result[i], in.begin_y, in.width_y, in.bins_y);
			
			for (int j = 0; j < out.bins_y; j++)
			{
				double x = out.begin_y + j * out.width_y;
				
				out.result[i][j] = a*x*x+b*x+c;
			}
		}
	}
	else
	{
		for (int i = 0; i < in.bins_y; i++)
		{
			
			double *chi_temp = new double[in.bins_x];
		
			for (int j = 0; j < in.bins_x; j++)
				chi_temp[j] = in.result[j][i];

			double a,b,c;
			parabola (a, b, c, chi_temp, in.begin_x, in.width_x, in.bins_x);
			
			for (int j = 0; j < out.bins_x; j++)
			{
				double x = out.begin_x + j * out.width_x;
				
				out.result[j][i] = a*x*x+b*x+c;
			}
			
			delete chi_temp;
		}
	}
}

double hdif (const hist2D& h1, const hist2D& h2)
{
	double res = 0;
	
	for (int i = 0; i < h1.bins_x; i++)
		for (int j = 0; j < h1.bins_y; j++)				
			res += h1.result[i][i] / h2.result[i][j];
		
	return res / h1.bins_x / h1.bins_y;
}

void find_min (const hist2D& h, double& ma, double& ds, double& minimum)
{
	minimum = 1e10;
	int x = 0, y = 0 ;
	
	for (int i = 0; i < h.bins_x; i++)
		for (int j = 0; j < h.bins_y; j++)
			if (h.result[i][j] < minimum)
			{
				minimum = h.result[i][j];
				x = i;
				y = j;
			}
			
	double a,b,c;
	
	parabola (a, b, c, h.result[x], h.begin_y, h.width_y, h.bins_y);
	
	ds = -b/2.0/a;
	
	double *h_temp = new double[h.bins_x];
		
	for (int i = 0; i < h.bins_x; i++)
			h_temp[i] = h.result[i][y];
	
	parabola (a, b, c, h_temp, h.begin_x, h.width_x, h.bins_x);
	
	ma = -b/2.0/a;
	
	delete h_temp;
}

void make_smooth (hist2D& h)
{
	double end_x = h.begin_x + h.bins_x * h.width_x;
	double end_y = h.begin_y + h.bins_y * h.width_y;
	
	hist2D chi_temp1 (h.name + "_s1", "", "", h.xlabel, h.ylabel, h.zlabel, h.bins_x, h.begin_x, end_x, h.bins_y, h.begin_y, end_y, 2, 0);	
	hist2D chi_temp2 (h.name + "_s2", "", "", h.xlabel, h.ylabel, h.zlabel, h.bins_x, h.begin_x, end_x, h.bins_y, h.begin_y, end_y, 2, 0);	
	
	make_smooth (h, chi_temp1, h.bins_x, h.begin_y, h.width_y, h.bins_y, 0);
	make_smooth (h, chi_temp2, h.bins_y, h.begin_x, h.width_x, h.bins_x, 1);
	
	//do
	//{
	//	make_smooth (chi_temp1, chi_temp2, h.bins_x, h.begin_y, h.width_y, h.bins_y, 0);
	//	make_smooth (chi_temp2, chi_temp1, h.bins_y, h.begin_x, h.width_x, h.bins_x, 1);
	//	
	//}
	//while (hdif (chi_temp2, chi_temp1) > 2);	
	
	double min_ma, min_ds, min;
	
	find_min (chi_temp1, min_ma, min_ds, min);
	
	string help = h.name + "_minimum.txt";
			
	ofstream file(help.c_str());
	
	file << min_ma << " " << min_ds << " " << min << endl;
	
	file.close();
	
	hist2D chi_ratio (h.name + "_test", "", "", h.xlabel, h.ylabel, h.zlabel, h.bins_x, h.begin_x, end_x, h.bins_y, h.begin_x, end_y, 2, 0);	
	
	for (int i = 0; i < h.bins_x; i++)
		for (int j = 0; j < h.bins_y; j++)
			chi_ratio.result[i][j] = h.result[i][j] / chi_temp1.result[i][j];
}
	
histogram :: histogram (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int typnormy, double norma)
				: name (name_of_hist), data (data_file), title (plot_title), xlabel (x_name), ylabel (y_name), bins_x (nof_bins), begin_x(begin), norm_type (typnormy), normalization (norma)
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
			if ((norm_type == 6 or norm_type == 7) and tnorm[i][j] != 0)
				result[j] += part_result[i][j] / tnorm[i][j];
			else if (norm_type == 5)
				result[j] += part_result[i][j] / xsec[i];
			else if (norm_type == 8)
				result[j] += part_result[i][j];
			else if (norm[i] != 0)
				result[j] += part_result[i][j] / norm[i] / width_x;

	if (norm_type == 1 or norm_type == 7)
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

void hist1D :: put (double x)
{
	int a = (x - begin_x) / width_x;
	
	if (a >= 0 and a < bins_x)
		result [a]++;
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
	string out3 = "analysis/" + name + "_points.txt";
	
	ofstream file1(out1.c_str());
	ofstream file2(out2.c_str());
	ofstream file3(out3.c_str());
		
	for (int i = 0; i < bins_x; i++)
		for (int j = 0; j < bins_y; j++)
			file3 << begin_x + width_x * i << " " << begin_y + width_y * j << " " << result[i][j] << endl;

	for (int i = 0; i < bins_x; i++)
	{	
		for (int j = 0; j < bins_y; j++)
		{
			if (max < result[i][j]) max = result[i][j];

			//if (result[i][j] != 0)
			{
				file1 << begin_x + width_x * i << " " << begin_y + width_y * j << " " << result[i][j] << endl;
				file1 << begin_x + width_x * i << " " << begin_y + width_y * (j + 1) << " " << result[i][j] << endl;
			}
		}
		
		file1 << endl;
		
		for (int j = 0; j < bins_y; j++)
		{
			//if (result[i][j] != 0)
			{
				file1 << begin_x + width_x * (i + 1) << " " << begin_y + width_y * j << " " << result[i][j] << endl;
				file1 << begin_x + width_x * (i + 1) << " " << begin_y + width_y * (j + 1) << " " << result[i][j] << endl;
			}
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

//	gnu << "set output '" << out3 << "'" << endl;
//	gnu << "plot '" << in2 << "' with lines lw 2 lt 1 lc rgb 'black' title 'NuWro (max = " << max << ")'";

	gnu << endl << endl;

	gnu << "set zlabel '" << zlabel << "' rotate left offset -0.5, 0" << endl;
	gnu << "set hidden3d" << endl;
	gnu << "set grid ztics" << endl;
	gnu << "set ticslevel 0" << endl;
	gnu << "set border 4095 lw 0.5" << endl;

//	gnu << "set output '" << out1 << "'" << endl;
//	gnu << "splot '" << in1 << "' with lines lw 2 lt 1 lc rgb 'black' title 'NuWro'";

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

void pattern :: startfsi ()
{			

	P.read("data/kaskada.txt");
	set_params ();

  // load the input data
  input.initialize( P );
  input.load_data();
	
	nucleus* nucl= make_nucleus(P);
	
	do
	{	
		for (int i = 0; i < events; i++)
		{
			event *e = new event();
			e->weight = 1;
			e->par = P;
			
			beam_uniform b(P);
			particle p0 = b.shoot();
			p0.r = start_point(nucl,P);
			e->out.push_back(p0);
			e->in.push_back(p0);
			
			k = new kaskada(P,*e,&input);
			k -> kaskadaevent();
			
			calculate(e);
				
			delete e;
			delete k;
			
			cout << name << ": " << 100*i/events << "%\r" << flush;
		}
				
		loop--;
	}
	while(loop > 0 and change_params());
	
	delete nucl;
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
		double cos=0;
		
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
			double cos = 0;

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
			double cos = 0;

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
	
	P.FSI_on = 0;
	
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
	
	P.FSI_on = 0;
	
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
		
	double mom = 0, cos = 0;
	
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

	P.FSI_on = 0;
	
	P.qel_nc_axial_mass = 1030;
	
	P.mec_kind = 1;
	
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
	double energy = e -> in[0].E();// / 1000.0;

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
	run_command ("mkdir -p analysis/mb_nce/temp/");
	
	zero ();
}

mb_nce :: ~mb_nce ()
{	
}

void mb_nce :: zero()
{
	for (int i = 0; i < 51; i++)
		res[i] = 0;
}

void mb_nce_he :: zero()
{
	for (int i = 0; i < 30; i++)
	{
		licznik[i] = 0;
		mianownik[i] = 0;
	}
}

void mb_nce :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
	
	N = new NuWro ();

	N -> set(P);

	if (mb_nce_mode == 3)
	{
		cout << "QEL on Hydrogen" << endl;
	
		run ();
		
		name = "qel_hydrogen.txt";
		save ();		
	}
	else if (mb_nce_mode == 2)
	{
		//QEL on Hydrogen

		for (int i = 0; i < mb_nce_nof; i++)
		{
			P.qel_nc_axial_mass = mb_nce_start + i*10;

			ff_configure (P);
				
			stringstream temp;
			string mas;

			temp << mb_nce_start + i*10;
			temp >> mas;
				
			name = "qelH_" + mas + ".txt";
				
			cout << "QEL on Hydrogen, Axial mass = " << mas << " MeV" << endl;
			
			run ();
			
			save ();
			
			zero ();
		}
	}
	else if (mb_nce_mode == 5)
	{
		cout << "RES on Carbon" << endl;

		P.read("data/target/C.txt");
	
		N -> refresh_target (P);
	
		P.dyn_qel_nc = 0;
		P.dyn_res_nc = 1;

		N -> refresh_dyn (P);
	
		singlepion(P);

		name = "res.txt";
			
		run ();
		save();
	}
	else if (mb_nce_mode == 6)
	{
		cout << "DIS on Carbon" << endl;

		P.read("data/target/C.txt");
	
		N -> refresh_target (P);
	
		P.dyn_qel_nc = 0;
		P.dyn_dis_nc = 1;
	
		N -> refresh_dyn (P);

		name = "dis.txt";

		run ();
		
		save();
	}
	else if (mb_nce_mode == 4)
	{
		P.read("data/target/C.txt");	
		N -> refresh_target (P);
						
		cout << "MEC on Carbon" << endl;
			
		P.dyn_qel_nc = 0;
		P.dyn_mec_nc = 1;	

		N -> refresh_dyn (P);

		run ();
		
		name = "mec.txt";
		save();
	}
	else if (mb_nce_mode == 1)
	{
		//QEL on Carbon

		P.read("data/target/C.txt");	
		N -> refresh_target (P);
		
		for (int i = 0; i < mb_nce_nof; i++)
		{
			P.qel_nc_axial_mass = mb_nce_start + i*10;

			ff_configure (P);
				
			stringstream temp;
			string mas;

			temp << mb_nce_start + i*10;
			temp >> mas;
				
			name = "qel_" + mas + ".txt";
							
			cout << "QEL on Carbon, Axial mass = " << mas << " MeV" << endl;
			
			run ();
			
			save ();
			
			zero ();
		}
	}
	else if (mb_nce_mode == 0)
		chi2 ();
				
	delete N;	
}

void mb_nce :: save()
{
	ofstream out (("analysis/mb_nce/temp/" + name).c_str());
	
	for (int i = 0; i < 51; i++)
		out << res[i] << endl;
	
	out.close();
}

void mb_nce :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/H.txt");
	
	P.qel_nc_axial_mass = 1030;
	
	P.mec_kind = 1;
	
	P.delta_s = 0;
			
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

	double total = 0;
			
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_proton or e -> post[i].pdg == pdg_neutron)
			total += e -> post[i].Ek();
				
	if (P.nucleus_p == 1)
	{
		int w;
		int bin = t2r (total, 0, w);
		if (bin >= 0)
			res [bin] += e -> weight * factor * 2.0 * wagi [0][w];
	}
	else if ( (e -> flag.res or e -> flag.dis) and (e -> fof (pdg_pi) + e -> fof (pdg_piP) + e -> fof (-pdg_piP) == 0))
	{
		int w;
		int bin = t2r (total, 4, w);
		if (bin >= 0)
			res [bin] += e -> weight * factor * 12.0 * wagi [4][w];
	}
	else if (e -> flag.mec)
	{
		int w;
		
		int help = 2;
		
		if (e -> in[1].pdg == pdg_neutron and e -> in[2].pdg == pdg_neutron)
			help = 3;

	    int bin = t2r (total, help, w);
	    if (bin >= 0)
			res [bin] += e -> weight * factor * 12.0 * wagi [help][w];
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
			res [bin] += e -> weight * factor * 12.0 * wagi [scenario][w];			
	}	
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

void mb_nce :: chi2 ()
{	
	mhist1D chi ("mb_nce/chi_mb_nce", "", "", "Axial Mass [MeV]", "{/Symbol x}^{2}", mb_nce_nof, mb_nce_start, mb_nce_start + 10*(mb_nce_nof - 1), 3, 0);

	chi.finalized = true;
	chi.cnames[0] = "without MEC";
	chi.cnames[1] = "with MEC";
	chi.cnames[2] = "with MEC + MAH = MAC";
	
	for (int i = 0; i < mb_nce_nof; i++)
	{				
		stringstream temp;
		string mas;

		temp << mb_nce_start + i*10;
		temp >> mas;
		
		name = "mb_nce/mb_nce_" + mas;
		
		mhist1D final (name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 3, 0);
		final.cnames[0] = "without MEC";
		final.cnames[1] = "with MEC";
		final.cnames[2] = "with MEC + MAH = MAC";
		final.finalized = true;	
		
		name = "mb_nce/mb_nce_" + mas + "_part_nomec";
		
		mhist1D part_nomec (name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 5, 0);
		part_nomec.cnames[0] = "total";
		part_nomec.cnames[1] = "hydrogen";
		part_nomec.cnames[2] = "qel";
		part_nomec.cnames[3] = "res/dis";
		part_nomec.cnames[4] = "background";
		part_nomec.finalized = true;	
					
		name = "mb_nce/mb_nce_" + mas + "_part_mec1";
		
		mhist1D part_mec1 (name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 6, 0);
		part_mec1.cnames[0] = "total";
		part_mec1.cnames[1] = "hydrogen";
		part_mec1.cnames[2] = "qel";
		part_mec1.cnames[3] = "res/dis";
		part_mec1.cnames[4] = "mec";
		part_mec1.cnames[5] = "background";
		part_mec1.finalized = true;	
		
		name = "mb_nce/mb_nce_" + mas + "_part_mec2";
		
		mhist1D part_mec2 (name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 6, 0);
		part_mec2.cnames[0] = "total";
		part_mec2.cnames[1] = "hydrogen";
		part_mec2.cnames[2] = "qel";
		part_mec2.cnames[3] = "res/dis";
		part_mec2.cnames[4] = "mec";
		part_mec2.cnames[5] = "background";
		part_mec2.finalized = true;	
					
		load (final, part_nomec, part_mec1, part_mec2, mas);

		chi.result[0][i] = calc_chi (final.result [0]);
		chi.result[1][i] = calc_chi (final.result [1]);
		chi.result[2][i] = calc_chi (final.result [2]);
	}
}

void mb_nce :: load (mhist1D& h, mhist1D& part_nomec, mhist1D& part_mec1, mhist1D& part_mec2, string Ma)
{
	double H[51] = {0};
	double RES[51] = {0};
	double MEC[51] = {0};
	double QEL[51] = {0};
	double HQEL[51] = {0};
	
	load (H, "qel_hydrogen.txt");
	load (RES, "resdis.txt");
	load (MEC, "mec.txt");
	load (QEL, "qel_" + Ma + ".txt");
	load (HQEL, "qelH_" + Ma + ".txt");
	
	for (int i = 0; i < 51; i++)
	{
		h.result[0][i] += H[i] + RES[i] + QEL[i] + mb_nce_bg[i];
		h.result[1][i] += H[i] + RES[i] + QEL[i] + MEC[i] + mb_nce_bg[i];
		h.result[2][i] += HQEL[i] + RES[i] + QEL[i] + + MEC[i] + mb_nce_bg[i];
		
		part_nomec.result[0][i] += H[i] + RES[i] + QEL[i] + mb_nce_bg[i];
		part_nomec.result[1][i] += H[i];
		part_nomec.result[2][i] += QEL[i];
		part_nomec.result[3][i] += RES[i];
		part_nomec.result[4][i] += mb_nce_bg[i];

		part_mec1.result[0][i] += H[i] + RES[i] + QEL[i] + MEC[i] + mb_nce_bg[i];
		part_mec1.result[1][i] += H[i];
		part_mec1.result[2][i] += QEL[i];
		part_mec1.result[3][i] += RES[i];
		part_mec1.result[4][i] += MEC[i];
		part_mec1.result[5][i] += mb_nce_bg[i];

		part_mec2.result[0][i] += HQEL[i] + RES[i] + QEL[i] + + MEC[i] + mb_nce_bg[i];
		part_mec2.result[1][i] += HQEL[i];
		part_mec2.result[2][i] += QEL[i];
		part_mec2.result[3][i] += RES[i];
		part_mec2.result[4][i] += MEC[i];
		part_mec2.result[5][i] += mb_nce_bg[i];
	}
}

void mb_nce :: load (double* x, string fname)
{
	ifstream in (("analysis/mb_nce/temp/" + fname).c_str());
	int i = 0;
		
	while (i < 51)
		in >> x[i++];
		
	in.close();
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

//MB NCE HE sample

mb_nce_he :: mb_nce_he () : pattern ("MB NCE HE", 100000000)
{	
	run_command ("mkdir -p analysis/mb_nce_he/");
	run_command ("mkdir -p analysis/mb_nce_he/temp");
			
	zero();
}

mb_nce_he :: ~mb_nce_he ()
{	

}

void mb_nce_he :: load (mhist1D& h, mhist1D& nomec, mhist1D& mec, string ds)
{
	double H1[30] = {0};
	double RES1[30] = {0};
	double MEC1[30] = {0};
	double QEL1[30] = {0};

	double H2[30] = {0};
	double RES2[30] = {0};
	double MEC2[30] = {0};
	double QEL2[30] = {0};
	
	load (H1, H2, "Hqel_" + ds + ".txt");
	load (RES1, RES2, "resdis.txt");
	load (MEC1, MEC2, "mec.txt");
	load (QEL1, QEL2, "qel_" + ds + ".txt");
		
	for (int i = 0; i < 30; i++)
	{
		h.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		h.result[1][i] += (H1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		
		nomec.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[1][i] += H1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[2][i] += QEL1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[3][i] += RES1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[4][i] += mb_nce_bg_he_sp[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		
		mec.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[1][i] += H1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[2][i] += QEL1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[3][i] += RES1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[4][i] += MEC1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[5][i] += mb_nce_bg_he_sp[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
	}
}

void mb_nce_he :: load (double* x, double *y, string fname)
{
	ifstream in (("analysis/mb_nce_he/temp/" + fname).c_str());
	int i = 0;
	
	while (i < 30)
	{
		in >> x[i];
		in >> y[i++];
	}
		
	in.close();
}

void mb_nce_he :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
	
	N = new NuWro ();

	N -> set(P);

	if (mb_nce_mode == 3)
	{
		cout << "QEL on Hydrogen" << endl;
	
		for (int i = 0; i < mb_nce_nof; i++)
		{
			P.delta_s = mb_nce_start + i*0.1;

			ff_configure (P);
				
			stringstream temp;
			string deltas;

			temp << mb_nce_start + i*0.1;
			temp >> deltas;
				
			name = "Hqel_" + deltas + ".txt";
				
			cout << "QEL on Hydrogen, Ds = " << deltas << endl;
			
			run ();
			
			save ();
			
			zero ();
		}
		
	}
	else if (mb_nce_mode == 5)
	{
		singlepion(P);

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
		
		name = "resdis.txt";
		save();
	}
	else if (mb_nce_mode == 4)
	{
		P.read("data/target/C.txt");	
		N -> refresh_target (P);
						
		cout << "MEC on Carbon" << endl;
			
		P.dyn_qel_nc = 0;
		P.dyn_mec_nc = 1;	

		N -> refresh_dyn (P);

		run ();
		
		name = "mec.txt";
		save();
	}
	else if (mb_nce_mode == 1)
	{
		//QEL on Carbon

		P.read("data/target/C.txt");	
		N -> refresh_target (P);
		
		P.qel_nc_axial_mass = mb_nce_he_ma;

		for (int i = 0; i < mb_nce_nof; i++)
		{
			P.delta_s = mb_nce_start + i*0.1;

			ff_configure (P);
				
			stringstream temp;
			string deltas;

			temp << mb_nce_start + i*0.1;
			temp >> deltas;
				
			name = "qel_" + deltas + ".txt";
				
			cout << "QEL on Carbon, Ds = " << deltas << " MeV" << endl;
			
			run ();
			
			save ();
			
			zero ();
		}
	}
	else if (mb_nce_mode == 0)
		chi2 ();
				
	delete N;
}

void mb_nce_he :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/H.txt");
	
	P.qel_nc_axial_mass = 1030;
	
	P.mec_kind = 1;
		
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

void mb_nce_he :: calculate (event *e)
{		
	using namespace PDG;
	
	const double factor = 2.95e45 / 14.0 / events / 18.0;

	double total = 0;
	double sp = 0;		
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_proton or e -> post[i].pdg == pdg_neutron)
		{
			total += e -> post[i].Ek();
			
			if (e -> post[i].pdg == pdg_proton and e -> post[i].primary and e -> post[i].p().z / e -> post[i].momentum() > 0.5)
				sp = e -> post[i].Ek();		
		}
				
	if (P.nucleus_p == 1)
	{
		int w_sp, w_all;
		int bin_sp = t2r (total, 0, w_sp, dystrR_he_sp);
		int bin_all = t2r (total, 0, w_all, dystrR_he_all);
			
		if (bin_sp >= 0 and e -> out[1].p().z / e -> out[1].momentum() > 0.5)
			licznik[bin_sp] += e -> weight * factor * 2.0 * wagi_he_sp [0][w_sp];

		if (bin_all >= 0)
			mianownik [bin_all] += e -> weight * factor * 2.0 * wagi_he_all [0][w_all];
	}
	else if ( (e -> flag.res or e -> flag.dis) and (e -> fof (pdg_pi) + e -> fof (pdg_piP) + e -> fof (-pdg_piP) == 0))
	{
		int w_sp, w_all;
		int bin_sp = t2r (total, 4, w_sp, dystrR_he_sp);
		int bin_all = t2r (total, 4, w_all, dystrR_he_all);
			
		if (bin_sp >= 0)
			licznik [bin_sp] += e -> weight * factor * 12.0 * wagi_he_sp [4][w_sp];

		if (bin_all >= 0)
			mianownik [bin_all] += e -> weight * factor * 12.0 * wagi_he_all [4][w_all];
	}
	else if (e -> flag.mec)
	{
		int w_sp, w_all;
		int bin_sp  = t2r (sp, 1, w_sp, dystrR_he_sp);
		
		int help = 2;
		
		if (e -> in[1].pdg == pdg_neutron and e -> in[2].pdg == pdg_neutron)
			help = 3;
		
		int bin_all = t2r (total, help, w_all, dystrR_he_all);
		
		if (bin_all >= 0)
		{
			if (sp > 0) licznik [bin_all] += e -> weight * factor * 12.0 * wagi_he_sp [help][w_sp];
			mianownik [bin_all] += e -> weight * factor * 12.0 * wagi_he_all [help][w_all];		
		}
	}
	else if (e -> flag.qel)
	{
		int scenario = 1;
			
		if (e -> in[1].pdg == pdg_neutron)
			scenario = 3;
		else if (e -> number_of_interactions ())
			scenario = 2;
			
		int w_sp, w_all;
						
		int bin_sp = t2r (total, scenario, w_sp, dystrR_he_sp);
		int bin_all = t2r (total, scenario, w_all, dystrR_he_all);
			
		if (bin_sp >= 0 and e -> out[1].p().z / e -> out[1].momentum() > 0.5)
			licznik [bin_sp] += e -> weight * factor * 12.0 * wagi_he_sp [scenario][w_sp];			
	
		if (bin_all >= 0)
			mianownik [bin_all] += e -> weight * factor * 12.0 * wagi_he_all [scenario][w_all];			
	}
}

int mb_nce_he :: t2r (double x, int s, int &w, const double R[5][30][30])
{
	if (x <= 0)
		return -1;
	
	if (x < 300)
		w = 0;
	else if (x > 900)
		w = 29;
	else
		w = (x - 300) * 28 / 600;
		
	double los = frandom00 ();
		
	return find_bin (s, w, los, R);
}

int mb_nce_he :: find_bin (int s, int col, double x, const double R[5][30][30])
{
	int a = 0;
	int b = 29;
	int c = (b - a) / 2 + a;
			
	if (R [s][0][col] == 1)
		return -1;
	
	if (x <= R [s][0][col])
		return 0;
	
	do
	{
		if (x > R [s][a][col] and x <= R [s][c][col])
			b = c;
		else
			a = c;
			
		c = (b - a) / 2 + a;
	}
	while (b - a > 1);

	return b;
}		

void mb_nce_he :: chi2 ()
{
	mhist1D chi ("mb_nce_he/chi_mb_nce", "", "", "{/Symbol d}s", "{/Symbol x}^{2}", mb_nce_nof, mb_nce_start, mb_nce_start + 0.1*mb_nce_nof, 2, 0);

	chi.finalized = true;
	chi.cnames[0] = "without MEC";
	chi.cnames[1] = "with MEC";
	
	for (int i = 0; i < mb_nce_nof; i++)
	{				
		stringstream temp;
		string deltas;
		
		double ds = mb_nce_start + i*0.1;
		
		if (ds > 0 and ds < 0.1)
			ds = 0;
		
		temp << ds;
		temp >> deltas;
		
		name = "mb_nce_he/mb_nce_" + deltas;
		
		mhist1D final (name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 2, 0);
		final.cnames[0] = "without MEC";
		final.cnames[1] = "with MEC";
		final.finalized = true;	
		
		name = "mb_nce_he/mb_nce_" + deltas + "_part_nomec";
		
		mhist1D nomec (name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 5, 0);	
		nomec.cnames[0] = "total";
		nomec.cnames[1] = "hydrogen";
		nomec.cnames[2] = "qel";
		nomec.cnames[3] = "res/dis";
		nomec.cnames[4] = "background";
		nomec.finalized = true;	

		name = "mb_nce_he/mb_nce_" + deltas + "_part_mec";
		
		mhist1D mec (name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 6, 0);	
		mec.cnames[0] = "total";
		mec.cnames[1] = "hydrogen";
		mec.cnames[2] = "qel";
		mec.cnames[3] = "res/dis";
		mec.cnames[4] = "mec";
		mec.cnames[5] = "background";
		mec.finalized = true;	
			
		load (final, nomec, mec, deltas);

		chi.result[0][i] = calc_chi (final.result [0]);		
		chi.result[1][i] = calc_chi (final.result [1]);
	}	
}

double mb_nce_he :: calc_chi (double *x)
{
	double res = 0;
	
	double dif [30] = {0};
		
	for (int i = 0; i < 30; i++)
		dif[i] = mb_nce_he_sp[i] / mb_nce_he_all[i] - x[i];
				
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 30; j++)
			res += dif[i] * Mrev_he [i][j] * dif[j];	

	return res;	
}

void mb_nce_he :: save()
{
	ofstream out (("analysis/mb_nce_he/temp/" + name).c_str());
	
	for (int i = 0; i < 30; i++)
		out << licznik[i] << " " << mianownik[i] << endl;
	
	out.close();
}

//MB NCEL both

mb_nce_both :: mb_nce_both () : pattern ("MB NCE BOTH", 1000000) //100000000)
{	
	run_command ("mkdir -p analysis/mb_nce_both/");
	run_command ("mkdir -p analysis/mb_nce_both/ma/");
	run_command ("mkdir -p analysis/mb_nce_both/ds/");
	run_command ("mkdir -p analysis/mb_nce_both/temp");
	run_command ("mkdir -p analysis/mb_nce_both/temp/ma");
	run_command ("mkdir -p analysis/mb_nce_both/temp/ds");
	
	//minimum = 0;
	
	zero();
}

mb_nce_both :: ~mb_nce_both ()
{	

}

void mb_nce_both :: zero()
{	
	for (int i = 0; i < 51; i++)
	{
		ma_res[i] = 0;
		
		for (int j = 0; j < 6; j++)
		{
			true_ma_res[j][i] = 0;
			rec_ma_res[j][i] = 0;
		}
		
		if (i < 30)
		{
			ds_licznik[i] = 0;
			ds_mianownik[i] = 0;
			
			for (int j = 0; j < 6; j++)
			{
				true_ratio[j][i] = 0;
				rec_licznik[j][i] = 0;
				rec_mianownik[j][i] = 0;
			}
		}
	}
}

void mb_nce_both :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/H.txt");
	
	P.qel_nc_axial_mass = mb_nce_ma;
	
	P.mec_kind = 1;
	
	P.delta_s = mb_nce_ds;
			
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

bool nonFile(string filename)
{
	fstream plik;
	plik.open(filename.c_str());
	
	if (plik.is_open())
	{
		plik.close();
		return false;
	}
	
	plik.close();
	return true;	
}

void mb_nce_both :: start ()
{	
	P.read("data/params.txt");
	
	set_params ();
	
	N = new NuWro ();

	N -> set(P);

	if (mb_nce_mode == 3)
	{
		P.FSI_on = 0;
		
		cout << "QEL on Hydrogen" << endl;
					
		stringstream ma_temp, ds_temp;
		
		string ma;
		string ds;

		ma_temp << mb_nce_ma;
		ma_temp >> ma;
		
		ds_temp << mb_nce_ds;
		ds_temp >> ds;
				
		name = "Hqel_" + ma + "_" + ds + ".txt";
				
		cout << "QEL on Hydrogen, Ds = " << ds << endl;
			
		run ();
		save ();			
	}
	else if (mb_nce_mode == 5)
	{
		cout << "RES on Carbon" << endl;

		P.read("data/target/C.txt");
	
		N -> refresh_target (P);
	
		P.dyn_qel_nc = 0;
		P.dyn_res_nc = 1;
		P.dyn_dis_nc = 1;
		
		N -> refresh_dyn (P);
	
		singlepion(P);

		name = "resdis.txt";

		P.dyn_res_nc = 1;
		P.dyn_dis_nc = 0;
			
		run ();

		P.dyn_res_nc = 1;
		P.dyn_dis_nc = 0;
		
		cout << "DIS on Carbon" << endl;
			
		run ();

		save();
	}
	else if (mb_nce_mode == 4)
	{
		P.read("data/target/C.txt");	
		N -> refresh_target (P);
						
		cout << "MEC on Carbon" << endl;
			
		P.dyn_qel_nc = 0;
		P.dyn_mec_nc = 1;	

		N -> refresh_dyn (P);

		run ();
		
		name = "mec.txt";
		save();
	}
	else if (mb_nce_mode == 1)
	{
		//QEL on Carbon

		P.read("data/target/C.txt");	
		N -> refresh_target (P);
					
		stringstream ma_temp, ds_temp;
		
		string ma;
		string ds;

		ma_temp << mb_nce_ma;
		ma_temp >> ma;
		
		ds_temp << mb_nce_ds;
		ds_temp >> ds;
				
		name = "qel_" + ma + "_" + ds + ".txt";
				
		cout << "QEL on Carbon, Ma = " << ma << " MeV, ds = " << ds << endl;
		
		string temp = "analysis/mb_nce_both/temp/ma/" + name;
		
		//if (nonFile(temp))
		{		
			run ();
			save ();
		}			
	}
	else if (mb_nce_mode == 0)
		chi2 ();
	
	if (N)						
	delete N;
	
	//cout << "Minimum = " << minimum << endl;
}

void mb_nce_both :: calculate (event *e)
{		
	using namespace PDG;
	
	const double factor = 2.95e45 / 14.0 / events / 18.0;

	double total = 0;
	double sp = 0;		
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_proton or e -> post[i].pdg == pdg_neutron)
		{
			total += e -> post[i].Ek();
			
			if (e -> post[i].pdg == pdg_proton and e -> post[i].Ek() > sp)// and e -> post[i].primary and e -> post[i].p().z / e -> post[i].momentum() > 0.5)
				sp = e -> post[i].Ek();		
		}
	
	static int delta_counter = 0;
			
	//sp = total;
		
	int w, w_sp, w_all, bin, bin_sp, bin_all, scenario = -1;
	double weight = 12.0; //carbon
				
	if (P.nucleus_p == 1)
	{
		scenario = 0;
		weight = 2.0;
	}
	else if ( (e -> flag.res or e -> flag.dis) and (e -> fof (pdg_pi) + e -> fof (pdg_piP) + e -> fof (-pdg_piP) == 0))
	{
		if (e -> flag.res)
			delta_counter++;
			
		scenario = 4;
	}
	else if (false) //e -> flag.res)
	{
		delta_counter++;
		
		if (delta_counter >= 5)
		{
			delta_counter -= 5;
			
			for (int i = 0; i < e -> post.size(); i++)
				if (e -> post[i].pdg == pdg_pi or e -> post[i].pdg == -pdg_piP or e -> post[i].pdg == pdg_piP)
					total += e -> post[i].E();
			
			scenario = 4;
		}
	}
	else if (e -> flag.mec)
	{
		scenario = 2;
		
		if (e -> in[1].pdg == pdg_neutron and e -> in[2].pdg == pdg_neutron)
			scenario = 3;		
	}
	else if (e -> flag.qel)
	{
		scenario = 1;
			
		if (e -> in[1].pdg == pdg_neutron)
			scenario = 3;
		else if (e -> number_of_interactions ())
			scenario = 2;
	}
	
	if (scenario != -1)
	{
	
		bin = t2r (total, scenario, w);
										
		if (e -> flag.mec)
			bin_sp  = t2r (sp, 1, w_sp, dystrR_he_sp);
		else	
			bin_sp = t2r (total, scenario, w_sp, dystrR_he_sp);
		
		bin_all = t2r (total, scenario, w_all, dystrR_he_all);
		
		int kategoria = scenario;
		
		if (e -> flag.mec)
			kategoria = 5;
				
		//if (bin_sp >= 0 and e -> out[1].p().z / e -> out[1].momentum() > 0.5 and (!e -> flag.mec or sp > 0))
		if (bin_sp >= 0 and (!e -> flag.mec or sp > 0))
		{
			ds_licznik [bin_sp] += e -> weight * factor * weight * wagi_he_sp [scenario][w_sp];
			rec_licznik [kategoria][bin_sp] += e -> weight * factor * weight * wagi_he_sp [scenario][w_sp];
		}
		
		if (bin_all >= 0)
		{
			ds_mianownik [bin_all] += e -> weight * factor * weight * wagi_he_all [scenario][w_all];
			rec_mianownik [kategoria][bin_all] += e -> weight * factor * weight * wagi_he_all [scenario][w_all];
			true_ratio [kategoria][w_all] += e -> weight * factor * weight;// * wagi_he_all [scenario][w_all];
		}
						
		if (bin >= 0)
		{
			ma_res [bin] += e -> weight * factor * weight * wagi [scenario][w];
			rec_ma_res [kategoria][bin] += e -> weight * factor * weight * wagi [scenario][w];
			true_ma_res [kategoria][w] += e -> weight * factor * weight;// * wagi [scenario][w];
		}
	}			
}

int mb_nce_both :: t2r (double x, int s, int &w)
{
	if (x <= 0)
		return -1;
	
	w = x / 18;
	
	if (w > 50)	w = 50;
	
	double los = frandom00 ();
		
	return find_bin (s, w, los);
}

int mb_nce_both :: find_bin (int s, int col, double x)
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

int mb_nce_both :: t2r (double x, int s, int &w, const double R[5][30][30])
{
	if (x <= 0)
		return -1;
	
	if (x < 300)
		w = 0;
	else if (x > 900)
		w = 29;
	else
		w = (x - 300) * 28 / 600;
		
	double los = frandom00 ();
		
	return find_bin (s, w, los, R);
}

int mb_nce_both :: find_bin (int s, int col, double x, const double R[5][30][30])
{
	int a = 0;
	int b = 29;
	int c = (b - a) / 2 + a;
			
	if (R [s][0][col] == 1)
		return -1;
	
	if (x <= R [s][0][col])
		return 0;
	
	do
	{
		if (x > R [s][a][col] and x <= R [s][c][col])
			b = c;
		else
			a = c;
			
		c = (b - a) / 2 + a;
	}
	while (b - a > 1);

	return b;
}

void mb_nce_both :: save()
{
	ofstream ma_out (("analysis/mb_nce_both/temp/ma/" + name).c_str());
	ofstream ds_out (("analysis/mb_nce_both/temp/ds/" + name).c_str());	

	for (int i = 0; i < 51; i++)
		ma_out << ma_res[i] << endl;
	
	for (int i = 0; i < 30; i++)
		ds_out << ds_licznik[i] << " " << ds_mianownik[i] << endl;
	
	ma_out.close();
	ds_out.close();
	
	ofstream true_ma_out (("analysis/mb_nce_both/temp/ma/" + name + "_true").c_str());
	ofstream rec_ma_out (("analysis/mb_nce_both/temp/ma/" + name + "_rec").c_str());
	ofstream true_ratio_out (("analysis/mb_nce_both/temp/ds/" + name + "_true").c_str());
	ofstream rec_ratio1_out (("analysis/mb_nce_both/temp/ds/" + name + "_rec1").c_str());
	ofstream rec_ratio2_out (("analysis/mb_nce_both/temp/ds/" + name + "_rec2").c_str());
	
	for (int j = 0; j < 51; j++)
	{
		true_ma_out << (j + 0.5) * 18 << " ";
		rec_ma_out << 40 + (j + 0.5) * 11.960784314 << " ";
		
		for (int i = 0; i < 6; i++)
		{
			true_ma_out << true_ma_res[i][j] << " ";
			rec_ma_out << rec_ma_res[i][j] << " ";
		}
		
		true_ma_out << endl;
		rec_ma_out << endl;
	}
	
	for (int j = 0; j < 30; j++)
	{
		true_ratio_out << 300 + (j - 0.5) * 21.428571429 << " ";
		rec_ratio1_out << 350 + (j + 0.5) * 7.5 << " ";
		rec_ratio2_out << 350 + (j + 0.5) * 7.5 << " ";
		
		for (int i = 0; i < 6; i++)
		{
			true_ratio_out << true_ratio[i][j] << " ";
			rec_ratio1_out << rec_licznik[i][j] << " ";
			rec_ratio2_out << rec_mianownik[i][j] << " ";
		}
		
		true_ratio_out << endl;
		rec_ratio1_out << endl;
		rec_ratio2_out << endl;
	}

	true_ma_out.close();
	rec_ma_out.close();
	rec_ratio1_out.close();
	rec_ratio2_out.close();
	true_ratio_out.close();
}

void mb_nce_both :: chi2 ()
{		
	int mb_nce_ma_bins = (mb_nce_ma_end - mb_nce_ma_start) / mb_nce_ma_step;
	int mb_nce_ds_bins = (mb_nce_ds_end - mb_nce_ds_start) / mb_nce_ds_step;
	
	hist2D chi ("mb_nce_both/chi_mb_nce", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
	hist2D chir ("mb_nce_both/chi_mb_nce_ratio", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
	hist2D chid ("mb_nce_both/chi_mb_nce_distr", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);

	chi.finalized = true;
	chir.finalized = true;
	chid.finalized = true;
	
	hist2D chi_mec ("mb_nce_both/chi_mb_nce_mec", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
	hist2D chi_mecr ("mb_nce_both/chi_mb_nce_mec_ratio", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
	hist2D chi_mecd("mb_nce_both/chi_mb_nce_mec_distr", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);

	chi_mec.finalized = true;
	chi_mecr.finalized = true;
	chi_mecd.finalized = true;
	
//	hist2D chi_mec_smooth1 ("mb_nce_both/chi_mb_nce_mec_smooth1", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
//	hist2D chi_mec_smooth2 ("mb_nce_both/chi_mb_nce_mec_smooth2", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
//	hist2D chi_mec_ratio ("mb_nce_both/chi_mb_nce_mec_ratio", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);

//	chi_mec_smooth1.finalized = true;
//	chi_mec_smooth2.finalized = true;
//	chi_mec_ratio.finalized = true;
	
//	hist2D chi_mec2 ("mb_nce_both/chi_mb_nce_mec2", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ds_bins, mb_nce_ds_start, mb_nce_ds_end, 2, 0);
//	chi_mec2.finalized = true;

//	hist2D chi_mec_smooth_expand ("mb_nce_both/chi_mb_nce_mec_smooth_expand", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, mb_nce_ma_bins/4, mb_nce_ds_start - 0.2, mb_nce_ds_end + 0.2, 2, 0);
//	hist2D chi_mec_smooth_expand1 ("mb_nce_both/chi_mb_nce_mec_smooth_expand1", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start, mb_nce_ma_end, 36, -0.8, 0.4, 2, 0);
//	hist2D chi_mec_smooth_expand2 ("mb_nce_both/chi_mb_nce_mec_smooth_expand2", "", "", "Axial Mass (MeV)", "{Symbol d}s", "{/Symbol x}^{2}", mb_nce_ma_bins, mb_nce_ma_start - 50, mb_nce_ma_end + 50, 36, -0.8, 0.4, 2, 0);
//	chi_mec_smooth_expand1.finalized = true;
//	chi_mec_smooth_expand2.finalized = true;
	
	for (int i = 0; i < mb_nce_ma_bins; i++)
	{	
		for (int j = 0; j < mb_nce_ds_bins; j++)
		{
						
			stringstream ds_temp, ma_temp;
			string ma, ds;
			
			double delta = mb_nce_ds_start + j*mb_nce_ds_step;
			double masa = mb_nce_ma_start + i*mb_nce_ma_step;
		
			if (delta > 0 and delta < 0.01)
				delta = 0;
		
			ds_temp << delta;
			ds_temp >> ds;
			
			ma_temp << masa;
			ma_temp >> ma;
		
			string ds_name = "mb_nce_both/ds/mb_nce_" + ma + "_" + ds;
			string ma_name = "mb_nce_both/ma/mb_nce_" + ma + "_" + ds;
			
			mhist1D ds_final (ds_name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 3, 0);
			ds_final.cnames[0] = "without MEC";
			ds_final.cnames[1] = "with MEC";
			ds_final.cnames[2] = "with MEC2";
			ds_final.finalized = true;	
		
			ds_name = "mb_nce_both/ds/mb_nce_" + ma + "_" + ds + "_part_nomec";		
			
			mhist1D ds_nomec (ds_name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 5, 0);	
			ds_nomec.cnames[0] = "total";
			ds_nomec.cnames[1] = "hydrogen";
			ds_nomec.cnames[2] = "qel";
			ds_nomec.cnames[3] = "res/dis";
			ds_nomec.cnames[4] = "background";
			ds_nomec.finalized = true;	

			ds_name = "mb_nce_both/ds/mb_nce_" + ma + "_" + ds + "_part_mec";
		
			mhist1D ds_mec (ds_name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 6, 0);	
			ds_mec.cnames[0] = "total";
			ds_mec.cnames[1] = "hydrogen";
			ds_mec.cnames[2] = "qel";
			ds_mec.cnames[3] = "res/dis";
			ds_mec.cnames[4] = "mec";
			ds_mec.cnames[5] = "background";
			ds_mec.finalized = true;	
			
			ds_name = "mb_nce_both/ds/mb_nce_" + ma + "_" + ds + "_part_mec2";
		
			mhist1D ds_mec2 (ds_name, "data/data/mb_ncel_ratio.txt", " ", "Reconstructed kinetic energy", "{/Symbol n}p -> {/Symbol n}p / {/Symbol n} N -> {/Symbol n}N", 30, 350, 800, 6, 0);	
			ds_mec2.cnames[0] = "total";
			ds_mec2.cnames[1] = "hydrogen";
			ds_mec2.cnames[2] = "qel";
			ds_mec2.cnames[3] = "res/dis";
			ds_mec2.cnames[4] = "mec";
			ds_mec2.cnames[5] = "background";
			ds_mec2.finalized = true;	
			
			mhist1D ma_final (ma_name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 3, 0);
			ma_final.cnames[0] = "without MEC";
			ma_final.cnames[1] = "with MEC";
			ma_final.cnames[2] = "with MEC2";
			ma_final.finalized = true;	
		
			ma_name = "mb_nce_both/ma/mb_nce_" + ma + "_" + ds + "_part_nomec";
		
			mhist1D ma_part_nomec (ma_name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 5, 0);
			ma_part_nomec.cnames[0] = "total";
			ma_part_nomec.cnames[1] = "hydrogen";
			ma_part_nomec.cnames[2] = "qel";
			ma_part_nomec.cnames[3] = "res/dis";
			ma_part_nomec.cnames[4] = "background";
			ma_part_nomec.finalized = true;	
					
			ma_name = "mb_nce_both/ma/mb_nce_" + ma + "_" + ds + "_part_mec1";
		
			mhist1D ma_part_mec1 (ma_name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 6, 0);
			ma_part_mec1.cnames[0] = "total";
			ma_part_mec1.cnames[1] = "hydrogen";
			ma_part_mec1.cnames[2] = "qel";
			ma_part_mec1.cnames[3] = "res/dis";
			ma_part_mec1.cnames[4] = "mec";
			ma_part_mec1.cnames[5] = "background";
			ma_part_mec1.finalized = true;	
		
			ma_name = "mb_nce_both/ma/mb_nce_" + ma + "_" + ds + "_part_mec2";
		
			mhist1D ma_part_mec2 (name, "data/data/mb_ncel.txt", "", "Reconstructed kinetic energy", "No. of events", 51, 40, 650, 6, 0);
			ma_part_mec2.cnames[0] = "total";
			ma_part_mec2.cnames[1] = "hydrogen";
			ma_part_mec2.cnames[2] = "qel";
			ma_part_mec2.cnames[3] = "res/dis";
			ma_part_mec2.cnames[4] = "mec";
			ma_part_mec2.cnames[5] = "background";
			ma_part_mec2.finalized = true;	
					
			load (ma_final, ma_part_nomec, ma_part_mec1, ma_part_mec2, ma, ds);
			load_ratio (ds_final, ds_nomec, ds_mec, ds_mec2, ma, ds);
			
			//if(calc_chi (ma_final.result [0]) + calc_chi_ratio (ds_final.result[0]) < 70)
				chi.result[i][j] = calc_chi (ma_final.result [0]) + calc_chi_ratio (ds_final.result[0]);
			//else
			//	chi.result[i][j] = 1e6;
				
			//if (calc_chi_ratio (ds_final.result[0]) < 40)
				chir.result[i][j] = calc_chi_ratio (ds_final.result[0]);
			//else
			//	chir.result[i][j] = 1e6;
			
			//if (calc_chi (ma_final.result [0]) < 40)
				chid.result[i][j] = calc_chi (ma_final.result [0]);
			//else
			//	chid.result[i][j] = 1e6;
				
			//if (calc_chi (ma_final.result [1]) + calc_chi_ratio (ds_final.result[1]) < 70)
				chi_mec.result[i][j] = calc_chi (ma_final.result [1]) + calc_chi_ratio (ds_final.result[1]);
			//else
			//	chi_mec.result[i][j] = 1e6;
			
			//if (calc_chi_ratio (ds_final.result[1]) < 40)
				chi_mecr.result[i][j] = calc_chi_ratio (ds_final.result[1]);
			//else
			//	chi_mecr.result[i][j] = 1e6;
			
			//if (calc_chi (ma_final.result [1]) < 40)
				chi_mecd.result[i][j] = calc_chi (ma_final.result [1]);
			//else
			//	chi_mecd.result[i][j] = 1e6;
			
			//chi_mec2.result[i][j] = calc_chi (ma_final.result [2]) + calc_chi_ratio (ds_final.result[2]);
		}
	}
	
	//make_smooth (chi);
	//make_smooth (chir);
	//make_smooth (chid);
	//make_smooth (chi_mec);
	//make_smooth (chi_mecr);
	//make_smooth (chi_mecd);
		
//	double min_ma, min_ds, min;
	
//	find_min (chi_mec_smooth1, min_ma, min_ds, min);
	
//	expand (chi_mec_smooth1, chi_mec_smooth_expand1, 0);
//	expand (chi_mec_smooth_expand1, chi_mec_smooth_expand2, 1);
		
	ofstream file("analysis/mb_nce_both/minimum.txt");
	double a,b,c;
	
	find_minimum (chi, a, b, c);
	file << chi.name << ": " << a << " " << b << " " << c << endl;
	find_minimum (chir, a, b, c);
	file << chir.name << ": " << a << " " << b << " " << c << endl;
	find_minimum (chid, a, b, c);
	file << chid.name << ": " << a << " " << b << " " << c << endl;
	find_minimum (chi_mec, a, b, c);
	file << chi_mec.name << ": " << a << " " << b << " " << c << endl;
	find_minimum (chi_mecr, a, b, c);
	file << chi_mecr.name << ": " << a << " " << b << " " << c << endl;
	find_minimum (chi_mecd, a, b, c);
	file << chi_mecd.name << ": " << a << " " << b << " " << c << endl;
	
	file.close();
	
	//for (int i = 0; i < mb_nce_ma_bins; i++)
	//	for (int j = 0; j < mb_nce_ds_bins; j++)
	//		chi_mec_ratio.result[i][j] = chi_mec.result[i][j] / chi_mec_smooth1.result[i][j];
	
	//for (int i = 0; i < chi_x; i++)
	//	for (int j = 0; j < chi_y; j++)
	//		temp_file << mb_nce_start + i * mb_nce_ma_step << " " << mb_nce_ds_start + j * mb_nce_ds_step << " " << chi2[i][j] << endl;
			
	//temp_file.close();	
}

void mb_nce_both :: load (mhist1D& h, mhist1D& part_nomec, mhist1D& part_mec1, mhist1D& part_mec2, string ma, string ds)
{
	double H[51] = {0};
	double RES[51] = {0};
	double MEC[51] = {0};
	double QEL[51] = {0};
	double HQEL[51] = {0};
	
	load (H, "Hqel_1030_" + ds + ".txt");
	//load (H, "qel_hydrogen.txt");
	load (RES, "resdis.txt");
	load (MEC, "mecfz.txt");
	load (QEL, "qel_" + ma + "_" + ds + ".txt");
	//load (HQEL, "Hqel_" + ma + "_" + ds + ".txt");
	
	for (int i = 0; i < 51; i++)
	{
		h.result[0][i] += H[i] + RES[i] + QEL[i] + mb_nce_bg[i];
		h.result[1][i] += H[i] + RES[i] + QEL[i] + MEC[i] + mb_nce_bg[i];
		h.result[2][i] += HQEL[i] + RES[i] + QEL[i] + + MEC[i] + mb_nce_bg[i];
		
		part_nomec.result[0][i] += H[i] + RES[i] + QEL[i] + mb_nce_bg[i];
		part_nomec.result[1][i] += H[i];
		part_nomec.result[2][i] += QEL[i];
		part_nomec.result[3][i] += RES[i];
		part_nomec.result[4][i] += mb_nce_bg[i];

		part_mec1.result[0][i] += H[i] + RES[i] + QEL[i] + MEC[i] + mb_nce_bg[i];
		part_mec1.result[1][i] += H[i];
		part_mec1.result[2][i] += QEL[i];
		part_mec1.result[3][i] += RES[i];
		part_mec1.result[4][i] += MEC[i];
		part_mec1.result[5][i] += mb_nce_bg[i];

		part_mec2.result[0][i] += HQEL[i] + RES[i] + QEL[i] + + MEC[i] + mb_nce_bg[i];
		part_mec2.result[1][i] += HQEL[i];
		part_mec2.result[2][i] += QEL[i];
		part_mec2.result[3][i] += RES[i];
		part_mec2.result[4][i] += MEC[i];
		part_mec2.result[5][i] += mb_nce_bg[i];
	}
}

void mb_nce_both :: load (double* x, string fname)
{
	ifstream in (("analysis/mb_nce_both/temp/ma/" + fname).c_str());
	int i = 0;
	
	if (in)
	{		
		while (i < 51)
			in >> x[i++];
	}
	else
		cout << "[Ma] Nie ma pliku: " << fname << endl;
	
	in.close();
}

void mb_nce_both :: load_ratio (mhist1D& h, mhist1D& nomec, mhist1D& mec, mhist1D& mec2, string ma, string ds)
{
	double H1[30] = {0};
	double RES1[30] = {0};
	double MEC1[30] = {0};
	double QEL1[30] = {0};
	double HQEL1[30] = {0};

	double H2[30] = {0};
	double RES2[30] = {0};
	double MEC2[30] = {0};
	double QEL2[30] = {0};
	double HQEL2[30] = {0};
	
	load (H1, H2, "Hqel_1030_" + ds + ".txt");
	//load (H1, H2, "qel_hydrogen.txt");
	load (RES1, RES2, "resdis.txt");
	load (MEC1, MEC2, "mecfz.txt");
	load (QEL1, QEL2, "qel_" + ma + "_" + ds + ".txt");
	//load (HQEL1, HQEL2, "Hqel_" + ma + "_" + ds + ".txt");
		
	for (int i = 0; i < 30; i++)
	{
		h.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		h.result[1][i] += (H1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		h.result[2][i] += (HQEL1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		
		nomec.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[1][i] += H1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[2][i] += QEL1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[3][i] += RES1[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		nomec.result[4][i] += mb_nce_bg_he_sp[i] / (H2[i] + RES2[i] + QEL2[i] + mb_nce_bg_he_all[i]);
		
		mec.result[0][i] += (H1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[1][i] += H1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[2][i] += QEL1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[3][i] += RES1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[4][i] += MEC1[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec.result[5][i] += mb_nce_bg_he_sp[i] / (H2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		
		mec2.result[0][i] += (HQEL1[i] + RES1[i] + QEL1[i] + MEC1[i] + mb_nce_bg_he_sp[i]) / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec2.result[1][i] += HQEL1[i] / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec2.result[2][i] += QEL1[i] / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec2.result[3][i] += RES1[i] / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec2.result[4][i] += MEC1[i] / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
		mec2.result[5][i] += mb_nce_bg_he_sp[i] / (HQEL2[i] + RES2[i] + QEL2[i] + MEC2[i] + mb_nce_bg_he_all[i]);
	}
}

void mb_nce_both :: load (double* x, double *y, string fname)
{
	ifstream in (("analysis/mb_nce_both/temp/ds/" + fname).c_str());
	int i = 0;
	
	if (in)
	{	
		while (i < 30)
		{
			in >> x[i];
			in >> y[i++];
		}
	}
	else
		cout << "[ds] Nie ma pliku: " << fname << endl;
	
	in.close();
}

double mb_nce_both :: calc_chi (double *x)
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

double mb_nce_both :: calc_chi_ratio (double *x)
{
	double res = 0;
	
	double dif [30] = {0};
		
	for (int i = 0; i < 30; i++)
		dif[i] = mb_nce_he_sp[i] / mb_nce_he_all[i] - x[i];
				
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 30; j++)
			res += dif[i] * Mrev_he [i][j] * dif[j];	

	return res;	
}

//end of mb ncel analysis		

void phil1a :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 1;
	P.nucleus_n = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1a :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 100)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil1b :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 0;
	P.nucleus_n = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1b :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 100)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil1c :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 0;
	P.nucleus_n = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1c :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 1)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil1d :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 1;
	P.nucleus_n = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1d :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 1)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil1e :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 0;
	P.nucleus_n = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1e :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 1)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil1f :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	
	P.nucleus_p = 0;
	P.nucleus_n = 1;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 1;
	P.dyn_dis_nc = 1;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil1f :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 10)
		h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phil2 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "600";
	
	P.nucleus_p = 1;
	P.nucleus_n = 0;
	
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 1;
	P.dyn_dis_cc = 1;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;	
}

void phil2 :: calculate (event *e)
{
	int pion = 100 * e -> nof (211) + 10 * e -> nof (-211) + e -> nof (111);
	
	if (pion == 100)
	{
		double Q2 = - e -> q2();
				
		double mom = 0;
		
		for (int i = 0; i < e -> out.size(); i++)
			if (e -> out[i].pdg == 211)
				mom = e -> out[i].momentum ();
				
		h2 -> put (Q2, e -> dyn, e -> weight);
		h1 -> put (mom, e -> dyn, e -> weight);
	}
}

void mec_cher :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/C.txt");
	
	P.qel_nc_axial_mass = 1030;
	
	P.mec_kind = 1;
			
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 1;	
}

void mec_cher :: calculate (event *e)
{
	double tk0 = 0;
	double tk1 = 0;
	
	for (int i = 0; i < e -> out.size(); i++)
		if (e -> out[i].nucleon() and e -> out[i].Ek () > tk0)
			tk0 = e -> out[i].Ek ();

	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].nucleon() and e -> post[i].Ek () > tk1)
			tk1 = e -> post[i].Ek ();
			
	h1 -> put (tk0, e->dyn, e->weight, 0);
	h1 -> put (tk1, e->dyn, e->weight, 1);
}

void mb_nce_test :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/C.txt");
					
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

void mb_nce_test2 :: set_params ()
{
	P.read("data/beam/newMB.txt");
	P.read("data/target/C.txt");
					
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
	
	P.mec_kind = 4;	
}

void mb_nce_test :: calculate (event *e)
{
	mianownik += e -> weight;
	
	if (e -> number_of_interactions() == 0)
		licznik += e -> weight;
}

void mb_nce_test2 :: calculate (event *e)
{
	double total = 0;
	double max = 0;
	
	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].pdg == pdg_proton or e -> post[i].pdg == pdg_neutron)
		{
			double E = e -> post[i].Ek();
			
			total += E;
			
			if (E > max)
				max = E;
		}
		
	h1 -> put (total, e -> dyn, e -> weight);
	h2 -> put (max, e -> dyn, e -> weight);
}

void tem_nieves_test :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "1000";
	P.read("data/target/C.txt");
					
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
	
	P.mec_kind = 3;	
}

void tem_nieves_test :: calculate (event *e)
{
	h1 -> put (e->out[0].p().z / e->out[0].momentum(), e->out[0].Ek(), e -> dyn, e -> weight);
}
void energy_test :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "800";
	P.read("data/target/C.txt");
					
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
	
	P.mec_kind = 1;	
	
	P.FSI_on = 1;
}

void niwg_nieves :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "500 1500";
	P.read("data/target/O.txt");
					
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
	
	P.mec_kind = 3;	
	
	P.FSI_on = 1;
}

void niwg_nieves :: calculate (event *e)
{
	double Ev = e->in[0].E();
	double mom = 0;
	
	for (int i = 0; i < e->post.size(); i++)
		if (e->post[i].pdg == 2212 and e->post[i].momentum() > mom)
			mom = e->post[i].momentum();
	
	h1norm -> put (Ev, e->dyn, 1);
			
	if (mom >= 1060)
		h1 -> put (Ev, e->dyn, 1);
		
	if (e->out[0].Ek() > 400)
	{
		h5norm -> put (Ev, e->dyn, 1);
			
		if (mom >= 1060)
			h5 -> put (Ev, e->dyn, 1);
	}
	else
	{
		h6norm -> put (Ev, e->dyn, 1);
			
		if (mom >= 1060)
			h6 -> put (Ev, e->dyn, 1);
	}
	
	if (Ev >= 500 and Ev <= 600)
		h2 -> put (mom, e->dyn, e->weight);
	else if (Ev >= 900 and Ev <= 1000)
		h3 -> put (mom, e->dyn, e->weight);
	else if (Ev >= 1400 and Ev <= 1500)
		h4 -> put (mom, e->dyn, e->weight);
}

void niwg_nieves2 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "1000";
	P.read("data/target/O.txt");
					
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
	
	P.mec_kind = 3;	
	
	P.FSI_on = 1;
}

void niwg_nieves2 :: calculate (event *e)
{
	double Ek = e->out[0].Ek();
	double cos = e->out[0].p().z/e->out[0].momentum();

	h1 -> put (cos, Ek, e->dyn, e->weight);

	double mom = 0;
	
	for (int i = 0; i < e->post.size(); i++)
		if (e->post[i].pdg == 2212 and e->post[i].momentum() > mom)
			mom = e->post[i].momentum();
			
	if (mom >= 1060)
		h2 -> put (cos, Ek, e->dyn, e->weight);
}

void energy_test :: calculate (event *e)
{	
	if (e -> fof(211) + e -> fof(-211) + e-> fof(111) == 0 and e-> number_of_interactions() == 0)
	{
		double energy0 = 0;
		double energy1 = 0;
	
		for (int i = 0; i < e -> out.size(); i++)
			energy0 += e->out[i].Ek();
			
		for (int i = 0; i < e -> post.size(); i++)
			energy1 += e->post[i].Ek();
		
		h1 -> put (energy0, e -> dyn, e -> weight, 0);
		h1 -> put (energy1, e -> dyn, e -> weight, 1);
	}
}

void bodek :: calculate (event *e)
{
	h2 -> put (e->in[1].momentum(), e->dyn, e->weight);
	
	using namespace PDG;
	
	double x = e->q0() + e->q2() / 2.0 / mass_proton;
	x /= 1000.0;
	int i;
	
	double Q2 = -e->q2() / 1000000.0;
	
	if (Q2 >= 0.05 and Q2 < 0.15)
		i = 0;
	else if (Q2 >= 0.15 and Q2 < 0.45)
		i = 1;
	else if (Q2 >= 0.45 and Q2 < 0.55)
		i = 2;
	else if (Q2 >= 0.55 and Q2 < 0.85)
		i = 3;
	else if (Q2 >= 0.85 and Q2 < 1.15)
		i = 4;
	else if (Q2 >= 1.15 and Q2 < 1.25)
		i = 5;
	else if (Q2 >= 1.25 and Q2 < 1.75)
		i = 6;
	else if (Q2 >= 1.75 and Q2 < 2.25)
		i = 7;
	else
		return;
		
	h1 -> put (x, e -> dyn, e -> weight, i);	
}

void bodek :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "10000";
	P.read("data/target/C.txt");
					
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
		
	P.FSI_on = 0;
	
	P.nucleus_target = 2; //1 - gfg, 2 - lfg
	P.sf_method = 1;
}

void phd1 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "1000";
	P.read("data/target/C.txt");
					
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	
	P.nucleus_target = 1; //1 - gfg, 2 - lfg
	P.sf_method = 0;
	P.sf_pb = 0;
}

void phd1 :: calculate (event *e)
{
	h1 -> put (e->out[1].momentum(), e->dyn, e->weight);
	h2 -> put (-e->q2()/1000000.0, e->dyn, e->weight);
}

void phd2 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "200 2000";
	P.read("data/target/C.txt");
					
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	P.pauli_blocking = 0;
	P.nucleus_target = 1; //1 - gfg, 2 - lfg
	P.sf_method = 0;
	P.sf_pb = 0;
}

void phd2 :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phd3 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "0 1500";
	P.read("data/target/C.txt");
	P.qel_cc_axial_mass = 1030;				
	P.dyn_qel_cc = 1;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 0;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	
	P.nucleus_target = 2; //1 - gfg, 2 - lfg
	P.sf_method = 0;
}

void phd3 :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	h1 -> put (e -> in[0].E(), e -> dyn, e -> weight);
}

void phd5 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "0 200000";

	P.read("data/target/C.txt");
	
	P.qel_cc_axial_mass = 1030;	
	
	P.coh_new = 1;
				
	P.dyn_qel_cc = 0;
	P.dyn_res_cc = 0;
	P.dyn_dis_cc = 0;
	P.dyn_coh_cc = 1;
	P.dyn_mec_cc = 0;

	P.dyn_qel_nc = 0;
	P.dyn_res_nc = 0;
	P.dyn_dis_nc = 0;
	P.dyn_coh_nc = 0;
	P.dyn_mec_nc = 0;
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	
	P.nucleus_target = 2; //1 - gfg, 2 - lfg
	P.sf_method = 0;
}

void phd5 :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E()/1000.0, e -> dyn);
	h1 -> put (e -> in[0].E()/1000.0, e -> dyn, e -> weight);
}

void phd4 :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "1000";
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

void phd4 :: calculate (event *e)
{
	double tk0 = 0;
	double tk1 = 0;
	double cos0 = 0;
	double cos1 = 0;
	
	for (int i = 0; i < e -> out.size(); i++)
		if (e -> out[i].nucleon() and e -> out[i].Ek () > tk0)
		{
			tk0 = e -> out[i].Ek ();
			cos0 = e -> out[i].p().z / e -> out[i].momentum();
		}

	for (int i = 0; i < e -> post.size(); i++)
		if (e -> post[i].nucleon() and e -> post[i].Ek () > tk1)
		{
			tk1 = e -> post[i].Ek ();
			cos1 = e -> post[i].p().z / e -> post[i].momentum();
		}
			
	h1 -> put (tk0, e->dyn, e->weight, 0);
	h1 -> put (tk1, e->dyn, e->weight, 1);
	h2 -> put (cos0, e->dyn, e->weight, 0);
	h2 -> put (cos1, e->dyn, e->weight, 1);
}

void phd6 :: set_params ()
{
	P.beam_particle = PDG::pdg_proton;
	P.beam_type = 0;
	P.beam_energy = intToStr(Tk[0]*1000 + PDG::mass_proton);
	P.read("data/target/Fe.txt");
		
	P.mec_kind = 1;
				
	P.formation_zone = "nofz";
	P.beam_placement = 1;
	P.first_step = 1;
}

void phd6 :: calculate (event *e)
{		
	if (e -> number_of_interactions() == 0)
		T[count]++;
}

bool phd6 :: change_params()
{
	count++;
	
	if (count < N)
	{
		P.beam_energy = intToStr(Tk[count]*1000 + PDG::mass_proton);
		return true;
	}
	else
		return false;
}

void phd7 :: set_params ()
{
	P.beam_particle = PDG::pdg_piP;
	P.beam_type = 0;
	P.beam_energy = intToStr(sqrt(mom[0]*mom[0]*1e6 + PDG::mass_piP*PDG::mass_piP));
	P.read("data/target/C.txt");
		
	P.mec_kind = 1;
				
	P.formation_zone = "trans";
	P.beam_placement = 1;
	P.first_step = 1;
}

void phd7 :: calculate (event *e)
{		
	if (e -> number_of_interactions() == 0)
		T[count]++;
}

bool phd7 :: change_params()
{
	count++;
	
	if (count < N)
	{
		P.beam_energy = intToStr(sqrt(mom[count]*mom[count]*1e6 + PDG::mass_piP*PDG::mass_piP));
		return true;
	}
	else
		return false;
}

void mec_bu :: set_params ()
{
	P.beam_particle = PDG::pdg_nu_mu;
	P.beam_type = 0;
	P.beam_energy = "0 2000";
	P.read("data/target/C.txt");
					
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
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	
	P.nucleus_target = 2; //1 - gfg, 2 - lfg
	P.sf_method = 0;
}

void mec_bu :: calculate (event *e)
{
	h1 -> histogram :: put (e -> in[0].E(), e -> dyn);
	h2 -> histogram :: put (e -> in[0].E(), e -> dyn);

	h1 -> put (e -> in[0].E()/1000.0, e -> dyn, e -> weight);

	if (e->out[1].p().z < 0)
		h2 -> put (e -> in[0].E()/1000.0, e -> dyn, e -> weight);
}

void mec_bu2 :: set_params ()
{
	P.read("data/beam/T2Knumu.txt");
	P.read("data/target/C.txt");
					
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
	
	P.mec_kind = 3;	
	
	P.FSI_on = 0;
	
	P.nucleus_target = 2; //1 - gfg, 2 - lfg
	P.sf_method = 0;
}

void mec_bu2 :: calculate (event *e)
{
	h1 -> put (e -> out[1].Ek()/1000.0, e -> dyn, e -> weight);

	if (e->out[1].p().z < 0)
		h2 -> put (e -> out[1].Ek()/1000.0, e -> dyn, e -> weight);
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
		case 34: wsk = new mb_nce_he; break;
		case 35: wsk = new phil1a; break;
		case 36: wsk = new phil1b; break;
		case 37: wsk = new phil1c; break;
		case 38: wsk = new phil1d; break;
		case 39: wsk = new phil1e; break;
		case 40: wsk = new phil1f; break;
		case 41: wsk = new phil2; break;
		case 42: wsk = new mec_cher; break;
		case 43: wsk = new mb_nce_both; break;
		case 44: wsk = new mb_nce_test; break;
		case 45: wsk = new mb_nce_test2; break;
		case 46: wsk = new tem_nieves_test; break;
		case 47: wsk = new energy_test; break;
		case 48: wsk = new niwg_nieves; break;
		case 49: wsk = new niwg_nieves2; break;
		case 50: wsk = new bodek; break;
		case 51: wsk = new phd1; break;
		case 52: wsk = new phd2; break;
		case 53: wsk = new mec_bu; break;
		case 54: wsk = new mec_bu2; break;
		case 55: wsk = new phd3; break;
		case 56: wsk = new phd4; break;
		case 57: wsk = new phd5; break;
		case 58: wsk = new phd6; break;
		case 59: wsk = new phd7; break;
		default: wsk = NULL;
	}
	
	return wsk;
}

int main (int argc, char **argv)
{
	set_dirs(argv[0]);
	init_genrand(time(NULL));
	
	if (argc > 1)
	{
		mb_nce_mode = atoi(argv[1]);
		
		mb_nce_start = atof(argv[2]);
		mb_nce_nof = atoi(argv[3]);
		
		if (argc == 5)
			mb_nce_he_ma = atoi(argv[4]);
			
		mb_nce_ma = atoi(argv[2]);
		mb_nce_ds = atof(argv[3]) / 100.0;
		
		if (argc == 8)
		{
			mb_nce_ma_step = atoi(argv[4]);
			mb_nce_ma_start = atoi(argv[2]);
			mb_nce_ma_end = atoi(argv[3]);
			mb_nce_ds_step = atof(argv[7]);
			mb_nce_ds_start = atof(argv[5]);
			mb_nce_ds_end = atof(argv[6]);
		}
	}
		
	run_command ("mkdir -p analysis/");
			
	for (int i = 0; i < nof_class; i++)
	{
		if (!active_class[i])
			continue;
			
		pattern *wsk = choose(i);
		
		if (wsk) 
		{
			if (wsk->fsi())
				wsk->startfsi();
			else
			wsk -> start();
			
			delete wsk;
		}
	}
			
	return 0;
}

