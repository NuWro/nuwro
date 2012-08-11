#include "event1.h"
#include "particle.h"

#include <TFile.h>
#include <TTree.h>

#include <string>

using namespace std;

void sim();
void sim_cc();
void re_sim();
void calc();
void chi2();
double calc_chi(string in);
void calcH(string in, double *res);
void calcC(string in, double res[5][51]);
void norm(double *tab, double x);
void q2_calc(string rootC);
void ratio_calc(string rootCnc);
void q2_calc(string in, double *res);
void ratio_calcH(string in, double *res);
void ratio_calcC(string in, double *res, double *reslike);

void run(string com)
{
	FILE *fp;
	fp = popen(com.c_str(), "w");
	pclose(fp);
}

double crosssection (string filename)
{
	string y;
	double x;
		
	vector<double> v(8);
	ifstream Input (filename.c_str());

	if (Input)
	{
		do
		{
			getline (Input, y);
			for (int j = 0; j < 8; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						Input>>y;
					}
					Input>>v[j];
				}
		}while (Input);
	}
	
	double result = 0;
	
	for (int i = 0; i < 8; i++)
	{
		result += v[i];
	}
	
	return result;
}

const string parQEL = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 2000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string parRES = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 2000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string parCC = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 2000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string H = "-p '@data/target/H.txt' ";
const string C = "-p '@data/target/C.txt' ";
const string fg = "-p 'sf_method = 0' ";
const string sf = "-p 'sf_method = 1' ";
const string ccma = "-p 'qel_cc_axial_mass = 1350' ";

string ncma(int val)
{
	stringstream temp;
	string x;
	
	temp << val;
	temp >> x;
	
	string ncma = "-p 'qel_nc_axial_mass = " + x + "' -p 'qel_s_axial_mass = " + x + "' ";
	return ncma;
}

const int events = 2000000;

const double flux = 5.22227e-10;
const double POT = 6.46165e20;
const double R = 610.6;
const double rho = 0.845;
const double deltaT = 18;
const double Na = 6.02214e23;
const double Nn = Na * rho * 4 * M_PI * R * R * R / 3;

void read_distr(string filename, double tab[5][51])
{
	ifstream Input (filename.c_str());
	
	if (Input)
	{
		do
		{
			for (int j = 0; j < 5; j++)
				for (int k = 0; k < 51; k++)
					Input >> tab[j][k];
		}while (Input);
	}
	else cout << "There is no tab file" << endl;
}

void read_rtab(string filename, double tab[51][51])
{
	ifstream Input (filename.c_str());
	
	if (Input)
	{
		do
		{
			for (int j = 0; j < 51; j++)
				for (int k = 0; k < 51; k++)
					Input >> tab[k][j];

		}while (Input);
	}
	else cout << "There is no tab file" << endl;
}

void read_Erec(string filename, double *tab)
{
	ifstream Input (filename.c_str());
	int i = 0;
	if (Input)
	{
		do
		{
			string help;
			Input >> help;
			Input >> tab[i];
			Input >> help;
			Input >> help;
			i++;
		}while (Input);
	}
	else cout << "There is no tab file" << endl;
}

void makeR(double *mu, double M[51][51], double R[51][51])
{
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			if (mu[i] != 0) R[i][j] = M[i][j] / mu[i];
			else R[i][j] = 0;
}

void trans (double *in, double T[][51], double *out)
{
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			out[i] += in[j]*T[j][i];
}

const double ibg[51] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.571555e-01, 1.364424e+00, 4.705236e+00, 1.758978e+01, 3.496415e+01, 6.816501e+01, 1.368511e+02, 2.066710e+02, 3.374507e+02, 5.005924e+02, 6.672170e+02, 8.545720e+02, 1.014480e+03, 1.165838e+03, 1.241060e+03, 1.300197e+03, 1.349314e+03, 1.328437e+03, 1.331230e+03, 1.305422e+03, 1.274116e+03, 1.217948e+03, 1.147513e+03, 1.085130e+03, 1.029805e+03, 1.012486e+03, 9.418484e+02, 8.868173e+02, 8.388830e+02, 7.663445e+02, 7.245175e+02, 6.952694e+02, 6.376685e+02, 6.139051e+02, 5.547155e+02, 5.363876e+02, 4.934474e+02, 4.632671e+02, 4.293046e+02, 4.191864e+02, 3.713735e+02, 3.634909e+02, 3.295256e+02, 3.028509e+02, 4.277227e+03};
const double bg[51] = {1.468132e+01,  1.751390e+02,  5.432570e+02,  6.452100e+02,  6.288804e+02,  5.888983e+02,  5.721696e+02,  5.599448e+02,  5.320706e+02,  5.089781e+02,  4.717286e+02,  4.501799e+02,  4.296122e+02,  3.620738e+02,  3.576632e+02,  3.276326e+02,  3.137620e+02,  2.857147e+02,  2.655224e+02,  2.647242e+02,  2.455104e+02,  2.285900e+02,  2.313580e+02,  2.262919e+02,  2.283181e+02,  2.124414e+02,  1.952493e+02,  1.962794e+02,  1.876680e+02,  2.057413e+02,  2.050090e+02,  1.850047e+02,  1.591487e+02,  1.803423e+02,  1.927868e+02,  2.144821e+02,  2.260172e+02,  2.518160e+02,  2.556633e+02,  2.555962e+02,  2.480172e+02,  2.373849e+02,  1.840215e+02,  1.800383e+02,  1.857559e+02,  2.038884e+02,  2.060577e+02,  2.033785e+02,  2.086383e+02,  2.217605e+02,  2.130483e+02};
const double data[51] = {4.000000e+01, 4.950000e+02, 1.902000e+03, 2.722000e+03, 2.908000e+03, 2.878000e+03, 2.934000e+03, 3.003000e+03, 2.892000e+03, 2.720000e+03, 2.651000e+03, 2.656000e+03, 2.633000e+03, 2.409000e+03, 2.364000e+03, 2.250000e+03, 2.104000e+03, 2.070000e+03, 1.929000e+03, 1.832000e+03, 1.762000e+03, 1.618000e+03, 1.546000e+03, 1.458000e+03, 1.371000e+03, 1.203000e+03, 1.174000e+03, 1.055000e+03, 9.930000e+02, 9.080000e+02, 8.320000e+02, 7.570000e+02, 6.840000e+02, 6.330000e+02, 6.520000e+02, 7.090000e+02, 7.010000e+02, 6.710000e+02, 5.890000e+02, 5.110000e+02, 4.890000e+02, 3.920000e+02, 3.590000e+02, 3.320000e+02, 3.450000e+02, 3.430000e+02, 3.560000e+02, 3.410000e+02, 3.150000e+02, 3.240000e+02, 3.360000e+02};

void true2rec(double Et[5][51], double *rec)
{
	double M[5][51][51];
	
	double mu[5][51];
	
	read_distr("mb_nce/tables/neut.txt", mu); //wczytuje rozkład energii true

	string tables[5] = {"mb_nce/tables/m1.txt", "mb_nce/tables/m2.txt", "mb_nce/tables/m3.txt", "mb_nce/tables/m4.txt", "mb_nce/tables/m5.txt"};
		
	for (int i = 0; i < 5; i++)
		read_rtab(tables[i], M[i]); //wczytuje tabele M 
	
	double R[5][51][51];
	
	for (int i = 0; i < 5; i++)
		makeR(mu[i], M[i], R[i]); //tworzy macierze R
	
	double help[5][51] = {{0}};
		
	for (int i = 0; i < 5; i++)
		trans(Et[i], R[i], help[i]); //przejscie true do rec dla kanałów 1 - 5
		
	for (int i = 0; i < 51; i++)
		rec[i] = bg[i];
	
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 5; j++)
			rec[i] += help[j][i];
}
