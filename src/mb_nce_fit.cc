#include "mb_nce_fit.h"

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

void read_distrHE(string filename, double tab[5][30])
{
	ifstream Input (filename.c_str());
	
	if (Input)
	{
		do
		{
			for (int j = 0; j < 5; j++)
				for (int k = 0; k < 30; k++)
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

void read_rtabHE(string filename, double tab[30][30])
{
	ifstream Input (filename.c_str());
	
	if (Input)
	{
		do
		{
			for (int j = 0; j < 30; j++)
				for (int k = 0; k < 30; k++)
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

void makeRHE(double *mu, double M[30][30], double R[30][30])
{
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 30; j++)
			if (mu[i] != 0) R[i][j] = M[i][j] / mu[i];
			else R[i][j] = 0;
}

void trans (double *in, double T[][51], double *out)
{
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			out[i] += in[j]*T[j][i];
}

void transHE (double *in, double T[][30], double *out)
{
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 30; j++)
			out[i] += in[j]*T[j][i];
}

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

void true2recp(double Et[5][30], double *rec)
{
	double M[5][30][30];
	
	double mu[5][30];
	
	read_distrHE("mb_nce/tables/neutHE.txt", mu); //wczytuje rozkład energii true

	string tables[5] = {"mb_nce/tables/rp1.txt", "mb_nce/tables/rp2.txt", "mb_nce/tables/rp3.txt", "mb_nce/tables/rp4.txt", "mb_nce/tables/rp5.txt"};
		
	for (int i = 0; i < 5; i++)
		read_rtabHE(tables[i], M[i]); //wczytuje tabele M 
	
	double R[5][30][30];
	
	for (int i = 0; i < 5; i++)
		makeRHE(mu[i], M[i], R[i]); //tworzy macierze R
	
	double help[5][30] = {{0}};
		
	for (int i = 0; i < 5; i++)
		transHE(Et[i], R[i], help[i]); //przejscie true do rec dla kanałów 1 - 5
		
	for (int i = 0; i < 30; i++)
		rec[i] = bgHEp[i];
	
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 5; j++)
			rec[i] += help[j][i];
}

void true2recN(double Et[5][30], double *rec)
{
	double M[5][30][30];
	
	double mu[5][30];
	
	read_distrHE("mb_nce/tables/neutHE.txt", mu); //wczytuje rozkład energii true

	string tables[5] = {"mb_nce/tables/rN1.txt", "mb_nce/tables/rN2.txt", "mb_nce/tables/rN3.txt", "mb_nce/tables/rN4.txt", "mb_nce/tables/rN5.txt"};
		
	for (int i = 0; i < 5; i++)
		read_rtabHE(tables[i], M[i]); //wczytuje tabele M 
	
	double R[5][30][30];
	
	for (int i = 0; i < 5; i++)
		makeRHE(mu[i], M[i], R[i]); //tworzy macierze R
	
	double help[5][30] = {{0}};
		
	for (int i = 0; i < 5; i++)
		transHE(Et[i], R[i], help[i]); //przejscie true do rec dla kanałów 1 - 5
		
	for (int i = 0; i < 30; i++)
		rec[i] = bgHEN[i];
	
	for (int i = 0; i < 30; i++)
		for (int j = 0; j < 5; j++)
			rec[i] += help[j][i];
}
