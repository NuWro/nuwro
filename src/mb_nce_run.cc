#include "mb_nce.h"
#include "TMinuit.h"

int classic()
{		
	//sim();
	//re_sim();
	//calc();
	//chi2();
	
	delta_s_sim(1350);
	delta_s_calc();
	
	//q2_calc("mbbeam_C_fg_QEL_1250.root");
	
	//sim_cc();
	//ratio_calc("mbbeam_C_fg_QEL_1250.root");
		
	return 1;
}

void chi2(double &chiq, double *par)
{
	fstream plik;
	plik.open("mb_nce/results/tm_chi2.txt", ios::out|ios::app);
	
	double ma = par[0];
	
	tm_sim(ma);
	tm_calc(ma);
	chiq = tm_chi(ma);
	
	plik << ma << " " << chiq << endl;
	
	plik.close();
}

void fituj(double &chiq, double *par, bool flag)
{
	chi2(chiq, par);
}

void funkcja(int &npar, double *gin, double &f, double *par, int iflag)
{
	double chiq=0;
	fituj(chiq, par, false);
	f = chiq;
}

int minuit()
{
	ofstream plik("mb_nce/results/tm_chi2.txt");
	plik.close();	
		
	TMinuit *gMinuit = new TMinuit();
	gMinuit->SetFCN(funkcja);

	gMinuit->DefineParameter(0, "Ma", 1300.0, 50.0, 800.0, 2000.0);

	//gMinuit->Migrad();
	gMinuit->mnsimp();
		
	double ma;
	double maerr;
	
	gMinuit->GetParameter(0,ma,maerr);
	
	ofstream out("mb_nce/results/tm_chi2_min.txt");
	
	out << "Ma = " << ma << " +/- " << maerr << endl;

	out.close();

	delete gMinuit;
	
	return 1;
}

int main(int argc, char **argv)
{
	run("mkdir -p mb_nce/");
	run("mkdir -p mb_nce/root/");
	run("mkdir -p mb_nce/results/");
	
	return classic();
	//return minuit();
}

