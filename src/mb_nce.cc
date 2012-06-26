#include "mb_nce.h"

//const string dir = "mb_nce/v0"; // v = 0
//const string dir = "mb_nce/v1"; // v
//const string dir = "mb_nce/v2"; // v -= 5
//const string dir = "mb_nce/v3"; // v + zmienne kf
//const string dir = "mb_nce/v4"; // v raz + praca wyjscia = 7
//const string dir = "mb_nce/v5"; // rav V raz + praca wyjscia = 0
//const string dir = "mb_nce/v6"; // czastka pamieta ef + zmienny potencjal + praca wyjscia = 0
//const string dir = "mb_nce/v7"; // -||- w = 7
//const string dir = "mb_nce/v8"; // czastka pamieta, srednia w kaskadzie, w = 7, ruch fermiego w kasadzie off
const string dir = "mb_nce/v9"; // -||- rf on
const string dir2 = "mb_nce/v9b"; //9*1.2

void sim()
{
	///run simulations for NC RES and DIS (used for background)
	
	string command1 = "./bin/nuwro -o '" + dir + "/root/mbbeam_C_fg_RESDIS.root' " + parRES + fg + ccma + ncma(1350) + C;
	//run(command1);
	
	string command2 = "./bin/nuwro -o '" + dir + "/root/mbbeam_H_fg_QEL_1030.root' " + parQEL + fg + ccma + ncma(1030) + H;
	//run(command2);
	
	///run simulations for NC EL for various value of Ma
	
	for (int i = 17; i <= 17; i++)
	{
		int ma = 1000 + i*50;
		
		stringstream temp;
		string x;
	
		temp << ma;
		temp >> x;
	
		string command = "./bin/nuwro -o '" + dir + "/root/mbbeam_C_fg_QEL_" + x + ".root' " + parQEL + fg + ccma + ncma(ma) + C;
		
		run(command);
	}
}

void sim_cc()
{
	///run simulations for CCQE (+RES+DIS)
	
	string command1 = "./bin/nuwro -o 'mb_nce/cc/root/mbbeam_C_fg.root' " + parCC + fg + ccma + C;
	run(command1);
	
	string command2 = "./bin/nuwro -o 'mb_nce/cc/root/mbbeam_H_fg.root' " + parCC + fg + ccma + H;
	run(command2);
}

void re_sim()
{
	///run simulations for NC RES and DIS (used for background)
	
	string command = "./bin/nuwro -o '" + dir + "/root/mbbeam_C_fg_RESDIS.root' -p 'kaskada_redo = 1'";
	run(command);
		
	///run simulations for NC EL for various value of Ma
	
	//for (int i = 0; i <= 65; i++)
	for (int i = 4; i <= 17; i++)
	{
		//int ma = 1000 + i*10;
		int ma = 1000 + i*50;
		
		stringstream temp;
		string x;
	
		temp << ma;
		temp >> x;
	
		string command = "./bin/nuwro -o '" + dir + "/root/mbbeam_C_fg_QEL_" + x + ".root' -p 'kaskada_redo = 1'";
		//run(command);
	}
	
}

void norm(double *tab, double x)
{	
	double factor = deltaT * Nn * POT * flux;

	for (int j = 0; j < 51; j++)
			tab[j] *= factor*x;
}

void calcH(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.qel == 1)
		{			
			double Tk = e1->out[1].Ek();
			int w = Tk/18;
			if (w > 50) w = 50;
			res[w]++;				
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void calcC(string in, double res[5][51])
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{		
		tt1->GetEntry(i);
		
		int piony = e1->fof(211) + e1->fof(-211) + e1->fof(111);
		int nukleony = e1->fof(2212) + e1->fof(2112);
		
		if (piony == 0 and nukleony != 0)
		{			
			double Tk = 0;
			
			for (int k = 0; k < e1->post.size(); k++)
				if (e1->post[k].pdg == 2112 or e1->post[k].pdg == 2212)
					Tk += e1->post[k].Ek();
						
			int a = Tk/18;
						
			if (a > 50)
				a = 50;
			else if (a < 0)
				continue;

			int b = 4;
			
			if (e1->flag.qel)
			{
				if (e1->in[1].pdg == 2112)
					b = 3;
				else if (e1->number_of_interactions() == 0)
					b = 1;
				else
					b = 2;
			}

			res[b][a]++;
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void calc()
{
	string rdC = dir + "/root/mbbeam_C_fg_RESDIS.root.fsi.root";
	string rdCtxt = "mb_nce/v0/root/mbbeam_C_fg_RESDIS.root.txt";

	string qelH = "mb_nce/v0/root/mbbeam_H_fg_QEL_1030.root";
	string qelHtxt = "mb_nce/v0/root/mbbeam_H_fg_QEL_1030.root.txt";

	double result[5][51] = {{0}};
	
	calcC(rdC, result);
	double resxsec = crosssection(rdCtxt);
	
	calcH(qelH, result[0]);
	double xsecH = crosssection(qelHtxt);

	for (int i = 0; i < 51; i++)
	{
		//result[4][i] = ibg[i];
		result[4][i] *= 1.2*resxsec/events/18.0;
		result[0][i] *= 1.2*xsecH/events/18.0;
	}
	
	norm(result[0], 2.0/14);
	norm(result[4], 12.0/14);
	
	for (int i = 4; i <= 15; i++)
	{
		int ma = 1000 + i*50;
		
		for (int i = 0; i < 51; i++)
		{
			result[1][i] = 0;
			result[2][i] = 0;
			result[3][i] = 0;
		}
		stringstream temp;
		string x;
	
		temp << ma;
		temp >> x;

		string qelC = dir + "/root/mbbeam_C_fg_QEL_" + x + ".root.fsi.root";
		string qelCtxt = "mb_nce/v0/root/mbbeam_C_fg_QEL_" + x + ".root.txt";
	
		string out = dir2 + "/results/fg_tk_true_" + x + ".txt";
		string out2 = dir2 + "/results/fg_tk_rec_" + x + ".txt";
	
		calcC(qelC, result);
		
		double xsecC = crosssection(qelCtxt);
		
		for (int i = 0; i < 51; i++)
		{
			result[1][i] *= 1.2*xsecC/events/18.0;
			result[2][i] *= 1.2*xsecC/events/18.0;
			result[3][i] *= 1.2*xsecC/events/18.0;
		}
		
		norm(result[1], 12.0/14);
		norm(result[2], 12.0/14);
		norm(result[3], 12.0/14);
		
		double result2[51];
		
		true2rec(result, result2);
		
		ofstream file(out.c_str());
		
		for (int i = 0; i < 51; i++)
		{
			file << i*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
			file << (i+1)*18 << " " << result[0][i] << " " << result[1][i] << " " << result[2][i] << " " << result[3][i] << " " << result[4][i] << endl;
		}
		
		file.close();

		ofstream file2(out2.c_str());
		
		for (int i = 0; i < 51; i++)
		{
			file2 << 40 + i*12 << " " << result2[i] << endl;
			file2 << 40 + (i+1)*12 << " " << result2[i] << endl;
		}
		
		file2.close();
	}
}

double calc_chi(string in)
{
	double M[51][51] = {{0}};
	double E[51] = {0};
		
	read_rtab("mb_nce/tables/Mrev.txt", M);
	read_Erec(in, E);
		
	double res = 0;
	
	double dif[51];
		
	for (int i = 0; i < 51; i++)
		dif[i] = data[i] - E[i];
				
	for (int i = 0; i < 51; i++)
		for (int j = 0; j < 51; j++)
			res += dif[i]*M[i][j]*dif[j];
	
	return res;
}

void chi2()
{
	string outchi = dir2 + "/results/chi2.txt";
	ofstream file(outchi.c_str());
	
	for (int i = 4; i <= 15; i++)
	{
		int ma = 1000 + i*50;
		
		stringstream temp;
		string x;
	
		temp << ma;
		temp >> x;
	
		string in = dir2 + "/results/fg_tk_rec_" + x + ".txt";
		
		double chi = calc_chi(in);
				
		file << ma-5 << " " << chi << endl;
		file << ma+5 << " " << chi << endl;
		cout << ma << " " << chi << endl;
	}
	
	file.close();

}

void q2_calc(string rootC)
{
	string rootH = "mb_nce/v0/root/mbbeam_H_fg_QEL_1050.root";
	rootC = dir + rootC;
	string txtH = rootH + ".txt";
	string txtC = rootC + ".txt";
	
	double H[50] = {0};
	double C[50] = {0};
	double CH2[50] = {0};
	
	q2_calc(rootH, H);
	q2_calc(rootC, C);
	
	double xsecH = crosssection(txtH);
	double xsecC = crosssection(txtC);
	
	for (int i = 0; i < 50; i++)
	{
		H[i] *= 2 * xsecH / events / 0.04;
		C[i] *= 12 * xsecC / events / 0.04;
		CH2[i] = (H[i] + C[i]) / 14;
	}
	string q2out = dir + "/results/q2_dist.txt";
	ofstream file(q2out.c_str());
	
	for (int i = 0; i < 50; i++)
		file << i*0.04 << " " << CH2[i]*1e40 << endl <<
				(i+1)*0.04 << " " << CH2[i]*1e40 << endl;
	
	file.close();
}

void q2_calc(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.qel == 1)
		{			
			double Q2 = -e1->q2()/1000000.0;
			int w = Q2/0.04;
			
			if (w < 50 and w >= 0)
				res[w]++;				
		}
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void ratio_calc(string rootCnc)
{
	string rootHnc = "mb_nce/v0/root/mbbeam_H_fg_QEL_1050.root";
	rootCnc = dir + rootCnc;
	
	string rootHcc = "mb_nce/cc/root/mbbeam_H_fg.root";
	string rootCcc = "mb_nce/cc/root/mbbeam_C_fg.root";
	
	string txtHnc = rootHnc + ".txt";
	string txtCnc = rootCnc + ".txt";
	string txtHcc = rootHcc + ".txt";
	string txtCcc = rootCcc + ".txt";
	
	double Hnc[50] = {0};
	double Cnc[50] = {0};
	double Cnclike[50] = {0};
	double CH2nc[50] = {0};
	double CH2nclike[50] = {0};

	double Hcc[50] = {0};
	double Ccc[50] = {0};
	double Ccclike[50] = {0};
	double CH2cc[50] = {0};
	double CH2cclike[50] = {0};
	
	ratio_calcH(rootHnc, Hnc);
	ratio_calcH(rootHcc, Hcc);
	
	ratio_calcC(rootCnc, Cnc, Cnclike);
	ratio_calcC(rootCcc, Ccc, Ccclike);
	
	double xsecHnc = crosssection(txtHnc);
	double xsecCnc = crosssection(txtCnc);

	double xsecHcc = crosssection(txtHcc);
	double xsecCcc = crosssection(txtCcc);
	
	for (int i = 0; i < 50; i++)
	{
		Hnc[i] *= 2 * xsecHnc / events / 0.04;
		Cnc[i] *= 12 * xsecCnc / events / 0.04;
		Cnclike[i] *= 12 * xsecCnc / events / 0.04;
		CH2nc[i] = (Hnc[i] + Cnc[i]) / 14;
		CH2nclike[i] = (Hnc[i] + Cnclike[i]) / 14;

		Hcc[i] *= 2 * xsecHcc / events / 0.04;
		Ccc[i] *= 12 * xsecCcc / events / 0.04;
		Ccclike[i] *= 12 * xsecCcc / events / 0.04;
		CH2cc[i] = (Hcc[i] + Ccc[i]) / 14;
		CH2cclike[i] = (Hcc[i] + Ccclike[i]) / 14;	
	}
	
	string rout = "mb_nce/cc/results/ratio.txt";
	string routlike = "mb_nce/cc/results/ratiolike.txt";
	
	ofstream file(rout.c_str());
	ofstream filelike(routlike.c_str());
	
	for (int i = 0; i < 50; i++)
	{
		file << i*0.04 << " " << CH2nc[i]/CH2cc[i] << (i+1)*0.04 << " " << CH2nc[i]/CH2cc[i] << endl;
		filelike << i*0.04 << " " << CH2nclike[i]/CH2cclike[i] << (i+1)*0.04 << " " << CH2nclike[i]/CH2cclike[i] << endl;
	}
	file.close();
	filelike.close();
}

void ratio_calcC(string in, double *res, double *reslike)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		bool nopi = (e1->nof(111) + e1->nof(211) + e1->nof(-211) == 0);
		bool nopilike = (e1->fof(111) + e1->fof(211) + e1->fof(-211) == 0);
		
		double Q2 = -e1->q2()/1000000.0;
		int w = Q2/0.04;
			
		if (w >= 50 or w < 0)
			continue;
		
		if (nopi) res[w]++;				
		if (nopilike) reslike[w]++;				
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

void ratio_calcH(string in, double *res)
{
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		bool nopi = (e1->nof(111) + e1->nof(211) + e1->nof(-211) == 0);
		
		double Q2 = -e1->q2()/1000000.0;
		int w = Q2/0.04;
			
		if (w >= 50 or w < 0)
			continue;
		
		if (nopi) res[w]++;				
		
		cout << in << ": " << 100*i/events << "%\r" << flush;
	}
	
	cout << in << ": done" << endl << endl;
		
	delete e1;
	delete tt1;
	delete tf1;
}

int main()
{	
	//sim();
	//calc();
	//chi2();
	//re_sim();
	
	q2_calc("/root/mbbeam_C_fg_QEL_1600.root.fsi.root");
	
	//sim_cc();
	//ratio_calc("/root/mbbeam_C_fg_QEL_1600.root.fsi.root");

	return 1;
}
