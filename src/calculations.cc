#include "fsi.h"
#include "data.h" 
#include "nucleusmaker.h"
#include "calculations.h"
#include "event1.h"
#include "dirs.h"
#include <TFile.h>
#include <TTree.h>

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
			for (int j = 0; j < 10; j++)
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

double crosssection (string filename, int dyn)
{
	string y;
	double x;
		
	vector<double> v(10);
	ifstream Input (filename.c_str());

	if (Input)
	{
		do
		{
			getline (Input, y);
			for (int j = 0; j < 10; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						Input>>y;
					}
					Input>>v[j];
				}
		}while (Input);
	}
	
	double result = v[dyn];
		
	return result;
}

void crosssection (string filename, double *xsec, int nc)
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
	
	xsec[4] = 0;
		
	for (int i = 0; i < 4; i++)
	{
		xsec[i] = v[2*i+nc];
		xsec[4] += v[2*i+nc];
	}
}

void calclog(string command)
{
	logfile<<day<<" "<<fullm<<" "<<year<<" "<<hour<<":"<<minute<<endl<<command<<endl<<endl;
}

void zero (double *tab, int bins)
{
	for (int i = 0; i < bins; i++) tab[i] = 0;
}

void zero (int *tab, int bins)
{
	for (int i = 0; i < bins; i++) tab[i] = 0;
}

void put (double value, double *source, double *target, double &rest, const int bins)
{
	if (value <= source[0]) target[0]++;
	else if (value > source[bins-1]) rest ++;
	else
	{
		for (int i = 0; i < bins-1; i++)
		{
			if (value > source[i] and value <= source[i+1])
			{
				target[i+1]++;
				break;
			}
		}
	}
}

void put (double value, double *source, double *target, double &rest, const int bins, double weight)
{
	if (value <= source[0]) target[0]+=weight;
	else if (value > source[bins-1]) rest+=weight;
	else
	{
		for (int i = 0; i < bins-1; i++)
		{
			if (value > source[i] and value <= source[i+1])
			{
				target[i+1]+=weight;
				break;
			}
		}
	}
}

void put (double value, const double *source, const double *sourceerr, double *target, const int bins)
{
	for (int i = 0; i < bins; i++)
	{
		if (value > (source[i] - sourceerr[i]) and value <= (source[i] + sourceerr[i]))
		{
			target[i]++;
			break;
		}
	}
}

void merge (double *tab1, double weight1, double *tab2, double weight2, double *result, const int bins)
{
	for (int i = 0; i < bins; i++)
	{
		result[i] = (tab1[i]*weight1 + tab2[i]*weight2)/(weight1 + weight2);
	}
}

double factor (double *target, const double *source, const int bins)
{
	double ss = 0;
	double st = 0;
	
	for (int i = 0; i < bins; i++)
	{
		ss += source[i];
		st += target[i];
	}
	
	return ss/st;
}

int calcK2K (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating K2K ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help = "root_files/K2K_NC_Oxygen_1m_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string Hfile = find_last("root_files/K2K_NC_Hydrogen_1m_*.root");
	string Ofile = find_last(help);
	
	if (noFile(Hfile) or noFile(Ofile) or noFile(Hfile+string(".txt")) or noFile(Ofile + string(".txt")) or noFile("root_files/K2K_CC_Hydrogen_0k.root.txt") or noFile("root_files/K2K_CC_Oxygen_0k.root.txt"))
	{
		logfile<<"K2K ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"K2K ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events     = 100000;
	const int bins       = 8;
	
	double H[bins]; zero(H, bins);
	double O[2][bins]; zero(O[0], bins); zero(O[1], bins);  //0 - before FSI, 1 - after FSi
	double vivi[2][4][bins];
	double mom[bins];
	
	for (int i = 0; i < bins; i++) mom[i] = K2Kx[i] + K2Kxerr[i];
	
	/*
	 * [0][0] - pi0 -> 0pi          [1][0] - 0pi -> pi0
	 * [0][1] - pi0 -> pi-          [1][1] - pi- -> pi0
	 * [0][2] - pi0 -> pi+          [1][2] - pi+ -> pi0
	 * [0][3] - pi0 -> more pis     [1][3] - more pis -> pi0
	 * 
	 */
	
	for (int i = 0; i < 4; i++) {zero(vivi[0][i], bins); zero(vivi[1][i], bins);}
	
	double restH = 0;
	double restO[2] = {0, 0};
	
	int counterH = 0;
	int counterO[2] = {0, 0};
	
	help = Hfile + string(".txt");
	double crossH = crosssection(help);
	help = Ofile + string(".txt");
	double crossO = crosssection(help);
	
	TFile *tf1 = new TFile(Hfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
				
		if (pion == 1)
		{
			counterH++;
			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg==111)
				{
					double val = e1->out[k].momentum();
					put(val, mom, H, restH, bins);
				}
			}
		}
		cout<<Hfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Hfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	TFile *tf2 = new TFile(Ofile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
		
		int pion    = 100*e2->nof(211) + 10*e2->nof(-211) + e2->nof(111);
		int pionfsi = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
				
		if (pion == 1)
		{
			counterO[0]++;
			
			for (int k = 0; k < e2->n(); k++)
			{
				if (e2->out[k].pdg==111)
				{
					double val = e2->out[k].momentum();
					put(val, mom, O[0], restO[0], bins);
					
					double rest;
					
					if (pionfsi == 0) put(val, mom, vivi[0][0], rest, bins);
					else if (pionfsi == 10) put(val, mom, vivi[0][1], rest, bins);
					else if (pionfsi == 100) put(val, mom, vivi[0][2], rest, bins);
					else if (pionfsi != 1) put(val, mom, vivi[0][3], rest, bins);
				}
			}
		}
		
		if (pionfsi == 1)
		{
			counterO[1]++;
			
			for (int k = 0; k < e2->f(); k++)
			{
				if (e2->post[k].pdg==111)
				{
					double val = e2->post[k].momentum();
					put(val, mom, O[1], restO[1], bins);
					
					double rest;
					
					if (pion == 0) put(val, mom, vivi[1][0], rest, bins);
					else if (pion == 10) put(val, mom, vivi[1][1], rest, bins);
					else if (pion == 100) put(val, mom, vivi[1][2], rest, bins);
					else if (pion != 1) put(val, mom, vivi[1][3], rest, bins);
				}
			}
		}
		
		cout<<Ofile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Ofile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;

	for (int i = 0; i < bins; i++)
	{
		O[0][i] *= crossO;
		O[1][i] *= crossO;
		H[i]    *= crossH;
		
		for (int k = 0; k < 4; k++)
		{
			vivi[0][k][i] *= crossO;
			vivi[1][k][i] *= crossO;
		}
	}
	
	restO[0] *= crossO;
	restO[1] *= crossO;
	restH    *= crossH;
	
	double H2O[2][bins];
	double H2Orest[2];
	
	H2Orest[0] = (2.0*restH + 16.0*restO[0])/(2.0 + 16.0);
	H2Orest[1] = (2.0*restH + 16.0*restO[1])/(2.0 + 16.0);	
		
	merge(H, 2.0, O[0], 16.0, H2O[0], bins);
	merge(H, 2.0, O[1], 16.0, H2O[1], bins);
	
	double factor0 = factor(H2O[0], K2Ky, bins);
	double factor1 = factor(H2O[1], K2Ky, bins);
	
	for (int i = 0; i < bins; i++)
	{
		H2O[0][i] *= factor0;
		H2O[1][i] *= factor1;
		
		for (int k = 0; k < 4; k++)
		{
			vivi[0][k][i] *= factor0;
			vivi[1][k][i] *= factor1;
		}
	}
	
	H2Orest[0] *= factor0;
	H2Orest[1] *= factor1;
		
	double crossHcc = crosssection("root_files/K2K_CC_Hydrogen_0k.root.txt");
	double crossOcc = crosssection("root_files/K2K_CC_Oxygen_0k.root.txt");
	
	double ratio[2];
	ratio[0] = (2.0*counterH*crossH + 16.0*counterO[0]*crossO)/(2.0*crossHcc + 16.0*crossOcc)/events;
	ratio[1] = (2.0*counterH*crossH + 16.0*counterO[1]*crossO)/(2.0*crossHcc + 16.0*crossOcc)/events;
	
	help = "results/K2K/";
	run(string("mkdir -p ") + help);
	help += string("k2k_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file(help.c_str());
	
	cout<<"NCpi0/CCall ratio before FSI: "<<ratio[0]<<endl<<endl<<"NCpi0/CCall ration after FSI: "<<ratio[1]<<endl<<endl;
	file<<"#NCpi0/CCall ratio before FSI: "<<ratio[0]<<endl<<endl<<"#NCpi0/CCall ration after FSI: "<<ratio[1]<<endl<<endl;

	cout<<"Momentum | K2K data | NuWro before FSI | NuWro after FSI ("<<fzname[fz]<<")"<<endl<<endl;
	file<<"#Momentum | NuWro before FSI | NuWro after FSI ("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		cout<<K2Kx[i]<<" "<<K2Ky[i]<<" "<<H2O[0][i]<<" "<<H2O[1][i]<<endl;
		file<<K2Kx[i]<<" "<<H2O[0][i]<<" "<<H2O[1][i]<<" "<<vivi[0][0][i]<<" "<<vivi[1][0][i]<<" "<<vivi[0][1][i] + vivi[0][2][i]<<" "<<vivi[1][1][i] + vivi[1][2][i]<<" "<<vivi[0][3][i]<<" "<<vivi[1][3][i];
		if (i < (bins - 2)) file<<endl;
		else if (i == (bins - 2)) file<<" "<<K2Kx[bins]-K2Kxerr[bins]<<" "<<H2Orest[0]<<" "<<H2Orest[1]<<endl;
		else if (i == (bins - 1)) file<<" "<<K2Kx[bins]+K2Kxerr[bins]<<" "<<H2Orest[0]<<" "<<H2Orest[1]<<endl;
	}
	
	cout<<K2Kx[bins]<<" "<<K2Ky[bins]<<" "<<H2Orest[0]<<" "<<H2Orest[1]<<endl;
	
	cout<<endl<<help<<" created."<<endl<<endl;
	help += string(" created from ") + Hfile + string(" and ")  + Ofile + string(" files");
	calclog(help);
	file.close();

	return 1;	
}

int calcMB (int fz, int xs, bool anti)
{	
	if (anti) cout<<endl<<endl<<"Calculating MiniBooNE antineutrino ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	else cout<<endl<<endl<<"Calculating MiniBooNE neutrino ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help;
	
	if (anti) help = "root_files/MB_anti_NC_Carbon_5m_";
	else help = "root_files/MB_NC_Carbon_5m_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string Hfile;
	
	if (anti) Hfile = find_last("root_files/MB_anti_NC_Hydrogen_5m_*.root");
	else Hfile = find_last("root_files/MB_NC_Hydrogen_5m_*.root");
	
	string Cfile = find_last(help);
	
	if (noFile(Hfile) or noFile(Cfile) or noFile(Hfile+string(".txt")) or noFile(Cfile + string(".txt")))
	{
		logfile<<"MiniBooNE ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"MiniBooNE ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events = 5000000;

	int ile;
	if (anti) ile = 10;
	else ile = 11;
	
	const int bins = ile;
	
	if (anti) ile = 10;
	else ile = 18;
	
	const int angbins = ile;
	
	double H[bins]; zero(H, bins);
	double C[2][bins]; zero(C[0], bins); zero(C[1], bins);  //0 - before FSI, 1 - after FSi
	double vivi[2][4][bins];
	double mom[bins];
	
	for (int i = 0; i < bins; i++)
	{
		if (anti) mom[i] = MBantiMomx[i] + MBantiMomxerr[i];
		else mom[i] = MBMomx[i] + MBMomxerr[i];
	}

	double Hang[angbins]; zero(Hang, angbins);
	double Cang[2][angbins]; zero(Cang[0], angbins); zero(Cang[1], angbins);  //0 - before FSI, 1 - after FSi
	double ang[angbins];
	
	for (int i = 0; i < angbins; i++)
	{
		if (anti) ang[i] = MBantiAnglex[i] + MBantiAnglexerr[i];
		else ang[i] = MBAnglex[i] + MBAnglexerr[i];
	}

	
	/*
	 * [0][0] - pi0 -> 0pi          [1][0] - 0pi -> pi0
	 * [0][1] - pi0 -> pi-          [1][1] - pi- -> pi0
	 * [0][2] - pi0 -> pi+          [1][2] - pi+ -> pi0
	 * [0][3] - pi0 -> more pis     [1][3] - more pis -> pi0
	 * 
	 */
	
	for (int i = 0; i < 4; i++) {zero(vivi[0][i], bins); zero(vivi[1][i], bins);}
		
	help = Hfile + string(".txt");
	double crossH = crosssection(help);
	help = Cfile + string(".txt");
	double crossC = crosssection(help);

	int counterH = 0;
	int counterC[2] = {0, 0};
		
	TFile *tf1 = new TFile(Hfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double rest;
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
				
		if (pion == 1)
		{
			counterH++;
			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg==111)
				{
					double val = e1->out[k].momentum();
					put(val, mom, H, rest, bins);
					val = e1->out[k].p().z/val;
					put(val, ang, Hang, rest, angbins);
					//int a = val/0.2 + 5;
					//if (a == 10) a--;
					//Hang[a]++;					
				}
			}
		}
		cout<<Hfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Hfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	TFile *tf2 = new TFile(Cfile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
		
		int pion    = 100*e2->nof(211) + 10*e2->nof(-211) + e2->nof(111);
		int pionfsi = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
				
		if (pion == 1)
		{
			counterC[0]++;
						
			for (int k = 0; k < e2->n(); k++)
			{
				if (e2->out[k].pdg==111)
				{
					double val = e2->out[k].momentum();
					put(val, mom, C[0], rest, bins);
					
					if (pionfsi == 0) put(val, mom, vivi[0][0], rest, bins);
					else if (pionfsi == 10) put(val, mom, vivi[0][1], rest, bins);
					else if (pionfsi == 100) put(val, mom, vivi[0][2], rest, bins);
					else if (pionfsi != 1) put(val, mom, vivi[0][3], rest, bins);
					
					val = e2->out[k].p().z/val;
					put(val, ang, Cang[0], rest, angbins);					
					//int a = val/0.2 + 5;
					//if (a == 10) a--;
					//Cang[0][a]++;	
				}
			}
		}
		
		if (pionfsi == 1)
		{
			counterC[1]++;
			
			for (int k = 0; k < e2->f(); k++)
			{
				if (e2->post[k].pdg==111)
				{
					double val = e2->post[k].momentum();
					put(val, mom, C[1], rest, bins);
					
					if (pion == 0) put(val, mom, vivi[1][0], rest, bins);
					else if (pion == 10) put(val, mom, vivi[1][1], rest, bins);
					else if (pion == 100) put(val, mom, vivi[1][2], rest, bins);
					else if (pion != 1) put(val, mom, vivi[1][3], rest, bins);
					
					val = e2->post[k].p().z/val;
					put(val, ang, Cang[1], rest, angbins);
					//int a = val/0.2 + 5;
					//if (a == 10) a--;
					//Cang[1][a]++;	
				}
			}
		}
		
		cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;	
	
	double cross[2];
	
	cross[0] = (2.0*counterH*crossH + 12.0*counterC[0]*crossC)/(2.0 + 12.0)/events;
	cross[1] = (2.0*counterH*crossH + 12.0*counterC[1]*crossC)/(2.0 + 12.0)/events;
	
	double CH2[2][bins];
	double CH2ang[2][angbins];
	
	for (int i = 0; i < bins; i++)
	{
		C[0][i] *= crossC;
		C[1][i] *= crossC;
		H[i]    *= crossH;
		
		for (int k = 0; k < 4; k++)
		{
			vivi[0][k][i] *= crossC;
			vivi[1][k][i] *= crossC;
		}
	}
	
	for (int i = 0; i < angbins; i++)
	{
		Cang[0][i] *= crossC;
		Cang[1][i] *= crossC;
		Hang[i]    *= crossH;
	}
		
	merge(H, 2.0, C[0], 12.0, CH2[0], bins);
	merge(H, 2.0, C[1], 12.0, CH2[1], bins);
	
	merge(Hang, 2.0, Cang[0], 12.0, CH2ang[0], angbins);
	merge(Hang, 2.0, Cang[1], 12.0, CH2ang[1], angbins);
	
	for (int i = 0; i < bins; i++)
	{
		double width;
		if (anti) width = 2*MBantiMomxerr[i]*events;
		else width = 2*MBMomxerr[i]*events;
				
		CH2[0][i] /= width;
		CH2[1][i] /= width;
				
		for (int k = 0; k < 4; k++)
		{
			vivi[0][k][i] /= width;
			vivi[1][k][i] /= width;
		}
	}	
	
	for (int i = 0; i < angbins; i++)
	{	
		double widthang;// = 0.2 * events;
		if (anti) widthang = 2.0*MBantiAnglexerr[i]*events;
		else widthang = 2.0*MBAnglexerr[i]*events;

		CH2ang[0][i] /= widthang;
		CH2ang[1][i] /= widthang;		
	}
	
	help = "results/MB/"; 
	run(string("mkdir -p ") + help);
	help += string("mb_");
	if (anti) help += string("anti_");
	help += fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file;
	file.open(help.c_str());
	
	file<<"#Cross section before FSI: "<<cross[0]<<endl<<"#Cross section after FSI: "<<cross[1]<<endl<<endl;
		
	file<<"#Momentum | NuWro before FSI | NuWro after FSI | vivi (0pi <-> 1otherpi <-> more pi) | Angle | MB data | NuWro before FSI | After FSI("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int i = 0; i < angbins; i++)
	{
		if (i < bins)
		{		
			if (anti) file<<MBantiMomx[i];
			else file<<MBMomx[i];
			file<<" "<<CH2[0][i]<<" "<<CH2[1][i]<<" "<<vivi[0][0][i]<<" "<<vivi[1][0][i]<<" "<<vivi[0][1][i] + vivi[0][2][i]<<" "<<vivi[1][1][i] + vivi[1][2][i]<<" "<<vivi[0][3][i]<<" "<<vivi[1][3][i]<<" ";
		}
		else file << "0 0 0 0 0 0 0 0 0 ";

		if (anti) file<<MBantiAnglex[i];
		else file<<MBAnglex[i];
		file<<" "<<CH2ang[0][i]<<" "<<CH2ang[1][i]<<endl;
	}
		
	cout<<endl<<help<<" created."<<endl<<endl;
	help += string(" created from ") + Hfile + string(" and ")  + Cfile + string(" files");
	calclog(help);
	file.close();
	
	double factor0;
	double factor1;
	
	double anglefac0;
	double anglefac1;
	
	if (anti)
	{
		factor0 = factor(CH2[0], MBantiMomy, bins);
		factor1 = factor(CH2[1], MBantiMomy, bins);
		
		anglefac0 = factor(CH2ang[0], MBantiAngley, angbins);
		anglefac1 = factor(CH2ang[1], MBantiAngley, angbins);
	}
	else
	{
		factor0 = factor(CH2[0], MBMomy, bins);
		factor1 = factor(CH2[1], MBMomy, bins);
		
		anglefac0 = factor(CH2ang[0], MBAngley, angbins);
		anglefac1 = factor(CH2ang[1], MBAngley, angbins);
	}
	
	for (int i = 0; i < bins; i++)
	{
		CH2[0][i] *= factor0;
		CH2[1][i] *= factor1;
		
		for (int k = 0; k < 4; k++)
		{
			vivi[0][k][i] *= factor0;
			vivi[1][k][i] *= factor1;
		}
	}
	
	for (int i = 0; i < angbins; i++)
	{
		CH2ang[0][i] *= factor0;
		CH2ang[1][i] *= factor1;
	}
	
	help = "results/MB/";
	run(string("mkdir -p ") + help);
	help += string("mb_normalized_");
	if (anti) help += string("anti_");
	help += fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	file.open(help.c_str());
	
	file<<"#Momentum | NuWro before FSI | NuWro after FSI | vivi (0pi <-> 1otherpi <-> more pi) | Angle | MB data | NuWro before FSI | After FSI("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int i = 0; i < angbins; i++)
	{
		if (i < bins)
		{		
			if (anti) file<<MBantiMomx[i];
			else file<<MBMomx[i];
			file<<" "<<CH2[0][i]<<" "<<CH2[1][i]<<" "<<vivi[0][0][i]<<" "<<vivi[1][0][i]<<" "<<vivi[0][1][i] + vivi[0][2][i]<<" "<<vivi[1][1][i] + vivi[1][2][i]<<" "<<vivi[0][3][i]<<" "<<vivi[1][3][i]<<" ";
		}
		else file << "0 0 0 0 0 0 0 0 0 ";

		if (anti) file<<MBantiAnglex[i];
		else file<<MBAnglex[i];
		file<<" "<<CH2ang[0][i]<<" "<<CH2ang[1][i]<<endl;
	}
		
	cout<<endl<<help<<" created."<<endl<<endl;
	help += string(" created from ") + Hfile + string(" and ")  + Cfile + string(" files");
	calclog(help);
			
	file.close();

	return 1;	
}

int calcMBback (int fz, int xs, bool anti)
{	
	if (anti) cout<<endl<<endl<<"Calculating MiniBooNE antineutrino ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	else cout<<endl<<endl<<"Calculating MiniBooNE neutrino ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help;
	
	if (anti) help = "root_files/MB_anti_Carbon_1m_";
	else help = "root_files/MB_Carbon_1m_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
		
	string Cfile = find_last(help);
	
	get_date();
	
	const int events = 1000000;

	const int backbins = 10;

	double Q2[backbins];
	double backpi[backbins]; zero(backpi, backbins);
	double backpr[backbins]; zero(backpr, backbins);
	double normpi[backbins]; zero(normpi, backbins);
	
	for (int i = 0; i < backbins; i++) Q2[i] = (i+1)*0.2;
		
	double rest;
	
	TFile *tf2 = new TFile(Cfile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
				
		double q2 = -e2->q2()/1000000.0;
		put(q2, Q2, normpi, rest, backbins);
		
		for (int k = 0; k < e2->f(); k++)
		{
			if (e2->post[k].pdg == 2212 && e2->post[k].p().z < 0) put(q2, Q2, backpr, rest, backbins);
			else if ((e2->post[k].pdg == 221 or e2->post[k].pdg == -211 or e2->post[k].pdg == 111) and e2->post[k].p().z < 0) put(q2, Q2, backpi, rest, backbins);
		}
				
		cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;
	
	for (int i = 0; i < backbins; i++)
	{
		if (normpi[i] != 0)
		{
			backpi[i] /= normpi[i];
			backpr[i] /= normpi[i];
		}
		else
		{
			backpi[i] = 0;
			backpr[i] = 0;
		}
		
		Q2[i] = (i+0.5)*0.2;
	}
	
	string bhelp = "results/MB/back_" + fzwork[fz] + string(".txt");
	if (anti) bhelp = "results/MB/back_anti_" + fzwork[fz] + string(".txt");
	
	ofstream wsteczne(bhelp.c_str());
		
	for (int i = 0; i < backbins; i++)
	{
		wsteczne << Q2[i] << " " << backpr[i] << " " << backpi[i] << endl;
	}
	
	wsteczne.close();	

	return 1;	
}

int calcMBCC (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating MiniBooNE CC ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Cfile = "root_files/MB_CC_Carbon_100k_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);
	
	string Hfile = find_last("root_files/MB_CC_Hydrogen_100k_*.root");
		
	if (noFile(Hfile) or noFile(Cfile) or noFile(Hfile+string(".txt")) or noFile(Cfile + string(".txt")))
	{
		logfile<<"MiniBooNE CC ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"MiniBooNE CC ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events = 100000;
	
	const int bins = 11;
	
	double H[bins]; zero(H, bins);
	double C[2][bins]; zero(C[0], bins); zero(C[1], bins);  //0 - before FSI, 1 - after FSi
	double mom[bins];
	
	for (int i = 0; i < bins; i++)
	{
		mom[i] = MBCCmom[i] + MBCCmomerr[i];
	}
	
	string help;	
	help = Hfile + string(".txt");
	double crossH = crosssection(help);
	help = Cfile + string(".txt");
	double crossC = crosssection(help);
		
	TFile *tf1 = new TFile(Hfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double rest;
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
				
		if (pion == 1)
		{
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg==111)
				{
					double val = e1->out[k].momentum()/1000.0;
					put(val, mom, H, rest, bins);
				}
			}
		}
		cout<<Hfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Hfile<<": done"<<endl<<endl;
	
	delete e1;
	delete tt1;
	delete tf1;
	
	TFile *tf2 = new TFile(Cfile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
		
		int pion    = 100*e2->nof(211) + 10*e2->nof(-211) + e2->nof(111);
		int pionfsi = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
				
		if (pion == 1)
		{
			for (int k = 0; k < e2->n(); k++)
			{
				if (e2->out[k].pdg==111)
				{
					double val = e2->out[k].momentum()/1000.0;
					put(val, mom, C[0], rest, bins);

				}
			}
		}
		
		if (pionfsi == 1)
		{
			for (int k = 0; k < e2->f(); k++)
			{
				if (e2->post[k].pdg==111)
				{
					double val = e2->post[k].momentum()/1000.0;
					put(val, mom, C[1], rest, bins);
				}
			}
		}
		
		cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;
	
	double CH2[2][bins];
	
	for (int i = 0; i < bins; i++)
	{
		C[0][i] *= crossC;
		C[1][i] *= crossC;
		H[i]    *= crossH;
	}
		
	merge(H, 2.0, C[0], 12.0, CH2[0], bins);
	merge(H, 2.0, C[1], 12.0, CH2[1], bins);
	
	for (int i = 0; i < bins; i++)
	{
		double width = 2.0*MBCCmomerr[i]*events/12.0;
		
		CH2[0][i] /= width;
		CH2[1][i] /= width;
	}	
	
	help = "results/MBCC/"; 
	run(string("mkdir -p ") + help);
	help += string("mb_cc_");
	help += fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file;
	file.open(help.c_str());
	
	file<<"#Momentum | NuWro before FSI | NuWro after FSI ("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{		
		file<<MBCCmom[i]<<" "<<CH2[0][i]<<" "<<CH2[1][i]<<endl;
	}
		
	cout<<endl<<help<<" created."<<endl<<endl;
	help += string(" created from ") + Hfile + string(" and ")  + Cfile + string(" files");
	calclog(help);
	file.close();
	
	double factor0 = factor(CH2[0], MBCCxsec, bins);
	double factor1 = factor(CH2[1], MBCCxsec, bins);

	for (int i = 0; i < bins; i++)
	{
		CH2[0][i] *= factor0;
		CH2[1][i] *= factor1;
	}
	
	help = "results/MBCC/"; 
	run(string("mkdir -p ") + help);
	help += string("mb_cc_normalized_");
	help += fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file2;
	file.open(help.c_str());
	
	file2<<"#Momentum | NuWro before FSI | NuWro after FSI ("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{		
		file<<MBCCmom[i]<<" "<<CH2[0][i]<<" "<<CH2[1][i]<<endl;
	}
		
	cout<<endl<<help<<" created."<<endl<<endl;
	help += string(" created from ") + Hfile + string(" and ")  + Cfile + string(" files");
	calclog(help);
	file2.close();
				
	return 1;	
}

int calcMBCCtotal (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating MiniBooNE CC total ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	get_date();
	
	const int events = 10000;

	const int bins = 14;
	
	string help = "results/MBCC/"; 
	run(string("mkdir -p ") + help);
	help += string("mb_cc_total_");
	help += fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file;
	file.open(help.c_str());

	file<<"#Neutrino energy | total xsec before FSI | total xsec after FSI ("<<fzname[fz]<<")"<<endl<<endl;
	
	for (int z = 0; z < bins; z++)
	{
		double en = MBCCnuen[z]*1000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string Cfile = "root_files/E" + energy + string("_MB_CC_Carbon_100k_");		
		Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");	
		Cfile = find_last(Cfile);		
		
		string Hfile = "root_files/E" + energy + string("_MB_CC_Hydrogen_100k_*.root");
		Hfile = find_last(Hfile);			
		
		double xsec[4]; //0 - hydrogen, 1 - carbon, 2 - carbon fsi
		zero(xsec, 3);
		
		help = Hfile + string(".txt");
		double crossH = crosssection(help);
		help = Cfile + string(".txt");
		double crossC = crosssection(help);
			
		TFile *tf1 = new TFile(Hfile.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
			
		tt1->SetBranchAddress("e",&e1);
		
		for (int i = 0; i < events; i++)
		{
			tt1->GetEntry(i);
			
			int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
								
			if (pion == 1) xsec[0]++;

			cout<<Hfile<<": "<<100*i/events<<"%\r"<<flush;
		}
		
		cout<<Hfile<<": done"<<endl<<endl;
			
		delete e1;
		delete tt1;
		delete tf1;
		
		TFile *tf2 = new TFile(Cfile.c_str());
		TTree *tt2 = (TTree*)tf2->Get("treeout");
		event *e2   = new event();
			
		tt2->SetBranchAddress("e",&e2);
		
		for (int i = 0; i < events; i++)
		{
			tt2->GetEntry(i);
			
			int pion    = 100*e2->nof(211) + 10*e2->nof(-211) + e2->nof(111);
			int pionfsi = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
				
			if (pion == 1) xsec[1]++;
			if (pionfsi == 1) xsec[2]++;
						
			cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
		}
		
		cout<<Cfile<<": done"<<endl<<endl;
			
		delete e2;
		delete tt2;
		delete tf2;		
				
		xsec[0] *= crossH/events;
		xsec[1] *= crossC/events;
		xsec[2] *= crossC/events;
		
		double result[2];
		
		result[0] = 2.0*xsec[0] + 12.0*xsec[1];
		result[1] = 2.0*xsec[0] + 12.0*xsec[2];
		
		file<<MBCCnuen[z]+MBCCnuenerr[z]<<" "<<result[0]<<" "<<result[1]<<endl;
		if (z==0) file<<MBCCnuen[z]-MBCCnuenerr[z]<<" "<<result[0]<<" "<<result[1]<<endl;
	}

	file.close();

	return 1;	
}

int calcSB (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating SciBooNE ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help = "root_files/MB_NC_Carbon_5m_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string Hfile = find_last("root_files/MB_NC_Hydrogen_5m_*.root");
	string Cfile = find_last(help);
	
	if (noFile(Hfile) or noFile(Cfile) or noFile(Hfile+string(".txt")) or noFile(Cfile + string(".txt")) or noFile("root_files/MB_CC_Hydrogen_0k.root.txt") or noFile("root_files/MB_CC_Carbon_0k.root.txt"))
	{
		logfile<<"SciBooNE ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"SciBooNE ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events     = 5000000;
	const int bins       = 9;
	const int binsang    = 10;
	
	double H[bins] = {0};
	double C[bins] = {0};
	
	double Hang[binsang] = {0};
	double Cang[binsang] = {0};
	
	int Hsinglepi0 = 0;
	int Hdoublepi0 = 0;
	int Hpi0pic = 0;

	int Csinglepi0 = 0;
	int Cdoublepi0 = 0;
	int Cpi0pic = 0;
	
	double restH = 0;
	double restC = 0;
	
	int counterH = 0;
	int counterC = 0;
	
	help = Hfile + string(".txt");
	double crossH = crosssection(help);
	help = Cfile + string(".txt");
	double crossC = crosssection(help);
		
	TFile *tf1 = new TFile(Hfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
				
		if (pion == 1)
		{
			counterH++;
			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg==111)
				{
					int a = e1->out[k].momentum()/80.0;
					
					if (a < 9) H[a]++;
					else restH++;
					
					int b = e1->out[k].p().z/e1->out[k].momentum()/0.2 + 5;
					if (b == 10) b--;
					
					Hang[b]++;					
				}
			}
			
			if (pion == 1) Hsinglepi0++;
			else if (pion == 2) Hdoublepi0++;
			else Hpi0pic++;
		}
		cout<<Hfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Hfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	TFile *tf2 = new TFile(Cfile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
		
		int pion = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
						
		if (pion == 1)
		{
			counterC++;
			
			for (int k = 0; k < e2->f(); k++)
			{
				if (e2->post[k].pdg==111)
				{
					int a = e2->post[k].momentum()/80.0;
					
					if (a < 9) C[a]++;
					else restC++;
					
					int b = e2->post[k].p().z/e2->post[k].momentum()/0.2 + 5;
					if (b == 10) b--;
					
					Cang[b]++;	
				}
			}
			
			if (pion == 1) Csinglepi0++;
			else if (pion == 2) Cdoublepi0++;
			else Cpi0pic++;
		}
		
		cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;

	for (int i = 0; i < bins; i++)
	{
		C[i] *= crossC;
		H[i] *= crossH;		
	}

	for (int i = 0; i < binsang; i++)
	{
		Cang[i] *= crossC;
		Hang[i] *= crossH;		
	}
	
	restC *= crossC;
	restH    *= crossH;
	
	double C8H8[bins];
	double C8H8ang[binsang];
	double C8H8rest;
	
	C8H8rest = (8.0*restH + 96.0*restC)/(8.0 + 96.0);	
		
	merge(H, 8.0, C, 96.0, C8H8, bins);
	merge(Hang, 8.0, Cang, 96.0, C8H8ang, binsang);
	
	double factormom = factor(C8H8, SBy, bins);
	double factorang = factor(C8H8ang, SByang, binsang);
	
	for (int i = 0; i < bins; i++)
	{
		//C8H8[i] *= factormom;
	}

	for (int i = 0; i < binsang; i++)
	{
		//C8H8ang[i] *= factorang;
	}
	
	//C8H8rest *= factormom;
		
	double crossHcc = crosssection("root_files/MB_CC_Hydrogen_0k.root.txt");
	double crossCcc = crosssection("root_files/MB_CC_Carbon_0k.root.txt");
	
	double ratio;
	ratio = (8.0*counterH*crossH + 96.0*counterC*crossC)/(8.0*crossHcc + 96.0*crossCcc)/events;
	
	double spi0 = (8.0*Hsinglepi0 + 96*Csinglepi0)/(8.0*counterH + 96.0*counterC);
	double dpi0 = (8.0*Hdoublepi0 + 96*Cdoublepi0)/(8.0*counterH + 96.0*counterC);
	double mixpi = (8.0*Hpi0pic + 96*Cpi0pic)/(8.0*counterH + 96.0*counterC);
	
	double sigma = (8.0*counterH*crossH + 96.0*counterC*crossC)/events;
	
	help = "results/SB/";
	run(string("mkdir -p ") + help);
	help += string("sb_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file(help.c_str());

	help = "results/SB/";
	run(string("mkdir -p ") + help);
	help += string("sb_ang_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream fileang(help.c_str());
	
	file<<"#NCpi0/CCall ratio: "<<ratio<<endl<<endl;
	file<<"#with singlepi0: "<< spi0 << " double pi0: " << dpi0 << " pi0+pic: " << mixpi << endl << endl;
	file<<"#sigma total: " << sigma << endl << endl;
	file<<"factorang: " << factorang << endl << endl;
	
	for (int i = 0; i < bins; i++)
	{
		file<<SBx[i]<<" "<<C8H8[i];
		if (i < (bins - 2)) file<<endl;
		else if (i == (bins - 2)) file<<" "<<SBx[bins]-SBxerr[bins]<<" "<<C8H8rest<<endl;
		else if (i == (bins - 1)) file<<" "<<SBx[bins]+SBxerr[bins]<<" "<<C8H8rest<<endl;
	}
	
	
	file.close();
	
	for (int i = 0; i < binsang; i++)
		fileang << SBxang[i] << " " << C8H8ang[i] << endl;
		
	fileang.close();

	return 1;	
}

int calcPNS (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Pion-Nucleus Scattering ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Cfile = "root_files/PionNucleus_Carbon_100k_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);

	string Ifile = "root_files/PionNucleus_Iron_100k_";
	Ifile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Ifile = find_last(Ifile);
	
	if (noFile(Cfile) or noFile(Ifile))
	{
		logfile<<"PionNucleus ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"PionNucleus ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events  = 100000;
	const double Emin = 0;
	const double Emax = 500;
	const int bins    = 10;
	
	double E[bins];
	double normC[bins]; zero(normC, bins);
	double normI[bins]; zero(normI, bins);
	
	for (int i = 0; i < bins; i++) E[i] = Emin + (i+1.0)*(Emax - Emin)/bins;

	double C[4][bins]; //0 - reaction, 1 - absorption, 2 - cex, 3 - inelastic
	double I[4][bins];
	
	for (int i = 0; i < 4; i ++)
	{
		zero(C[i], bins);
		zero(I[i], bins);
	}
	
	TFile *tf1 = new TFile(Cfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
	
		int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		double Ek = e1->out[0].Ek();
		double rest;
		
		put(Ek, E, normC, rest, bins);
		
		if (e1->number_of_interactions() != 0)
		{
			put(Ek, E, C[0], rest, bins);
			
			if (pion == 0) put(Ek, E, C[1], rest, bins);
			else if (pion == 1) put(Ek, E, C[2], rest, bins);
			else put (Ek, E, C[3], rest, bins);
		}
		
		cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
	}
		
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;

	TFile *tf2 = new TFile(Ifile.c_str());
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event *e2   = new event();
		
	tt2->SetBranchAddress("e",&e2);
	
	for (int i = 0; i < events; i++)
	{
		tt2->GetEntry(i);
	
		int pion = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
		double Ek = e2->out[0].Ek();
		double rest;

		put(Ek, E, normI, rest, bins);
		
		if (e2->number_of_interactions() != 0)
		{
			put(Ek, E, I[0], rest, bins);
			
			if (pion == 0) put(Ek, E, I[1], rest, bins);
			else if (pion == 1) put(Ek, E, I[2], rest, bins);
			else put (Ek, E, I[3], rest, bins);
		}
		
		cout<<Ifile<<": "<<100*i/events<<"%\r"<<flush;
	}
		
	cout<<Ifile<<": done"<<endl<<endl;
		
	delete e2;
	delete tt2;
	delete tf2;
	
	params pC;
	params pI;
	
	pC.nucleus_p = 6;
	pC.nucleus_n = 6;
	pI.nucleus_p = 26;
	pI.nucleus_n = 30;
	pC.nucleus_model = 1;
	pI.nucleus_model = 1;
	
	nucleus *cnucleus = make_nucleus (pC);
    nucleus & jadroC = *cnucleus;

	nucleus *inucleus = make_nucleus (pI);
    nucleus & jadroI = *inucleus;

	double radiusC = jadroC.radius();
	double radiusI = jadroI.radius();
 
    double factorC = 10.0 * M_PI * radiusC * radiusC / fermi / fermi;
	double factorI = 10.0 * M_PI * radiusI * radiusI / fermi / fermi;
	
	for (int i = 0; i < bins; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if (normC[i] != 0) C[k][i] *= factorC/normC[i];
			else C[k][i] = 0;
			if (normI[i] != 0) I[k][i] *= factorI/normI[i];
			else I[k][i] = 0;
		}
		
		E[i] = Emin + (i+0.5)*(Emax - Emin)/bins;
	}

	string help = "results/PNS/";
	run(string("mkdir -p ") + help);
	
	string chelp = help + string("pns_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string ihelp = help + string("pns_iron_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileC(chelp.c_str());
	ofstream fileI(ihelp.c_str());
	
	cout<<"Carbon: "<<endl<<endl<<"Kinetic energy | Reaction | Absorption | CEX | Inelastic"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<E[i]<<" "<<C[0][i]<<" "<<C[1][i]<<" "<<C[2][i]<<" "<<C[3][i]<<endl;
	cout<<endl<<"Iron: "<<endl<<endl<<"Kinetic energy | Reaction | Absorption | CEX | Inelastic"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<E[i]<<" "<<I[0][i]<<" "<<I[1][i]<<" "<<I[2][i]<<" "<<I[3][i]<<endl;
	
	fileC<<"#Kinetic energy | Reaction | Absorption | CEX | Inelastic"<<endl<<endl;
	fileI<<"#Kinetic energy | Reaction | Absorption | CEX | Inelastic"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		fileC<<E[i]<<" "<<C[0][i]<<" "<<C[1][i]<<" "<<C[2][i]<<" "<<C[3][i]<<endl;
		fileI<<E[i]<<" "<<I[0][i]<<" "<<I[1][i]<<" "<<I[2][i]<<" "<<I[3][i]<<endl;
	}
	
	cout<<endl<<chelp<<" created."<<endl<<endl;
	cout<<endl<<ihelp<<" created."<<endl<<endl;
	
	chelp += string(" created from ") + Cfile;
	ihelp += string(" created from ") + Ifile;
	
	calclog(chelp);
	calclog(ihelp);
	
	fileC.close();
	fileI.close();
	
	return 1;
}

void calcTrans (string filename, double *x, double *y, double *norm, const int bins, int events, int argument)
{
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();
		
	tt1->SetBranchAddress("e",&e);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		double arg = e->out[0].Ek();
		double rest;
		
		if (argument == 1) arg = e->out[0].Ek()*2.0*938.0;//*pow(10.0, -6.0);
		else if (argument == 2) arg = e->out[0].momentum();
		
		put(arg, x, norm, rest, bins);
		
		if (e->number_of_interactions() == 0) put(arg, x, y, rest, bins);
		/*else if (e->number_of_interactions() == 1 and (e->nod[0] == 1 or e->nod[4] == 1))
		{
			for (int k = 0; k < e->f(); k++)
			{
				if (e->post[k].pdg = e->out[0].pdg)
				{
					double cos = e->post[k].p()*e->out[0].p()/e->post[k].momentum()/e->out[0].momentum();
					if (cos < 0.95 and cos > -0.95) put(arg, x, y, rest, bins);
				}
			}
		}*/
	
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
		
	cout<<filename<<": done"<<endl<<endl;
		
	delete e;
	delete tt1;
	delete tf1;
	
}

int calcPrThe (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Proton Transparency for high energy ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Cfile = "root_files/ProtonTransparency_he_Carbon_100k_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);

	string Ifile = "root_files/ProtonTransparency_he_Iron_100k_";
	Ifile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Ifile = find_last(Ifile);
	
	if (noFile(Cfile) or noFile(Ifile))
	{
		logfile<<"ProtonTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"ProtonTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events  = 100000;
	const double Q2min = 500000;
	const double Q2max = 7000000;
	const int bins    = 20;
	
	double Q2[bins];
	double normC[bins]; zero(normC, bins);
	double normI[bins]; zero(normI, bins);
	
	for (int i = 0; i < bins; i++) Q2[i] = Q2min + (i+1.0)*(Q2max - Q2min)/bins;

	double C[bins]; zero(C, bins);
	double I[bins]; zero(I, bins);

	calcTrans(Cfile, Q2, C, normC, bins, events, 1);
	calcTrans(Ifile, Q2, I, normI, bins, events, 1);
	
	for (int i = 0; i < bins; i++)
	{
		if (normC[i] != 0) C[i] /= normC[i];
		else C[i] = 0;
		if (normI[i] != 0) I[i] /= normI[i];
		else I[i] = 0;

		Q2[i] = Q2min + (i+0.5)*(Q2max - Q2min)/bins;
		Q2[i] /= 1000000.0;
	}
	
	string help = "results/PrTrans/";
	run(string("mkdir -p ") + help);
	
	string chelp = help + string("prtrans_he_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string ihelp = help + string("prtrans_he_iron_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileC(chelp.c_str());
	ofstream fileI(ihelp.c_str());
	
	cout<<"Carbon: "<<endl<<endl<<"Q2 [GeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Q2[i]<<" "<<C[i]<<endl;
	cout<<endl<<"Iron: "<<endl<<endl<<"Q2 [GeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Q2[i]<<" "<<I[i]<<endl;
	
	fileC<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	fileI<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		fileC<<Q2[i]<<" "<<C[i]<<endl;
		fileI<<Q2[i]<<" "<<I[i]<<endl;
	}
	
	cout<<endl<<chelp<<" created."<<endl<<endl;
	cout<<endl<<ihelp<<" created."<<endl<<endl;
	
	chelp += string(" created from ") + Cfile;
	ihelp += string(" created from ") + Ifile;
	
	calclog(chelp);
	calclog(ihelp);
	
	fileC.close();
	fileI.close();
	
	return 1;
}

int calcPiThe (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Pion Transparency for high energy ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Cfile = "root_files/PionTransparency_he_Carbon_100k_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);

	string Afile = "root_files/PionTransparency_he_Aluminium_100k_";
	Afile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Afile = find_last(Afile);
	
	string Mfile = "root_files/PionTransparency_he_Copper_100k_";
	Mfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Mfile = find_last(Mfile);
		
	if (noFile(Cfile) or noFile(Afile) or noFile(Mfile))
	{
		logfile<<"PionTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"PionTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events  = 100000;
	const int bins    = 6;
	
	double ppi[bins];
	double normC[bins]; zero(normC, bins);
	double normA[bins]; zero(normA, bins);
	double normM[bins]; zero(normM, bins);
	
	ppi[0] = 1000.0*(PTppi[0] - (PTppi[1] - PTppi[0])/2.0);
	for (int i = 0; i < bins - 1; i++) 
        ppi[i+1] = 1000.0*(PTppi[i] + (PTppi[i+1] - PTppi[i])/2.0);
	ppi[bins-1] = 1000.0*(PTppi[bins-2] + (PTppi[bins-2] - PTppi[bins-3])/2.0);
	
	double C[bins]; zero(C, bins);
	double A[bins]; zero(A, bins);
	double M[bins]; zero(M, bins);

	calcTrans(Cfile, ppi, C, normC, bins, events, 2);
	calcTrans(Afile, ppi, A, normA, bins, events, 2);
	calcTrans(Mfile, ppi, M, normM, bins, events, 2);
	
	for (int i = 0; i < bins; i++)
	{
		if (normC[i] != 0) C[i] /= normC[i];
		else C[i] = 0;
		if (normA[i] != 0) A[i] /= normA[i];
		else A[i] = 0;
		if (normM[i] != 0) M[i] /= normM[i];
		else M[i] = 0;
	}
	
	string help = "results/PiTrans/";
	run(string("mkdir -p ") + help);
	
	string chelp = help + string("pitrans_he_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string ahelp = help + string("pitrans_he_aluminium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string mhelp = help + string("pitrans_he_copper_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileC(chelp.c_str());
	ofstream fileA(ahelp.c_str());
	ofstream fileM(mhelp.c_str());
		
	cout<<"Carbon: "<<endl<<endl<<"Q2 [GeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins-1; i++) cout<<PTQ2[i]<<" "<<C[i+1]<<endl;
	cout<<endl<<"Aluminium: "<<endl<<endl<<"Q2 [GeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins-1; i++) cout<<PTQ2[i]<<" "<<A[i+1]<<endl;
	cout<<endl<<"Copper: "<<endl<<endl<<"Q2 [GeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins-1; i++) cout<<PTQ2[i]<<" "<<M[i+1]<<endl;
		
	fileC<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	fileA<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	fileM<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	
	for (int i = 0; i < bins-1; i++)
	{
		fileC<<PTQ2[i]<<" "<<C[i+1]<<endl;
		fileA<<PTQ2[i]<<" "<<A[i+1]<<endl;
		fileM<<PTQ2[i]<<" "<<M[i+1]<<endl;
	}
	
	cout<<endl<<chelp<<" created."<<endl<<endl;
	cout<<endl<<ahelp<<" created."<<endl<<endl;
	cout<<endl<<mhelp<<" created."<<endl<<endl;
	
	chelp += string(" created from ") + Cfile;
	ahelp += string(" created from ") + Afile;
	mhelp += string(" created from ") + Mfile;
	
	calclog(chelp);
	calclog(ahelp);
	calclog(mhelp);
		
	fileC.close();
	fileA.close();
	fileM.close();
	
	return 1;
}

void calcTrans2 (string filename, double *x, double *y, double *norm, const int bins, int events)
{
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();
		
	tt1->SetBranchAddress("e",&e);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		double rest;
		
		if (e->nof(211) + e->nof(-211) + e->nof(111) == 1 and e->nof(2212) + e->nof(2112) == 1)
		{
			double mom = 0;
			
			for (int k = 0; k < e->n(); k++)
			{
				if (e->out[k].pdg == 211 or e->out[k].pdg == -211 or e->out[k].pdg == 111)
				{
					mom = e->out[k].momentum()/1000.0;
					break;
				}
			}
			
			if (mom > 2.5)
			{
				put(mom, x, norm, rest, bins);
		
				if (e->number_of_pion_elastic() + e->number_of_pion_ce() + e->number_of_pion_spp() + e->number_of_pion_dpp() + e->number_of_pion_abs() == 0) put(mom, x, y, rest, bins);
			}
		}
		
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
		
	cout<<filename<<": done"<<endl<<endl;
		
	delete e;
	delete tt1;
	delete tf1;	
}

int calcPiTrans(int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating Pion Transparency for high energy v2("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Cfile = "root_files/PiTrans_he_Carbon_1m_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);

	string Afile = "root_files/PiTrans_he_Aluminium_1m_";
	Afile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Afile = find_last(Afile);
	
	string Mfile = "root_files/PiTrans_he_Copper_1m_";
	Mfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Mfile = find_last(Mfile);
		
	//if (noFile(Cfile) or noFile(Afile) or noFile(Mfile))
	//{
	//	logfile<<"PionTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
	//	cout<<"PionTransparency ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
	//	return 1;
	//}
		
	get_date();
	
	const int events  = 1000000;
	const int bins    = 10; //6;
	
	double ppi[bins];
	double normC[bins]; zero(normC, bins);
	double normA[bins]; zero(normA, bins);
	double normM[bins]; zero(normM, bins);
		
	double C[bins]; zero(C, bins);
	double A[bins]; zero(A, bins);
	double M[bins]; zero(M, bins);

//	ppi[0] = 1000.0*(PTppi[0] - (PTppi[1] - PTppi[0])/2.0);
//	for (int i = 0; i < bins - 1; i++) ppi[i+1] = 1000.0*(PTppi[i] + (PTppi[i+1] - PTppi[i])/2.0);
//	ppi[bins-1] = 1000.0*(PTppi[bins-2] + (PTppi[bins-2] - PTppi[bins-3])/2.0);
	
	for (int i = 0; i < bins; i++) ppi[i] = (i+1)*0.2+2.5;
	
	calcTrans2(Cfile, ppi, C, normC, bins, events);
	calcTrans2(Afile, ppi, A, normA, bins, events);
	//calcTrans2(Mfile, ppi, M, normM, bins, events);
	
	for (int i = 0; i < bins; i++)
	{
		if (normC[i] != 0) C[i] /= normC[i];
		else C[i] = 0;
		if (normA[i] != 0) A[i] /= normA[i];
		else A[i] = 0;
		if (normM[i] != 0) M[i] /= normM[i];
		else M[i] = 0;
	}
	
	string help = "results/PiTrans/";
	run(string("mkdir -p ") + help);
	
	string chelp = help + string("pitrans_he_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string ahelp = help + string("pitrans_he_aluminium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string mhelp = help + string("pitrans_he_copper_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileC(chelp.c_str());
	ofstream fileA(ahelp.c_str());
	ofstream fileM(mhelp.c_str());
			
	fileC<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	fileA<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	fileM<<"#Q2 [GeV] | Transparency"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		ppi[i] = (i+0.5)*0.2 + 2.5;		
		fileC<<ppi[i]<<" "<<C[i]<<endl;
		fileA<<ppi[i]<<" "<<A[i]<<endl;
		fileM<<ppi[i]<<" "<<M[i]<<endl;
	}
	
	cout<<endl<<chelp<<" created."<<endl<<endl;
	cout<<endl<<ahelp<<" created."<<endl<<endl;
	cout<<endl<<mhelp<<" created."<<endl<<endl;
	
	chelp += string(" created from ") + Cfile;
	ahelp += string(" created from ") + Afile;
	mhelp += string(" created from ") + Mfile;
	
	calclog(chelp);
	calclog(ahelp);
	calclog(mhelp);
		
	fileC.close();
	fileA.close();
	fileM.close();
	
	return 1;
}

double calcPiTrans3b(string filename)
{
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e   = new event();
		
	tt1->SetBranchAddress("e",&e);
	
	double counter = 0;
	
	for (int i = 0; i < 100000; i++)
	{
		tt1->GetEntry(i);
		
		if (e->number_of_interactions() == 0) counter++;
	
		cout<<filename<<": "<<100*i/100000<<"%\r"<<flush;
	}
		
	cout<<filename<<": done"<<endl<<endl;
		
	delete e;
	delete tt1;
	delete tf1;

	return counter/100000.0;
}

int calcPiTrans3(int fz)
{
	double C[5];
	double A[5];
	double M[5];
	
	for (int i = 0; i < 5; i++)
	{
		int energy;
		
		energy = PTppi[i]*1000.0;
		energy = energy*energy + 140.0*140.0;
		energy = sqrt(energy);
		
		stringstream temp;
		string en;
		
		temp << energy;
		temp >> en;
		
		string rf = "root_files/PiTrans_e" + en + string("_carbon_100k_") + fzwork[fz] + rot;
		C[i] = calcPiTrans3b(rf);
		
		rf = "root_files/PiTrans_e" + en + string("_aluminium_100k_") + fzwork[fz] + rot;
		//A[i] = calcPiTrans3b(rf);

		rf = "root_files/PiTrans_e" + en + string("_copper_100k_") + fzwork[fz] + rot;
		M[i] = calcPiTrans3b(rf);
	}
	
	string name = "results/car_" + fzwork[fz] + ".txt";
	
	ofstream fileC(name.c_str());

	name = "results/alum_" + fzwork[fz] + ".txt";
	ofstream fileA(name.c_str());

	name = "results/cop_" + fzwork[fz] + ".txt";
	ofstream fileM(name.c_str());
		
	for (int i = 0; i < 5; i++)
	{
		fileC << PTppi[i] << " " << C[i] << endl;
		fileA << PTppi[i] << " " << A[i] << endl;
		fileM << PTppi[i] << " " << M[i] << endl;
	}
	
	fileC.close();
	fileA.close();
	fileM.close();
	
	return 1;
}

int calcPrTle (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Proton Transparency for low energy ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Lfile = "root_files/ProtonTransparency_le_Lithium_100k_";
	Lfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Lfile = find_last(Lfile);

	string Cfile = "root_files/ProtonTransparency_le_Carbon_100k_";
	Cfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);

	string Afile = "root_files/ProtonTransparency_le_Aluminium_100k_";
	Afile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Afile = find_last(Afile);
	
	if (noFile(Lfile) or noFile(Cfile) or noFile(Afile))
	{
		logfile<<"ProtonTransparency low energy ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"ProtonTransparency low energy ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events  = 100000;
	const double Ekmin = 100;
	const double Ekmax = 300;
	const int bins    = 10;
	
	double Ek[bins];
	double normL[bins]; zero(normL, bins);
	double normC[bins]; zero(normC, bins);
	double normA[bins]; zero(normA, bins);
	
	for (int i = 0; i < bins; i++) Ek[i] = Ekmin + (i+1.0)*(Ekmax - Ekmin)/bins;

	double L[bins]; zero(L, bins);
	double C[bins]; zero(C, bins);
	double A[bins]; zero(A, bins);

	calcTrans(Lfile, Ek, L, normL, bins, events, 0);
	calcTrans(Cfile, Ek, C, normC, bins, events, 0);
	calcTrans(Afile, Ek, A, normA, bins, events, 0);
	
	for (int i = 0; i < bins; i++)
	{
		if (normL[i] != 0) L[i] /= normL[i];
		else L[i] = 0;
		if (normC[i] != 0) C[i] /= normC[i];
		else C[i] = 0;
		if (normA[i] != 0) A[i] /= normA[i];
		else A[i] = 0;

		Ek[i] = Ekmin + (i+0.5)*(Ekmax - Ekmin)/bins;
	}
	
	string help = "results/PrTrans/";
	run(string("mkdir -p ") + help);
	
	string lhelp = help + string("prtrans_le_lithium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string chelp = help + string("prtrans_le_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string ahelp = help + string("prtrans_le_aluminium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileL(lhelp.c_str());
	ofstream fileC(chelp.c_str());
	ofstream fileA(ahelp.c_str());
	
	cout<<"Lithium: "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<L[i]<<endl;
	cout<<endl<<"Carbon: "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<C[i]<<endl;
	cout<<endl<<"Aluminium: "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<A[i]<<endl;
		
	fileL<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileC<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileA<<"#Ek [GeV] | Transparency"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		fileL<<Ek[i]<<" "<<L[i]<<endl;
		fileC<<Ek[i]<<" "<<C[i]<<endl;
		fileA<<Ek[i]<<" "<<A[i]<<endl;
	}
	
	cout<<endl<<lhelp<<" created."<<endl<<endl;
	cout<<endl<<chelp<<" created."<<endl<<endl;
	cout<<endl<<ahelp<<" created."<<endl<<endl;
	
	lhelp += string(" created from ") + Lfile;
	chelp += string(" created from ") + Cfile;
	ahelp += string(" created from ") + Afile;
	
	calclog(lhelp);
	calclog(chelp);
	calclog(ahelp);
	
	fileL.close();
	fileC.close();
	fileA.close();
	
	return 1;
}

int calcPiTle (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Pion Transparency for low energy ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string Lpfile = "root_files/PionTransparency_le_pip_Lithium_100k_";
	Lpfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Lpfile = find_last(Lpfile);

	string L0file = "root_files/PionTransparency_le_pi0_Lithium_100k_";
	L0file += fzwork[fz] + sep + xsec[xs] + string("*.root");
	L0file = find_last(L0file);

	string Cpfile = "root_files/PionTransparency_le_pip_Carbon_100k_";
	Cpfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cpfile = find_last(Cpfile);

	string C0file = "root_files/PionTransparency_le_pi0_Carbon_100k_";
	C0file += fzwork[fz] + sep + xsec[xs] + string("*.root");
	C0file = find_last(C0file);

	string Apfile = "root_files/PionTransparency_le_pip_Aluminium_100k_";
	Apfile += fzwork[fz] + sep + xsec[xs] + string("*.root");
	Apfile = find_last(Apfile);

	string A0file = "root_files/PionTransparency_le_pi0_Aluminium_100k_";
	A0file += fzwork[fz] + sep + xsec[xs] + string("*.root");
	A0file = find_last(A0file);

	if (noFile(Lpfile) or noFile(Cpfile) or noFile(Apfile) or noFile(L0file) or noFile(C0file) or noFile(A0file))
	{
		logfile<<"PionTransparency low energy ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"PionTransparency low energy ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events  = 100000;
	const double Ekmin = 0;
	const double Ekmax = 500;
	const int bins    = 10;
	
	double Ek[bins];
	double normLp[bins]; zero(normLp, bins);
	double normCp[bins]; zero(normCp, bins);
	double normAp[bins]; zero(normAp, bins);
	double normL0[bins]; zero(normL0, bins);
	double normC0[bins]; zero(normC0, bins);
	double normA0[bins]; zero(normA0, bins);
	
	for (int i = 0; i < bins; i++) Ek[i] = Ekmin + (i+1.0)*(Ekmax - Ekmin)/bins;

	double Lp[bins]; zero(Lp, bins);
	double Cp[bins]; zero(Cp, bins);
	double Ap[bins]; zero(Ap, bins);
	double L0[bins]; zero(L0, bins);
	double C0[bins]; zero(C0, bins);
	double A0[bins]; zero(A0, bins);

	calcTrans(Lpfile, Ek, Lp, normLp, bins, events, 0);
	calcTrans(Cpfile, Ek, Cp, normCp, bins, events, 0);
	calcTrans(Apfile, Ek, Ap, normAp, bins, events, 0);
	calcTrans(L0file, Ek, L0, normL0, bins, events, 0);
	calcTrans(C0file, Ek, C0, normC0, bins, events, 0);
	calcTrans(A0file, Ek, A0, normA0, bins, events, 0);
	
	for (int i = 0; i < bins; i++)
	{
		if (normLp[i] != 0) Lp[i] /= normLp[i];
		else Lp[i] = 0;
		if (normCp[i] != 0) Cp[i] /= normCp[i];
		else Cp[i] = 0;
		if (normAp[i] != 0) Ap[i] /= normAp[i];
		else Ap[i] = 0;
		if (normL0[i] != 0) L0[i] /= normL0[i];
		else L0[i] = 0;
		if (normC0[i] != 0) C0[i] /= normC0[i];
		else C0[i] = 0;
		if (normA0[i] != 0) A0[i] /= normA0[i];
		else A0[i] = 0;

		Ek[i] = Ekmin + (i+0.5)*(Ekmax - Ekmin)/bins;
	}
	
	string help = "results/PiTrans/";
	run(string("mkdir -p ") + help);
	
	string lphelp = help + string("pitrans_le_pip_lithium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string cphelp = help + string("pitrans_le_pip_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string aphelp = help + string("pitrans_le_pip_aluminium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string l0help = help + string("pitrans_le_pi0_lithium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string c0help = help + string("pitrans_le_pi0_carbon_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string a0help = help + string("pitrans_le_pi0_aluminium_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");

	ofstream fileLp(lphelp.c_str());
	ofstream fileCp(cphelp.c_str());
	ofstream fileAp(aphelp.c_str());
	ofstream fileL0(l0help.c_str());
	ofstream fileC0(c0help.c_str());
	ofstream fileA0(a0help.c_str());
		
	cout<<"Lithium (charged pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<Lp[i]<<endl;
	cout<<endl<<"Carbon (charged pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<Cp[i]<<endl;
	cout<<endl<<"Aluminium (charged pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<Ap[i]<<endl;

	cout<<endl<<"Lithium (neutral pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<L0[i]<<endl;
	cout<<endl<<"Carbon (neutral pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<C0[i]<<endl;
	cout<<endl<<"Aluminium (neutral pions): "<<endl<<endl<<"Ek [MeV] | Transparency"<<endl<<endl;
	for (int i = 0; i < bins; i++) cout<<Ek[i]<<" "<<A0[i]<<endl;
		
	fileLp<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileCp<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileAp<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileL0<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileC0<<"#Ek [GeV] | Transparency"<<endl<<endl;
	fileA0<<"#Ek [GeV] | Transparency"<<endl<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		fileLp<<Ek[i]<<" "<<Lp[i]<<endl;
		fileCp<<Ek[i]<<" "<<Cp[i]<<endl;
		fileAp<<Ek[i]<<" "<<Ap[i]<<endl;
		fileL0<<Ek[i]<<" "<<L0[i]<<endl;
		fileC0<<Ek[i]<<" "<<C0[i]<<endl;
		fileA0<<Ek[i]<<" "<<A0[i]<<endl;
	}
	
	cout<<endl<<lphelp<<" created."<<endl<<endl;
	cout<<endl<<cphelp<<" created."<<endl<<endl;
	cout<<endl<<aphelp<<" created."<<endl<<endl;
	cout<<endl<<l0help<<" created."<<endl<<endl;
	cout<<endl<<c0help<<" created."<<endl<<endl;
	cout<<endl<<a0help<<" created."<<endl<<endl;
	
	lphelp += string(" created from ") + Lpfile;
	cphelp += string(" created from ") + Cpfile;
	aphelp += string(" created from ") + Apfile;
	l0help += string(" created from ") + L0file;
	c0help += string(" created from ") + C0file;
	a0help += string(" created from ") + A0file;
	
	calclog(lphelp);
	calclog(cphelp);
	calclog(aphelp);
	calclog(l0help);
	calclog(c0help);
	calclog(a0help);
	
	fileLp.close();
	fileCp.close();
	fileAp.close();
	fileL0.close();
	fileC0.close();
	fileA0.close();
	
	return 1;
}

/*void calcFZ ()
{	
	ofstream pqfile("tmp/pion_qel.txt");
	ofstream pifile("tmp/pion_inel.txt");
	ofstream nqfile("tmp/nucleon_qel.txt");
	ofstream nifile("tmp/nucleon_inel.txt");
	
	pqfile<<"#momentum";
	pifile<<"#momentum";
	nqfile<<"#momentum";
	nifile<<"#momentum";
	
	for (int i = 0; i < nof_fz; i++)
	{
		pqfile<<" | "<<fzname[i];
		pifile<<" | "<<fzname[i];
		nqfile<<" | "<<fzname[i];
		nifile<<" | "<<fzname[i];
	}
	
	pqfile<<endl;
	pifile<<endl;
	nqfile<<endl;
	nifile<<endl;
	
	particle pion;
	particle nucleon;
	
	pion.set_mass(mass_pi);
	pion.set_pi();
	nucleon.set_mass(mass_proton);
	
	const int bins = 50;
	const int events = 1000;
	
	double fzpq[nof_fz][bins];
	double fzpi[nof_fz][bins];
	double fznq[nof_fz][bins];
	double fzni[nof_fz][bins];

	double momentum[bins];
	
	for (int i = 0; i < nof_fz; i++)
	{
		zero(fzpq[i], bins);
		zero(fzpi[i], bins);
		zero(fznq[i], bins);
		zero(fzni[i], bins);
	}
	
	vec dir = rand_dir();
	vect z4(0,0,0,0);
	params p;
	p.first_step = 1;
	
	for (int k = 0; k < bins; k++)
	{
		momentum[k] = k*100;
			
		pion.set_momentum(100.0*k*dir);
		nucleon.set_momentum(100.0*k*dir);
		
		for (int j = 0; j < nof_fz; j++)
		{
			//if (strcmp(fzwork[j].c_str(), "delta") == 0) j++;
						
			p.formation_zone = fzwork[j];
			
			for (int i = 0; i < events; i++)
			{
				fzpq[j][k] += formation1(pion, p, z4, 1, z4, 0, 2)/fermi;
				fzpi[j][k] += formation1(pion, p, z4, 0, z4, 0, 2)/fermi;
				fznq[j][k] += formation1(nucleon, p, z4, 1, z4, 0, 2)/fermi;
				fzni[j][k] += formation1(nucleon, p, z4, 0, z4, 0, 2)/fermi;
			}
			
			fzpq[j][k] /= events;
			fzpi[j][k] /= events; 
			fznq[j][k] /= events;
			fzni[j][k] /= events; 
		}
	}
	
	for (int k = 0; k < bins; k++)
	{
		pqfile<<momentum[k]<<" ";
		pifile<<momentum[k]<<" ";
		nqfile<<momentum[k]<<" ";
		nifile<<momentum[k]<<" ";
		
		for (int i = 0; i < nof_fz; i++)
		{
			//if (strcmp(fzwork[i].c_str(), "delta") == 0) i++;

			pqfile<<fzpq[i][k]<<" ";
			pifile<<fzpi[i][k]<<" ";
			nqfile<<fznq[i][k]<<" ";
			nifile<<fzni[i][k]<<" ";
		}
		
		pqfile<<endl;
		pifile<<endl;
		nqfile<<endl;
		nifile<<endl;
	}
	
	pqfile.close();
	pifile.close();
	nqfile.close();
	nifile.close();
}*/
/*
int viviNomad3(int fz, int xs)
{
	string help = "root_files/Nomad_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);

	const int events = 1000000;
	const int bins = 10;
	
	double Q2[bins];
	for (int i = 0; i < bins; i++) Q2[i] = (i + 1) * 10.0;
	
	double vivi[7][bins] = {{0}};
	double vqel[5][bins] = {{0}};
	double vspp[5][bins] = {{0}};
	double vdpp[5][bins] = {{0}};
	double vtpp[5][bins] = {{0}};
	double vnpp[5][bins] = {{0}};
	double vother[8][bins] = {{0}};
	
	double norm[bins] = {0};
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);

		if (-e1->q2()/1000000.0 < Q2[bins-1])
		{
			int which = 0;
					
			for (int l = 0; l < bins-1; l++)
			{
				if (-e1->q2()/1000000.0 > Q2[l] and -e1->q2()/1000000.0 <= Q2[l+1])
				{
					which = l;
					break;
				}
			}

			norm[which]++;
			
			for (int k = 0; k < e1->post.size(); k++)
			{							
				if (e1->post[k].pdg == -211 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 350 && e1->post[k].momentum() < 800)
				{
					int ile = e1->number_of_interactions();
					int qel = e1->post[k].his_pqel + e1->post[k].his_pcex;
					int npp = e1->post[k].his_nspp + e1->post[k].his_ndpp;				
					int spp = e1->post[k].his_pspp;				
					int dpp = e1->post[k].his_pdpp;				
					int tpp = e1->post[k].his_ptpp;				

					vivi[6][which]++;
					
					if (ile != 0)
					{
						if (npp + spp + dpp + tpp == 0)
						{
							vivi[0][which]++;
							if (qel <= 3) vqel[qel][which]++;
							else vqel[4][which]++;
						}
						else if (npp + dpp + tpp == 0)
						{
							vivi[1][which]++;
							if (qel <= 3) vspp[qel][which]++;
							else vspp[4][which]++;
						}
						else if (npp + spp + tpp == 0)
						{
							vivi[2][which]++;
							if (qel <= 3) vdpp[qel][which]++;
							else vdpp[4][which]++;
						}
						else if (npp + spp + dpp == 0)
						{
							vivi[3][which]++;
							if (qel <= 3) vtpp[qel][which]++;
							else vtpp[4][which]++;
						}
						else if (spp + dpp + tpp == 0)
						{
							vivi[4][which]++;
							if (qel <= 3) vnpp[qel][which]++;
							else vnpp[4][which]++;
						}
						else
						{
							vivi[5][which]++;
							vother[7][which]++;
							
							if (tpp + npp == 0) vother[0][which]++;
							else if (dpp + npp == 0) vother[1][which]++;
							else if (spp + npp == 0) vother[2][which]++;
							else if (dpp + tpp == 0) vother[3][which]++;
							else if (spp + tpp == 0) vother[4][which]++;
							else if (spp + dpp == 0) vother[5][which]++;
							else vother[6][which]++;
						}						
					}
				}
			}
		}
		
		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream fvivi("results/nomvivi/vivi.txt");
	ofstream fother("results/nomvivi/other.txt");
	ofstream fqel("results/nomvivi/qel.txt");
	ofstream fspp("results/nomvivi/spp.txt");
	ofstream fdpp("results/nomvivi/dpp.txt");
	ofstream ftpp("results/nomvivi/tpp.txt");
	ofstream fnpp("results/nomvivi/npp.txt");
	
	for (int i = 0; i < bins; i++)
	{
		fvivi << (i + 0.5) * 10.0;
		fother << (i + 0.5) * 10.0;
		fqel << (i + 0.5) * 10.0;
		fspp << (i + 0.5) * 10.0;
		fdpp << (i + 0.5) * 10.0;
		ftpp << (i + 0.5) * 10.0;
		fnpp << (i + 0.5) * 10.0;
		
		for (int k = 0; k < 7; k++) fvivi << " " << vivi[k][i]/norm[i];
		for (int k = 0; k < 8; k++) fother << " " << vother[k][i]/norm[i];
		for (int k = 0; k < 5; k++) fqel << " " << vqel[k][i]/norm[i];
		for (int k = 0; k < 5; k++) fspp << " " << vspp[k][i]/norm[i];
		for (int k = 0; k < 5; k++) fdpp << " " << vdpp[k][i]/norm[i];
		for (int k = 0; k < 5; k++) ftpp << " " << vtpp[k][i]/norm[i];
		for (int k = 0; k < 5; k++) fnpp << " " << vnpp[k][i]/norm[i];

		fvivi << endl;
		fother << endl;
		fqel << endl;
		fspp << endl;
		fdpp << endl;
		ftpp << endl;
		fnpp << endl;
	}
	
	fvivi.close();
	fother.close();
	fqel.close();
	fspp.close();
	fdpp.close();
	ftpp.close();
	fnpp.close();
	
	return viviNomad4(fz, xs);
}
*/

/*int viviNomad2(int fz, int xs)
{
	string help = "root_files/Nomad_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);

	const int events = 1000000;
	const int bins = 50;
	
	double Q2[bins];
	for (int i = 0; i < bins; i++) Q2[i] = (i + 1) * 2.0;
	
	double nqel[bins] = {0};
	double nspp[bins] = {0};
	double ndpp[bins] = {0};
	double pqel[bins] = {0};
	double pcex[bins] = {0};
	double pspp[bins] = {0};
	double pdpp[bins] = {0};
	double ptpp[bins] = {0};
	
	double norm[bins] = {0};
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		for (int k = 0; k < e1->post.size(); k++)
		{
			if ((e1->post[k].pdg == 211 or e1->post[k].pdg == -211 or e1->post[k].pdg == 111) and e1->post[k].p().z < 0 and e1->post[k].mother_pdg < 2500)
			{
				int which = 0;
				
				for (int l = 0; l < bins-1; l++)
				{
					if (-e1->q2()/1000000.0 > Q2[l] and -e1->q2()/1000000.0 <= Q2[l+1])
					{
						which = l;
						break;
					}
				}
				
				if (-e1->q2()/1000000.0 > Q2[bins-1]) break;
				
				norm[which]++;
				
				nqel[which] += e1->post[k].his_nqel;
				nspp[which] += e1->post[k].his_nspp;
				ndpp[which] += e1->post[k].his_ndpp;				
				pqel[which] += e1->post[k].his_pqel;				
				pcex[which] += e1->post[k].his_pcex;				
				pspp[which] += e1->post[k].his_pspp;				
				pdpp[which] += e1->post[k].his_pdpp;				
				ptpp[which] += e1->post[k].his_ptpp;				
			}
		}
		
		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream plik("results/nomvivi/nofc.txt");
	
	for (int i = 0; i < bins; i++)
	{
		plik << (i + 0.5) * 2.0 << " " << nqel[i]/norm[i] << " " << nspp[i]/norm[i] << " " << ndpp[i]/norm[i] << " " << pqel[i]/norm[i] << " " << pcex[i]/norm[i] << " " << pspp[i]/norm[i] << " " << pdpp[i]/norm[i] << " " << ptpp[i]/norm[i] << endl;
	}
	
	return 1;
}	
*/

/*int viviNomad_new(int fz, int xs)
{
	string help = "root_files/Nomad_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);

	const int events = 5000000;
	
	double x[7] = {0};
	double sum = 0;
	
	
	 //~ * 0 - in primary vertex
	 //~ * 1 - in single pion production
	 //~ * 2 - in double pion production
	 //~ * 3 - in triple pion production
	 //~ * 4 - in nucleon pion production
	 //~ * 5 - in nucleon double pp
	 //~ * 6 - other there was more than one pion production in FSI
	 
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		for (int k = 0; k < e1->post.size(); k++)
		{
			if (e1->post[k].pdg == -211 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 350 && e1->post[k].momentum() < 800)
			{
				int spp = e1->post[k].his_pspp;
				int dpp = e1->post[k].his_pdpp;
				int tpp = e1->post[k].his_ptpp;
				int nspp = e1->post[k].his_nspp;
				int ndpp = e1->post[k].his_ndpp;
				
				//int el = e1->post[k].his_pqel;				
				//int cx = e1->post[k].his_pcex;
				
				if (spp + dpp + tpp + nspp + ndpp == 0)	x[0]++;
				else if (spp == 1 and dpp + tpp + nspp + ndpp == 0) x[1]++;		
				else if (dpp == 1 and spp + tpp + nspp + ndpp == 0) x[2]++;		
				else if (tpp == 1 and spp + dpp + nspp + ndpp == 0) x[3]++;		
				else if (nspp == 1 and spp + dpp + tpp + ndpp == 0) x[4]++;		
				else if (ndpp == 1 and spp + dpp + tpp + nspp == 0) x[5]++;		
				else if (spp + dpp + tpp + nspp + ndpp > 1) x[6]++;
				
				sum++;
			}
		}

		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream plik("nomad_vivi.txt");
	
	for (int i = 0; i < 7; i++)
	{
		plik << 100*x[i]/sum << " ";
	}
	
	plik.close();
	return 1;
}
*/
int viviNomad4(int fz, int xs)
{
	string help = "root_files/Nomad_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);
	
	const int events = 5000000;
	const int bins = 10;
	
	double mom2[bins];
	for (int i = 0; i < bins; i++) mom2[i] = (i + 1.0) * 0.1;
	
	double distr[bins] = {0};
	double norma[bins] = {0};
	double norm = 0;
	double rest;
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);

		if (e1->in[0].pdg == 14)
		{
			norm++;

			for (int k = 0; k < e1->post.size(); k++)
			{
				if (e1->post[k].pdg == -211)
				{
					if (e1->post[k].p().z < 0)
					{
						double P = e1->post[k].momentum()/1000.0;
						double E = e1->post[k].E()/1000.0;
						P *= P;
						
						put(P, mom2, distr, rest, bins);//, E/P);
					}
				}
			}
		}				
		
		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream plik("results/rozklad.txt");
		
	for (int i = 0; i < bins; i++)
	{
		mom2[i] = (i + 0.5)*0.1;
		
		double P = sqrt(mom2[i]);
		double E = sqrt(mom2[i] + 0.0196);
		
		plik << mom2[i] << " " << E/P*distr[i]/norm/0.1 << endl;
	}
	
	plik.close();
	
	return 1;
}
	

int calcNomad (int fz, int xs)
{	
	//return viviNomad3(fz, xs);
	
	cout<<endl<<endl<<"Calculating Nomad ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help = "root_files/Nomad_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);
	
	if (noFile(file))
	{
		logfile<<"Nomad ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"Nomad ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events     = 5000000;
	const int prbins     = 8;
	const int pibins     = 6;
	
	double pr[prbins]; zero(pr, prbins);
	double pi[pibins]; zero(pi, pibins);
	
	double pivivi[7][pibins];
	for (int z = 0; z < 7; z++) zero(pivivi[z], pibins);
	
	double onlyqel[pibins] = {0};
	
	double normpr[prbins]; zero(normpr, prbins);
	double normpi[pibins]; zero(normpi, pibins);
	double Q2pr[prbins];
	double Q2pi[pibins];
	
	int counterPr[4] = {0, 0, 0, 0};
	int counterPi[4] = {0, 0, 0, 0};
	int countpr = 0;
	int countpi = 0;
	
	for (int i = 0; i < prbins; i++) Q2pr[i] = NOMADprotonsQ2[i] + NOMADprotonsQ2err[i];
	for (int i = 0; i < pibins; i++) Q2pi[i] = NOMADpionsQ2[i] + NOMADpionsQ2err[i];
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int ile = e1->f();
		
		double q = -e1->q2()/1000000.0;
		double rest = 0;
		
		put(q, Q2pr, normpr, rest, prbins);
		put(q, Q2pi, normpi, rest, pibins);
		
		for (int k = 0; k < ile; k++)
		{	
			if (e1->post[k].pdg == 2212 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 370 && e1->post[k].momentum() < 700)
			{
				put(q, Q2pr, pr, rest, prbins);
				countpr++;
			}
			else if (e1->post[k].pdg == -211 && e1->post[k].p().z < 0 && e1->post[k].momentum() > 350 && e1->post[k].momentum() < 800)
			{
				put(q, Q2pi, pi, rest, pibins);
				countpi++;
				
				int all = e1->number_of_interactions();
				int qel = e1->number_of_pion_elastic();
				int cex = e1->number_of_pion_ce();
				int spp = e1->number_of_pion_spp();
				int dpp = e1->number_of_pion_dpp();
				int tpp = e1->number_of_pion_tpp();
				int nspp = e1->number_of_nucleon_spp();
				int ndpp = e1->number_of_nucleon_dpp();
				
				if (all != 0)
				{
					if (qel == 0) put(q, Q2pi, pivivi[0], rest, pibins);
					if (cex == 0) put(q, Q2pi, pivivi[1], rest, pibins);
					if (spp == 0) put(q, Q2pi, pivivi[2], rest, pibins);
					if (dpp == 0) put(q, Q2pi, pivivi[3], rest, pibins);
					if (tpp == 0) put(q, Q2pi, pivivi[4], rest, pibins);
					if (nspp == 0) put(q, Q2pi, pivivi[5], rest, pibins);
					if (ndpp == 0)put(q, Q2pi, pivivi[6], rest, pibins);
					if (spp+dpp+tpp+nspp+ndpp == 0) put(q, Q2pi, onlyqel, rest, pibins);
				}
			}		
		}
		
		switch (countpr)
		{
			case 0: counterPr[0]++; break;
			case 1: counterPr[1]++; break;
			case 2: counterPr[2]++; break;
			case 3: counterPr[3]++; break;
			default: break;
		}
		
		switch (countpi)
		{
			case 0: counterPi[0]++; break;
			case 1: counterPi[1]++; break;
			case 2: counterPi[2]++; break;
			case 3: counterPi[3]++; break;
			default: break;
		}
		
		countpi = 0;
		countpr = 0;

		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int i = 0; i < prbins; i++)
	{
		if (normpr[i] != 0) pr[i] /= normpr[i];
		else pr[i] = 0;
	}

	for (int i = 0; i < pibins; i++)
	{
		if (normpi[i] != 0)
		{
			onlyqel[i] /= normpi[i];
			pi[i] /= normpi[i];
			for (int j = 0; j < 7; j++) pivivi[j][i] /= normpi[i];
		}
	}
	
	for (int i = 0; i < 4; i++)
	{
		counterPr[i] *= (double)944019.0/events;
		counterPi[i] *= (double)944019.0/events;
	}
	
	help = "results/Nomad/";
	run(string("mkdir -p ") + help);
	
	string prhelp = help + string("nomad_protons_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	string pihelp = help + string("nomad_pions_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	
	ofstream prfile(prhelp.c_str());
	ofstream pifile(pihelp.c_str());

	prfile<<"#number of events with 0,1,2,3 backwards protons: "<<counterPr[0]<<" "<<counterPr[1]<<" "<<counterPr[2]<<" "<<counterPr[3]<<endl;
	pifile<<"#number of events with 0,1,2,3 backwards pions: "<<counterPi[0]<<" "<<counterPi[1]<<" "<<counterPi[2]<<" "<<counterPi[3]<<endl;
	
	cout<<"Backwards protons: "<<endl;
	cout<<"Q2[GeV] | NOMAD data | NuWro ("<<fzname[fz]<<")"<<endl<<endl;
	prfile<<"#Q2[GeV] | NuWro ("<<fzname[fz]<<")"<<endl<<endl;
	pifile<<"#Q2[GeV] | NuWro ("<<fzname[fz]<<")"<<endl<<endl;
		
	for (int i = 0; i < prbins; i++)
	{
		cout<<NOMADprotonsQ2[i]<<" "<<NOMADprotonsBack[i]<<" "<<pr[i]<<endl;
		prfile<<NOMADprotonsQ2[i]<<" "<<pr[i]<<endl;
	}
	
	cout<<endl<<"Backwards pions: "<<endl;
	cout<<"Q2[GeV] | NOMAD data | NuWro ("<<fzname[fz]<<")"<<endl<<endl;

	for (int i = 0; i < pibins; i++)
	{
		cout<<NOMADpionsQ2[i]<<" "<<NOMADpionsBack[i]<<" "<<pi[i]<<endl;
		pifile<<NOMADpionsQ2[i]<<" "<<pi[i]<<endl;
	}
	
	cout<<endl<<prhelp<<" created."<<endl<<endl;
	cout<<endl<<pihelp<<" created."<<endl<<endl;

	help = prhelp + string(" and ") + pihelp + string(" created from ") + file;
	calclog(help);
	prfile.close();
	pifile.close();
	
	string wynik;
	
	if (fz == 0) wynik = "results/nomvivi/nofz.txt";
	else wynik = "results/nomvivi/fz.txt";
	
	ofstream ftest(wynik.c_str());
	
	for (int i = 0; i < pibins; i++)
	{
		ftest << NOMADpionsQ2[i] << " " << pi[i];
		
		for (int j = 0; j < 7; j++)
		{
			ftest << " " << pivivi[j][i];
		}
		
		ftest << " " << onlyqel[i] << endl;
	}
	
	return 1; //viviNomad(fz, xs);
}

int calcAtmNu (int fz, int xs)
{
	cout<<endl<<endl<<"Calculating Atmospheric Neutrinos ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
	
	string help = "root_files/AtmNu_100k_";
	help += fzwork[fz] + sep + xsec[xs] + string("*.root");
	
	string file = find_last(help);
	
	if (noFile(file))
	{
		logfile<<"Atmospheric Neutrinos ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		cout<<"Atmospheric Neutrinos ("<<fzname[fz]<<"): there are no root files, make simulation first!"<<endl<<endl;
		return 1;
	}
		
	get_date();
	
	const int events     = 100000;
	const int bins = 8;
	double count[2][bins];
	for (int i = 0; i < bins; i++) {count[0][i] = 0; count[1][i] = 0;}
	
	/* 0 - 0pi
	 * 1 - pi+
	 * 2 - pi0
	 * 3 - pi-pi+
	 * 4 - 2pi0
	 * 5 - pi+ n*pi0 n>0
	 * 6 - 2pi+ n*pi0 n>=0
	 * 7 - pi-
	 */
	
	string names[8] = {"0pi", "pi+", "pi0", "pi-pi+", "2pi0", "pi+(>0*pi0)", "2pi+(>=0*pi0)", "pi-"};
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);

	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pions = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		int pionsfsi = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		switch (pions){
			case   0: count[0][0]++; break;
			case 100: count[0][1]++; break;
			case   1: count[0][2]++; break;
			case 110: count[0][3]++; break;
			case   2: count[0][4]++; break;
			case  10: count[0][7]++; break;
			default: break;
		}
		
		if (e1->nof(211) == 1 && e1->nof(111) > 0) count[0][5]++;
		else if (e1->nof(211) == 2 && e1->nof(111) >= 0) count[0][6]++;
		
		switch (pionsfsi){
			case   0: count[1][0]++; break;
			case 100: count[1][1]++; break;
			case   1: count[1][2]++; break;
			case 110: count[1][3]++; break;
			case   2: count[1][4]++; break;
			case  10: count[1][7]++; break;
			default: break;
		}
		
		if (e1->fof(211) == 1 && e1->fof(111) > 0) count[1][5]++;
		else if (e1->fof(211) == 2 && e1->fof(111) >= 0) count[1][6]++;
		
		cout<<file<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<file<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;

	for (int i = 0; i < bins; i++)
	{
		count[0][i]/=events;
		count[1][i]/=events;
	}
	
	help = "results/AtmNu/"; 
	run(string("mkdir -p ") + help);
	
	help = help + string("atmnu_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	
	ofstream afile(help.c_str());
		
	for (int i = 0; i < bins; i++)
	{
		cout<<names[i]<<" "<<count[0][i]<<" "<<count[1][i]<<endl;
		afile<<names[i]<<" "<<count[0][i]<<" "<<count[1][i]<<endl;
	}
	
	cout<<endl<<help<<" created."<<endl<<endl;
	
	help = help + string(" created from ") + file;
	calclog(help);
	afile.close();
	return 1;	
}

void xsectable (string name, double *en, double xsec[][5], int bins)
{
	string tex = name + "_cc_total.tex";
	string dvi = name + "_cc_total.dvi";
	string ps  = name + "_cc_total.ps";
	string pdf = name + "_cc_total.pdf";
	string txt = name + "_cc_total.txt";
	
	string ev = "GeV";
	string title;
	
	if (strcmp(name.c_str(), "sb") == 0) {title = "SciBooNE CC total cross section on $C_8H_8$"; ev = "MeV";}
	if (strcmp(name.c_str(), "nomad") == 0) title = "NOMAD CC total cross section on carbon";
	if (strcmp(name.c_str(), "minos") == 0) title = "MINOS CC total cross section on iron";
	if (strcmp(name.c_str(), "minos_anti") == 0) title = "MINOS CC total cross section on iron (anti-neutrino)";
	 
	ofstream texfile(tex.c_str());
	ofstream txtfile(txt.c_str());
	

	texfile<<fixed;
	texfile<<setprecision (2);
	texfile<<"\\documentclass[titlepage]{article}"<<endl<<endl<<"\\usepackage[margin = 0in, tmargin=0.5in, portrait]{geometry}"<<endl<<"\\pagestyle{empty}"<<endl;
	texfile<<"\\usepackage{array}"<<endl<<"\\newcolumntype{j}{m{50pt}}"<<"\\newcolumntype{R}{>{\\hfill\\arraybackslash}m{2cm}<{}}"<<endl<<"\\usepackage{multirow}"<<endl;
	
	texfile<<endl<<"\\begin{document}"<<endl;

	texfile<<"\\begin{table}[!ht]"<<endl;
	texfile<<"\\begin{center}"<<endl;
	texfile<<"\\begin{tabular}{l||R|R|R|R|R}"<<endl;
	texfile<<"\\multirow{2}{*}{Energy ["<<ev<<"]} & \\multicolumn{5}{c}{Cross section [$10^{-38}cm^{2}/nucleon$]} \\tabularnewline\\cline{2-6}"<<endl;
	texfile<<" & QE & RES & DIS & COH & Total \\\\"<<endl;
	texfile<<"\\hline \\hline"<<endl;
	
	txtfile<<"#Energy ["<<ev<<"] | QE | RES | DIS | COH | Total"<<endl;
	
	for (int i = 0; i < bins; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
		
		texfile << energy << " & " << xsec[i][0]*1e38 <<  " & " << xsec[i][1]*1e38 <<  " & " << xsec[i][2]*1e38 <<  " & " << xsec[i][3]*1e38 <<  " & " << xsec[i][4]*1e38 << " \\\\" << endl;
		txtfile << energy << " " << xsec[i][0]*1e38 <<  " " << xsec[i][1]*1e38 <<  " " << xsec[i][2]*1e38 <<  " " << xsec[i][3]*1e38 <<  " " << xsec[i][4]*1e38 << endl;		
	}

	texfile<<"\\end{tabular}"<<endl<<"\\caption{"<<title<<"}"<<"\\end{center}"<<endl<<"\\end{table}"<<endl<<"\\end{document}"<<endl;
	texfile.close();
	txtfile.close();
	
	run(string("latex ") + tex);
	run(string("dvips ") + dvi);
	run(string("ps2pdf ") + ps);
	
	string dir = "results/" + name;
	string texdir = "results/" + name + "/tex/";
	
	run(string("mkdir -p ") + dir);
	run(string("mkdir -p ") + texdir);
	run(string("mv ") + pdf + string(" ") + dir); 
	run(string("mv ") + txt + string(" ") + dir); 
	run(string("mv ") + name + string("_cc_total*") + string (" ") + texdir);
}

int calcSBCCtotal()
{
	double en[6] = {380, 620, 870, 1110, 1430, 2470};
	double Hxsec[6][5];
	double Cxsec[6][5];
	double xsec[6][5];	

	for (int i = 0; i < 6; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
		
		string Hfile = "root_files/E" + energy + string("_SB_CC_Hydrogen_0k_*.txt");
		Hfile = find_last(Hfile);
		
		string Cfile = "root_files/E" + energy + string("_SB_CC_Carbon_0k_*.txt");
		Cfile = find_last(Cfile);
		
		crosssection(Hfile, Hxsec[i], 0);
		crosssection(Cfile, Cxsec[i], 0);
		
		for (int k = 0; k < 5; k++)
		{
			xsec[i][k] = (8.0*Hxsec[i][k] + 8.0*12.0*Cxsec[i][k])/(8.0+8.0*12.0);
		}	
	}
	
	xsectable("sb", en, xsec, 6);
	
	return 1;
}

int calcNOMADCCtotal()
{
	double en[30] = {4.6, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 16.2, 18.7, 21.2, 23.7, 26.2, 28.7, 32.3, 37.3, 42.4, 47.4, 54.6, 64.7, 74.8, 84.8, 94.8, 107, 122, 136.9, 165.9, 228.3};
	double xsec[30][5];	

	for (int i = 0; i < 30; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i]*1000.0;
		temp >> energy;
		
		string file = "root_files/E" + energy + string("_NOMAD_CC_Carbon_0k_*.txt");
		file = find_last(file);
		
		crosssection(file, xsec[i], 0);
	}
	
	xsectable("nomad", en, xsec, 30);
	
	return 1;
}

int calcMINOSCCtotal()
{
	double en[13] = {3.48, 4.45, 5.89, 7.97, 10.45, 13.43, 16.42, 19.87, 23.88, 27.89, 32.81, 38.87, 45.77};
	double enbar[11] = {6.07, 7.99, 10.43, 13.42, 16.41, 19.82, 23.82, 27.84, 32.72, 38.74, 45.61};
	double xsec[13][5];	
	double xsecbar[11][5];

	for (int i = 0; i < 13; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i]*1000.0;
		temp >> energy;
		
		string file = "root_files/E" + energy + string("_MINOS_CC_Iron_0k_*.txt");
		file = find_last(file);
		
		crosssection(file, xsec[i], 0);
		
		if (i < 11)
		{
			stringstream tempbar;
			string energybar;
		
			tempbar << enbar[i]*1000.0;
			tempbar >> energybar;
			
			string filebar = "root_files/E" + energybar + string("_MINOS_CC_Iron_0k_anti*.txt");
			filebar = find_last(filebar);
			
			crosssection(filebar, xsecbar[i], 0);
		}
	}
	
	xsectable("minos", en, xsec, 13);
	xsectable("minos_anti", enbar, xsecbar, 11);
	
	return 1;
}

int calcsfg()
{
	ofstream res("results/sfg.txt");
	
	for (int i = 5; i <= 30; i++)
	{
		stringstream temp;
		string energy;
		
		temp << i*50;
		temp >> energy;
		
		string file  = "root_files/E" + energy + string("_carbon_fg.root.txt");
		string file2 = "root_files/E" + energy + string("_carbon_sf.root.txt");
		
		double ratio = crosssection(file2)/crosssection(file);
		
		res << i*50 << " " << ratio << endl;
	}
	
	res.close();
	
	return 1;
}

int calcbg()
{	
	cout<<endl<<endl<<"Calculating backg: "<<endl<<endl;
	
	string file = "NIWG/tech/root_files/numu_oxygen_fg.root";
	
	const int events     = 100000;
	
	double ccqe = 0;
	double ccpi = 0;
	double ncpi = 0;
	
	TFile *tf1 = new TFile(file.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion  = e1->nof(211) + e1->nof(-211) + e1->nof(111);
		int piona = e1->fof(211) + e1->fof(-211) + e1->fof(111);
		
		if (e1->flag.cc and e1->flag.qel and piona == 0) ccqe++;
		else if (e1->flag.nc and e1->fof(211) == 0 and e1->fof(-211) == 0 and e1->fof(111) == 1) ncpi++;
		else if (e1->flag.cc and pion == 1 and piona == 0) ccpi++;
		
	}		
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream plik("results/background.txt");
	
	plik << "# ccqe / ccpi / ncpi "<<endl<<endl;
	plik << ccqe/events*100.0 << " " << ccpi/events*100.0 << " " << ncpi/events*100.0;
	
	plik.close();
	return 1;	
}

int calcMBCCrat (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating MB CC ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
		
	get_date();
	
	const int events     = 10000;
	const int bins       = 13;
		
	string help = "results/MBratio/";
	run(string("mkdir -p ") + help);
	help += string("mb_ratio_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file(help.c_str());
	
	file << "#neutrino energy | ccqe-like (1c) | ccpi+-like (1b) | ccpi+-like/ccqe-like (1a) | ccpi0 (5)| ccqe-like/total (6-1) | ccpi+-like/total (6-2) | ccqe-like/ccqe-true (7) | ccpi+-like/ccpi+-true" << endl;

	for (int i = 0; i < bins; i++)
	{

		double Hqe = 0;
		double Hpi = 0;
		double Cqe = 0;
		double Cpi = 0;
		double Htot = 0;
		double Ctot = 0;
		double Hpi0 = 0;
		double Cpi0 = 0;
		double Hqetrue = 0;
		double Cqetrue = 0;
		double Hpitrue = 0;
		double Cpitrue = 0;
		
		double res = 0;
		
		double en = mbraten[i]*1000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string Hfile = "root_files/E" + energy + string("_MB_CC_Hydrogen_100k_*.root");
		string Cfile = "root_files/E" + energy + string("_MB_CC_Carbon_100k_") + fzwork[fz] + sep + xsec[xs] + string("*.root");
		
		Hfile = find_last(Hfile);
		Cfile = find_last(Cfile);
				
		double crossH = crosssection(Hfile + string(".txt"));
		double crossC = crosssection(Cfile + string(".txt"));
		
		double total = (12.0*crossC + 2.0*crossH);
		
		TFile *tf1 = new TFile(Hfile.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
			
		tt1->SetBranchAddress("e",&e1);
		
		for (int k = 0; k < events; k++)
		{
			tt1->GetEntry(k);
			
			if (e1->flag.qel) Hqetrue++;
						
			int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);

			if (pion == 0) Hqe++;
			else if (pion == 100) Hpi++;
			else if (pion == 1) Hpi0++;

			cout<<Hfile<<": "<<100*k/events<<"%\r"<<flush;
		}
		
		cout<<Hfile<<": done"<<endl<<endl;
			
		delete e1;
		delete tt1;
		delete tf1;
		
		TFile *tf2 = new TFile(Cfile.c_str());
		TTree *tt2 = (TTree*)tf2->Get("treeout");
		event *e2   = new event();
			
		tt2->SetBranchAddress("e",&e2);
		
		for (int k = 0; k < events; k++)
		{
			tt2->GetEntry(k);
			
			if (e2->flag.qel) Cqetrue++;
			
			int pion = 100*e2->fof(211) + 10*e2->fof(-211) + e2->fof(111);
			int pion0 = 100*e2->nof(211) + 10*e2->nof(-211) + e2->nof(111);
					
			if (pion == 0) Cqe++;
			else if (pion == 100) Cpi++;
			else if (pion == 1) Cpi0++;
			
			if (pion0 == 100) Cpitrue++;
			
			cout<<Cfile<<": "<<100*i/events<<"%\r"<<flush;
		}
		
		cout<<Cfile<<": done"<<endl<<endl;
			
		delete e2;
		delete tt2;
		delete tf2;

		Cqe  *= 12.0*crossC; 
		Hqe  *=  2.0*crossH;
		Cpi  *= 12.0*crossC;
		Hpi  *=  2.0*crossH;
		Cpi0 *= 12.0*crossC;
		Hpi0 *=  2.0*crossH;
		Cqetrue *= 12.0*crossC;
		Hqetrue *=  2.0*crossH;
		Hpitrue = Hpi;
		Cpitrue *= 12.0*crossC;
		
		double qe  = (Cqe + Hqe)/events;
		double pi  = (Cpi + Hpi)/events;
		double pi0 = (Cpi0 + Hpi0)/events;
		double qetrue = (Cqetrue + Hqetrue)/events;
		double pitrue = (Cpitrue + Hpitrue)/events;
								
		file << mbraten[i] << " " << qe << " " << pi << " " << pi/qe << " " << pi0 << " " << qe/total << " " << pi/total << " " << qe/qetrue << " " << pi/pitrue << " " << qetrue << " " << pitrue << endl;
	}	
	
	cout<<endl<<help<<" created."<<endl<<endl;
	file.close();

	return 1;	
}

int calcnuintTr (int fz, int xs)
{	
	cout<<endl<<endl<<"Calculating NUINT11 pr trans ("<<fzname[fz]<<", "<<xsec[xs]<<"): "<<endl<<endl;
		
	get_date();
	
	const int events     = 1000000;
	const int bins       = 35;
	
	double momentum[bins];
	for (int i = 0; i < bins; i++) momentum[i] = (i+6)*50.0;	
	
	double some[bins]; zero(some, bins);
	double pip[bins]; zero(pip, bins);
	double pim[bins]; zero(pim, bins);
	double pi0[bins]; zero(pi0, bins);
	double mpi[bins]; zero(mpi, bins);
	double norm[bins]; zero(norm, bins);
	double rest;
	
	string help = "results/NUINT/";
	run(string("mkdir -p ") + help);
	help += string("proton_") + fzwork[fz] + sep + xsec[xs] + sep + date + string(".txt");
	ofstream file(help.c_str());

	string Cfile = "root_files/Proton_1000_2000_Carbon_1m_" + fzwork[fz] + sep + xsec[xs] + string("*.root");
	Cfile = find_last(Cfile);
						
	TFile *tf1 = new TFile(Cfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);
		
		double mom = e1->out[0].momentum();
		put(mom, momentum, norm, rest, bins);
		
		int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		if (e1->number_of_interactions() != 0)
		{
			put(mom, momentum, some, rest, bins);
		
			if (pion == 100) put(mom, momentum, pip, rest, bins);
			else if (pion == 10) put(mom, momentum, pim, rest, bins);
			else if (pion == 1) put(mom, momentum, pi0, rest, bins);
			else if (pion != 0) put(mom, momentum, mpi, rest, bins);
		}
		
		cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
	}
		
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
	
	file<<"#proton momentum / there was interaction / there was one pi / there was one pi+ / there was one pi- / there was one pi0 / there was more than one pi"<<endl;
		
	for (int i = 1; i < bins; i++)
	{
		if (norm[i] != 0)
		{
			some[i] /= norm[i];
			pip[i]  /= norm[i];
			pim[i]  /= norm[i];
			pi0[i]  /= norm[i];
			mpi[i]  /= norm[i];
		}
		
		file << (i+5.5)*50.0 << " " << some[i] << " " << pip[i] + pim[i] + pi0[i] << " " << pip[i] << " " << pim[i] << " " << pi0[i] << " " << mpi[i] << endl;
	}
		
	cout<<endl<<help<<" created."<<endl<<endl;
	file.close();

	return 1;	
}

int roman_calc ()
{		
	const int events     = 100000;
	const int bins       = 21;
	
	double en[bins] = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
	double N1[bins]; zero(N1, bins);
	double N2[bins]; zero(N2, bins);
	double N3[bins]; zero(N3, bins);
	
	for (int i = 0; i < bins; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
		
		string Cfile = "root_files/E" + energy + "_carbon_100k.root";
							
		TFile *tf1 = new TFile(Cfile.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
		
		tt1->SetBranchAddress("e",&e1);

		for (int k = 0; k < events; k++)
		{
			tt1->GetEntry(k);
		
			int pion   = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
			int proton = e1->fof(2212);
			int slowpr = 0;
			
			for (int z = 0; z < e1->f(); z++)
			{
				if (e1->post[z].pdg == 2212 and e1->post[z].Ek() > 100) slowpr++;
			}
			
			if (pion == 0)
			{
				N1[i]++;
				if (proton == 0) N3[i]++;
				if (slowpr == 0) N2[i]++;
			}		
			cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
		}
		
		cout<<Cfile<<": done"<<endl<<endl;
			
		delete e1;
		delete tt1;
		delete tf1;
	}
	
	ofstream file("results/nuwro_fsi2.txt");
	file<<"#energy | N1 | N2 | N3 | N2/N1 | N3/N1"<<endl;
		
	for (int i = 0; i < bins; i++)
	{
		file << en[i] << " " << N1[i] << " " << N2[i] << " " << N3[i] << " " << N2[i]/N1[i] << " " << N3[i]/N1[i] << endl;
	}
		
	cout<<endl<<"results/nuwro_fsi2.txt created."<<endl<<endl;
	file.close();

	return 1;	
}

int roman_3b ()
{		
	const int events     = 100000;
	const int bins       = 100;
	
	double en[3] = {500, 1000, 2000};
	double mom[bins];
	
	for (int i = 0; i < 3; i++)
	{
		double b[bins]; zero(b, bins);
		
		double cp1[bins]; zero(cp1, bins);
		double cp2[bins]; zero(cp2, bins);
		
		double co1[bins]; zero(co1, bins);
		double co2[bins]; zero(co2, bins);
		
		double rest;

		for (int l = 0; l < bins; l++) mom[l] = (l + 1)*0.025;

		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
		
		string Cfile = "root_files/E" + energy + "_carbon_100k.root";
							
		TFile *tf1 = new TFile(Cfile.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
		
		tt1->SetBranchAddress("e",&e1);

		for (int k = 0; k < events; k++)
		{
			tt1->GetEntry(k);
			
			int pion   = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
			
			if (e1->nof(2212) != 0 and e1->fof(2212) != 0)
			{
				double momentum = 0;
				
				for (int j = 0; j < e1->f(); j++)
				{
					double x = e1->post[j].momentum()/1000.0;
					if (e1->post[j].pdg == 2212 and x > momentum) momentum = x;
				}
				
				put(momentum, mom, cp2, rest, bins);
			
				momentum = 0;
				
				for (int j = 0; j < e1->n(); j++)
				{
					double x = e1->out[j].momentum()/1000.0;
					if (e1->out[j].pdg == 2212 and x > momentum) momentum = x;
				}
				
				put(momentum, mom, co2, rest, bins);
				
				if (pion == 0)
				{
					momentum = 0;
				
					for (int j = 0; j < e1->f(); j++)
					{
						double x = e1->post[j].momentum()/1000.0;
						if (e1->post[j].pdg == 2212 and x > momentum) momentum = x;
					}
					
					put(momentum, mom, cp1, rest, bins);
					put(momentum, mom, b, rest, bins);
					
					momentum = 0;
					
					for (int j = 0; j < e1->n(); j++)
					{
						double x = e1->out[j].momentum()/1000.0;
						if (e1->out[j].pdg == 2212 and x > momentum) momentum = x;
					}
					
					put(momentum, mom, co1, rest, bins);
				}			
			}

			cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
		}
		
		cout<<Cfile<<": done"<<endl<<endl;
			
		delete e1;
		delete tt1;
		delete tf1;
		
		string txt1 = "results/nuwro_3b_" + energy + ".txt";
		string txt2 = "results/nuwro_3c_all_" + energy + ".txt";
		string txt3 = "results/nuwro_3c_nopi_" + energy + ".txt";
		
		ofstream file1(txt1.c_str());
		ofstream file2(txt2.c_str());
		ofstream file3(txt3.c_str());

		file1<<"#momentum [GeV] | number of events"<<endl;
		file2<<"#momentum [GeV] | interaction point | final state"<<endl;
		file3<<"#momentum [GeV] | interaction point | final state"<<endl;
		
		for (int l = 0; l < bins; l++) mom[l] = (l + 0.5)*0.025;
			
		for (int l = 0; l < bins; l++)
		{
			file1 << mom[l] << " " << b[l] << endl;
			file2 << mom[l] << " " << co2[l] << " " << cp2[l] << endl;
			file3 << mom[l] << " " << co1[l] << " " << cp1[l] << endl;
		}

		file1.close();
		file2.close();
		file3.close();
	}
	
	return 1;	
}

int roman_3d ()
{		
	const int events     = 100000;
	const int bins       = 100;
	
	double energy[bins];
	double nofeve1[bins]; zero(nofeve1, bins);
	double nofeve2[bins]; zero(nofeve2, bins);

	double rest;
	
	int counter = 0;

	for (int l = 0; l < bins; l++) energy[l] = (l + 1)*0.012;
	
	double nof_pr1[8]; zero(nof_pr1, 8);
	double nof_pr2[8]; zero(nof_pr2, 8);

	string Cfile = "root_files/E1000_carbon_100k.root";
							
	TFile *tf1 = new TFile(Cfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);

		int pion   = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		if (e1->flag.qel and pion == 0)
		{
			counter++;
			
			int ile = e1->fof(2212);
			nof_pr1[ile]++;
			
			ile = 0;
			
			double energia1 = e1->out[0].energy()/1000.0;
			double energia2	= energia1;		
			for (int a = 0; a < e1->f(); a++)
			{
				if (e1->post[a].pdg == 2212)
				{
					energia1 += e1->post[a].Ek()/1000.0;
					energia2 += e1->post[a].Ek()/1000.0;
				}
				else if (e1->post[a].pdg == 2112) energia2 += e1->post[a].Ek()/1000.0;
				 
				if (e1->post[a].pdg == 2212 and e1->post[a].Ek() > 100) ile++;
			}
			
			put(energia1, energy, nofeve1, rest, bins);
			put(energia2, energy, nofeve2, rest, bins);
			
			nof_pr2[ile]++;
		}
			
		cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
	}
		
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
		
	string txt1 = "results/nuwro_3d.txt";
	string txt2 = "results/nuwro_3e.txt";
		
	ofstream file1(txt1.c_str());
	ofstream file2(txt2.c_str());

	file1<<"#number of protons | number of events | number of events (for Tk > 100)"<<endl;
	file2<<"#energy | number of events (muon+protons) | number of events (muon+protons+neutrons)"<<endl;
	
	for (int i = 0; i < 8; i++)
	{
		nof_pr1[i] *= 1000.0/counter;
		nof_pr2[i] *= 1000.0/counter;
	}
	
	for (int i = 0; i < bins; i++)
	{
		nofeve1[i] *= 1000.0/counter;
		nofeve2[i] *= 1000.0/counter;
	}
		
	for (int l = 0; l < bins; l++) energy[l] = (l + 0.5)*0.012;
	
	for (int i = 0; i < 8; i++) file1<<i<<" "<<nof_pr1[i]<<" "<<nof_pr2[i]<<endl;
	
	for (int i = 0; i < bins; i++) file2<<energy[i]<<" "<<nofeve1[i]<<" "<<nofeve2[i]<<endl;

	file1.close();
	file2.close();

	return 1;	
}

int roman_7b ()
{		
	const int events     = 100000;
	const int bins       = 100;
	
	double momentum[bins];
	double angle[bins];

	double dist1[bins][bins]; 
	double dist2[bins][bins];
	
	int count1 = 0;
	int count2 = 0;
	
	for (int i = 0; i < bins; i++)
	{
		momentum[i] = (i+1)*10.0;
		angle[i]    = (i+1)*0.02 - 1;
		
		zero(dist1[i], bins);
		zero(dist2[i], bins);
	}

	double rest;
	
	string Cfile = "root_files/E1000_carbon_100k.root";
							
	TFile *tf1 = new TFile(Cfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);
		
		double mom = e1->out[0].momentum();
		double ang = e1->out[0].p().z/mom;
		
		int wmom = bins - 1;
		int wang = bins - 1;
				
		for (int l = 0; l < bins-1; l++)
		{
			if (mom > momentum[l] and mom <= momentum[l+1])
			{
				wmom = l;
				break;
			}
		}
		
		for (int l = 0; l < bins-1; l++)
		{
			if (ang > angle[l] and ang <= angle[l+1])
			{
				wang = l;
				break;
			}
		}
		
		if (e1->flag.qel)
		{
			dist1[wmom][wang]++;
			count1++;
		}
		
		if (100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111) == 0)
		{
			dist2[wmom][wang]++;		
			count2++;
		}
			
		cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
	}
		
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
		
	string txt1 = "results/nuwro_7d1.txt";
	string txt2 = "results/nuwro_7d2.txt";
		
	ofstream file1(txt1.c_str());
	ofstream file2(txt2.c_str());
	
	for (int i = 0; i < bins; i++)
	{
		momentum[i] = (i+0.5)*10.0;
		angle[i]    = (i+0.5)*0.02 - 1;
	}
	
	for (int i = 0; i < bins; i++)
	{
		for (int k = 0; k < bins; k++)
		{
			file1 << momentum[i] << " " << angle[k] << " " << dist1[i][k] << endl;
			file2 << momentum[i] << " " << angle[k] << " " << dist2[i][k] << endl;
		}
		
		file1<<endl;
		file2<<endl;
	}

	file1.close();
	file2.close();

	return 1;	
}

int roman_extra ()
{		
	const int events     = 100000;
	const int bins       = 150;
	
	double energy[bins];
	double nofeve1[bins]; zero(nofeve1, bins);
	double nofeve2[bins]; zero(nofeve2, bins);

	double rest;
	
	int counter1 = 0;
	int counter2 = 0;

	for (int l = 0; l < bins; l++) energy[l] = (l + 1)*0.001 + 0.9;

	string Cfile = "root_files/E1000_carbon_100k.root";
							
	TFile *tf1 = new TFile(Cfile.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);

		int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		int nucl = 100*e1->fof(2212) + e1->fof(2112);
		
		if (pion == 0)
		{
			double energia1 = e1->out[0].energy()/1000.0;
			double energia2	= energia1;
	
			if (nucl == 100)
			{
				counter1++;
						
				for (int a = 0; a < e1->f(); a++)
				{
					if (e1->post[a].pdg == 2212) energia1 += e1->post[a].Ek()/1000.0;
				}

				put(energia1, energy, nofeve1, rest, bins);

			}
			else if (nucl == 200)
			{
				counter2++;
						
				for (int a = 0; a < e1->f(); a++)
				{
					if (e1->post[a].pdg == 2212) energia2 += e1->post[a].Ek()/1000.0;
				}

				put(energia2, energy, nofeve2, rest, bins);
			}
		}
			
		cout<<Cfile<<": "<<100*k/events<<"%\r"<<flush;
	}
		
	cout<<Cfile<<": done"<<endl<<endl;
		
	delete e1;
	delete tt1;
	delete tf1;
		
	string txt1 = "results/other/nuwro_extra.txt";
		
	ofstream file1(txt1.c_str());

	file1<<"#energy | mup | mupp "<<endl;
	
	
	for (int i = 0; i < bins; i++)
	{
		nofeve1[i] *= 1000.0/counter1;
		nofeve2[i] *= 1000.0/counter2;
	}
		
	for (int l = 0; l < bins; l++) energy[l] = (l + 0.5)*0.001 + 0.9;
	
	for (int i = 0; i < bins; i++) file1<<energy[i]<<" "<<nofeve1[i]<<" "<<nofeve2[i]<<endl;

	file1.close();

	return 1;	
}

/*void oset_count (string filename, double *n, double *en, double *en2, int events)
{
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double norm1 = 0;
	double norm2 = 0;
	
	double rest = 0;
	
	double energy[25];
	
	for (int i = 0; i < 25; i++) energy[i] = (i+1)*10.0;

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);
		
		if (e1->nopp and !e1->abs)
		{
			if (e1->nofi != 0) norm1++;
			
			double tk;
			
			for (int z = 0; z < e1->f(); z++) if (e1->post[z].pdg == 111 or e1->post[z].pdg == 211 or e1->post[z].pdg == -211) tk = e1->post[z].Ek();
			
			if (e1->nofi != 0) put (tk, energy, en2, rest, 25);
			
			switch (e1->nofi)
			{
				case 1: n[0]++; break;
				case 2: n[1]++; break;
				case 3: n[2]++; break;
				case 4: n[3]++; break;
				default: break;
			}
		}
		else if (e1->nopp and e1->abs)
		{
			norm2++;
			
			put(e1->pabsen, energy, en, rest, 25);
			
			switch (e1->nofi)
			{
				case 0: n[4]++; break;
				case 1: n[5]++; break;
				case 2: n[6]++; break;
				case 3: n[7]++; break;
				case 4: n[8]++; break;
				default: break;
			}
		}
		
		cout<<filename<<": "<<100*k/events<<"%\r"<<flush;
	}
		
	cout<<filename<<": done"<<endl<<endl;
			
	delete e1;
	delete tt1;
	delete tf1;

	for (int i = 0; i < 4; i++) n[i] /= norm1;
	for (int k = 4; k < 9; k++) n[k] /= norm2;
}
*/

/*void oset_impac(string filename, double *qel, double *abs, int events)
{
	double bin[28];
	for (int i = 0; i < 28; i++) bin[i] = (i+1)*0.25;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double norm[28]; zero(norm, 28);
	double rest = 0;

	for (int k = 0; k < events; k++)
	{
		tt1->GetEntry(k);
		
		if (e1->nopp)
		{
			double x = e1->out[0].r.x;
			double y = e1->out[0].r.y;

			double b = sqrt(x*x + y*y)/fermi;

			put(b, bin, norm, rest, 28);	
		
			if (e1->abs) put(b, bin, abs, rest, 28);
			else if (e1->nofi != 0) put(b, bin, qel, rest, 28);			
		}
	}
	
	for (int i = 0; i < 28; i++)
	{
		qel[i] /= norm[i];
		abs[i] /= norm[i];
	}
}
*/

/*int oset_calc ()
{		
	const int events     = 1000000;

	double n1[9]; zero(n1, 9);
	double n2[9]; zero(n2, 9);
	double x[28];
	double qel[28]; zero(qel, 28);
	double abs[28]; zero(abs, 28);
	double en1[25]; zero(en1, 25);
	double en2[25]; zero(en2, 25);
	double en3[25]; zero(en3, 25);
	double en4[25]; zero(en4, 25);
	double energy[25];
	
	for (int i = 0; i < 25; i++) energy[i] = (i+0.5)*10.0; 
	for (int i = 0; i < 28; i++) x[i] = (i+0.5)*0.25;
	
	string file1 = "root_files/pionE85_1m.root";
	string file2 = "root_files/pionE245_1m.root";
	string file3 = "root_files/pionE180_1m.root";
	
	ofstream out1("results/oset_85.txt");
	ofstream out2("results/oset_245.txt");
	
	oset_count(file1, n1, en1, en3, events);
	oset_count(file2, n2, en2, en4, events);
	oset_impac(file3, qel, abs, events);
		
	out1<<fixed;
	out1<<setprecision(2);
	out1<<"Probability that quasielastic scattering proceeds through 1 | 2 | 3 | 4 collisions:"<<endl<<endl;						
	out1<<n1[0]<<" "<<n1[1]<<" "<<n1[2]<<" "<<n1[3]<<endl<<endl;
	out1<<"Probability that absorption occurs after 0 | 1 | 2 | 3 | 4 quasielastic collisions:"<<endl<<endl;						
	out1<<n1[4]<<" "<<n1[5]<<" "<<n1[6]<<" "<<n1[7]<<" "<<n1[8];
	out1.close();

	out2<<fixed;
	out2<<setprecision(2);
	out2<<"Probability that quasielastic scattering proceeds through 1 | 2 | 3 | 4 collisions:"<<endl<<endl;						
	out2<<n2[0]<<" "<<n2[1]<<" "<<n2[2]<<" "<<n2[3]<<endl<<endl;
	out2<<"Probability that absorption occurs after 0 | 1 | 2 | 3 | 4 quasielastic collisions:"<<endl<<endl;						
	out2<<n2[4]<<" "<<n2[5]<<" "<<n2[6]<<" "<<n2[7]<<" "<<n2[8];
	out2.close();
		
	cout<<endl<<"results/oset*.txt created."<<endl<<endl;
	
	ofstream texfile("osetcomp.tex");
	
	texfile<<fixed;
	texfile<<setprecision (2);
	texfile<<"\\documentclass[titlepage]{article}"<<endl<<endl<<"\\usepackage[margin = 0in, tmargin=0.5in, portrait]{geometry}"<<endl<<"\\pagestyle{empty}"<<endl;
	texfile<<"\\usepackage{array}"<<endl<<"\\newcolumntype{j}{m{50pt}}"<<"\\newcolumntype{R}{>{\\hfill\\arraybackslash}m{2cm}<{}}"<<endl<<"\\usepackage{multirow}"<<endl;
	
	texfile<<endl<<"\\begin{document}"<<endl;

	texfile<<"\\begin{table}[!ht]"<<endl;
	texfile<<"\\begin{center}"<<endl;
	texfile<<"\\begin{tabular}{lcccccc}"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<" $T_k [MeV]$ & & & n = 1 & n = 2 & n = 3 & n = 4 \\\\ \\hline"<<endl;
	texfile<<"\\multirow{2}{*}{85} & Oset & & "<<pq1[0]<<" & "<<pq1[1]<<" & "<<pq1[2]<<" & "<<pq1[3]<<" \\\\"<<endl;
	texfile<<" & NuWro & & "<<n1[0]<<" & "<<n1[1]<<" & "<<n1[2]<<" & "<<n1[3]<<" \\\\"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<"\\multirow{2}{*}{245} & Oset & & "<<pq2[0]<<" & "<<pq2[1]<<" & "<<pq2[2]<<" & "<<pq2[3]<<" \\\\"<<endl;
	texfile<<" & NuWro & & "<<n2[0]<<" & "<<n2[1]<<" & "<<n2[2]<<" & "<<n2[3]<<" \\\\"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<"\\end{tabular}"<<endl;
	texfile<<"\\caption{Probability that the quasielastic scattering proceeds through $n$ collisions}"<<endl;
	texfile<<"\\end{center}"<<endl<<"\\end{table}"<<endl;
	
	texfile<<"\\begin{table}[!ht]"<<endl;
	texfile<<"\\begin{center}"<<endl;
	texfile<<"\\begin{tabular}{lccccccc}"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<" $T_k [MeV]$ & & & n = 0 & n = 1 & n = 2 & n = 3 & n = 4 \\\\ \\hline"<<endl;
	texfile<<"\\multirow{2}{*}{85} & Oset & & "<<pa1[0]<<" & "<<pa1[1]<<" & "<<pa1[2]<<" & "<<pa1[3]<<" & "<<pa1[4]<<" \\\\"<<endl;
	texfile<<" & NuWro & & "<<n1[4]<<" & "<<n1[5]<<" & "<<n1[6]<<" & "<<n1[7]<<" & "<<n1[8]<<" \\\\"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<"\\multirow{2}{*}{245} & Oset & & "<<pa2[0]<<" & "<<pa2[1]<<" & "<<pa2[2]<<" & "<<pa2[3]<<" & "<<pa2[4]<<" \\\\"<<endl;
	texfile<<" & NuWro & & "<<n2[4]<<" & "<<n2[5]<<" & "<<n2[6]<<" & "<<n2[7]<<" & "<<n2[8]<<" \\\\"<<endl;
	texfile<<"\\hline"<<endl;
	texfile<<"\\end{tabular}"<<endl;
	texfile<<"\\caption{Probability that pion absorption occurs after $n$ quasielastic collision}"<<endl;
	texfile<<"\\end{center}"<<endl<<"\\end{table}"<<endl<<"\\end{document}"<<endl;
	texfile.close();
	
	run(string("latex osetcomp.tex"));
	run(string("dvips osetcomp.dvi"));
	run(string("ps2pdf osetcomp.ps"));
	run(string("mv osetcomp.pdf results"));
	run(string("rm osetcomp.*"));
	
	ofstream plik("results/osetcompb.txt");
	
	for (int i = 0; i < 28; i++)
	{
		plik << x[i] << " " << qel[i] << " " << abs[i] <<endl;
	}
	
	plik.close();
	
	ofstream plik2("results/enl.txt");
	ofstream plik3("results/enlqe.txt");
	
	for (int i = 0; i < 25; i++)
	{
		plik2 << energy[i] << " " << en1[i] << " " << en2[i] << endl;
		plik3 << energy[i] << " " << en3[i] << " " << en4[i] << endl;
	}
	
	plik2.close();
	plik3.close();

	return 1;	
}
*/
void kt()
{
	TFile *tf1 = new TFile("root_files/kaskada_new.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
			
	tt1->SetBranchAddress("e",&e1);
	
	double E[60];
	
	for (int i = 0; i < 60; i++) E[i] = (i + 1) * 50.0;
	
	double res[7][60];
	
	for (int i = 0; i < 7; i++) zero(res[i], 60);
	
	double norm[60]; zero(norm, 60);
	
	double help = 0;
	
	double cos[20];
	for (int i = 0; i < 20; i++) cos[i] = (i + 1)*0.1 - 1.0;
	
	double cres[20]; zero(cres, 20);
	
	for (int i = 0; i < 1000000; i++)
	{
		tt1->GetEntry(i);
		
					int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
			int nofp = e1->fof(211) + e1->fof(-211) + e1->fof(111);
			double energy = e1->out[0].E();
		
		if (e1->number_of_interactions() != 0)
		{
			
			switch (pion)
			{
				case   0: put(energy, E, res[0], help, 60); break;
				case   1: put(energy, E, res[3], help, 60); break;
				case  10: put(energy, E, res[2], help, 60); break;
				case 100: put(energy, E, res[1], help, 60); break;
				default: break;
			}
			
			switch (nofp)
			{
				case 0: break;
				case 1: break;
				case 2: put(energy, E, res[4], help, 60); break;
				case 3: put(energy, E, res[5], help, 60); break;
				default: put(energy, E, res[6], help, 60); break;
			}
			
			if (nofp == 1)
			{
				double theta = 0;
				
				for (int a = 0; a < e1->f(); a++)
				{
					if (e1->post[a].pdg == 211 or e1->post[a].pdg == 111 or e1->post[a].pdg == -211)
					{
						theta = e1->post[a].p().z/e1->post[a].momentum();
						break;
					}
				}
				
				put(theta, cos, cres, help, 20);
			}
		}
		
		put(energy, E, norm, help, 60);
		
		cout<<100*i/1000000<<"%\r"<<flush;
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream plik("results/kaskada_new.txt");
	
	for (int i = 0; i < 60; i++)
	{
		if (norm[i] != 0) for (int k = 0; k < 7; k++) res[k][i] /= norm[i];
		
		plik << (i + 0.5) * 50.0;
		
		for (int k = 0; k < 7; k++) plik << " " << res[k][i];
		
		plik << endl;
	}
	
	plik.close();
	
	ofstream plik2("results/angle_new.txt");
	
	for (int i = 0; i < 20; i++)
	{
		plik2 << (i + 1.0) * 0.1 - 1.0 << " " << cres[i] << endl;
	}
	
	plik2.close();
}

void t2k_anal_calc()
{
	string names[4] = {"carbon_fg.root", "carbon_sf.root", "oxygen_fg.root", "oxygen_sf.root"};
	
	const int events = 1000000;
	const int bins   = 30;
	const int abins  = 20;
	
	double energy[bins];
	double momentum[bins];
	
	for (int i = 0; i < bins; i++)
	{
		energy[i]   = (i + 1.0) * 50.0;
		momentum[i] = (i + 1.0) * 50.0;
	}
	
	double angle[abins];
	
	for (int i = 0; i < abins; i++) angle[i] = (i + 1.0) * 0.1 - 1.0; 
	
	//0-6 ccqe, ccpi+, ccpi-, ccpi0, ncpi+, ncpi-, ncpi0 // 0-1 before fsi, after fsi // 0-3 see names 
	
	double fen[7][2][4][bins] = {{0}, {0}, {0}}; 

	double fmom[7][2][4][bins] = {{0}, {0}, {0}}; 

	double fang[7][2][4][abins] = {{0}, {0}, {0}}; 
	
	double normen[4][bins]   = {{0}};
	double normmom[4][bins]  = {{0}};
	double normang[4][abins] = {{0}};
	
	double help = 0;

	for (int k = 0; k < 4; k++)
	{
		string rootfile = "root_files/T2K_" + names[k];
		string rtxtfile = rootfile + ".txt";
		
		TFile *tf1 = new TFile(rootfile.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
		
		tt1->SetBranchAddress("e",&e1);
	
		for (int i = 0; i < events; i++)
		{
			tt1->GetEntry(i);
			
			double E = e1->in[0].E();
			double P = e1->out[0].momentum();
			double A = e1->out[0].p().z/P;
			
			put(E, energy, normen[k], help, bins);
			put(P, momentum, normmom[k], help, bins);
			put(A, angle, normang[k], help, abins);
			
			int pion0 = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
			int pion1 = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
			int nucl0 = 10*e1->nof(pdg_proton) + e1->nof(pdg_neutron);
			int nucl1 = 10*e1->fof(pdg_proton) + e1->fof(pdg_neutron);
			
			if (e1->flag.cc)
			{
				if (e1->flag.qel)
				{
					put(E, energy, fen[0][0][k], help, bins);
					put(P, momentum, fmom[0][0][k], help, bins);
					put(A, angle, fang[0][0][k], help, abins);
				}
				if (pion1 == 0 and nucl1 == 10)
				{
					put(E, energy, fen[0][1][k], help, bins);
					put(P, momentum, fmom[0][1][k], help, bins);
					put(A, angle, fang[0][1][k], help, abins);					
				}
				if (pion0 == 100)
				{
					put(E, energy, fen[1][0][k], help, bins);
					put(P, momentum, fmom[1][0][k], help, bins);
					put(A, angle, fang[1][0][k], help, abins);
				}
				if (pion1 == 100)
				{
					put(E, energy, fen[1][1][k], help, bins);
					put(P, momentum, fmom[1][1][k], help, bins);
					put(A, angle, fang[1][1][k], help, abins);					
				}
				if (pion0 == 10)
				{
					put(E, energy, fen[2][0][k], help, bins);
					put(P, momentum, fmom[2][0][k], help, bins);
					put(A, angle, fang[2][0][k], help, abins);
				}
				if (pion1 == 10)
				{
					put(E, energy, fen[2][1][k], help, bins);
					put(P, momentum, fmom[2][1][k], help, bins);
					put(A, angle, fang[2][1][k], help, abins);
				}
				if (pion0 == 1)
				{
					put(E, energy, fen[3][0][k], help, bins);
					put(P, momentum, fmom[3][0][k], help, bins);
					put(A, angle, fang[3][0][k], help, abins);
				}
				if (pion1 == 1)
				{
					put(E, energy, fen[3][1][k], help, bins);
					put(P, momentum, fmom[3][1][k], help, bins);
					put(A, angle, fang[3][1][k], help, abins);
				}
			}
			else
			{
				if (pion0 == 100)
				{
					put(E, energy, fen[4][0][k], help, bins);
					put(P, momentum, fmom[4][0][k], help, bins);
					put(A, angle, fang[4][0][k], help, abins);					
				}
				if (pion1 == 100)
				{
					put(E, energy, fen[4][1][k], help, bins);
					put(P, momentum, fmom[4][1][k], help, bins);
					put(A, angle, fang[4][1][k], help, abins);					
				}
				if (pion0 == 10)
				{
					put(E, energy, fen[5][0][k], help, bins);
					put(P, momentum, fmom[5][0][k], help, bins);
					put(A, angle, fang[5][0][k], help, abins);					
				}
				if (pion1 == 10)
				{
					put(E, energy, fen[5][1][k], help, bins);
					put(P, momentum, fmom[5][1][k], help, bins);
					put(A, angle, fang[5][1][k], help, abins);					
				}
				if (pion0 == 1)
				{
					put(E, energy, fen[6][0][k], help, bins);
					put(P, momentum, fmom[6][0][k], help, bins);
					put(A, angle, fang[6][0][k], help, abins);					
				}
				if (pion1 == 1)
				{
					put(E, energy, fen[6][1][k], help, bins);
					put(P, momentum, fmom[6][1][k], help, bins);
					put(A, angle, fang[6][1][k], help, abins);					
				}
			}

			cout<<rootfile<<": "<<100*i/events<<"%\r"<<flush;
		}
	
		cout<<rootfile<<": done"<<endl<<endl;
		
		double cross = crosssection(rtxtfile);
		
		for (int i = 0; i < bins; i++)
		{
			if (normen[k][i] != 0)
			{
				for (int l = 0; l < 2; l++)
				{
					for (int m = 0; m < 7; m++)
					{
						fen[m][l][k][i] *= cross/normen[k][i];
					}
				}
			}
			if (normmom[k][i] != 0)
			{
				for (int l = 0; l < 2; l++)
				{
					for (int m = 0; m < 7; m++)
					{
						fmom[m][l][k][i] *= cross/normmom[k][i];
					}
				}
			}
		}			

		for (int i = 0; i < abins; i++)
		{
			if (normang[k][i] != 0)
			{
				for (int l = 0; l < 2; l++)
				{
					for (int m = 0; m < 7; m++)
					{
						fang[m][l][k][i] *= cross/normang[k][i];
					}
				}
			}
		}
		
		delete e1;
		delete tt1;
		delete tf1;			
	}
	
	for (int i = 0; i < bins; i++)
	{
		energy[i]   = (i + 0.5) * 50.0;
		momentum[i] = (i + 0.5) * 50.0;
	}
		
	for (int i = 0; i < abins; i++) angle[i] = (i + 0.5) * 0.1 - 1.0; 
	
	run("mkdir -p results/t2k_anal/");
	
	string f1[7] = {"results/t2k_anal/CCQE_", "results/t2k_anal/CCPiP_", "results/t2k_anal/CCPiM_", "results/t2k_anal/CCPi0_", "results/t2k_anal/NCPiP_", "results/t2k_anal/NCPiM_", "results/t2k_anal/NCPi0_"};
	string f2[2] = {"before_fsi_", "after_fsi_"};
	string f3[3] = {"neutrino_energy_", "lepton_momentum_", "lepton_angle_"};
	string f4[2] = {"carbon.txt", "oxygen.txt"};
	
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 2; l++)
				{
					string fname = f1[i] + f2[j] + f3[k] + f4[l];
					ofstream plik(fname.c_str());
					
					switch (k)
					{
						case 0: plik << "#neutrino energy | "; break;
						case 1: plik << "#lepton momentum | "; break;
						case 2: plik << "#lepton angle | "; break;
						default: break;
					}
					
					plik << " fermi gas | spectral function | fg/sf " << endl;
					
					int ile = bins;
					if (k == 2) ile = abins;
					
					int tar = 0;
					if (l == 1) tar = 2;
					
					for (int b = 0; b < ile; b++)
					{
						switch (k)
						{
							case 0: plik << energy[b] << " " << fen[i][j][tar][b] << " " << fen[i][j][tar+1][b] << " " << fen[i][j][tar][b]/fen[i][j][tar+1][b] << endl; break;
							case 1: plik << momentum[b] << " " << fmom[i][j][tar][b] << " " << fmom[i][j][tar+1][b] << " " << fmom[i][j][tar][b]/fmom[i][j][tar+1][b] << endl; break;
							case 2: plik << angle[b] << " " << fang[i][j][tar][b] << " " << fang[i][j][tar+1][b] << " " << fang[i][j][tar][b]/fang[i][j][tar+1][b] << endl; break;
							default: break;
						}
					}
					
					plik.close();
				}
			}
		}
	}
}

void hayato_calc()
{
	const int events = 25000;
	
	ofstream file("results/oxygen.txt");
	file << "#neutrino energy [GeV] | xsec on proton (fg) | xsec on neutron (fg) | xsec on proton (sf) | xsec on neutron (sf)" << endl;
	
	for (int i = 0; i < 111; i++)
	{
		int k = i;
		double en = k*25.0 + 250.0;
		
		if (i >= 30)
		{
			k = i - 30;
			en = k*50.0 + 1000.0;
		}
		if (i >= 60)
		{
			k = i - 60;
			en = k*250.0 + 2500;
		}
		if (i >= 90)
		{
			k = i - 90;
			en = k*1000.0 + 10000;
		}
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string fg = "root_files/hayato/E" + energy + "_oxygen_fg.root";
		string sf = "root_files/hayato/E" + energy + "_oxygen_sf.root";
		
		string fgtxt = "root_files/hayato/E" + energy + "_oxygen_fg.root.txt";
		string sftxt = "root_files/hayato/E" + energy + "_oxygen_sf.root.txt";
		
		int proton_fg = 0;
		int neutron_fg = 0;

		int proton_sf = 0;
		int neutron_sf = 0;
		
		TFile *tf1 = new TFile(fg.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
		
		tt1->SetBranchAddress("e",&e1);
	
		for (int a = 0; a < events; a++)
		{
			tt1->GetEntry(a);
			
			if (e1->in[1].pdg == pdg_proton) proton_fg++;
			else neutron_fg++;
			
			cout<<fg<<": "<<100*a/events<<"%\r"<<flush;
		}
	
		cout<<fg<<": done"<<endl<<endl;
		
		delete e1;
		delete tt1;
		delete tf1;
		
		TFile *tf2 = new TFile(sf.c_str());
		TTree *tt2 = (TTree*)tf2->Get("treeout");
		event *e2   = new event();
		
		tt2->SetBranchAddress("e",&e2);
	
		for (int a = 0; a < events; a++)
		{
			tt2->GetEntry(a);
			
			if (e2->in[1].pdg == pdg_proton) proton_sf++;
			else neutron_sf++;
			
			cout<<sf<<": "<<100*a/events<<"%\r"<<flush;
		}
	
		cout<<sf<<": done"<<endl<<endl;
		
		delete e2;
		delete tt2;
		delete tf2;
		
		double cross_fg = 2.0*crosssection(fgtxt);
		double cross_sf = 2.0*crosssection(sftxt);
		
		file << energy << " " << proton_fg * cross_fg / events << " " << neutron_fg * cross_fg / events << " " << proton_sf * cross_sf / events << " " << neutron_sf * cross_sf / events << endl;
	}
	
	file.close();
	
}

/*void dens_test_calc()
{
	TFile *tf1 = new TFile("root_files//E40k_carbon.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
	
	tt1->SetBranchAddress("e",&e1);
	
	int nofp[10] = {0};
	int nofn[10] = {0};
	
	double dens[10][50] = {{0}};
	
	double r[10] = {0};
	
	int norm[10] = {0};
	
	for (int i = 0; i < 25000; i++)
	{		
		tt1->GetEntry(i);
		
		int n = e1->number_of_interactions();
						
		if (n < 10 and n >= 0)
		{
			norm[n]++;
			nofp[n] += e1->protons_hist;
			nofn[n] += e1->neutrons_hist;
			
			r[n] += e1->radius_hist;
			
			for (int l = 0; l < 50; l++) dens[n][l] += e1->density_hist[l];
		}
		
		cout<<0.001*i<<"%\r"<<flush;
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	ofstream p1("results/dens/test.txt");
	for (int i = 0; i < 10; i++)
	{
		p1 << i << " " << r[i]/norm[i] << " " << (double)nofp[i]/norm[i] << " " << (double)nofn[i]/norm[i] << endl;
				
		stringstream temp;
		string nofi;

		temp << i;
		temp >> nofi;
		
		string help = "results/dens/d" + nofi + ".txt";
		ofstream p2(help.c_str());
		for (int k = 0; k < 50; k++) p2 << (k + 1.0) * 0.2 << " " << dens[i][k]/norm[i] << endl;
		p2.close();
	}
	p1.close();	
}
*/
/*void ptcalc()
{
	const int events = 100000;
	const int bins = 20;
	
	double Q2[bins];
	for (int i = 0; i < bins; i++)
		Q2[i] = (i + 1.0) * 0.25;
	
	double mom[bins];
	for (int i = 0; i < bins; i++)
		mom[i] = (i + 1.0) * 100.0;
	
	double mdis[10][bins] = {{0}};
	double mres[10][bins] = {{0}};

	double qdis[10][bins] = {{0}};
	double qres[10][bins] = {{0}};

	double mdisnorm[10][bins] = {{0}};
	double mresnorm[10][bins] = {{0}};

	double qdisnorm[10][bins] = {{0}};
	double qresnorm[10][bins] = {{0}};
	
	double edis[10] = {0};
	double eres[10] = {0};
	
	double eresnorm[10] = {0};
	double edisnorm[10] = {0};	
	
	double rest = 0;
	
	string out[10];
	
	for (int i = 0; i < 10; i++)
	{
		double en = (i+1.0)*500.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string in  = "root_files/ptE" + energy + rot;
		out[i] = "results/E" + energy + ".txt";
		
		TFile *tf1 = new TFile(in.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
	
		tt1->SetBranchAddress("e",&e1);
		
		for (int z = 0; z < events; z++)
		{		
			tt1->GetEntry(z);

			double przekaz = -e1->q2()/1000000.0;
			
			for (int k = 0; k < e1->f(); k++)
			{
				if ((e1->post[k].pdg == 211 or e1->post[k].pdg == -211 or e1->post[k].pdg == 111) and e1->post[k].mother_pdg == 0)
				{
					double ped = e1->post[k].momentum();
					double fz = e1->post[k].fozo();
					
					if (e1->flag.res)
					{
						put(ped, mom, mresnorm[i], rest, bins);
						put(ped, mom, mres[i], rest, bins, fz);
						
						put(przekaz, Q2, qresnorm[i], rest, bins);
						put(przekaz, Q2, qres[i], rest, bins, fz);
						
						eres[i] += fz;
						eresnorm[i]++;
					}
					else if (e1->flag.dis)
					{
						put(ped, mom, mdisnorm[i], rest, bins);
						put(ped, mom, mdis[i], rest, bins, fz);
						
						put(przekaz, Q2, qdisnorm[i], rest, bins);
						put(przekaz, Q2, qdis[i], rest, bins, fz);
						
						edis[i] += fz;
						edisnorm[i]++;
					}					
				}
			}
		
			cout << in << " " << 100.0*z/events << "%\r" << flush;			
		}
		
		cout << endl;
		
		delete e1;
		delete tt1;
		delete tf1;
	}
	
	ofstream wynik("results/fe.txt");
	
	for (int i = 0; i < 10; i ++)
	{
		for (int k = 0; k < bins; k++)
		{
			mres[i][k] /= mresnorm[i][k];
			qres[i][k] /= qresnorm[i][k];
			mdis[i][k] /= mdisnorm[i][k];
			qdis[i][k] /= qdisnorm[i][k];
		}
		
		wynik << (i + 1.0) * 500.0 << " " << eres[i]/eresnorm[i] << " " << edis[i]/edisnorm[i] << endl;
	}
	
	wynik.close();
	
	for (int i = 0; i < 10; i++)
	{
		ofstream plik(out[i].c_str());
		plik << "#momentum | res | dis | Q2 | res | dis" << endl;
		
		for (int k = 0; k < bins; k++)
			plik << (k + 0.5) * 100.0 << " " << mres[i][k] << " " << mdis[i][k] << " " << (k + 0.5) * 0.25 << " " << qres[i][k] << " " << qdis[i][k] << endl;
			
		plik.close();
	}
}
*/
/*void ang_calc()
{
	const int events = 100000;
	const int bins = 18;
	
	double angle[bins];
	for (int i = 0; i < bins; i++)
		angle[i] = (i + 1.0) * 10;
	
	double ce[10][bins] = {{0}};
	
	double rest = 0;
	
	string out[10];
	
	for (int i = 0; i < 10; i++)
	{
		double en = (i+1.0)*100.0+100.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string in  = "root_files/iiE" + energy + rot;
		out[i] = "results/angle/iiE" + energy + ".txt";
		
		TFile *tf1 = new TFile(in.c_str());
		TTree *tt1 = (TTree*)tf1->Get("treeout");
		event *e1   = new event();
	
		tt1->SetBranchAddress("e",&e1);
		
		for (int z = 0; z < events; z++)
		{		
			tt1->GetEntry(z);
			
			for (int k = 0; k < e1->f(); k++)
			{				
				if (e1->post[k].pdg == 211 and e1->post[k].mother_proc == 20)
				{
					double kat = e1->post[k].p()*e1->post[k].mother_momentum/e1->post[k].momentum()/e1->post[k].mother_momentum.length();
					kat = acos(kat)*180.0/M_PI; 
					put(kat, angle, ce[i], rest, bins);
				}
			}
		
			cout << in << " " << 100.0*z/events << "%\r" << flush;			
		}
		
		cout << endl;
		
		delete e1;
		delete tt1;
		delete tf1;
	}
		
	for (int i = 0; i < 10; i++)
	{
		ofstream plik(out[i].c_str());
		
		for (int k = 0; k < bins; k++)
			plik << (k + 0.5) * 10 << " " << ce[i][k] << endl;
			
		plik.close();
	}
}
*/
void make_dist(string filename, double res[100][20])
{	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < 500000; i++)
	{
		tt1->GetEntry(i);
		
		int a = e1->out[0].momentum()/50;
		if (a == 100) a--;
		
		int b = e1->out[0].p().z/e1->out[0].momentum()/0.1 + 10;
		if (b == 20) b--;
		
		if (a < 100 and b < 20)
			res[a][b]++;
					
		cout << filename << ": "<< i/5000 << "\r" << flush; 
	}
	
	cout << endl;
	
	delete e1;
	delete tt1;
	delete tf1;
}

void kendall_calc(string pdg, string p)
{
	string fname = "ccqe/xsec/xsec_" + pdg + "_" + p + ".txt";
	
	ofstream f(fname.c_str());
	
	for (int i = 3; i < 23; i++)
	{
		int en = (i+1)*50;
		if (i > 18) en = (i-18)*1000;
		
		int en2 = (i+2)*50;
		if (i > 18) en2 = (i-17)*1000;
		
		f << "#energy | FG | SF | FG/SF " << endl;
		
		f << (en2 - en)/2.0 + en << " ";

		stringstream temp;
		stringstream temp2;
		string energy;
		string energy2;

		temp << en;
		temp >> energy;
		temp2 << en2;
		temp2 >> energy2;
		
		string file1 = "ccqe/E" + energy + "_" + energy2 + "_" + p + "_" + p + "_" + pdg + "_FG.root";
		string file2 = "ccqe/E" + energy + "_" + energy2 + "_" + p + "_" + p + "_" + pdg + "_SF.root";
		string file1txt = file1 + ".txt";
		string file2txt = file2 + ".txt";
		
		double cross1 = crosssection(file1txt);
		double cross2 = crosssection(file2txt);
		
		f << cross1 << " " << cross2 << " " << cross1/cross2 << " " << (cross1 - cross2)/cross1 << endl;
		
		cross1 /= 500000.0;
		cross2 /= 500000.0;
		
		double res1[100][20] = {{0}};
		double res2[100][20] = {{0}};
		
		make_dist(file1, res1);
		make_dist(file2, res2);
		
		string txtname = "ccqe/results/ccqe_E" + energy + "_" + energy2 + "_PDG_" + pdg + "_TARGET_" + p + ".txt";
		
		ofstream plik(txtname.c_str());
		
		for (int a = 0; a < 100; a++)
		{
			for (int b = 0; b < 20; b ++)
			{
				plik << (a+1)*50 << " " << b*0.1 - 0.9 << " " << res1[a][b] << " " << res2[a][b] << " " << res1[a][b]*cross1 << " " << res2[a][b]*cross2 << endl;
			}
		}
		
		plik.close();
	}
	
	f.close();
}

void xsec_calc()
{
	ofstream file("towork/xsec.txt");
	file << "#energy | qel cc | qel nc | res cc | res nc | dis cc | dis nc | cctotal | nctotal" << endl;
	
	for (int i = 0; i < 23; i++)
	{
		int en = (i+1)*100;
		if (i > 9) en = 1000 + (i-9)*500;
		if (i > 17) en = 5000 + (i-17)*1000;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
//		string fproton = "towork/E" + energy + "_proton.root.txt";
//		string fneutron = "towork/E" + energy + "_neutron.root.txt";
		string tmec = "towork/E" + energy + "_mec.root.txt";
		
//		double proton[8];
//		double neutron[8];
//		double res[8];
		double npnh;
		
		double energia = en/1000.0;
		file << energia;
		
		//~ for (int i = 0; i < 10; i++)
		//~ {
			//~ int a = i;
			//~ 
			//~ if (i == 6 or i == 7)
				//~ continue;
			//~ else if (i > 7)
				//~ a = i - 2;
			//~ 
			//~ proton[i] = 1e38*crosssection(fproton, a);
			//~ neutron[i] = 1e38*crosssection(fneutron, a);
			//~ res[i] = (proton[i] + neutron[i])/2.0;
			//~ file << " " << res[i]/energia;
		//~ }
		
		npnh = 1e38*crosssection(tmec,8);
		
		//file << " " << (res[0] + res[2] + res[4] + res[6])/energia << " " << (res[1] + res[3] + res[5] + res[7])/energia << endl;
		file << " " << npnh/energia << endl;
	}
	
	file.close();
}

void xsec_calc2()
{
	ofstream file("towork/xsec_anti.txt");
	file << "#energy | qel cc | qel nc | res cc | res nc | dis cc | dis nc | cctotal | nctotal" << endl;
	
	for (int i = 0; i < 23; i++)
	{
		int en = (i+1)*100;
		if (i > 9) en = 1000 + (i-9)*500;
		if (i > 17) en = 5000 + (i-17)*1000;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string fproton = "towork/E" + energy + "_proton_anti.root.txt";
		string fneutron = "towork/E" + energy + "_neutron_anti.root.txt";
		
		double proton[6];
		double neutron[6];
		double res[6];
		
		double energia = en/1000.0;
		file << energia;
		
		for (int i = 0; i < 6; i++)
		{
			proton[i] = 1e38*crosssection(fproton, i);
			neutron[i] = 1e38*crosssection(fneutron, i);
			res[i] = (proton[i] + neutron[i])/2.0;
			file << " " << res[i]/energia;
		}
		
		file << " " << (res[0] + res[2] + res[4])/energia << " " << (res[1] + res[3] + res[5])/energia << endl;
	}
	
	file.close();
}

void xsec_calc4()
{
	ofstream file("towork/xsec_1350.txt");
	file << "#energy | qel cc | qel nc | res cc | res nc | dis cc | dis nc | cctotal | nctotal" << endl;
	
	for (int i = 0; i < 23; i++)
	{
		int en = (i+1)*100;
		if (i > 9) en = 1000 + (i-9)*500;
		if (i > 17) en = 5000 + (i-17)*1000;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string fproton = "towork/E" + energy + "_proton_1350.root.txt";
		string fneutron = "towork/E" + energy + "_neutron_1350.root.txt";
		
		double proton[6];
		double neutron[6];
		double res[6];
		
		double energia = en/1000.0;
		file << energia;
		
		for (int i = 0; i < 6; i++)
		{
			proton[i] = 1e38*crosssection(fproton, i);
			neutron[i] = 1e38*crosssection(fneutron, i);
			res[i] = (proton[i] + neutron[i])/2.0;
			file << " " << res[i]/energia;
		}
		
		file << " " << (res[0] + res[2] + res[4])/energia << " " << (res[1] + res[3] + res[5])/energia << endl;
	}
	
	file.close();
}

void xsec_calc3()
{
	ofstream file("towork/xsec_coh.txt");
	file << "#energy | NC coh" << endl;
		
	for (int i = 0; i < 23; i++)
	{
		int en = (i+1)*100;
		if (i > 9) en = 1000 + (i-9)*500;
		if (i > 17) en = 5000 + (i-17)*1000;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string cohin = "towork/E" + energy + "_carbon_coh.root.txt";
				
		double energia = en/1000.0;
		file << energia << " " << 1e40*crosssection(cohin) << endl;		
	}
	
	file.close();
}

void ccqe_calc(string filename, double tab[][20], double *q2)
{
	string txt = filename + ".txt";
	double xsec = crosssection(txt);
	
	int events = 5000000;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		double tk = e1->out[0].Ek();
		double cos = e1->out[0].p().z/e1->out[0].momentum();
		
		int a = tk/50;
		int b = cos/0.1 + 10;
		if (b == 20) b--;
		
		tab[a][b]++;
		
		double q = -e1->q2()/1000000.0;
		int c = q/0.05;
		if (c == 20) c--;
		
		q2[c]++;
		
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<filename<<": done"<<endl<<endl;
	
	
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int a = 0; a < 20; a ++)
		for (int b = 0; b < 20; b++)
			tab[a][b] *= xsec/events/100/0.1;
			
	for (int c = 0; c < 20; c++)
		q2[c] *= xsec/events/0.05;
}

void res_calc(string filename, double *mom, double tab[][20])
{
	string txt = filename + ".txt";
	double xsec = crosssection(txt);
	
	int events = 5000000;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		double tk = e1->out[0].Ek();
		double cos = e1->out[0].p().z/e1->out[0].momentum();
		
		int a = tk/50;
		int b = cos/0.1 + 10;
		if (b == 20) b--;
		
		tab[a][b]++;
		
		int piony = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		
		if (piony == 100)
		{
			for (int x = 0; x < e1->n(); x++)
			{
				if (e1->out[x].pdg == 211)
				{
					int y = e1->out[x].momentum()/1000/0.05;
					mom[y]++;
				}
			}
		}
		
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<filename<<": done"<<endl<<endl;

	
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int a = 0; a < 20; a ++)
		for (int b = 0; b < 20; b++)
			tab[a][b] *= xsec/events/100/0.1;

	for (int y = 0; y < 20; y++)
		mom[y] *= xsec/events/0.05;
}

void calc_new_distr(string filename, double tab[][20])
{
	string txt = filename + ".txt";
	double xsec = crosssection(txt);
	
	int events = 5000000;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		double tk = e1->out[0].momentum();//Ek();
		double cos = e1->out[0].p().z/e1->out[0].momentum();
		
		int a = tk/50;
		int b = cos/0.1 + 10;
		if (b == 20) b--;
		
		tab[a][b]++;
				
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<filename<<": done"<<endl<<endl;

	
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int a = 0; a < 20; a ++)
		for (int b = 0; b < 20; b++)
			tab[a][b] *= xsec/events/100/0.1;
}

void prepare_distr(double tab[][20], string out)
{
	ofstream plik(out.c_str());
	
	double max = 0;
	
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 20; j++)
			if (max < tab[i][j])
				max = tab[i][j];
				
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			double scale = sqrt(tab[i][j]/max);
			double shiftx = 0.05 - scale*0.05;
			double shifty = 0.1 - scale*0.1;			
			
			double x1 = i*0.05 + shiftx/2;
			double x2 = (i+1)*0.05 - shiftx/2;
			
			double y1 = j*0.1 - 1.0 + shifty/2;
			double y2 = (j+1)*0.1 - 1.0 - shifty/2;
			
			plik << x1 << " " << y1 << endl;
			plik << x2 << " " << y1 << endl;
			plik << x2 << " " << y2 << endl;
			plik << x1 << " " << y2 << endl;
			plik << x1 << " " << y1 << endl;
			
			plik << endl;
		}
	}

	plik.close();
}

void new_distr()
{
	string fg = "towork/E1_ccqe_carbon_fg.root";
	string sf = "towork/E1_ccqe_carbon_sf.root";
	string res = "towork/E1_carbon.root";
	
	double qfg[20][20] = {{0}};
	double qsf[20][20] = {{0}};
	double lep[20][20] = {{0}};
	
	calc_new_distr(fg, qfg);
	calc_new_distr(sf, qsf);
	calc_new_distr(res, lep);
	
	prepare_distr(qfg, "towork/fg_box_mom.txt");
	prepare_distr(qsf, "towork/sf_box_mom.txt");
	prepare_distr(lep, "towork/res_box_mom.txt");
}	

void towork_calc()
{
	string fg = "towork/E1_ccqe_carbon_fg.root";
	string sf = "towork/E1_ccqe_carbon_sf.root";
	
	double fgt[20][20] = {{0}};
	double sft[20][20] = {{0}};
	
	double qfg[20] = {0};
	double qsf[20] = {0};
	
	ccqe_calc(fg, fgt, qfg);
	ccqe_calc(sf, sft, qsf);
	
	ofstream fgr("towork/fg.txt");
	ofstream sfr("towork/sf.txt");
	
	for (int a = 0; a < 20; a++)
	{
		for (int b = 0; b < 20; b++)
		{
			fgr << a*0.05 << " " << b*0.1 - 1.0 << " " << fgt[a][b] << endl;
			fgr << a*0.05 << " " << b*0.1 - 0.9 << " " << fgt[a][b] << endl;
			sfr << a*0.05 << " " << b*0.1 - 1.0 << " " << sft[a][b] << endl;
			sfr << a*0.05 << " " << b*0.1 - 0.9 << " " << sft[a][b] << endl;

		}

		fgr << endl;
		sfr << endl;
					
		for (int b = 0; b < 20; b++)
		{
			fgr << a*0.05 + 0.05 << " " << b*0.1 - 1.0 << " " << fgt[a][b] << endl;
			fgr << a*0.05 + 0.05 << " " << b*0.1 - 0.9 << " " << fgt[a][b] << endl;
			sfr << a*0.05 + 0.05 << " " << b*0.1 - 1.0 << " " << sft[a][b] << endl;
			sfr << a*0.05 + 0.05 << " " << b*0.1 - 0.9 << " " << sft[a][b] << endl;
		}
		
		fgr << endl;
		sfr << endl;
	}
	
	fgr.close();
	sfr.close();
	
	ofstream qfgr("towork/qfg.txt");
	ofstream qsfr("towork/qsf.txt");
	
	for (int c = 0; c < 20; c++)
	{
		qfgr << 0.025 + c*0.05 << " " << qfg[c] << endl;
		qsfr << 0.025 + c*0.05 << " " << qsf[c] << endl;
	}
	
	qfgr.close();
	qsfr.close();
	
	string res = "towork/E1_carbon.root";
	
	double pimom[20] = {0};
	double lep[20][20] = {{0}};
	
	res_calc(res, pimom, lep);
	
	ofstream resl("towork/reslepton.txt");
	
	for (int a = 0; a < 20; a++)
	{
		for (int b = 0; b < 20; b++)
		{
			resl << a*0.05 << " " << b*0.1 - 1.0 << " " << lep[a][b] << endl;
			resl << a*0.05 << " " << b*0.1 - 0.9 << " " << lep[a][b] << endl;
		}
		
		resl << endl;
		
		for (int b = 0; b < 20; b++)
		{
			resl << a*0.05 + 0.05 << " " << b*0.1 - 1.0 << " " << lep[a][b] << endl;
			resl << a*0.05 + 0.05 << " " << b*0.1 - 0.9 << " " << lep[a][b] << endl;
		}
		
		resl << endl;
	}
	
	resl.close();
	
	ofstream pi("towork/pimom.txt");
	
	for (int c = 0; c < 20; c++)
		pi << 0.025 + c*0.05 << " " << pimom[c] << endl;
	
	pi.close();
}

void ccpip_js_calcH(string filename, double *res)
{
	int events = 2000000;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		
		if (pion == 100)
		{
			double Tk = 0;
			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg == 211)
				{
					Tk = e1->out[k].Ek();
					break;
				}
			}
			
			int a = 0;
			
			if (Tk > 50) a = (Tk - 50)/25 + 1;
			
			if (a < 15) res[a]++;
			
		}
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<filename<<": done"<<endl<<endl;

	delete e1;
	delete tt1;
	delete tf1;

}

void ccpip_js_calcC(string filename, double *res)
{
	int events = 2000000;
	
	TFile *tf1 = new TFile(filename.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		if (pion == 100)
		{
			double Tk = 0;
			
			for (int k = 0; k < e1->f(); k++)
			{
				if (e1->post[k].pdg == 211)
				{
					Tk = e1->post[k].Ek();
					break;
				}
			}
			
			int a = 0;
			
			if (Tk > 50) a = (Tk - 50)/25 + 1;
			
			if (a < 15) res[a]++;
			
		}
		cout<<filename<<": "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<filename<<": done"<<endl<<endl;

	delete e1;
	delete tt1;
	delete tf1;

}
	
void ccpip_js_calc()
{
	double H[15] = {0};
	double C[15] = {0};
	
	double xsecH = crosssection("ccpip/H.root.txt");
	double xsecC = crosssection("ccpip/C.root.txt");
	
	ccpip_js_calcH("ccpip/H.root", H);
	ccpip_js_calcC("ccpip/C.root", C);
	
	double res[15] = {0};
	
	for (int i = 0; i < 15; i++)
		res[i] = (2.08*xsecH*H[i] + 12.0*xsecC*C[i])*1e41/2000000;
		
	ofstream file("ccpip/ccpip.txt");
	
	file << "0 " << res[0]/50 << endl;
	file << "50 " << res[0]/50 << endl;
	
	for (int i = 0; i < 14; i++)
		file << 50 + i*25 << " " << res[i+1]/25 << endl << 50 + (i+1)*25 << " " << res[i+1]/25 << endl;
	
	file.close();

	ofstream file2("ccpip/ccpip2.txt");
	
	file2 << "25 " << res[0]/50 << endl;
	
	for (int i = 0; i < 14; i++)
		file2 << 62.5 + i*25.0 << " " << res[i+1]/25 << endl;
		
	file2.close();
}	

void hayato_calc0812()
{	
	ofstream file("results/numu_oxygen_total.txt");
	file << "#neutrino energy [MeV] | total xsec per nucleon with dipole FF | total xsec per nucleon with BBBA05, hep-ex/0602017 BBBA05 for Q2<18 GeV" << endl;
	
	for (int en = 50; en <= 1250; en += 50)
	{
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
		
		string dip = "root_files/hayato/E" + energy + "_oxygen_dip.root.txt";
		string bbb = "root_files/hayato/E" + energy + "_oxygen_bbb.root.txt";
		
		double cross_dip = crosssection(dip);
		double cross_bbb = crosssection(bbb);
		
		file << energy << " " << cross_dip << " " << cross_bbb << endl;
	}
	
	file.close();
	
}

void test_calc()
{
	int events = 100000;
	double res[20] = {0};
	
	TFile *tf1 = new TFile("test.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		if (pion == 100)
		{
			double Tk = 0;
			
			for (int k = 0; k < e1->f(); k++)
			{
				if (e1->post[k].pdg == 211)
				{
					Tk = e1->post[k].Ek();
					break;
				}
			}
			
			int a = Tk / 50;
			if (a < 20) res[a]++;
			
		}
		cout<<"test: "<<100*i/events<<"%\r"<<flush;
	}
	
	cout<<"test: done"<<endl<<endl;

	delete e1;
	delete tt1;
	delete tf1;
	
	double xsec = crosssection("test.root.txt");
	
	cout << "Total: " << xsec << endl << endl;
	
	for (int i = 0; i < 20; i++)
		cout << 50 * (i + 0.5) << " " << res[i] * xsec / events / 50.0 << endl;

}
