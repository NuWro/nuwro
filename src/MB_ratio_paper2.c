#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "event1.h"
#include "particle.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include "TBrowser.h"
#include "TGraph2D.h"

ofstream File ("Paper_MB_Ratio_CC1PiPlusOverCCQE.dat");

vector<double> przekrojeH(8);
vector<double> przekrojeC(8);

void crosssections (vector<double> &v, ifstream& Input)
{
string y;
double x;
if (Input)
{
do
{
getline (Input, y);

for (int j=0; j<8; j++)
{
for (int k=0; k<3; k++)
{Input>>y;}
Input>>v[j];
}
} while (Input);
}
}

int events=100000;

double PI = 4*atan(1);

double Ccounter1PiPlusIso=0;
double CcounterQEIso=0;
int Ccounter1PiPlusFSI=0;
int CcounterQEFSI=0;

int Hcounter1PiPlusIso=0;
int HcounterQEIso=0;
int Hcounter1PiPlusFSI=0;
int HcounterQEFSI=0;

int energia[13]={450,550,650,750,850,950,1050,1150,1250,1350,1500,1700,2100};

double xx[13]={450/1e3,550/1e3,650/1e3,750/1e3,850/1e3,950/1e3,1050/1e3,1150/1e3,1250/1e3,1350/1e3,1500/1e3,1700/1e3,2100/1e3};

double xxer[13]={50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 50/1e3, 100/1e3, 100/1e3, 300/1e3};
double yy1[13]={0.036, 0.1, 0.191, 0.278, 0.371, 0.465, 0.551, 0.607, 0.677, 0.7, 0.777, 0.904, 1.022};
double yy1er[13]={0.005, 0.011, 0.019, 0.028, 0.04, 0.053, 0.066, 0.077, 0.091, 0.097, 0.109, 0.137, 0.161};
double yy2[13]={0.045, 0.13, 0.258, 0.381, 0.52, 0.656, 0.784, 0.855, 0.957, 0.985, 1.073, 1.233, 1.318};
double yy2er[13]={0.008, 0.018, 0.033, 0.047, 0.064, 0.082, 0.1, 0.114, 0.132, 0.141, 0.157, 0.207, 0.247};

double Tina [17] = { 517/1e3, 559/1e3, 636/1e3, 669/1e3, 741/1e3, 777/1e3,  849/1e3, 951/1e3, 1054/1e3, 1150/1e3, 1250/1e3, 1351/1e3, 1485/1e3, 1521/1e3, 1702/1e3, 1779/1e3, 1896/1e3 };
double Leitner [17] = {0.086, 0.116, 0.174,  0.2, 0.257, 0.282, 0.338, 0.413, 0.473, 0.512, 0.551, 0.59, 0.646, 0.658, 0.697,  0.717, 0.741};

double wynikFSI[13];
double wynikIso[13];

TGraphErrors *MBRatioFSI = new TGraphErrors (13, xx, yy1, xxer, yy1er);
TGraphErrors *MBRatioIso = new TGraphErrors (13, xx, yy2, xxer, yy2er);
TGraph *GiBUU = new TGraph (17, Tina, Leitner);

char crossC[130];
char crossH[130];
char eventsC[130];
char eventsH[130];

void MB_ratio_paper()
{

File<<"# NuWro: CC1Pi+/CCQE ratio "<<endl;
File<<"# Energy ... CH2 with FSI ... carbon without FSI"<<endl;

for (int m=0; m<13; m++)//should be 13
{double energy = xx[m];//now we have to substitute...

Ccounter1PiPlusIso=0;
CcounterQEIso=0;
Ccounter1PiPlusFSI=0;
CcounterQEFSI=0;

Hcounter1PiPlusIso=0;
HcounterQEIso=0;
Hcounter1PiPlusFSI=0;
HcounterQEFSI=0;

sprintf(crossH, "MB_E=%d_Hydrogen_100Kilo.root.txt" , energia[m]);
sprintf(crossC, "MB_E=%d_Carbon_100Kilo.root.txt" , energia[m]);
sprintf(eventsH, "MB_E=%d_Hydrogen_100Kilo.root" , energia[m]);
sprintf(eventsC, "MB_E=%d_Carbon_100Kilo.root" , energia[m]);


ifstream InputH (crossH);
ifstream InputC (crossC);

crosssections(przekrojeH, InputH);
crosssections(przekrojeC, InputC);

double crossC = przekrojeC[0] + przekrojeC[2] + przekrojeC[4] + przekrojeC[6];
double crossH = przekrojeH[2] + przekrojeH[4];

cout<<crossC<<"   "<<crossH<<endl;

	TFile *tf1 = new TFile(eventsC);
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event * e= new event();
		
	tt1->SetBranchAddress("e",&e);
	
	for( int i=0; i < events; i++ )
	{
	tt1->GetEntry(i);

	if (e->nof(211)==0 && e->nof(-211)==0 && e->nof(111)==0)
	CcounterQEIso++;

	if (e->nof(211)==1 && e->nof(-211)==0 && e->nof(111)==0)
	Ccounter1PiPlusIso++;
	
	if (e->fof(211)==0 && e->fof(-211)==0 && e->fof(111)==0)
	CcounterQEFSI++;

	if (e->fof(211)==1 && e->fof(-211)==0 && e->fof(111)==0)
	Ccounter1PiPlusFSI++;

	}
	delete e;

	TFile *tf2 = new TFile(eventsH);
	TTree *tt2 = (TTree*)tf2->Get("treeout");
	event * e= new event();
		
	tt2->SetBranchAddress("e",&e);
	
	for( int i=0; i < events; i++ )
	{
	tt2->GetEntry(i);

		if (e->nof(211)==1 && e->nof(-211)==0 && e->nof(111)==0 )//selection of single Pi+ 
	Hcounter1PiPlusFSI++;
	
	if (e->nof(211)==0 && e->nof(-211)==0 && e->nof(111)==0 )//selection of single Pi+ 
	HcounterQEFSI++;

	}
	
	delete e;

double ratioFull=( Ccounter1PiPlusFSI*crossC*12  +  Hcounter1PiPlusFSI*crossH*2 )/( CcounterQEFSI*crossC*12  +  HcounterQEFSI*crossH*2 );
double ratioIso=Ccounter1PiPlusIso/CcounterQEIso;

wynikFSI[m]=ratioFull;
wynikIso[m]=ratioIso;

File<<xx[m]<<"   "<<ratioFull<<"   "<<ratioIso<<endl;
cout<<xx[m]<<"   "<<ratioFull<<"   "<<ratioIso<<endl;

}//the end of loop in energy


/*
ifstream Result ("MB_Ratio_CC1PiPlusOverCCQE.dat");

string yy; double nada;
if (Result)
{
do
{
getline (Result, yy);
getline (Result, yy);

for (int jj=0; jj<3; jj++)
{
Result>>nada; 
Result>>wynikFSI[jj]; 
Result>>wynikIso[jj];
} 
} while (Result);
}
*/
TGraph *NuWroFSI = new TGraph (13, xx, wynikFSI);
TGraph *NuWroIso = new TGraph (13, xx, wynikIso);

/*
c1 = new TCanvas("c1","Test",200,10,700,500);

TCanvas *anl1 = new TCanvas("anl1", "anl",26,26,613,826);
   anl1->Range(-0.2509977,-0.1717173,1.254176,1.076966);

   anl1->SetBorderSize(2);
   anl1->SetTickx();
   anl1->SetTicky();
   anl1->SetLeftMargin(0.1303828);
   anl1->SetRightMargin(0.06937799);
   anl1->SetTopMargin(0.0612855);
   anl1->SetBottomMargin(0.1375187);
   anl1->SetFrameBorderMode(0);
   anl1->SetFrameBorderMode(0);


MBRatioIso->SetTitle();
MBRatioIso->SetMarkerColor(4);
MBRatioIso->SetMarkerStyle(21);
MBRatioIso->GetXaxis()->SetTitle("Neutrino energy [MeV]");
MBRatioIso->GetXaxis()->CenterTitle();
MBRatioIso->Draw("AP");

NuWroIso->SetLineWidth(3);
NuWroIso->SetLineColor(2);
NuWroIso->Draw();

leg = new TLegend(0.5,0.15,0.89,0.3);
leg->AddEntry(MBRatioIso,"MiniBooNE data","P");
leg->AddEntry(NuWroIso,"NuWro","L");
leg->Draw();

leg2 = new TLegend (0.7, 0.9, 0.9, 1.0);
leg2->SetHeader("June 3, 2009"); 
leg2->Draw();

c1->Update();
c1->Print("Ratio_Iso_June3.eps");
*/
TCanvas *c2 = new TCanvas("c2","Test2",200,10,700,500);

//TCanvas *anl1 = new TCanvas("anl1", "anl",26,26,613,826);

   c2->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   c2->SetFillColor(0);
   c2->SetBorderSize(2);
   c2->SetTickx();
   c2->SetTicky();
   c2->SetLeftMargin(0.1303828);
   c2->SetRightMargin(0.06937799);
   c2->SetTopMargin(0.0612855);
   c2->SetBottomMargin(0.1375187);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);

MBRatioFSI->SetTitle();
MBRatioFSI->SetMarkerColor(1);
MBRatioFSI->SetMarkerStyle(21);
MBRatioFSI->Draw("AP");

NuWroFSI->SetLineWidth(6);
NuWroFSI->SetLineStyle(2);
NuWroFSI->SetLineColor(1);
MBRatioFSI->GetXaxis()->SetTitle("Neutrino energy [GeV]");
MBRatioFSI->GetXaxis()->CenterTitle();
MBRatioFSI->GetYaxis()->SetTitle("Ratio CC1#pi^{+}/CCQE on CH_{2}");
MBRatioFSI->GetYaxis()->CenterTitle();
NuWroFSI->Draw();

leg = new TLegend(0.3,0.15,0.89,0.3, NULL,"brNDC");
   leg->SetTextFont(22);
   leg->SetTextSize(0.05);
   leg->SetLineColor(0);
   leg->SetLineStyle(2);
   leg->SetLineWidth(3);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

leg->AddEntry(MBRatioFSI,"MiniBooNE data","P");
leg->AddEntry(NuWroFSI,"NuWro: dipole fit, M_{A}^{QE}=1.03 GeV","L");
leg->Draw();

c2->Update();
c2->Print("Paper_Ratio_CH2_October12.eps");


}
//koniec void
