#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include "jednostki.h"
#include "data.h"

using namespace std;

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

int main()
{
	ofstream file("ccqe_test2_numu_sf.txt");
	//ofstream file("ccqe_test2_numu_fg.txt");
	
	file << "#muon neutrino energy = 550 - 600, target = carbon, cos(theta) = 0.9 - 1.0, Spectral Function" << endl << endl;
	//file << "#muon neutrino energy = 550 - 600, target = carbon, cos(theta) = 0.9 - 1.0, Fermi Gas (binding energy = 27)" << endl << endl;
	
	TFile *tf1 = new TFile("ccqe/E550_600_6_6_14_SF.root");
	//TFile *tf1 = new TFile("ccqe/E550_600_6_6_14_FG.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	int count = 0;
	
	for (int i = 0; i < 500000; i++)
	{
		tt1->GetEntry(i);
		
		double mom = e1->out[0].momentum();
		double cos = e1->out[0].p().z/mom;
				
		if (cos >= 0.9)
		{
			//if ((mom >= 425 and mom < 475) or (mom >= 575 and mom < 625))
			if (mom >= 575 and mom < 625)
			{
				count++;
								
				file << "LEPTON MOMENTUM: " << mom << endl;
				file << "COSINUS(THETA): " << cos << endl << endl;
				
				particle p1 = e1->in[0];
				file << "Neutrino: (" << p1.E() << ", " << p1.p().x << ", " << p1.p().y << ", " << p1.p().z << ") " << endl;
				particle p2 = e1->in[1];
				file << "Initial nucleon: (" << p2.E() << ", " << p2.p().x << ", " << p2.p().y << ", " << p2.p().z << ") " << endl;
				particle p3 = e1->out[0];
				file << "Lepton: (" << p3.E() << ", " << p3.p().x << ", " << p3.p().y << ", " << p3.p().z << ") " << endl;
				particle p4 = e1->out[1];
				file << "Final nucleon: (" << p4.E() << ", " << p4.p().x << ", " << p4.p().y << ", " << p4.p().z << ") " << endl << endl;

				file << " final energy - inital energy = " << p3.E() + p4.E() - p2.E() - p1.E() + 27 << endl;
				
				file << " final momentum - inital momentum = (" << p4.p().x + p3.p().x - p2.p().x - p1.p().x << ", " << p4.p().y + p3.p().y - p2.p().y - p1.p().y << ", " << p4.p().z + p3.p().z - p2.p().z - p1.p().z << ")" << endl;
				
				file << endl << "------------------------------------------------------------" << endl << endl;
			}
		}
		
		if (count == 9) break;
			
		cout << 100*i/500000 << "%\r" << flush;
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	return 0;
}

/*
int main()
{
	double oxygen[36];
	double carbon[36];
	
	for (int i = 0; i <= 35; i++)
	{
		double energy = 250 + i*50;
		
		stringstream temp;
		string x;

		temp << energy;
		temp >> x;
	
		string cfile = "ccqe/numu_carbon_" + x + ".txt";
		string ofile = "ccqe/numu_oxygen_" + x + ".txt";
		
		carbon[i] = crosssection(cfile);
		oxygen[i] = crosssection(ofile);
	}
	
	ofstream file("ccqe.txt");
	
	for (int i = 0; i <= 35; i++)
		file << 250 + i*50 << " " << carbon[i] << " " << oxygen[i] << endl;
		
	file.close();
}

/*void put (double value, double *source, double *target, double &rest, const int bins)
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
int main(int argc, char** argv)
{
    string in = argv[1];
    string out = argv[2];
    string out2 = argv[3];
    string out3 = argv[4];
    string out4 = argv[5];
    string out5 = argv[6];
    
    const int events = 1000000;
		
	int pi[3][5] = {0};
	string przedzialy[3][5] = {
	{"przod - przod", "przod - lekki tyl", "przod - mocny tyl", "przod - abs", "przod - inny"}, 
	{"lekki tyl - przod", "lekki tyl - lekki tyl", "lekki tyl - mocny tyl", "lekki tyl - abs", "lekki tyl - inny"},
	{"mocny tyl - przod", "mocny tyl - lekki tyl", "mocny tyl - mocny tyl", "mocny tyl - abs", "mocny tyl - inny"}
	};
	
	int mom[3][20] = {{0}};
    int q2[2][10] = {{0}};
    double norm[10] = {0};
    
    int momentum[2][25] = {{0}};
    
    int pt[2][25] = {{0}};
    
    TFile *tf1 = new TFile (in.c_str());
    TTree *tt1 = (TTree*)tf1->Get("treeout");
    event *e1 = new event();
    
    tt1->SetBranchAddress("e",&e1);
    
    for (int i = 0; i < events; i++)
    {
		tt1->GetEntry(i);
		
		int pion0 = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		int pion1 = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);
		
		int przed;
		int po;
		
		int q = -e1->q2()/1000000/0.2;
		
		if (q < 10)
			norm[q]++;
				
		if (pion0 == 1)
		{			
			for (int k = 0; k < e1->n(); k++)
			{
				if (e1->out[k].pdg == 111)
				{
					double cos = e1->out[k].p().z/e1->out[k].momentum();
					
					int a = e1->out[k].momentum()/50;
					if (a < 20)
						mom[przed][a]++;
					
					int momprzed = e1->out[k].momentum()/20;
					
					if (cos >= 0) przed = 0;
					else if (cos >= -0.5) przed = 1;
					else przed = 2;
					
					if (cos < 0 and q < 10) q2[0][q]++;
					if (cos < 0 and momprzed < 25) momentum[0][momprzed]++;
					
					double ped = e1->out[k].p().y*e1->out[k].p().y + e1->out[k].p().x*e1->out[k].p().x;
					
					int b = sqrt(ped)/20;
					if (b < 25) pt[0][b]++;
				}
			}	
		}
		
		if (pion1 == 1)
		{
			for (int k = 0; k < e1->f(); k++)
			{
				if (e1->post[k].pdg == 111)
				{
					double cos = e1->post[k].p().z/e1->post[k].momentum();
					
					if (cos >= 0) po = 0;
					else if (cos >= -0.5) po = 1;
					else po = 2;
					
					if (cos < 0 and q < 10) q2[1][q]++;
					
					int mompo = e1->post[k].momentum()/20;
					
					if (cos < 0 and mompo < 25) momentum[1][mompo]++;
					
					double ped = e1->post[k].p().y*e1->post[k].p().y + e1->post[k].p().x*e1->post[k].p().x;
					
					int b = sqrt(ped)/20;
					if (b < 25) pt[1][b]++;
				}				
			}
		}
		else if (pion1 == 0) po = 3;
		else po = 4;
			
		if (pion0 == 1)
			pi[przed][po]++;
			
				
		cout << 100*i/events << "%\r" << flush;
	
    }
    
    delete tt1;
    delete tf1;
    delete e1;
    
    ofstream plik (out.c_str());
    
    plik << setprecision(3);
    
    for (int i = 0; i < 3; i++)
	{
		for (int k = 0; k < 5; k++)
			plik << przedzialy[i][k] << ": " << pi[i][k] << " (" << 100.0*pi[i][k]/(pi[i][0]+pi[i][1]+pi[i][2]+pi[i][3]+pi[i][4]) << " %)" << endl;
			
		plik << endl;
	}
    
    plik.close();
    
    ofstream plik2 (out2.c_str());
    
    for (int i = 0; i < 20; i++)
		plik2 << 50*i + 25 << " " << mom[0][i] << " " << mom[1][i] << " " << mom[2][i] << endl;
		
	plik2.close();
	
	ofstream plik3 (out3.c_str());
	
	for (int i = 0; i < 10; i++)
		plik3 << 0.2*i + 0.1 << " " << q2[0][i]/norm[i] << " " << q2[1][i]/norm[i] << endl;
		
	plik3.close();
	
	ofstream plik4 (out4.c_str());
	
	double sum[2] = {0};
	double sumpt[2] = {0};
	
	for (int i = 0; i < 25; i++)
		{
			sum[0] += momentum[0][i];
			sum[1] += momentum[1][i];
			sumpt[0] += pt[0][i];
			sumpt[1] += pt[1][i];
		}
		
	for (int i = 0; i < 25; i++)
		plik4 << 20*i + 10 << " " << (double)momentum[0][i]/sum[0] << " " << (double)momentum[1][i]/sum[1] << endl;
	
	plik4.close();
	
	ofstream plik5(out5.c_str());
	
	for (int i = 0; i < 25; i++)
		plik5 << 20*i + 10 << " " << (double)pt[0][i]/sumpt[0] << " " << (double)pt[1][i]/sumpt[1] << endl;
	
    cout << endl << "Done" << endl;
    
    return 1;
}

/*
int main(int argc,  char** argv)
{
	string model = argv[1];
	string in = "nomad/root/" + model + ".root";
	string out1 = "nomad/main/" + model + ".txt";
	string out2 = "nomad/distr/" + model + ".txt";	
		
	const int events = 5000000;
	
	const int bins1 = 6;
	double backpi[bins1] = {0};
	double norm1[bins1] = {0};
	double Q2[bins1];
	
	int counterPi[4] = {0, 0, 0, 0};
	int countpi = 0;
	
	for (int i = 0; i < bins1; i++) Q2[i] = NOMADpionsQ2[i] + NOMADpionsQ2err[i];
	
	const int bins2 = 10;
	
	double mom2[bins2];
	for (int i = 0; i < bins2; i++) mom2[i] = (i + 1.0) * 0.1;
	
	double distr[bins2] = {0};
	double norm2 = 0;
	
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		double q = -e1->q2()/1000000.0;
		double rest = 0;
		
		put(q, Q2, norm1, rest, bins1);
		
		if (e1->in[0].pdg == 14)
			norm2++;
		
		for (int k = 0; k < e1->f(); k++)
		{	
			if (e1->post[k].pdg == -211 && e1->post[k].p().z < 0)			
			{
				if (e1->post[k].momentum() > 350 && e1->post[k].momentum() < 800)
				{
					put(q, Q2, backpi, rest, bins1);		
					countpi++;
				}
				
				if (e1->in[0].pdg == 14)
				{
					double P = e1->post[k].momentum()/1000.0;
					double E = e1->post[k].E()/1000.0;
					P *= P;
						
					put(P, mom2, distr, rest, bins2);
				}
			}
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

		cout << 100*i/events << "%\r" << flush;
	}
	
	delete e1;
	delete tt1;
	delete tf1;
	
	for (int i = 0; i < 4; i++)
		counterPi[i] *= (double)944019.0/events;
	
	ofstream plik(out1.c_str());
	
	plik<<"#number of events with 0,1,2,3 backwards pions: "<<counterPi[0]<<" "<<counterPi[1]<<" "<<counterPi[2]<<" "<<counterPi[3]<<endl;	
	
	for (int i = 0; i < bins1; i++)
		plik<<NOMADpionsQ2[i]<<" "<<backpi[i]/norm1[i]<<endl;
	
	plik.close();
	
	ofstream plik2(out2.c_str());
	
	for (int i = 0; i < bins2; i++)
	{
		mom2[i] = (i + 0.5)*0.1;
		
		double P = sqrt(mom2[i]);
		double E = sqrt(mom2[i] + 0.0196);
		
		plik2 << mom2[i] << " " << E/P*distr[i]/norm2/0.1 << endl;
	}
	
	plik2.close();
	
	return 1;
}	*/
/*
int main(int argc,  char** argv)
{
	string in = argv[1];
	string out = argv[2];
		
	const int events = 1000000;
	
	double pi0[4] = {0};
	double nopi0[4] = {0}; 
		
	TFile *tf1 = new TFile(in.c_str());
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	for (int i = 0; i < events; i++)
	{
		tt1->GetEntry(i);
		
		int pion0 = 100*e1->nof(211) + 10*e1->nof(-211) + e1->nof(111);
		int pion1 = 100*e1->fof(211) + 10*e1->fof(-211) + e1->fof(111);

		if (pion0 == 1)
		{
			switch (pion1)
			{
				case   1: pi0[0]++; break;
				case   0: pi0[1]++; break;
				case  10: pi0[2]++; break;
				case 100: pi0[2]++; break;
				default: pi0[3]++; break;
			}
		}
		
		if (pion1 == 1)
		{
			switch (pion0)
			{
				case   1: nopi0[0]++; break;
				case   0: nopi0[1]++; break;
				case  10: nopi0[2]++; break;
				case 100: nopi0[2]++; break;
				default: nopi0[3]++; break;
			}
		}
	}
	
	ofstream plik(out.c_str());
			
	double sum0 = pi0[0] + pi0[1] + pi0[2] + pi0[3];		
	double sum1 = nopi0[0] + nopi0[1] + nopi0[2] + nopi0[3];		
			
	plik << "There was 1pi0 in primary vertex: " << endl << endl;
	plik << 100*pi0[0]/sum0 << "% events had 1pi0 after fsi" << endl;		
	plik << 100*pi0[1]/sum0 << "% events had no pi after fsi" << endl;		
	plik << 100*pi0[2]/sum0 << "% events had charged pi after fsi" << endl;		
	plik << 100*pi0[3]/sum0 << "% events had more pi after fsi" << endl;		
	
	plik << endl << endl;

	plik << "There was 1pi0 after fsi: " << endl << endl;
	plik << 100*nopi0[0]/sum1 << "% events had 1pi0 before fsi" << endl;		
	plik << 100*nopi0[1]/sum1 << "% events had no pi before fsi" << endl;		
	plik << 100*nopi0[2]/sum1 << "% events had charged pi before fsi" << endl;		
	plik << 100*nopi0[3]/sum1 << "% events had more pi before fsi" << endl;		
				
	plik.close();
	delete e1;
	delete tt1;
	delete tf1;
	return 1;	
}	
/*
int main()
{
	TFile *tf1 = new TFile("fz1.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double pi[100] = {0};
	double nu[100] = {0};
	int pin[100] = {0};
	int nun[100] = {0}; 
		
	for (int i = 0; i < 5000000; i++)
	{
		tt1->GetEntry(i);
		
		int ile = e1->post.size();
		
		for (int k = 0; k < ile; k++)
		{
			double momentum = e1->post[k].momentum();
			int pdg = e1->post[k].pdg; 
			double fz = e1->post[k].fz;
			
			int a = momentum/20;
						
			if (a < 100 and fz > 0 and fz < 5)
			{
				if (pdg == 2212 or pdg == 2112)
				{
					nun[a]++;
					nu[a] += fz;
				}
				else if (pdg == 211 or pdg == -211 or pdg == 111)
				{
					pin[a]++;
					pi[a] += fz;
				}
			}
		}
		
		cout << 100*i/5000000 << "%\r" << flush;
	}
	
	ofstream fzplik("fz1.txt");
	
	for (int o = 0; o < 100; o++)
	    fzplik << (o + 0.5)*20.0 << " " << pi[o]/pin[o] << " " << nu[o]/nun[o] << endl;
			
	fzplik.close();
	delete e1;
	delete tt1;
	delete tf1;
	return 1;
}

/*
	params pC, pI, pt;
	
	pC.nucleus_p = 6;
	pC.nucleus_n = 6;
	pC.nucleus_model = 1;
	
	pI.nucleus_p = 26;
	pI.nucleus_n = 30;
	pI.nucleus_model = 1;

	pt.nucleus_p = 26;
	pt.nucleus_n = 30;
	pt.nucleus_model = 1;
	
	nucleus *pnucleusC = make_nucleus (pC);
    nucleus & jadroC = *pnucleusC;
    
    nucleus *pnucleusI = make_nucleus (pI);
    nucleus & jadroI = *pnucleusI;
    
    nucleus *pnucleust = make_nucleus (pt);
    nucleus &jadrot = *pnucleust;
    
	double radiusC = jadroC.radius();
	double radiusI = jadroI.radius();
	double radiust = jadrot.radius();
	
	cout<<endl<<radiust/fermi<<endl;

	ofstream plikC("carbon-test.txt");
	ofstream plikI("iron-test.txt");
	
	ofstream fileC("dens_carbon.txt");
	ofstream fileI("dens_iron.txt");
	
	plikC << "radius: "<< radiusC/fermi << endl<<"------------------------------------------------------------"<< endl<< endl;
	plikI << "radius: "<< radiusI/fermi << endl<<"------------------------------------------------------------"<< endl<< endl;
	
	double r = -0.2*fermi;
	
	const double energy[8] = {11, 30, 49, 85, 128, 184, 250, 350};
	
	do
	{
		r += 0.2*fermi;
		
		double densC = jadroC.density(0, r);
		
		int xsec = 2;
		
		plikC << "r: " << r/fermi << "          rho/rho0: " << densC*fermi*fermi*fermi/0.17 << endl << "------------------------------" << endl << endl;
		
		plikC << "Ek     sii     sij     sabs" <<endl<<endl;
		
		fileC << r/fermi << " " << densC/0.17*fermi3 << " " << densC*radiusC*radiusC*fermi << endl;
				
		for (int k = 0; k < 8; k++)
		{
			double Ek = energy[k];
			
			piondata pdC (pion_params (Ek, xsec, densC), xsec);
			double siiC = pdC.sigma (0, Ek, xsec, densC);
			double sijC = pdC.sigma (1, Ek, xsec, densC);
			double sabsC = pdC.sigma (2, Ek, xsec, densC);

			plikC << Ek << " " << siiC/millibarn << " " << sijC/millibarn << " " << sabsC/millibarn << endl;

		}
		
		plikC << endl << endl;
		
	}while(r < radiusC);

	r = -0.2*fermi;

	do
	{
		r += 0.2*fermi;
		
		double denst = jadrot.density(0, r)*fermi3;
		
		cout << r/fermi << " " << denst << endl;
		
	}while(r < radiust);

	r = -0.2*fermi;
	
	ofstream file("density.txt");
	
	double factor = 0.5/jadroI.density(0,0);
	
	do
	{
		r += 0.2*fermi;
				
		double densI = jadroI.density(0, r);
		
		file<<r/fermi<<" "<<densI*factor<<endl;
		
		int xsec = 3;
		
		plikI << "r: " << r/fermi << "          rho/rho0: " << densI*fermi3/0.17 << endl << "------------------------------" << endl << endl;
		
		plikI << "Ek     sii     sij     sabs" <<endl<<endl;
		
		fileI << r/fermi << " " << densI*fermi3/0.17 << " " << densI *radiusI *radiusI*fermi <<endl;
		
		for (int k = 0; k < 8; k++)
		{
			double Ek = energy[k];
			
			piondata pdI (pion_params (Ek, xsec, densI), xsec);
			double siiI = pdI.sigma (0, Ek, xsec, densI);
			double sijI = pdI.sigma (1, Ek, xsec, densI);
			double sabsI = pdI.sigma (2, Ek, xsec, densI);
			
			plikI << Ek << " " << siiI/millibarn << " " << sijI/millibarn << " " << sabsI/millibarn << endl;
		}
		
		plikI << endl << endl;
		
	}while(r < radiusI);
	
	plikC.close();
	plikI.close();
	file.close();
	fileI.close();
	fileC.close();
	
	ifstream inf("metro-test.txt");
	
	double ra;
	double radius[21] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
	int rad[21];
	
	for (int k = 0; k < 21; k++) rad[k] = 0;
	
	if (inf)
	{
		do
		{
			inf >> ra;
			
			for (int k = 0; k < 20; k ++)
			{
				if (ra > radius[k] && ra <= radius[k+1]) rad[k]++;
			}

		}while (inf);
	}
	
	ofstream res("metro-res.txt");
	
	for (int k = 0; k < 20; k++)
	{
		res << "r: " << radius[k] << " - " << radius[k+1] << "     number of abs: " << rad[k] <<endl;
	}
	
	res.close();
	
	ofstream comp("oset_comp.txt");
	
	r = -0.05*fermi;
	
	do
	{
		r += 0.05*fermi;
		
		double d = jadroI.density(0, r);
			
		piondata pd (pion_params (165.0, 1, d), 1);
		double sii = pd.sigma (0, 165.0, 1, d); // millibarn/10.0;
		double sij = pd.sigma (1, 165.0, 1, d); //millibarn/10.0;
		double sabs = pd.sigma (2, 165.0, 1, d); //millibarn/10.0;
		
		double d2 = d/(140.0*140.0*140.0); 
		
		//d *= fermi3;
		
		//double Pabs = 1.0 - exp(-sabs*d);
		
		//double Pqe = 1.0 - exp(-(sii + sij)*d);
			
		double Pabs = 0.5*sabs*d*fermi;
		double Pqe = 0.5*(sii+sij)*d*fermi;	
				
		comp << r/fermi << " " << Pabs << " " << Pqe << " " << d2 << endl;
		
	}while(r/fermi < 7.0);
	
	r = -0.01*fermi;
	double A = 0;
	const double Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862;
	double promien;
	
	do
	{
		r += 0.01*fermi;
		
		double d = jadroI.density(0, r)*fermi3;
		promien = r/fermi;
		
		A += d*promien*promien*4.0*Pi*0.01;
				
	}while(promien < 20.0);
	
	cout << "A = " << A << endl;
	
	comp.close();
	return 1;
}
*/
