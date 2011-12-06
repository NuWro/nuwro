#include "vivisection.h"

int main()
{
	TFile *tf1 = new TFile("root_files/PiTrans_he_Carbon_1m_nofz_Oset_2011.05.04.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	const int bins = 30;
	
	double r[bins];
	double ile[bins]; zero(ile, bins);
	double rest;
	
	for (int i = 0; i < bins; i++) r[i] = (i + 1)*0.15;
	
	for (int i = 0; i < 1000000; i++)
	{
		tt1->GetEntry(i);
		
		if (e1->flag.dis and e1->nof(111) + e1->nof(211) + e1->nof(-211) == 1) put(e1->out[0].r.length()/fermi, r, ile, rest, bins);
	}

	int norma = 0;

	for (int i = 0; i < bins; i++) r[i] = (i + 0.5)*0.15;
	for (int i = 0; i < bins; i++) norma += ile[i];
	
	ofstream plik("newr.txt");
	
	for (int i = 0; i < bins; i++) plik << r[i] << " " << 1000.0*ile[i]/norma << endl;
	
	plik.close();
	
	delete e1;
	delete tt1;
	delete tf1;
	return 1;
}	
/*		
	TFile *tf1 = new TFile("fztest.root");
	TTree *tt1 = (TTree*)tf1->Get("treeout");
	event *e1   = new event();
		
	tt1->SetBranchAddress("e",&e1);
	
	double mom[50];
	double pi[50];
	double pi2[50];
	double nu[50];
	int pin[50];
	int pin2[50];
	int nun[50]; 
	
	for (int z = 0; z < 50; z++)
	{
		mom[z] = (z+1)*100.0;
		pi[z] = 0;
		pi2[z] = 0;
		nu[z] = 0;
		pin[z] = 0;
		pin2[z] = 0;
		nun[z] = 0;
	}	
	
	for (int i = 0; i < 1000000; i++)
	{
		tt1->GetEntry(i);
		
		int ile = e1->post.size();
		int pion = e1->nof(211)+e1->nof(-211)+e1->nof(111);
		
		for (int k = 0; k < ile; k++)
		{
			if (e1->flag.res)
			{
			double momentum = e1->post[k].momentum();
			int pdg = e1->post[k].pdg; 
			double fz = e1->post[k].fz/fermi;
			
			int ktory = 0;
			
			for (int a = 1; a < 50; a++)
			{
				if (momentum >= mom[a-1] and momentum < mom[a])
				{
					ktory = a;
					break;
				}
			}
		
			if (pdg == 211 or pdg == -211 or pdg == 111)
			{
				pin[ktory]++;
				pi[ktory] += fz;
				if (pion == 1)
				{
					pi2[ktory] += fz;
					pin2[ktory]++;
				}
			}
			else if (pdg == 2112 or pdg == 2212)
			{
				nun[ktory]++;
				nu[ktory] += fz;
			}
		}
		}
	}
	
	ofstream fzplik("fztest.txt");
	
	for (int o = 0; o < 50; o++)
	{
		if (nun[o] != 0) nu[o] /= nun[o];
		if (pin[o] != 0) {pi[o] /= pin[o]; pi2[o] /= pin2[o];}
	    
	    fzplik << (o + 0.5)*100.0 << " " << nu[o] << " " << pi[o] << " " << pi2[o] << endl;
	}
			
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
