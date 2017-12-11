#include <iostream>
#include <string>
//#include "event1.h" // uncomment when using root 6 (with nuwro/src in include path)

using namespace std;

#define M_PI 3.14159265358979323846

void saveHist (TH1* h, const char* filename)
{
   ofstream file(filename);

   for (int i = 1; i <= h->GetNbinsX(); i++)
      file << h->GetBinLowEdge(i) + h->GetBinWidth(i)/2 << " " << h->GetBinContent(i) << "\n";

   file.close();
}

double getSigma (string filename, int dyn = -1)
{
    static const unsigned int nDyn = 10;
    static const unsigned int nCol = 4;

    vector <double> xsec(nDyn);

    string dummy;

    ifstream input(filename.c_str());

    getline (input, dummy);

    for (unsigned int i = 0; i < nDyn; i++)
    	for (unsigned int j = 0; j < nCol; j++) input >> xsec[i];

    if (dyn >= 0 && dyn < nDyn) return xsec[dyn];

    double total = 0;

    for (unsigned int k = 0; k < nDyn; k++) total += xsec[k];

    return total;
}

// normalize to cross section per nucleon
void normalize(TH1* h, const double factor)
{
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        h->SetBinContent(i, h->GetBinContent(i) * factor / h->GetBinWidth(i));
    }
}

void qel(const string prefix)
{
    string root_file = prefix + ".root";

    TFile *file = new TFile(root_file.c_str());

    TTree *tree = (TTree*)file->Get("treeout");
    event *e    = new event();

    tree->SetBranchAddress("e", &e);

    TH1D* h_Q2 = new TH1D("hQ2", "Q2", 50, 0, 2);
    TH1D* h_Tk = new TH1D("hTk", "Lepton kinetic energy", 50, 0, 2);
    TH1D* h_ang = new TH1D("hang", "Scattering angle", 50, -1, 1);
    
    const unsigned int nEvents = tree->GetEntries();

    for (unsigned int i = 0; i < nEvents; i++)
    {
        tree->GetEntry(i);

        const double Q2 = -e->q2() / 1000000.0;
        const double Tk = e->out[0].Ek() / 1000.0;
        const double ang = e->out[0].p().z / e->out[0].momentum();

        h_Q2 -> Fill(Q2, e->weight);
        h_Tk -> Fill(Tk, e->weight);
        h_ang -> Fill(ang, e->weight);
    }
    
    saveHist(h_Q2, (prefix + "_Q2.txt").c_str());
    saveHist(h_Tk, (prefix + "_Tk.txt").c_str());
    saveHist(h_ang, (prefix + "_ang.txt").c_str());
}
