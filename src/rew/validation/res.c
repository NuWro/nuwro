#include <iostream>
#include <string>
//#include "event1.h" // uncomment when using root 6 (with nuwro/src in include path)

using namespace std;

#define M_PI 3.14159265358979323846

void saveHist(TH1* h, const char* filename) {
  ofstream file(filename);

  for (int i = 1; i <= h->GetNbinsX(); i++)
    file << h->GetBinLowEdge(i) + h->GetBinWidth(i) / 2 << " " << h->GetBinContent(i) << "\n";

  file.close();
}

double getSigma(string filename, int dyn = -1) {
  static const unsigned int nDyn = 10;
  static const unsigned int nCol = 4;

  vector<double> xsec(nDyn);

  string dummy;

  ifstream input(filename.c_str());

  getline(input, dummy);

  for (unsigned int i = 0; i < nDyn; i++)
    for (unsigned int j = 0; j < nCol; j++) input >> xsec[i];

  if (dyn >= 0 && dyn < nDyn) return xsec[dyn];

  double total = 0;

  for (unsigned int k = 0; k < nDyn; k++) total += xsec[k];

  return total;
}

// normalize to cross section per nucleon
void normalize(TH1* h, const double factor) {
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    h->SetBinContent(i, h->GetBinContent(i) * factor / h->GetBinWidth(i));
  }
}

void res(const string prefix) {
  string root_file = prefix + ".root";

  TFile* file = new TFile(root_file.c_str());

  TTree* tree = (TTree*)file->Get("treeout");
  event* e = new event();

  tree->SetBranchAddress("e", &e);

  TH1D* h_pip = new TH1D("h_pip", "pip", 50, 0, 2);
  TH1D* h_pi0 = new TH1D("h_pi0", "pi0", 50, 0, 2);
  TH1D* h_pim = new TH1D("h_pim", "pim", 50, 0, 2);

  const unsigned int nEvents = tree->GetEntries();

  double sum = 0;

  for (unsigned int i = 0; i < nEvents; i++) {
    tree->GetEntry(i);

    if (e->nof(211) + e->nof(111) + e->nof(-211) != 1) continue;
    // if (e->W() > 1210) continue;

    for (unsigned int j = 0; j < e->out.size(); ++j)
      if (e->out[j].pion()) switch (e->out[j].pdg) {
          case 211: h_pip->Fill(e->out[j].momentum() / 1000.0, e->weight); break;
          case 111: h_pi0->Fill(e->out[j].momentum() / 1000.0, e->weight); break;
          case -211: h_pim->Fill(e->out[j].momentum() / 1000.0, e->weight); break;
        }

    sum += e->weight;
    cout << 100*i/nEvents << "%\r" << flush;
  }

  saveHist(h_pip, (prefix + "_pip.txt").c_str());
  saveHist(h_pi0, (prefix + "_pi0.txt").c_str());
  saveHist(h_pim, (prefix + "_pim.txt").c_str());

  cout << "Sum of weights = " << sum << "\n\n";
}
