#include <iostream>
#include "../event1.h"
#include "../nucleus.h"
#include "Reweighters.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

void SetupSPP(params&);  // defined in rewRES.cc

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cerr << "[INFO] Usage: "
            "reweight_to <nuwro_output.root> -o <results_filename.root> -p par1 val1 -p par2 val2 ...\n";
    exit(-1);
  }

  // setup input tree
  TFile* f1 = new TFile(argv[1]);
  if (f1->IsZombie()) {
    std::cout << "Error opening file" << std::endl;
    exit(-1);
  }

  // load events tree
  TTree* t1 = dynamic_cast<TTree*>(f1->Get("treeout"));
  if (!t1) {
    std::cerr << "[ERROR]: Couldn't find TTree (\"treeout\") in file " << argv[1] << "." << std::endl;
    exit(-1);
  }

  // set up event pointer
  event* e = new event;
  t1->SetBranchAddress("e", &e);

  vector<RewParam*> args;    // the list of parameters to reweight
  vector<double> vals;       // the list of new parameters values
  int nargs = 0;             // the number of parameters to reweight
  char* outname = NULL;      // the full NuWro output file with new weight
  char* weightsname = NULL;  // the output ROOT file with weight only

  for (int i = 2; i < argc; i++) {
    if (string(argv[i]) == "-o") {
      outname = argv[++i];
      continue;
    } else if (string(argv[i]) == "-w") {
      weightsname = argv[++i];
      continue;
    } else if (string(argv[i]) == "-p") {
      RewParam& p = rew(argv[++i]);
      if (p.name == "") {
        cerr << "[Error] parameter \"" << argv[i] << "\" can not be reweighted. Try one of:\n";
        rew.list(cerr);
        exit(1);
      } else {
        args.push_back(&p);
        stringstream s(argv[++i]);
        double x;
        if (i < argc && (s >> x)) {
          vals.push_back(x);
          REW(p.engine).active = true;
          nargs++;
        } else {
          cerr << "[Error] parameter \"" << argv[i - 1] << "\" must have numeric value (not \"" << argv[i] << "\").\n";
          exit(1);
        }
      }
    } else {
      cerr << "[Error] -p or -o expected instead of\"" << argv[i] << "\".\n";
      cerr << "[INFO] Usage: "
              "reweight_to <nuwro_output.root> -o <results_filename.root> -p par1 val1 -p par2 val2 ...\n";
      exit(1);
    }
  }

  if (nargs == 0 or (outname == NULL && weightsname == NULL)) {
    cerr << "[INFO] Usage: "
            "reweight_to <nuwro_output.root> (-o <weighted_events.root> | -w <weights_only.root>) -p par1 val1 -p par2 "
            "val2 ...\n";
    exit(1);
  }

  /// create output file
  double weight;

  TFile* f2 = NULL;
  TTree* t2 = NULL;
  TFile* f3 = NULL;
  TTree* t3 = NULL;

  if (weightsname) {
    f2 = new TFile((string(outname) + ".weights").c_str(), "recreate");
    t2 = new TTree("weights", "Tree of weights");
    t2->Branch("weight", &weight, "weight/D");
  }

  if (outname) {
    f3 = new TFile(outname, "recreate");
    t3 = new TTree("treeout", "Tree of events");
    t3->Branch("e", "event", &e);
  }

  // Calculate and save weights for each event
  int n = t1->GetEntries();

  for (int ie = 0; ie < n; ie++) {
    t1->GetEntry(ie);
    REW.init(e->par);
    nucleus t(e->par);
    if (ie == 0) SetupSPP(e->par);

    double nominal = REW.weight(*e, e->par, t);

    for (int j = 0; j < nargs; j++) args[j]->set(vals[j]);

    ff_configure(e->par);

    weight = REW.weight(*e, e->par, t) / nominal;

    cout << weight << " ";  // endl;
    if (t2) t2->Fill();

    if (t3) {
      e->weight *= weight;
      t3->Fill();
    }
  }

  f1->Close();
  delete f1;
  if (f2) {
    f2->Write();
    f2->Close();
  }
  if (f3) {
    f3->Write();
    f3->Close();
  }

  delete e;
  delete f2;
  delete f3;
  cout << "Output file: \"" << outname << "\"" << endl;
  return 0;
}
