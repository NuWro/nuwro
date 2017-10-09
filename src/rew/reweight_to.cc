#include <iostream>
#include "../event1.h"
#include "../nucleus.h"
#include "Reweighters.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

void SetupSPP(params&);  // defined in rewRES.cc

void info(const string msg) { std::cout << "\033[32m[INFO] " << msg << "\033[0m\n"; }

void usage() {
  // print usage and exit
  info(
      "Usage: "
      "reweight_to <nuwro_output.root> "
      "-p par1 val1 -p par2 val2 ..."
      "[-o <weighted_events.root>] "
      "[--no_weights] [--no_events]");
  exit(1);
}

void error(const string msg) {
  std::cerr << "\033[31m[ERROR] " << msg << "\033[0m\n";
  usage();
}

void not_reweightable(const string par) {
  // print available parameters to reweight and exit
  error("Parameter \"" + par + "\" can not be reweighted. Try one of:\n");
  rew.list(cerr);
  exit(1);
}

enum {
  both = 0,  // create both outputs (default)
  events,    // create only file with weighted events
  weights    // create only file with weights
} mode;

int main(int argc, char* argv[]) {
  // at least input and output files and on parameter to reweight is required
  if (argc < 7) usage();

  // setup input tree
  TFile* f1 = new TFile(argv[1]);
  if (f1->IsZombie()) {
    std::cout << "[ERROR] Couldn't open input file: " << argv[1] << std::endl;
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
  char* outname = NULL;      // the full NuWro output file with new weight
  char* weightsname = NULL;  // the output ROOT file with weights only

  for (int i = 2; i < argc; i++) {
    // parse command line arguments
    if (string(argv[i]) == "-o")
      outname = argv[++i];
    else if (string(argv[i]) == "--no_events")
      mode = weights;
    else if (string(argv[i]) == "--no_weights")
      mode = events;
    else if (string(argv[i]) == "-p") {
      // check if parameter is given after -p
      if (not argv[i + 1]) error("Parameter name is required after -p");
      // check it parameter value is given
      if (not argv[i + 2]) error("Parameter value is required after -p parameter_name");

      // create parameter to reweight
      RewParam& p = rew(argv[++i]);

      // check if it is reweightable (empty string -> not on RewParams list)
      if (p.name == "") not_reweightable(argv[i]);

      // add a paramter to the reweighters list
      args.push_back(&p);

      // exit if given parameter value is NaN
      try {
        vals.push_back(stod(argv[++i]));
      } catch (std::exception& e) {
        error("Parameter value must be a number");
      }

      // parameter added successfully
      REW(p.engine).active = true;

    } else
      error("Unexpected flag");
  }

  // obviously, there must be something to reweight
  if (args.size() == 0) error("Provide at least one parameter to reweight");

  // if output filename is not given - use input filename with .reweighted suffix
  if (outname == NULL) outname = &(string(argv[1]) + ".reweighted")[0];

  // use .weights suffix for a file with weights only
  weightsname = &(string(outname) + ".weights")[0];

  /// create output files
  TFile* f2 = NULL;
  TTree* t2 = NULL;
  TFile* f3 = NULL;
  TTree* t3 = NULL;

  double weight;  // weights holder

  if (mode != events) {
    f2 = new TFile(weightsname, "recreate");
    t2 = new TTree("weights", "Tree of weights");
    t2->Branch("weight", &weight, "weight/D");
  }

  if (mode != weights) {
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

    for (int j = 0; j < args.size(); j++) args[j]->set(vals[j]);

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
