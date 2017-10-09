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
  both = 0,    // create both outputs (default)
  no_weights,  // create only file with weighted events
  no_events    // create only file with weights
} mode;

int main(int argc, char* argv[]) {
  // at least input and one parameter to reweight is required
  if (argc < 5) usage();

  // setup input tree
  TFile* tfile_input = new TFile(argv[1]);
  if (tfile_input->IsZombie()) {
    std::cout << "[ERROR] Couldn't open input file: " << argv[1] << std::endl;
    exit(-1);
  }

  // load events tree
  TTree* ttree_input = dynamic_cast<TTree*>(tfile_input->Get("treeout"));
  if (!ttree_input) {
    std::cerr << "[ERROR]: Couldn't find TTree (\"treeout\") in file " << argv[1] << "." << std::endl;
    exit(-1);
  }

  // set up event pointer
  event* e = new event;
  ttree_input->SetBranchAddress("e", &e);

  vector<RewParam*> args;      // the list of parameters to reweight
  vector<double> vals;         // the list of new parameters values
  string out_reweighted = "";  // the full NuWro output file with new weight
  string out_weights = "";     // the output ROOT file with weights only

  for (int i = 2; i < argc; i++) {
    // parse command line arguments
    if (string(argv[i]) == "-o")
      out_reweighted = argv[++i];
    else if (string(argv[i]) == "--no_events")
      mode = no_events;
    else if (string(argv[i]) == "--no_weights")
      mode = no_weights;
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
  if (out_reweighted == "") out_reweighted = string(argv[1]) + ".reweighted";

  // use .weights suffix for a file with weights only
  out_weights = out_reweighted + ".weights";

  // prepare outputs
  TFile* tfile_weights = NULL;  // just weights
  TTree* ttree_weights = NULL;
  TFile* tfile_reweighted = NULL;  // reweighted events
  TTree* ttree_reweighted = NULL;

  double weight;  // weights holder

  if (mode != no_weights) {
    tfile_weights = new TFile(out_weights.c_str(), "recreate");
    ttree_weights = new TTree("weights", "Tree of weights");
    ttree_weights->Branch("weight", &weight, "weight/D");
  }

  if (mode != no_events) {
    tfile_reweighted = new TFile(out_reweighted.c_str(), "recreate");
    ttree_reweighted = new TTree("treeout", "Tree of events");
    ttree_reweighted->Branch("e", "event", &e);
  }

  // Calculate and save weights for each event
  int n_events = ttree_input->GetEntries();

  for (int ie = 0; ie < n_events; ie++) {
    // get i-th event from the ttree
    ttree_input->GetEntry(ie);

    // initialize reweighters with event parameters
    REW.init(e->par);
    nucleus t(e->par);

    // call single pion production setup for first event
    // TODO: disable if pion production is not considered
    if (ie == 0) SetupSPP(e->par);

    // calculate cross section for nominal parameters
    double nominal = REW.weight(*e, e->par, t);

    // change all parameters set up to reweight
    for (int j = 0; j < args.size(); j++) args[j]->set(vals[j]);

    // update form factors
    ff_configure(e->par);

    // calculate cross section with new parameters
    weight = REW.weight(*e, e->par, t) / nominal;

    if (ttree_weights) ttree_weights->Fill();

    if (ttree_reweighted) {
      e->weight *= weight;
      ttree_reweighted->Fill();
    }

    cout << "Reweighting \"" << argv[1] << "\": " << 100 * ie / n_events << "%\r" << flush;
  }

  tfile_input->Close();

  if (tfile_weights) {
    tfile_weights->Write();
    tfile_weights->Close();
  }

  if (tfile_reweighted) {
    tfile_reweighted->Write();
    tfile_reweighted->Close();
  }

  delete e;
  delete tfile_input;
  delete tfile_weights;
  delete tfile_reweighted;

  if (mode != no_events) cout << "Reweighted events saved in: " << out_reweighted << endl;

  if (mode != no_weights) cout << "Weights saved in: " << out_weights << endl;

  return 0;
}
