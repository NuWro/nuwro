#include "kaskada7.h"
#include "beam.h"
#include "beam_uniform.h"
#include <TFile.h>
#include <TTree.h>
#include "args.h"
#include "nucleusmaker.h"
#include "dirs.h"
#include "input_data.h"
#include "output.h"

void raport(double i, double n, const char* text, int precision)
{
  static int prev=-1;
  int proc=precision*i/n;
  if(proc!=prev)
  {
    prev=proc;
    cerr.precision(3);
    cerr<< "        ";
    cerr<< showpoint << proc*100.0/precision << " % of ";
    cerr<< text << '\r' << flush;
  }
  cerr.precision();
}

int main(int argc,  char** argv)
{
  cout << R"(
  ____________________________________________________________________________
 |                                                                            |
 |                                                                            |
 |                                      `.``      `.       .-.        `       |
 |                                  `-/+os+s  ./ohmN:.    sNNNy:`   .---.+`   |
 |    |\ |     |  |  _  _           :oooyysy: +oodMMd-`  .MMMMM/-   -----s/   |
 |    | \| |_| |/\| |  (_)           `.`oyy+d`   `mMMo-   yMMMyo`   .----d.   |
 |             __        __   __        .yyyoo    :MMN:.  :MMho.    `---h:    |
 |              _) /|   /  \ (__\        :yyoh-    sMMh-`.mMho-     ---h:     |
 |             /__  | . \__/  __/         oyy+h    `mMM+-mMho-     ---h:      |
 |                                        .yyys+    -MMNmMh+-     ---h:       |
 |                                         :yyod.    sMMMh+-     ---h:        |
 |   Wrocław Neutrino Event Generator       oyy+h    `mMh-+.    ---h:         |
 |   https://github.com/NuWro/nuwro         .yyss/  .s/y--oh   .--y-          |
 |                                           :yy+d`.sy+----d/ .--y-           |
 |   J. T. Sobczyk et al.                     oyy++syoy`.--:s.--y-            |
 |   Institute of Theoretical Physics         .yyssyoy.  ------y-             |
 |   University of Wrocław                     :yyyoh.   `----y-              |
 |   Poland                                     osoh.     .-:y-               |
 |                                              `-:.       .:-                |
 |                                                                            |
 |____________________________________________________________________________|
             )";
  cout << R"(
  ____________________________________________________________________________
 |                                                                            |
 |              __                                                            |
 |             /    _   _  _  _   _|  _   |\/|  _   _|  _                     |
 |             \__ (_| _) (_ (_| (_| (-   |  | (_) (_| (-                     |
 |                                                                            |
 |                                                                            |
 |   Hadrons are introduced directly to the cascade model. Only selected      |
 |   parameters are active. The incident particle starting point follows      |
 |   the beam_placement parameter as                                          |
 |   (0) nucleus center                                                       |
 |   (1) random nucleon's position:           transparency mode               |
 |   (2) just under the surface of the nucleus: scattering mode               |
 |                                                                            |
 |____________________________________________________________________________|
             )" << endl;

  frame_top("Simulation parameters");

  // initialize the simulation
  set_dirs(argv[0]);
  args a("kaskada","kaskada.txt","kaskada.root");
  a.read (argc, argv);
  params p;
  p.read (a.input);
  p.read (a.params, "command line");
  p.list (cout);
  p.list (string(a.output)+".par");
  frandom_init(p.random_seed);

  frame_bottom();

  frame_top("Initialize the simulation");

  // prepare the root output
  event *e = new event;
  TFile *f = new TFile(a.output,"recreate");
  TTree *t2= new TTree("treeout","Tree of events");
  t2->Branch("e","event",&e);   // tree1 has only one branch (branch of events)

  // make the nucleus and beam
  nucleus* nucl;
  beam_uniform* b;
  try
  {
    cout << "     -> Building the target nuclei..." << endl;
    nucl = make_nucleus(p);

    cout << "     -> Creating the beam..." << endl;
    b = new beam_uniform(p);
    b->check_energy();
  }
  catch(const char* ex)
  {
    cout << "Exception: " << ex << endl;
    return 1;
  }

  // load the input data
  input_data input;
  try
  {
    cout << "     -> Loading external physical data..." << endl;
    input.initialize(p);
    input.load_data();
  }
  catch(char const* ex)
  {
    cout << "Exception: " << ex << endl;
    return 1;
  }

  frame_bottom();

  frame_top("Run cascade events");

  for(int i=0;i<p.number_of_events;i++)
  {
    e = new event;
    e->weight = 1;
    e->par = p;
    particle p0=b->shoot();
    p0.r=start_point(nucl,p);
    e->out.push_back(p0);
    particle pi = nucl->get_nucleon(p0.r);
    e->in.push_back(pi);
    e->in.push_back(pi);
    kaskada k(p,*e,&input);
    k.kaskadaevent();
    t2->Fill();
    delete e;
    raport(i+1,p.number_of_events,"cascade events ready...",1000);
  }
  f->Write();

  cout << "        100. % of cascade events ready..." << endl;

  frame_bottom();

  frame_top("Finalize the simulation");
  cout << "     " << "-> Generated the output file: \"" << a.output << "\"" << endl;
  frame_bottom();

  delete nucl;
  delete b;
  delete t2;
  delete f;

  genrand_write_state();

  return 0;
}
