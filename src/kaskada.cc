#include "kaskada7.h"
#include "beam.h"
#include "beam_uniform.h"
#include <TFile.h>
#include <TTree.h>
#include "args.h"
#include "nucleusmaker.h"
#include "dirs.h"
#include "input_data.h"

int main(int argc,  char** argv)
{
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
  
  // load the input data
  input_data input;
  try
  {
    input.initialize( p );
    input.load_data();
  }
  catch(char const* ex)
  {
    cout << "Exception: " << ex << endl;
    return 1;
  }

  // prepare the root output
  event *e = new event;
  TFile *f = new TFile(a.output,"recreate");
  TTree *t2= new TTree("treeout","Tree of events");
  t2->Branch("e","event",&e);   // tree1 has only one branch (branch of events)

  // make the nucleus and beam
  try
  {
    nucleus* nucl= make_nucleus(p);
    beam_uniform b(p);
    b.check_energy();
  }
  catch(const char* ex)
  {
    cout << "Exception: " << ex << endl;
    return 1;
  }

  beam_uniform b(p);                // WTF?
  nucleus* nucl= make_nucleus(p);   // WTF?

  for(int i=0;i<p.number_of_events;i++)
  {
    cout<<"event="<<i<<" begin      \r";
    e=new event;
    e->weight = 1;
    e->par = p;
    particle p0=b.shoot();
    p0.r=start_point(nucl,p);
    e->out.push_back(p0);
    e->flag.qel = 0;
    e->flag.res = 0;
    particle pi = nucl->get_nucleon(p0.r);
    e->in.push_back(pi);
    e->in.push_back(pi);
    //e->in[1].p4().t = 0;
    kaskada k(p,*e,&input);
    k.kaskadaevent();
    t2->Fill();
    delete e;
    cout<<"event "<<i<<": completed.\r";
  }

  f->Write();
  delete nucl;
  delete t2;
  delete f;
  genrand_write_state();
  return 0;
}
