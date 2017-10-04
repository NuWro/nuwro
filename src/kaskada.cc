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
  set_dirs(argv[0]);
  args a("kaskada","kaskada.txt","kaskada.root");
  a.read(argc,argv);
  params p;
  p.read(a.input);
  p.read(a.params,"Command line");
  p.list();
  frandom_init(p.random_seed);
  
  try
  {
    input_data input_test( p );
    input_test.initialize();
    input_test.load_data();
    data_container *test = input_test.get_data_container();
    // for(int i=0;i<20000000;i++)
    //   test->set_input_point(i);
    test->set_input_point(800);
    cout << test->get_value(2) << "\n";
  }
  catch( char const* ex )
  {
    cout << ex << endl;
    return 1;
  }

  event *e=new event;
  TFile *f= new TFile(a.output,"recreate");
  TTree * t2=new TTree("treeout","Tree of events");
  t2->Branch("e","event",&e); ///< tree1 has only one branch (branch of events)

  try
  {
	  nucleus* nucl= make_nucleus(p);
	  beam_uniform b(p);
	  b.check_energy();         
  }
  catch(const char* w)
  {
	  cout<<endl<<"Exception:     "<<w<<endl<<endl;
	  return 1;
  }

  beam_uniform b(p);
  nucleus* nucl= make_nucleus(p); 

  for(int i=0;i<p.number_of_events;i++)
  {  cout<<"event="<<i<<" begin      \r"; 
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
	 kaskada k(p,*e);
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
