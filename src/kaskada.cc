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
  a.read(argc, argv);
  params p;
  p.read(a.input);
  p.read(a.params, "command line");
  p.list(cout);
  p.list(string(a.output)+".par");
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

  // process external events
  // format: weight,dyn,PDG,E,px,py,pz,PDG,...
  std::ifstream events_file;
  stringstream name_sstream;
  if(p.kaskada_events)
  {
    name_sstream << p.kaskada_events_file;
    try
    {
      cout << "     -> Preparing the events file..." << endl;
      events_file.open((name_sstream.str()).c_str(), std::ios::in);
      if(!events_file.is_open())
        throw std::runtime_error("unable to load the events file");
    }
    catch(const std::exception& ex)
    {
      cout << "Exception: " << ex.what() << endl;
      return 1;
    }
  }

  frame_bottom();

  frame_top("Run cascade events");

  for(int i=0;i<p.number_of_events;i++)
  {
    e = new event;
    e->par = p;

    // process external events
    if(p.kaskada_events)
    {
      try
      {
        std::string line, field;
        getline(events_file, line);
        if(line.empty())
          throw std::runtime_error("unequal number of events and lines or empty lines somewhere");
        stringstream line_sstream(line);

        // read the cross section
        getline(line_sstream, field, ',');
        e->weight = std::stod(field);

        // read the interaction channel code
        getline(line_sstream, field, ',');
        int code = std::stoi(field);
        // everything except coherent spp and neutrino-lepton scattering
        if((code < 0 || code > 5) && (code < 8 || code > 10) && (code < 20 || code > 21))
        {
          throw std::runtime_error("unknown interaction channel code at event " + std::to_string(i));
        }
        e->dyn = code;

        // read pdg + four-momentum per each particle
        std::vector<particle> event_particles;
        int lepton_idx = -1;
        while(std::getline(line_sstream, field, ','))
        {
          // read 
          int pdg = std::stoi(field);
          particle new_particle(pdg, PDG::mass(pdg));

          // note if it is the lepton
          if(new_particle.lepton())
          {
            if(lepton_idx < 0)
            {
              lepton_idx = event_particles.size();
            }
            else
            {
              throw std::runtime_error("more than one lepton in the final state at event "  + std::to_string(i));
            }
          }

          // read four-momentum
          std::vector<double> fourmomentum;
          for(int j=0; j<4; j++)
          {
            if(std::getline(line_sstream, field, ','))
              fourmomentum.push_back(std::stod(field));
            else
              throw std::runtime_error("events file format is corrupted at event " + std::to_string(i));
          }
          new_particle.set_momentum(vec(fourmomentum[1],fourmomentum[2],fourmomentum[3]));

          // check the correctness of the four-momentum
          if(std::fabs(fourmomentum[0] - new_particle.E()) > 1.e-1)
            throw std::runtime_error("particles specified incorrectly (off-shell?) at event "  + std::to_string(i));

          // push the particle to the particles vector
          event_particles.push_back(new_particle);
        }

        // make sure there was a lepton in the event after all
        if(lepton_idx < 0)
        {
          throw std::runtime_error("no lepton in the final state at event "  + std::to_string(i));
        }

        // generate the point of interaction
        vec position = nucl->get_random_r()*rand_dir();

        // isolate the lepton and put to out
        particle lepton = event_particles[lepton_idx];
        lepton.r = position;
        e->out.push_back(lepton);
        event_particles.erase(event_particles.begin() + lepton_idx);

        // put the hadrons to out
        int initial_charge = 0;
        for(int j = 0; j<event_particles.size(); j++)
        {
          initial_charge += event_particles[j].charge();
          event_particles[j].r = position;
          e->out.push_back(event_particles[j]);
        }

        // dummy incoming neutrino
        int neutrino_pdg;
        if(code < 20 && code % 2 == 0) // cc
        {
          if(lepton.pdg > 0)
          {
            neutrino_pdg = lepton.pdg+1;
            initial_charge--;
          }
          else
          {
            neutrino_pdg = lepton.pdg-1;
            initial_charge++;
          }
        }
        else // nc or el
        {
          neutrino_pdg = lepton.pdg;
        }
        particle neutrino(neutrino_pdg, PDG::mass(neutrino_pdg));
        e->in.push_back(neutrino);

        // dummy nucleon target(s)
        if(e->dyn < 8 || e->dyn > 9) // not MEC, so only one target nucleon
        {
          if(initial_charge > 1)
            throw std::runtime_error("illegal initial charge at event "  + std::to_string(i));
          particle nucleon;
          nucleon.r = position;
          if(initial_charge)
            nucleon.set_proton();
          else
            nucleon.set_neutron();
          e->in.push_back(nucleon);
        }
        else                         // MEC, we need two initial nucleons
        {
          if(initial_charge > 2)
            throw std::runtime_error("illegal initial charge at event "  + std::to_string(i));
          particle nucleon1, nucleon2;
          nucleon1.r = position;
          nucleon2.r = position;
          if(initial_charge == 2)
          {
            nucleon1.set_proton();
            nucleon2.set_proton();
          }
          else if(initial_charge == 1)
          {
            nucleon1.set_proton();
            nucleon2.set_neutron();
          }
          else
          {
            nucleon1.set_neutron();
            nucleon2.set_neutron();
          }
          e->in.push_back(nucleon1);
          e->in.push_back(nucleon2);
        }
      }
      catch(const std::exception& ex)
      {
        cout << "Exception: " << ex.what() << endl;
        return 1;
      }
    }
    // or use the standard kaskada functionality
    else
    {
      e->weight = 1;
      particle p0=b->shoot();
      p0.r=start_point(nucl,p);
      e->out.push_back(p0);
      particle pi = nucl->get_nucleon(p0.r);
      e->in.push_back(pi);
      e->in.push_back(pi);
    }

    // Run the cascade
    kaskada k(p,*e,&input);
    k.kaskadaevent();
    t2->Fill();
    delete e;
    raport(i+1,p.number_of_events,"cascade events ready...",1000);
  }
  f->Write();

  // close the file behind yourself
  if(p.kaskada_events)
    events_file.close();

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
