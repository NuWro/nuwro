#include "event1.h"
#include "params.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "generatormt.h"
#include "vect.h"
#include <queue>
#include "pdg.h"
#include "nucleus.h"
#include <TROOT.h>
#include <TTree.h>
#include "beam.h"
//#include "Metropolis.h"
#include "Interaction.h"
#include "proctable.h"
#include "nucleusmaker.h"
#include "fsi.h"

using namespace std;
using namespace PDG;

#define LOG(x) cout<< x <<endl;
#define LOG2(x,y) cout<< x << y <<endl;
#define ERR(x) cerr<< x <<endl;
#define ERR2(x,y) cerr<< x << y <<endl;

class kaskada
{
	private:
		queue < particle > parts;
		params par;
		event *e;
		nucleus *nucl;
		particle *p;
		Interaction *I;
		interaction_parameters X;
		double max_step;
		double radius;

		void prepare_particles(); 						///take all nucleons and pions from input vector, others copy directly to the output vector
		interaction_parameters prepare_interaction();	///set density and cross sections
		bool make_interaction(); 					    ///generate a kinematic
		bool finalize_interaction(); 					///copy new particle to a queue, apply formation zone if it is on 
		bool leave_nucleus(); 							///move particle to the output particles vector and reduce momentum for nucleons
		bool move_particle();				///propagate particle and change nucleons energy according to nuclear potential

		bool check (particle & p1, particle & p2, particle *spect, int n, particle p[],int k);  ///check if charge is conserved
		bool check2 (particle & p1, particle & p2, particle *spect, int n, particle p[],int k); ///check if fourmomentum is conserved
		
		void clean ();
		
		//void procinfo(particle p1, particle p2, int n, particle p[]); /// (empty function?)

	public:
		kaskada(params &p, event &e1);
		~kaskada();
		int kaskadaevent();
};
