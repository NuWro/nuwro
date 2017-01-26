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

using namespace std;
using namespace PDG;

#define LOG(x) cout<< x <<endl;
#define LOG2(x,y) cout<< x << y <<endl;
#define ERR(x) cerr<< x <<endl;
#define ERR2(x,y) cerr<< x << y <<endl;

//! Semi-classical cascade model of the final state interactions.
/*! The model runs for nucleons and pions from the primary vertex. The probability
		of interaction between the particle and the residual nucleus is integrated using
		a Monte Carlo method. In each step, the free path is calculated (for a given
		density and total cross section). If the distance is lower than the maximal step,
		it propagates and interacts. If not, there is no interaction and the particle
		is propagated with the maximal step. All of the interaction products are also put
		in the cascade. */

class kaskada
{
		queue < particle > parts;											//!< Queue for the particles in the cascade.
		params par;																		//!< Params of the simulation.
		event *e;																			//!< Current event.
		nucleus *nucl;																//!< Nucleus for the use of the cascade.
		particle *p;																	//!< Current particle in the cascade.
		Interaction *I;																//!< Interaction in the cascade.
		interaction_parameters X;											//!< Keeps the information from prepare to make interaction.
		double max_step;															//!< Maximal step of propagation.
		double radius;																//!< Radius of the nucleus.

	public:
		kaskada(params &p, event &e1);								//!< The default constructor.
																									/*!< Takes the params file and the current event.
																											 Generates a new nucleus for the cascade. */
		~kaskada();																		//!< The default destructor.
		int kaskadaevent();														//!< Runs the cascade.

	private:
		void prepare_particles(); 										//!< Handles the particles from the input (out) vector.
																									/*!< Nucleons and pions are prepared and added to the queue as off-shell particles.
																											 Other particles are copied directly to the output vector (post) */
		interaction_parameters prepare_interaction();	//!< Calculates the free path.
																									/*!< The free path depends on the density and the total cross section. 
																											 The density is set for a current position.
																											 The cross section is set according to the kinetic energy and particle type. */
		bool move_particle();													//!< Propagates the particle.
																									/*!< Particle is propagated by no more than max_step.
																											 If its kinetic energy is lower than binding it remains jailed. */
		bool leave_nucleus();													//!< Handles the particles propagated outside the nucleus.
																									/*!< If the particle is not jailed: it escapes, returns on-shell
																											 and is added to the output vector (post). */
		bool make_interaction(); 											//!< Generates kinematics.
																									/*!< The interaction is rejected if the chosen kinematics violates Pauli blocking. */
		bool finalize_interaction(); 									//!< Copies new particles to a queue.
																									/*!< If on, the formation zone is applied. */
		void clean ();																//!< Cleans after the cascade.
																									/*!< Clears the queue and remembers the residual nucleus. */

		bool check  (particle & p1, particle & p2,
								 particle *spect, int n, particle p[], int k);	//!< Checks if the charge is conserved
		bool check2 (particle & p1, particle & p2,
								 particle *spect, int n, particle p[], int k);	//!< Checks if the fourmomentum is conserved
		
		//void procinfo(particle p1, particle p2, int n, particle p[]); // (empty function?)
};
