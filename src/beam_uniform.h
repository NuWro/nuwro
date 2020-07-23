#ifndef _beam_uniform_h_
#define _beam_uniform_h_
#include "EnergyProfile.h"
#include "params.h"
#include "particle.h"
#include "beam.h"


/// The beam shoots with identical particles (whose energy is defined by some 
/// energy profile but the direction of momentum is identical)
/// beam: this class provides the input of the initial beam particle 
/// into the generator.The beam particle is defined by providing 
/// its energy, direction and its PDG code to the constructor
/// of the beam class 

class beam_uniform : public beam
{
protected:
   EnergyProfile g;
   vec dir;
   int pdg;
   double mass;
public:   
   beam_uniform(NUWRO::params& p):g(p.beam_energy),dir(p.beam_direction.dir()),pdg(p.beam_particle), mass(PDG::mass(p.beam_particle))
   {
   }
   /// get next particle form the beam
   virtual particle shoot(bool dis=0)
   { particle p(pdg,mass);
     p.set_momentum(dir);
     p.set_energy(g.shoot(dis));
     return p;
   }

   double disratio()
   {
	   return g.disratio();
   }
   
  
   void check_energy()
   {
	   if (g.minE() < mass) throw "energy can't be lower than particle mass!";	   
   }
};
#endif // _beam_uniform_h_
