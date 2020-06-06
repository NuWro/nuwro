#ifndef _beam_singleroothist_h_
#define _beam_singleroothist_h_
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

class beam_singleroothist : public beam
{
protected:
   EnergyProfile g;
   vec dir;
   int pdg;
   double mass;
public:
 beam_singleroothist(NSNWRO::params& p):g(p.beam_inputroot,p.beam_inputroot_flux),dir(p.beam_direction.dir()),pdg(p.beam_particle), mass(PDG::mass(p.beam_particle))
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

   int GetPDG() const { return pdg; }
   
   EnergyProfile const & EProf() const { return g; }
};
#endif // _beam_singleroothist_h_
