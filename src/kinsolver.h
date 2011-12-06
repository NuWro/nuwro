#ifndef _kinsover_
#define _kinsover_
#include "particle.h"

class kinsolver
{private:
   double _E_lab;
   vec _cms_p_dir;
   vec _cmsspeed;
   particle  _nu;
   particle  _N0;
   particle & _lepton;
   particle & _N1;
    
 public:
 
   kinsolver(double E_lab, 
              vec cms_p_dir,
	      vec cmsspeed,                    
	      particle nu, particle N0,        // in particles
	      particle &lepton, particle &N1   // out particles
	   );
	      
  double dE(double pnew);
  double findmomentum();
};
#endif
