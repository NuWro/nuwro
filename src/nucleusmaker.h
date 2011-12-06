#ifndef _nucleusmaker_h_
#define _nucleusmaker_h_
#include "nucleus.h"
#include "flatnucleus.h"
#include "anynucleus.h"
#include "params.h"
#include "generatormt.h"

inline nucleus * make_nucleus(params& p)
{ 
  switch(p.nucleus_model)
  {
    case 0: return new flatnucleus(p);
    case 1: return new anynucleus(p);
    default: throw "unknown nucleus model";
  }
}

inline vec start_point(nucleus *nucl,params &p)
{
  switch(p.beam_placement)
  {
    //  nucleus center
    case 0:  return vec(0,0,0);

    // transparency mode: interaction starts at random nucleon 
    case 1: return nucl->get_random_r()*rand_dir();

    // pion or nucleon scattering mode - the particle starts jusst
    // under the surface of the nucleus 
    case 2: 
   
	  double x, y, z;
	  do
	    {
	      x = 2 * frandom () - 1;
	      y = 2 * frandom () - 1;
	    }
	  while (x * x + y * y > 1);
	  
	  z = -sqrt (1 - x * x - y * y);
	  
	  return 0.999* vec (x, y, z) * nucl->radius ();
		
    // wrong parameter 
    default: throw "invalid beam_placement";
  }
} 


#endif
