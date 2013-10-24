
#ifndef __BEAM_H_
#define __BEAM_H_

#include "particle.h"
#include "params.h"

/// The beam shoots with identical particles (whose energy is defined by some 
/// energy profile but the direction of momentum is identical)
/// beam: this class provides the input of the initial beam particle 
/// into the generator.The beam particle is defined by providing 
/// its energy, direction and its PDG code to the constructor
/// of the beam class 

/// beam interface
class beam
{
public:
	virtual particle shoot(bool dis=0) = 0;
	virtual double nu_per_POT(){ return 0;}
	virtual ~beam(){}
};


/// create new beam object
beam * create_beam(params &p);


#endif	/// __BEAM_H_
