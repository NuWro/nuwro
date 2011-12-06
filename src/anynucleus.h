#ifndef _anynucleus_h_
#define _anynucleus_h_

#include "nucleus.h"
#include "params.h"

//______________________________________________________________________
/// this class implements the target nucleus
/// it should take into account of the its radius, charge ,
/// proton and neutron density as well as local fermi momentum and Pauli blocking
/// the number of nucleons can change as a result of scattering and the density 
/// and Fermi momentum should be consistently updated - this is still NOT implemented

/// nucleus concrete class 
/// very accurate density profiles implemented
/// density dependent local Fermi momentum
/// shrinks as nucleons escape 
//______________________________________________________________________

class anynucleus:public nucleus
{
    
private:

  double calculate_radius();                   ///< needed to initialize r;

public:

  anynucleus (params & par);                   ///< Construct nucleus from params given in the input file 
  double V ()  { return _V;  }                 ///< Potential
  double radius ()  {   return _r;  }          ///< Radius
  double density (double r);                   ///< calculates nuclear density of proton. neutrons at radius r from the center
  double get_random_r ();                      ///< random distance from the center weighted with nuclear density

};
#endif
