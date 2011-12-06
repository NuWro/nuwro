#ifndef _flatnucleus_h_
#define _flatnucleus_h_

#include "nucleus.h"
#include "particle.h"
#include "pdg.h"
#include "params.h"
#include "generatormt.h"
#include <algorithm>

using namespace std;

/// this class implements the target nucleus
/// it should take into account of the its radius, charge ,
/// proton and neutron density as well as local fermi momentum and Pauli blocking
/// the number of nucleons can change as a result of scattering and the density 
/// and Fermi momentum should be consistently updated - this is still NOT implemented

struct flatnucleus :public nucleus
{
private:

  double Volume;

public:

  flatnucleus (params & par);                 ///< Construct nucleus from params given in the input file 
  double radius () {return _r;}                ///< Radius  
  double V(){return _V;}                      ///< Potential
  double density (double r);                  ///< calculates nuclear density of proton. neutrons at radius r from the center
  double get_random_r ();                     ///< random distance from the center waighted with nuclear density       

};
////////////////////////////////////////////////////////////////////////
//            Implementation
////////////////////////////////////////////////////////////////////////

inline  flatnucleus::flatnucleus (params & par): nucleus(par)
  {
    using namespace PDG;
    if (p + n == 1)
      Eb = kf = _r = local_kf = 0;
    _V = sqrt (mass_proton * mass_proton + kf * kf) - mass_proton + 7 * MeV;
    
    double A = n + p;
    double Ap = pow (A, 1.0 / 3);
      _r = Ap *1.2* fermi;		//????
    Volume = 4.0 / 3 * M_PI * _r * _r * _r;
  }


///////////////////////////////////////////////////////////////////////////////
/// calculates nuclear density of proton. neutrons at radius r from the center
///////////////////////////////////////////////////////////////////////////////
inline  double flatnucleus::density (double rx )
  { if(rx>radius()) return 0;
    return (pr + nr) / Volume;
  }

///////////////////////////////////////////////////////////////////////////////
// random distance from the center waighted with nuclear density
///////////////////////////////////////////////////////////////////////////////
inline  double flatnucleus::get_random_r ()
  { double r1= _r * frandom_sqr();	// rozklad kwadratowy
    return  r1;
  }


#endif
