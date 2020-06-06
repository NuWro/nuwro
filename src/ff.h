#ifndef _FF_H_
#define _FF_H_
#include <utility>
#include "params.h"

using namespace std;
//_______________________________________________________________
//
// Form Factor Library Access
//_______________________________________________________________

/// Form Factor configuration
void ff_configure(NSNWRO::params& p);

/// Calculate F1, F2 form factors
/// usage: list(F1,F2)=f12(q2,  kind);
/// q2:   fourmomentum transfer squared
/// kind: 0 - cc, 
///       1 - nc on proton, 
///       2 - nc on neutron
pair <double,double> f12(double q2, int kind);

///Calculate Fa, Fp form factors
/// usage: list(Fa,Fp)=fap(q2,  kind);
/// q2:   fourmomentum transfer squared
/// kind: 0 - cc, 
///       1 - nc on proton, 
///       2 - nc on neutron
pair <double,double> fap(double q2, int kind);


/// Pair Assigment helper structure
struct list
{ double &a;
  double &b;	
  list(double&a0,double&b0):a(a0),b(b0){}
  void operator=(pair<double,double> ab)
        {a=ab.first;b=ab.second;}
};

#endif
