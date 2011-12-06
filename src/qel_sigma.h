/////////////////////////////////////////////////////////////////////////////////////
// caa informacja o kwazielastycznym rozpraszaniu neutrin na nukleonach 
// na podstawie Llewelyn-Smith 
// (C) C. Juszczak 2001
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _QEL_h_
#define _QEL_h_

#include "jednostki.h"
#include "masses.h"
#include "calga.h"
#include "params.h"
#include "ff.h"

//extern double qelmlep;



/// semielastic neutrino nucleon scattering 
class QEL
{ 
  bool anty;
  double Enu;          //Eneria neutrina w układzie spoczywającej tarczy  
  bool use_lambda;    // ??????
  double m, mm, M, MM; // lepton mass and hadrom mass (sqared)
  
public:
   QEL(params& p);
   QEL(double lepton_mass0, 
	   double hadron_mass = M12,
	   bool anty0=0,
	   bool use_lambda2 = true
	   );
// 0 - cc //  1 - nc proton // 2- nc neutron
   double sigma (double q2, int kind);            // cross section formula
   double qq (double factor);                  // granice cakowania
//   double total_sigma(double E);
   void set_energy(double E){Enu=E;}

}; //QEL


 
 
#endif
