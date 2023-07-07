#ifndef _dis2res_h_
#define _dis2res_h_
//#include "delta.h"



void funkcja_spp_0 ();

double funkcja_spp_3 (double W, int nukleon_in, int pion, bool current);

double funkcja_spp_2 (double W, int nukleon_in, bool current);

int meson_out_ (double W, int lepton_in, int nukleon_in, bool current);

int nukleon_out_ (double W, int lepton_in, int nukleon_in, int meson_out,
		  bool current);
double funkcja_dis (double W, double W_min, double W_max, double alfa);


// KN: The following functions are commented out until anyone proves their usefulness.
//     Please do not leave garbage in the code.

// double kombinacja (int FFset, double E, double W, double nu, int lepton_in,
// 		   int nukleon_in, int nukleon_out, int meson_out,
// 		   bool current, double alfa, double W_min, double W_max);
// //nowa kombinacja dis i res z produkcja stanow koncowych

// double cr_dis_res (int FFset, double E, double W, double nu, int lepton_in,
// 		   int nukleon_in, int meson_out, int nukleon_out,
// 		   bool current);

#endif
