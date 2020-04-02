#ifndef _resevent_hybrid_h_
#define _resevent_hybrid_h_

#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"

// Cuts in Q2 and W for the hybrid model
const double  Q2max_hybrid  =3.901*GeV2;
const double  Q2min_hybrid  =0.001*GeV2;
const double  Q2spc_hybrid  =0.100*GeV2;
const int     Q2bin_hybrid  =(Q2max_hybrid-Q2min_hybrid)/Q2spc_hybrid+1;
const double  Q2max_2_hybrid=0.201*GeV2;
const double  Q2min_2_hybrid=0.001*GeV2;
const double  Q2spc_2_hybrid=0.010*GeV2;
const int     Q2bin_2_hybrid=(Q2max_2_hybrid-Q2min_2_hybrid)/Q2spc_2_hybrid+1;
const double   Wmax_hybrid  =1500*MeV;
const double   Wmin_hybrid  =1075*MeV;
const double   Wspc_hybrid  =   5*MeV;
const int      Wbin_hybrid  =(Wmax_hybrid-Wmin_hybrid)/Wspc_hybrid+1;
const double   Wmax_2_hybrid=1100*MeV;
const double   Wmin_2_hybrid=1075*MeV;
const double   Wspc_2_hybrid=   1*MeV;
const int      Wbin_2_hybrid=(Wmax_2_hybrid-Wmin_2_hybrid)/Wspc_2_hybrid+1;
const double cthmax_hybrid  =     1;
const double cthmin_hybrid  =    -1;
const double cthspc_hybrid  = 2./19;
const int    cthbin_hybrid  =(cthmax_hybrid-cthmin_hybrid)/cthspc_hybrid+1;

// Generate RES event
void resevent_hybrid(params& p, event& e, bool cc);

// Double-differential cross section in CMS without LAB factors
// - in terms of tabularized hadronic responses
double hybrid_dsdQ2dW_tab(res_kinematics* kin, int params[4], vect final_pion=NULL);

// Three-fold differential cross section in CMS without LAB factors
// - in terms of tabularized hadronic responses
double hybrid_dsdQ2dWdcth_tab(res_kinematics* kin, int params[4], vect final_pion);
// - in terms of the ABCDE decomposition
double hybrid_dsdQ2dWdcth    (res_kinematics* kin, int params[4], vect final_pion);

// Four-fold differential cross section in CMS without LAB factors
//  - in terms of the ABCDE decomposition
double hybrid_dsdQ2dWdOm(res_kinematics* kin, int params[4], vect final_pion);

// Sample cos_theta_pi^* using ABCDE decomposition
// - using ABCDE decomposition, fitting a polynomial, inverting using bisection
double hybrid_sample_costh  (double Enu, double Q2, double W, int params[4]);

// Sample phi_theta^* using ABCDE decomposition
// - fitting a polynomial, inverting using bisection
double hybrid_sample_phi  (double Enu, double Q2, double W, int params[4], double costh_rnd);
// - integrating analytically, inverting using bisection
double hybrid_sample_phi_2(double Enu, double Q2, double W, int params[4], double costh_rnd);
// - integrating analytically, inverting using Newton's method
double hybrid_sample_phi_3(double Enu, double Q2, double W, int params[4], double costh_rnd);

// Fit polynomial to given points
void hybrid_poly_fit(const int N, double* xpts, double* ypts, double* coeffs);

// Cumulative distribuant of a polynomial
double hybrid_poly_dist(const int N, double* coeffs, double x_min, double x);

// Random variable from plynomial using bisection
double hybrid_poly_rnd(const int N, double* coeffs, double x_min, double x_max, double epsilon);

// Random variable directly from ABCDE decomposition
// - using bisection
double hybrid_dcmp_rnd  (double (*ABCDE)[5], double x_min, double x_max, double epsilon);
// - using Newton's method
double hybrid_dcmp_rnd_2(double (*ABCDE)[5], double x_min, double x_max, double epsilon);

// Modify final hadron directions as specified in the Adler frame (k, kp in CMS)
vec hybrid_dir_from_adler(double costh, double phi, vect k, vect kp);

// Function to be called in finishevent to resample the hadronic angles
// - resample the whole direction
void resevent_dir_hybrid(event& e);
// - resample only the phi angle
void resevent_phi_hybrid(event& e);

// Interpolations, axes xyz
double    linear_interp(double f0,   double f1,   double xd);
double  bilinear_interp(double f00,  double f10,  double f01,  double f11,  double xd, double yd);
double trilinear_interp(double f000, double f010, double f100, double f110,
                        double f001, double f011, double f101, double f111, double xd, double yd, double zd);

#endif