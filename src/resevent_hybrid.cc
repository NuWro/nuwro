/*
This function calculates RES events from the Hybrid model
*/

#include "resevent_hybrid.h"
#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"
#include "particle.h"
#include "vect.h"
#include "dis/LeptonMass.h"
#include "hybrid_RES.h"
#include "hybrid/hybrid_gateway.h"
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TVectorD.h"
#include "TH1D.h"


void resevent_hybrid(params &p, event &e, bool cc) {      // free nucleon only!
  e.weight = 0;           // if kinematically forbidden

  res_kinematics kin(e);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(1500)) return;

  const double our_W_val  = 1230;
  const double our_W_wid  =    1;
  if( fabs(kin.W-our_W_val) > our_W_wid ) return;

  const double our_Q2_val = 100000;
  const double our_Q2_wid =    100;
  if( fabs((-kin.q*kin.q)-our_Q2_val) > our_Q2_wid ) return;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  // final state particles
  particle final_pion, final_nucleon;
  particle final_pion2, final_nucleon2; // needed for choosing a decay channel

  // selection of the final state

  // specify the total electric charge of the pion-nucleon system
  int j = kin.neutrino.pdg < 0;
  int k = not cc;
  int l = kin.target.pdg != PDG::pdg_proton;
  int final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  double xsec_pip=0, xsec_pi0=0, xsec_pim=0, xsec_inclusive=0; // cross sections

  // choose a random direction in CMS
  vec kierunek = rand_dir();

  // specify the params needed for ABCDE (note strange order!)
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (kin.neutrino.pdg > 0));  // helicity
  // params[3] is the target nucleon, params[1] is the decay channel

  // cross section function
  double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dW;
  //double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dWdcth;
  //double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dWdOm;

  switch (final_charge) { // calculate the cross section with only the CMS variables
    case  2:  // pi+ + proton (nu_11)
      {params[3] = 1; params[1] = 1;
       final_pion.set_pdg_and_mass( PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pip = hybrid_xsec(&kin, params, final_pion);}
      xsec_inclusive = xsec_pip;
      if ( not (xsec_inclusive > 0) ) return;
      break;
    case  1:  // pi+ + neutron (nu_22) or pi0 + proton (nu_21)
      {params[3] = 2; params[1] = 2;
       final_pion.set_pdg_and_mass( PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pip = hybrid_xsec(&kin, params, final_pion);}
      {params[3] = 2; params[1] = 1;
       final_pion2.set_pdg_and_mass( PDG::pdg_pi ); final_nucleon2.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg, final_nucleon2, final_pion2, kierunek);
       xsec_pi0 = hybrid_xsec(&kin, params, final_pion2);}
      xsec_inclusive = xsec_pip + xsec_pi0;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pip / xsec_inclusive < frandom() ) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 2; params[1] = 1;
      }
      else // make sure the params are okay for "1"
      {
        params[3] = 2; params[1] = 2;
      }
      break;
    case  0:  // pi0 + neutron (anu_11) or pi- + proton (anu_12)
      {params[3] = 1; params[1] = 1;
       final_pion.set_pdg_and_mass( PDG::pdg_pi ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pi0 = hybrid_xsec(&kin, params, final_pion);}
      {params[3] = 1; params[1] = 2;
       final_pion.set_pdg_and_mass( -PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg, final_nucleon2, final_pion2, kierunek);
       xsec_pim = hybrid_xsec(&kin, params, final_pion2);}
      xsec_inclusive = xsec_pi0 + xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pi0 / xsec_inclusive < frandom() ) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 1; params[1] = 2;
      }
      else // make sure the params are okay for "1"
      {
        params[3] = 1; params[1] = 1;
      }
      break;
    case -1:  // pi- + neutron (anu_22)
      {params[3] = 2; params[1] = 2;
       final_pion.set_pdg_and_mass( -PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pim = hybrid_xsec(&kin, params, final_pion);}
      xsec_inclusive = xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      break;
    default:
      cerr << "[WARNING]: Reaction charge out of range\n";
  };

  // Omega_pi^* was chosen in hadronic CMS

  // // Choose cos_theta^* for dsdQ2dW, in Adler frame
  // double costh_rnd = hybrid_sample_costh(kin.neutrino.E(), -kin.q*kin.q, kin.W, params);

  // // Choose phi^* for dsdQ2dWdcosth, in Adler frame
  // double phi_rnd = hybrid_sample_phi(kin.neutrino.E(), -kin.q*kin.q, kin.W, params, costh_rnd);

  // // Modify final hadron directions as specified in the Adler frame
  // vect nu  = kin.neutrino; nu.boost(-kin.hadron_speed); 
  // vect lep = kin.lepton;  lep.boost(-kin.hadron_speed);
  // vec dir_rnd = hybrid_dir_from_adler(costh_rnd, phi_rnd, nu, lep);

  // // Recalculate the hadronic kinematics, dir_rnd is the new direction of pion
  // double momentum = final_nucleon.momentum();
  // final_nucleon = vect(final_nucleon.E(), -momentum * dir_rnd.x, -momentum * dir_rnd.y, -momentum * dir_rnd.z);
  // final_pion = vect(final_pion.E(), momentum * dir_rnd.x, momentum * dir_rnd.y, momentum * dir_rnd.z);

  // set event weight
  e.weight = xsec_inclusive;

  // the cross section needs a factor (initial lepton momentum)^-2 in LAB
  // the cross section needs a jacobian: dw = dQ2/2M
  e.weight /= e.in[0].E() * e.in[0].E() / 2 / res_kinematics::avg_nucleon_mass / kin.jacobian;
  // coupling
  e.weight *= G*G*cos2thetac/2;
  // units
  e.weight /= cm2;

  // boost back to LAB frame
  final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
  final_nucleon.p4() = final_nucleon.boost(kin.target.v());

  final_pion.p4() = final_pion.boost(kin.hadron_speed);
  final_pion.p4() = final_pion.boost(kin.target.v());

  // save final state hadrons
  // warning: the order is essential
  e.out.push_back(final_pion);
  e.out.push_back(final_nucleon);

  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
}

double hybrid_dsdQ2dW(res_kinematics *kin, int params[4], vect final_pion)
{
  double *hybrid_grid;

  if(kin->neutrino.pdg > 0)
  {
    switch (params[3]*10 + params[1])
    {
      case 11:
        hybrid_grid = hybrid_grid_11;
        break;

      case 21:
        hybrid_grid = hybrid_grid_21;
        break;

      case 22:
        hybrid_grid = hybrid_grid_22;
        break;

      default:
        cerr << "[WARNING]: no hybrid inclusive grid!\n";
    }
  }
  else
  {
    switch (params[3]*10 + params[1])
    {
      case 11:  // anu_11 has the tables of nu_21
        hybrid_grid = hybrid_grid_21;
        break;

      case 12:  // anu_12 has the tables of nu_22
        hybrid_grid = hybrid_grid_22;
        break;

      case 22:  // anu_22 has the tables of nu_11
        hybrid_grid = hybrid_grid_11;
        break;

      default:
        cerr << "[WARNING]: no hybrid inclusive grid!\n";
    }
  }

  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // we build leptonic tensor in CMS with q along z, see JES paper App. A
  vect kl_inc_lab = kin->neutrino.p4();   //
  vect kl_lab     = kin->lepton.p4();     //  They are not in LAB, but in target rest frame!
  vect q_lab      = kin->q;               //

  double v = q_lab.length() / (q_lab[0] + res_kinematics::avg_nucleon_mass);
  double g = 1 / sqrt(1 - v*v);
  double c = (kl_inc_lab.length() - kl_lab[3]) / q_lab.length();
  double s = sqrt(pow(kl_lab.length(),2) - pow(kl_lab[3],2)) / q_lab.length();

  vect kl_inc (g*(kl_inc_lab[0] - v*kl_inc_lab.length()*c), kl_inc_lab.length()*s,
               0, g*(-v*kl_inc_lab[0] + kl_inc_lab.length()*c));
  vect kl     (g*(kl_lab[0] - v*(kl_inc_lab.length()*c - q_lab.length())), kl_inc_lab.length()*s,
               0, g*(-v*kl_lab[0] + kl_inc_lab.length()*c - q_lab.length()));

  double kl_inc_dot_kl = kl_inc * kl;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // interpolating the nuclear tensor elements
  double w00=0;
  double w03=0;
  double w33=0;
  double w11=0;
  double w12=0;

  //if( Q2 < Q2min_hybrid ) Q2 = Q2min_hybrid;
  //if( Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  if( Q2 >= Q2min_hybrid && Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  {
    int m=int((Q2-Q2min_hybrid)/Q2spc_hybrid);
    int n=int((W-Wmin_hybrid)/Wspc_hybrid);
    int pos=(n*Q2bin_hybrid+m)*5;

    int a=pos;
    int b=a+5;
    int c=pos+Q2bin_hybrid*5;
    int d=c+5;

    double H1=W-Wmin_hybrid-n*Wspc_hybrid;
    double H2=Q2-Q2min_hybrid-m*Q2spc_hybrid;

    w00=(H2*(H1*hybrid_grid[d]+(Wspc_hybrid-H1)*hybrid_grid[b])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c]+(Wspc_hybrid-H1)*hybrid_grid[a])/Wspc_hybrid)/Q2spc_hybrid;
    w03=(H2*(H1*hybrid_grid[d+1]+(Wspc_hybrid-H1)*hybrid_grid[b+1])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+1]+(Wspc_hybrid-H1)*hybrid_grid[a+1])/Wspc_hybrid)/Q2spc_hybrid;
    w33=(H2*(H1*hybrid_grid[d+2]+(Wspc_hybrid-H1)*hybrid_grid[b+2])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+2]+(Wspc_hybrid-H1)*hybrid_grid[a+2])/Wspc_hybrid)/Q2spc_hybrid;
    w11=(H2*(H1*hybrid_grid[d+3]+(Wspc_hybrid-H1)*hybrid_grid[b+3])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+3]+(Wspc_hybrid-H1)*hybrid_grid[a+3])/Wspc_hybrid)/Q2spc_hybrid;
    w12=(H2*(H1*hybrid_grid[d+4]+(Wspc_hybrid-H1)*hybrid_grid[b+4])/Wspc_hybrid+(Q2spc_hybrid-H2)*(H1*hybrid_grid[c+4]+(Wspc_hybrid-H1)*hybrid_grid[a+4])/Wspc_hybrid)/Q2spc_hybrid;
  }

  double l00 = (2.*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  double l03 = (-kl_inc[0]*kl[3] + -kl[0]*kl_inc[3]);
  double l33 = (2*kl_inc[3]*kl[3] + kl_inc_dot_kl);
  double l11 = (2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  double l22 = (2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  double l12 = (kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);

  result = l00*w00 + 2*l03*w03 + l33*w33 + 0.5*(l11+l22)*w11 - 2*l12*w12;

  return result;
}

double hybrid_dsdQ2dWdcth(res_kinematics* kin, int params[4], vect final_pion)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result = 2*Pi*ABCDE[0][0];
  result *= pion_momentum / pow(2*Pi,4);
  result *= 2; // Phase space!

  return result;
}

double hybrid_dsdQ2dWdOm(res_kinematics* kin, int params[4], vect final_pion)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result = ABCDE[0][0] + ABCDE[0][1]*cos(pion_phi) + ABCDE[0][2]*cos(2*pion_phi)
                       + ABCDE[0][3]*sin(pion_phi) + ABCDE[0][4]*sin(2*pion_phi);
  result *= pion_momentum / pow(2*Pi,4);
  result *= 4 * Pi; // Phase space!

  return result;
}

double hybrid_sample_costh(double Enu, double Q2, double W, int params[4])
{
  // Choosing cos_th^* from dsdQ2dWdcosth
  double costh_rnd;

  // Specify points for the polynomial interpolation
  const int costh_pts = 3;             // 3, 5, 7, 9, ...
  double costh[costh_pts];             // Points of interpolation
  for( int i = 0; i < costh_pts; i++ ) // Fill costh with evenly spaced points
    costh[i] = -cos( Pi / (costh_pts-1) * i );

  // Get the A function, \propto dsdQ2dWdcosth
  double ABCDE[costh_pts][5]; double ds[costh_pts];
  hybrid_ABCDE(Enu, Q2, W, costh, costh_pts, params, ABCDE);
  for( int i = 0; i < costh_pts; i++ )
    ds[i] = ABCDE[i][0];

  // Fit a polynomial to given number of points in ds(costh)
  double poly_coeffs[costh_pts];
  hybrid_poly_fit(costh_pts, costh, ds, poly_coeffs);

  // Normalize the coefficients to obtain a probability density
  double norm = hybrid_poly_dist(costh_pts, poly_coeffs, costh[0], costh[costh_pts-1]);
  for( int i = 0; i < costh_pts; i++ )
    poly_coeffs[i] /= norm;

  // Choose a value for costh_rnd
  costh_rnd = hybrid_poly_rnd(costh_pts, poly_coeffs, costh[0], costh[costh_pts-1], 0.001);

  return costh_rnd;
}

double hybrid_sample_phi(double Enu, double Q2, double W, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Specify points for the polynomial interpolation
  const int phi_pts = 7;             // 3, 5, 7, 9, ...
  double phi[phi_pts];               // Points of interpolation
  for( int i = 0; i < phi_pts; i++ ) // Fill phi with evenly spaced points
    phi[i] = 2*Pi / (phi_pts-1) * i - Pi;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5]; double ds[phi_pts];
  hybrid_ABCDE(Enu, Q2, W, &costh_rnd, 1, params, ABCDE);
  for( int i = 0; i < phi_pts; i++ )
    ds[i] = ABCDE[0][0] + ABCDE[0][1]*cos(phi[i]) + ABCDE[0][2]*cos(2*phi[i])
                        + ABCDE[0][3]*sin(phi[i]) + ABCDE[0][4]*sin(2*phi[i]);

  // Fit a polynomial to given number of points in ds(phi)
  double poly_coeffs[phi_pts];
  hybrid_poly_fit(phi_pts, phi, ds, poly_coeffs);

  // Normalize the coefficients to obtain a probability density
  double norm = hybrid_poly_dist(phi_pts, poly_coeffs, phi[0], phi[phi_pts-1]);
  for( int i = 0; i < phi_pts; i++ )
    poly_coeffs[i] /= norm;

  // Choose a value for phi_rnd
  phi_rnd = hybrid_poly_rnd(phi_pts, poly_coeffs, phi[0], phi[phi_pts-1], 0.001);

  return phi_rnd;
}

void hybrid_poly_fit(const int N, double* xpts, double* ypts, double* coeffs)
{
  // We perform a polynomial interpolation
  // Inspired by https://en.wikipedia.org/wiki/Polynomial_interpolation

  // Create matrix of x powers
  TMatrixD A(N,N);
  for( int i = 0; i < N; i++ ) // rows
  {
    for( int j = 0; j < N; j++ ) // columns backwards
    {
      A[i][j] = pow(xpts[i],(N-1-j));
    }
  }

  // Perform an LU decomposition
  TDecompLU LU(A);

  // Create a vector of y points
  TVectorD b(N);
  for( int i = 0; i < N; i++ )
    b[i] = ypts[i];

  // Solve the equation
  LU.Solve(b);
  for( int i = 0; i < N; i ++ )
    coeffs[i] = b[i];
}

double hybrid_poly_dist(const int N, double* coeffs, double x_min, double x)
{
  double result = 0.;
  for( int i = 0; i < N; i++ )
  {
    result += coeffs[i]/(N-i) * (pow(x,N-i)-pow(x_min,N-i));
  }
  return result;
}

double hybrid_poly_rnd(const int N, double* coeffs, double x_min, double x_max, double epsilon)
{
  double a = x_min, b = x_max;
  double x, fx;
  double y = frandom();

  do
  {
    x = (b-a)/2 + a;
    fx = hybrid_poly_dist(N, coeffs, x_min, x);
    if( fx > y )
      b = x;
    else
      a = x;
  }
  while( fabs(fx - y) > epsilon );

  return x;
}

vec hybrid_dir_from_adler(double costh, double phi, vect k, vect kp)
{
  vec Zast = k - kp;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,Zast);
  Zast.normalize(); Yast.normalize(); Xast.normalize();

  vec kierunek = costh * Zast + sqrt(1 - costh*costh) * (cos(phi) * Xast + sin(phi) * Yast);

  return kierunek;
}

void resevent_dir_hybrid(event& e)
{
  // get all 4-vectors from the event and boost them to the N-rest frame
  vect target   = e.in[1];  
  target.t     -= get_binding_energy (e.par, target);
  vect neutrino = e.in[0];  neutrino.boost(-target.v());
  vect lepton   = e.out[0]; lepton.boost  (-target.v());
  vect pion     = e.out[1]; pion.boost    (-target.v());
  vect nucleon  = e.out[2]; nucleon.boost (-target.v());

  // calculate hadron_speed and boost to CMS
  vect q = neutrino - lepton;
  double Mef = min(sqrt(target * target), res_kinematics::avg_nucleon_mass);
  double q3  = sqrt(pow(Mef + q.t,2) - e.W()*e.W());
  vec hadron_speed = q / sqrt(e.W()*e.W() + q3 * q3);
  neutrino.boost (-hadron_speed);
  lepton.boost   (-hadron_speed);
  pion.boost     (-hadron_speed);
  nucleon.boost  (-hadron_speed);

  // generate params for the Ghent code
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (e.out[0].pdg > 0));
  if( e.in[0].pdg > 0 ) // neutrino
  {
    if( e.in[1].pdg == PDG::pdg_proton )
    {
      params[3] = 1; params[1] = 1;
    }
    else
    {
      params[3] = 2;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }
  else                   // antineutrino
  {
    if( e.in[1].pdg == PDG::pdg_neutron )
    {
      params[3] = 2; params[1] = 2;
    }
    else
    {
      params[3] = 1;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }

  // Choose cos_theta^* for dsdQ2dW, in Adler frame
  double costh_rnd = hybrid_sample_costh(neutrino.t, -e.q2(), e.W(), params);

  // Choose phi^* for dsdQ2dWdcosth, in Adler frame
  double phi_rnd = hybrid_sample_phi(neutrino.t, -e.q2(), e.W(), params, costh_rnd);

  // Modify final hadron directions as specified in the Adler frame
  vec dir_rnd = hybrid_dir_from_adler(costh_rnd, phi_rnd, neutrino, lepton);

  // Recalculate the hadronic kinematics, dir_rnd is the new direction of pion
  double momentum = nucleon.length();
  nucleon = vect(nucleon.t, -momentum * dir_rnd.x, -momentum * dir_rnd.y, -momentum * dir_rnd.z);
  pion = vect(pion.t, momentum * dir_rnd.x, momentum * dir_rnd.y, momentum * dir_rnd.z);

  // Boost hadrons back
  nucleon.boost(hadron_speed);
  nucleon.boost(target.v());
  pion.boost   (hadron_speed);
  pion.boost   (target.v());

  // Correct the particles in the out vector
  e.out[1].p4() = pion;
  e.out[2].p4() = nucleon;
}
