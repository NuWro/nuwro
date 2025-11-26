// Vector and Axial Form Factor Calculations for lepton-Nucleus Interactions
// This code implements various models for calculating vector and axial vector form factors,
#include "ff.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include "jednostki.h"
#include "params.h"
#include "pdg.h"
#include "event1.h"
#include "rew/rewparams.h"
#include <vector>
#include "event1.h"
#include "kinematics.h"
#include "nucleus.h"
#include "sfevent.h"
#include "sf/GConstants.h"
using namespace std;
constexpr inline double pow2(double x) { return x * x; }

static double mva_errorBar = 0.0;//new JS
static double deut_errorBar = 0.0;//new JS
// strange =0 nie strange FF  strange =1 old implementation (recover old bahaviour)  strange =2 new implementation (uses strange axial mass != nc axial mass)
struct FF
{
  double Q2;
  double GEp, GEn;
  double GMp, GMn;
  FF() : Q2(0), GEp(0), GEn(0), GMp(0), GMn(0) {}
  inline pair<double, double> f12(int kind);
};
//_________________________________________________________
/// Pointer to current vector form factors model
static FF (*FFfromq2)(const double) = 0;

/// Pointer to current axial form factors model
static double (*Axialfromq2)(const double, const double) = 0;
// _________________________________________________________
double axialcorr(int axialFF, double q2);

/// Functions calculating form factors
FF DipoleFF(const double q2);   // 1. dipole electric form factor G_E^V
FF bba03_FF(const double q2);   // 2. hep-ex/0308005 BBA-2003 for Q2<6 GeV
FF bbba05_FF(const double q2);  // 3. arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV
FF bbba07_FF(const double q2);
double GetA(const double q2, const double* x);
FF JLab_FF(const double q2);  // 4. PHYSICAL REVIEW C, VOLUME 65, 051001(R)
                              // PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
FF kg_FF(const double q2);    // 5. K. Graczyk ...
FF npar_FF(const double q2);  // 6. nowa (1990:) parametryzacja JS z qelcc
/// FA Functions
double dipole_FA(const double q2, const double ma); // Standard dipole
double deuterium_FA(const double qSq); //axial form factor from neutrino-deuteron scattering data, Meyer et al., Phys. Rev. D 93, 113015 (2016)
double bbba07_FA(const double q2, const double ma); // BBBA07 Form Factors
double comp2_FA(const double q2, const double ma); // 2 Component Model
double comp3_FA(const double q2, const double ma); // 3 Component Model
double zexp_FA(const double q2, const double ma); // Z-expansion Model
double MINERvA_FA(const double q2, const double sigma ); // axial form factor determined by MINERvA, Nature 614, 48 (2023)//changed JS
double LQCD_FA(const double q2, const double bandSwitch); // Average LQCD form factor obtained from A. Meyer via private communication
// IMPLEMENTATION
/// Calculate F1,F2
pair<double, double> FF::f12(int kind)
{
  double Ge = 0, Gm = 0, f1 = 0, f2 = 0, F1s = 0, F2s = 0;
  const double tau = Q2 / (4 * Mass2);
  //C Thorpe added Dec 2018
  //needed by hyperon production channels
  double f1p=0,f1n=0,f2p=0,f2n=0;
  switch (kind)
  {
    case 0:
    case 6:  // cc and mec qel part
      Ge = GEp - GEn;
      Gm = GMp - GMn;
      break;
    case 1:
    case 7:  // nc proton and mec qel part
      Ge = 0.5 * (GEp - GEn) - 2 * sin2thetaW * GEp;  //-0.5*GEs;
      Gm = 0.5 * (GMp - GMn) - 2 * sin2thetaW * GMp;  //-0.5*GMs;
      break;
    case 2:
    case 8:  // nc neutron and mec qel part
      Ge = 0.5 * (GEn - GEp) - 2 * sin2thetaW * GEn;  //-0.5*GEs;
      Gm = 0.5 * (GMn - GMp) - 2 * sin2thetaW * GMn;  //-0.5*GMs;
      break;
    case 3:  // cc mec
      Ge = GEp - GEn;
      Gm = sqrt(1 + 6.0 * Q2 / 1e6 * exp(-Q2 / 0.35 / 1e6)) * (GMp - GMn);
      break;
    case 4:                                           // nc proton mec
      Ge = 0.5 * (GEp - GEn) - 2 * sin2thetaW * GEp;  //-0.5*GEs;
      Gm =
          sqrt(1 + 6.0 * Q2 / 1e6 * exp(-Q2 / 0.35 / 1e6)) * 0.5 * (GMp - GMn) -
          2 * sin2thetaW * GMp;  //-0.5*GMs;
      break;
    case 5:                                           // nc neutron mec
      Ge = 0.5 * (GEn - GEp) - 2 * sin2thetaW * GEn;  //-0.5*GEs;
      Gm = sqrt(1 + 6.0 * Q2 / 1e6 * exp(-Q2 / 0.35 / 1e6)) *
           (0.5 * (GMn - GMp) - 2 * sin2thetaW * GMn);  //-0.5*GMs;
      break;
    case 10: //elastic ep scattering
        Ge=GEp;
        Gm=GMp;
        break;
    case 11: //elastic en scattering
        Ge=GEn;
        Gm=GMn;
        break;
    //hyperon channels: 12,13,14
    case 12:
    case 13:
    case 14:
      f1p = (1/(1+tau))*(GEp + tau*GMp);
      f1n = (1/(1+tau))*(GEn + tau*GMn);
      f2p = (1/(1+tau))*(GMp - GEp);
      f2n = (1/(1+tau))*(GMn - GEn);
      break;
  }
  f1 = (Ge + tau * Gm) / (1 + tau);
  f2 = (Gm - Ge) / (1 + tau);
  if ((kind == 1 or kind == 2) and
      strangeEM)  // strangeness in F1, F2 (only for kind!=0 i.e. nc)
    switch (strangeEM)
    {
      case 1:
      {
        double mian = (1 + tau) * pow2(1 + Q2 / MV2);
        f1 -= 0.5 * (0.53 / mian);
        f2 -= 0.5 * (-0.40 / mian);
        break;
      }
      case 2:  // more complex version from KG
      {
        double mian = (1 + tau) * pow2(1 + Q2 / MV2);
        ;
        double mus = -0.39;
        double f1s = 0.49;
        f1 -= 0.5 * mus / mian;
        f2 -= 0.5 * f1s * Q2 / GeV2 / mian;
        break;
      }
      default:
        break;  // no strange correction
    }
  //hyperon channels
  //Lambda zero
  if(kind == 12)
  {
    f1 = (-1)*pow(1.5,0.5)*f1p;
    f2 = (-1)*pow(1.5,0.5)*f2p;
    //SU(3) symmetry breaking correction
    if(sym_break == true)
    {
      f1 *= 0.976;
    }
  }
  //Sigma zero
  if(kind == 13)
  {
    f1 = (-1)*(f1p+2*f1n)/(pow(2,0.5));
    f2 = (-1)*(f2p+2*f2n)/(pow(2,0.5));
    //SU(3) symmetry breaking correction
    if(sym_break == true)
    {
      f1 *= 0.975;
    }
  }
  //Sigma minus
  if(kind==14)
  {
    f1 = (-1)*(f1p+2*f1n);
    f2 = (-1)*(f2p+2*f2n);
    if(sym_break == true)
    {
      f1 *= 0.975;
    }
  }
  return pair<double, double>(f1, f2);
}
// dipole electric form factor G_E^V
FF DipoleFF(const double q2)
{
  double a = 1.0 - q2 / MV2;
  double a2 = a * a;
  double tau = -q2 / (4 * Mass2);
  FF ff;
  ff.Q2 = -q2;
  ff.GEp = 1.0 / a2;
  //C Thorpe: Updated dipole FFs
  ff.GEn = (-1)*mu_n*tau/(1+Dipole_Lambda*tau)/a2;
  ff.GMp = mu_p / a2;
  ff.GMn = mu_n / a2;
  return ff;
};
FF bba03_FF(const double q2)
{
  const double Q2 = -q2 / GeV2;       // Q2 in GeV2
  const double tau = -q2 / (4 * Mass2);  // must be dimensionless
  FF ff;
  ff.Q2 = -q2;
  // hep-ex/0308005 BBA-2003 for Q2<6 GeV2
  ff.GEp =
      1.0 /
      (1.0 +
       Q2 * (3.253 +
             Q2 * (1.422 +
                   Q2 * (0.08582 +
                         Q2 * (0.3318 + Q2 * (-0.09371 + Q2 * 0.01076))))));
  ff.GEn = -mu_n * 0.942 * tau / (1 + 4.61 * tau) / pow2(1 - q2 / MV2);
  ff.GMp =
      mu_p / (1.0 +
              Q2 * (3.104 +
                    Q2 * (1.428 +
                          Q2 * (0.1112 +
                                Q2 * (-0.006981 +
                                      Q2 * (0.0003705 + Q2 * -0.7063e-5))))));
  ff.GMn =
      mu_n /
      (1.0 +
       Q2 * (3.043 +
             Q2 * (0.8548 + Q2 * (0.6806 + Q2 * (-0.1287 + Q2 * 0.008912)))));
  return ff;
}
FF bbba05_FForig(const double q2)
{
  const double tau = -q2 / (4 * Mass2);
  FF ff;
  ff.Q2 = -q2;
  // arXiv:hep-ex/0602017 BBBA05 for Q2<18 GeV
  ff.GEp =
      (1.0 - tau * 0.0578) / (1.0 + tau * (11.1 + tau * (13.6 + tau * 33.0)));
  ff.GEn = tau * (1.25 + tau * 1.30) /
           (1.0 + tau * (-9.86 + tau * (305.0 + tau * (-758.0 + tau * 802.0))));
  ff.GMp = mu_p * (1.0 + tau * 0.15) /
           (1.0 + tau * (11.1 + tau * (19.6 + tau * 7.54)));
  ff.GMn = mu_n * (1.0 + tau * 1.81) /
           (1.0 + tau * (14.1 + tau * (20.7 + tau * 68.7)));
  return ff;
}
double GetA(const double q2, const double* x)
{
  double eps_int[7] = {0., 1. / 6., 1. / 3., 1. / 2., 2. / 3., 5. / 6., 1.};
  double Q2 = -q2;
  double tau = -q2 / 4.0 / Mass2;
  double eps = 1.0 / (1 + sqrt(1 + (1. / tau)));
  double A = 0.0;
  // Sum up A
  for (int j = 0; j < 7; j++) {
    double mod = x[j];
    for (int k = 0; k < 7; k++) {
      if (k == j) continue;
      mod *= (eps - eps_int[k]) / (eps_int[j] - eps_int[k]);
    }
    A += mod;
  }
  return A;
};
FF bbba07_FF(const double q2)
{
  FF ff;
  ff.Q2 = -q2;
  // Set Variables
  double Q2 = -q2;
  double tau = -q2 / 4.0 / Mass2;
  double eps = 2.0 / (1 + sqrt(1 + (1. / tau)));
  // Calculate Lagrange
  double AEp = GetA(q2, p_AEp);
  double GK_Ep =
      (1. - 0.24 * tau) / (1.0 + tau * (10.98 + tau * (12.82 + tau * (21.97))));
  double AMp = GetA(q2, p_AMp);
  double GK_Mp = (1. + 0.1717 * tau) /
                 (1.0 + tau * (11.26 + tau * (19.32 + tau * (8.33))));
  double AEn = GetA(q2, p_AEn);
  //  double Ep = 0.0;
  double AMn = GetA(q2, p_AMn);
  //  double Mp = 0.0;
  // Define p form factors
  ff.GEp = AEp * GK_Ep;
  ff.GMp = AMp * GK_Mp * mu_p;
  // Define n form factors
  // a = 1.7
  // b = 3.3
  ff.GEn = AEn * ff.GEp * ((1.7 * tau) / (1 + 3.3 * tau));
  ff.GMn = AMn * ff.GMp * mu_n / mu_p;
  return ff;
};
FF JLab_FF(const double q2)
{
  const double Q2 = -q2 / GeV2;
  const double Q = sqrt(Q2);
  const double tau = -q2 / 4 / Mass2;
  FF ff;
  ff.Q2 = -q2;
  // PHYSICAL REVIEW C, VOLUME 65, 051001(R)
  ff.GEp = (1.0 - 0.13 * (Q2 - 0.04)) /
           (1.0 + (0.116 + (0.241 + 0.345 * Q2) * Q2) * Q +
            (2.874 + 1.006 * Q2) * Q2);
  ff.GMp = mu_p / (1.0 + (0.116 + (0.241 + 0.345 * Q2) * Q2) * Q +
                   (2.874 + 1.006 * Q2) * Q2);
  // PHYSICAL REVIEW C, VOLUME 51, 409 (1995)
  ff.GEn = -1.25 * mu_n * tau / (1 + 18.3 * tau) / (1 + Q2 * 1.0e+6 / MV2) /
           (1 + Q2 * 1.0e+6 / MV2);
  ff.GMn = mu_n / (1.0 - (1.74 + 7.63 * Q2) * Q + (9.29 + 4.63 * Q2) * Q2);
  return ff;
};
FF npar_FF(const double q2) {
  const double Q2 = -q2 / GeV2;
  const double Q = sqrt(Q2);
  const double tau = -q2 / 4 / Mass2;
  FF ff;
  ff.Q2 = -q2;
  ff.GEp = 1 / (1 + Q * (0.62 + Q * (0.68 + Q * (2.80 + Q * 0.83))));
  ff.GMp =
      mu_p / (1 + Q * (0.35 + Q * (2.44 + Q * (0.50 + Q * (1.04 + Q * 0.34)))));
  ff.GMn = mu_n / (1 + Q * (-0.74 + Q * (9.29 + Q * (-7.63 + Q * 4.63))));
  const double a = 1.25;
  const double b = 10.4;
  double GD = 1 / pow2(1 - q2 / ( MV2 ));
  ff.GEn = -a * mu_n * tau * GD / (1 + b * tau);
  bool strange = false;  // comment out next 2 lines
  if (strange) ff.GEp = ff.GMp / mu_p * (1 - 0.14 * (-q2 / GeV2 - 0.3));
  return ff;
};
static double funkcja1(double Q2, double* w) {
  double q2 = Q2 / GeV / GeV;

  return w[4] / (1 + exp(-(q2 * w[0] + w[1]))) +
         w[5] / (1 + exp(-(q2 * w[2] + w[3]))) + w[6];
}
static double funkcja2(double Q2, double* w) {
  double q2 = Q2 / GeV / GeV;

  return w[6] / (1 + exp(-(q2 * w[0] + w[1]))) +
         w[7] / (1 + exp(-(q2 * w[2] + w[3]))) +
         w[8] / (1 + exp(-(q2 * w[4] + w[5]))) + w[9];
}
FF kg_FF(const double q2) {
  static double tab_gen[7] = {10.19704,  2.36812,  -1.144266, -4.274101,
                              0.8149924, 2.985524, -0.7864434};
  static double tab_gmn[10] = {3.19646,    2.565681, 6.441526,  -2.004055,
                               -0.2972361, 3.606737, -3.135199, 0.299523,
                               1.261638,   2.64747};
  static double tab_gep[10] = {3.930227,   0.1108384, -5.325479, -2.846154,
                               -0.2071328, 0.8742101, 0.4283194, 2.568322,
                               2.577635,   -1.185632};
  static double tab_gmp[10] = {-2.862682,  -1.560675, 2.321148, 0.1283189,
                               -0.2803566, 2.794296,  1.726774, 0.861083,
                               0.4184286,  -0.1526676};
  static double tab_axial[10] = {-26.10885, 1.823041,  -8.391283, -7.737312,
                                 15.27646,  0.3992788, -1.350184, -0.2021121,
                                 -2.870517, 3.879841};
  static double MA_nnff = 1.015 * GeV;
  double Q2 = -q2;
  double y = pow2(1.0 + Q2 / MV2);
  FF ff;
  ff.Q2 = Q2;
  ff.GEp = funkcja2(Q2, tab_gep) / y;
  ff.GEn = funkcja1(Q2, tab_gen);
  ff.GMp = funkcja2(Q2, tab_gmp) * mu_p / y;
  ff.GMn = funkcja2(Q2, tab_gmn) * mu_n / y;
  return ff;
}
FF bbba05_FF(const double q2)
{
  FF ff;
  ff.Q2 = -q2;
  double Q2 = -q2;
  double tau = -q2 / (4.0 * Mass2);
  ff.GEp = (1.0 - tau * 0.0578) / (1.0 + tau * (11.1 + tau * (13.6 + tau * 33.0)));
  ff.GEn = tau * (1.25 + tau * 1.30) / (1.0 + tau * (-9.86 + tau * (305.0 + tau * (-758.0 + tau * 802.0))));
  ff.GMp = mu_p * (1.0 + tau * 0.15) / (1.0 + tau * (11.1 + tau * (19.6 + tau * 7.54)));
  ff.GMn = mu_n * (1.0 + tau * 1.81) / (1.0 + tau * (14.1 + tau * (20.7 + tau * 68.7)));
  return ff;
}
// Calculate vector form factors
pair<double, double> f12(double q2, int kind)
{

  FF ff = FFfromq2(q2);
  ff.Q2 = -q2;

  return ff.f12(kind);
}

//____________________________________
// AXIAL FORM FACTOR FUNCTIONS
//____________________________________

/// BBBA07 Dipole
double bbba07_FA(double q2, double ma)
{
  double FA = dipole_FA(q2, ma);
  double AAx = GetA(q2, p_AAx);
  return FA* AAx;
}
//____________________________________
/// 2 Component Model
double comp2_FA(double q2, double ma)
{
  (void) ma; // MA Ignored for this function
  double ma_axl = 1.230;
  double gterm = 1.0 / pow2(1.0 + axial_ff_gamma * q2);
  double aterm = (1.0 - axial_ff_alpha + \
		  (axial_ff_alpha * (ma_axl * ma_axl) / (ma_axl*ma_axl + q2)));
  return gA * gterm * aterm;
}
//____________________________________
/// 3 Component model
double comp3_FA(double q2, double ma)
{
  (void) ma; // MA Ignored for this function
  double comp2_term = comp2_FA(q2, ma);
  double exp_term = gA * sqrt( axial_ff_theta) * axial_ff_beta * \
    exp(axial_ff_theta + axial_ff_beta * q2);
  return comp2_term + exp_term;
}
//____________________________________
/// Z Expansion Model
static const int kmax = 15;
static double zexp_aterms[15];
static int zexp_nterms;
static bool zexp_q4limit;
static double zexp_tc;
static double zexp_t0;
double zexp_GetZ(double q2)
{
  // T Cut
  double num = sqrt(zexp_tc - q2) - sqrt(zexp_tc - zexp_t0);
  double den = sqrt(zexp_tc - q2) + sqrt(zexp_tc - zexp_t0);
  return num/den;
}
void PrintZExpTerms(bool showFA)
{
  cout << " ZEXP State! " << endl;
  cout << " ------------------" << endl;
  cout << " T0 = " << zexp_t0 << endl;
  cout << " TC = " << zexp_tc << endl;
  int ncount = zexp_nterms;
  if (zexp_q4limit > 0) ncount += 4;
  for (int i = 0; i <= ncount; i++)
  {
    cout << "ZEXP A" << i << " = " << zexp_aterms[i] << endl;
  }
  if (showFA)
  {
    cout << " FA Values " << endl;
    cout << " FAZ(0.00) = " << zexp_FA(0.00,0.0) << endl;
    cout << " FAZ(0.25) = " << zexp_FA(0.25,0.0) << endl;
    cout << " FAZ(0.50) = " << zexp_FA(0.50,0.0) << endl;
    cout << " FAZ(0.75) = " << zexp_FA(0.75,0.0) << endl;
    cout << " FAZ(1.00) = " << zexp_FA(1.00,0.0) << endl;
    cout << " FAZ(1.50) = " << zexp_FA(1.50,0.0) << endl;
    cout << " FAZ(2.00) = " << zexp_FA(2.00,0.0) << endl;
    cout << " FAZ(3.00) = " << zexp_FA(3.00,0.0) << endl;
  }
}
double zexp_FA(double q2, double ma)
{
  (void) ma; // MA Ignored for this function
  // Read Params
  q2 = -fabs(q2);
  // Calculate z
  double z = zexp_GetZ(q2);
  double FA = 0.0;
  int ncount = zexp_nterms;
  if (zexp_q4limit > 0) ncount += 4;
  for (int i = 0; i <= ncount; i++)
  {
    FA += pow(z,i) * zexp_aterms[i];
  }
  return FA;
}
void zexp_applysumrules()
{
  //  PrintZExpTerms(false);
  // The Code below is from private correspondence
  // with Aaron Meyer on the calculation of sum Rules.
  // - P. Stowell
  // Gives the Q^-4 format at high Q^2
  double k0 = (double)zexp_nterms;
  double z0 = zexp_GetZ(0.0);
  double k1 = (double)zexp_nterms+1;
  double z1 = pow(z0, (int)k1);
  double k2 = (double)zexp_nterms+2;
  double z2 = pow(z0, (int)k2);
  double k3 = (double)zexp_nterms+3;
  double z3 = pow(z0, (int)k3);
  double k4 = (double)zexp_nterms+4;
  double z4 = pow(z0, (int)k4);
  // Get Delta (z shifts through terms)
  double del =  6.0
    - 1.0 * k4 * k3 * k2 * z1
    + 3.0 * k4 * k3 * z2 * k1
    - 3.0 * k4 * z3 * k2 * k1
    + 1.0 * z4 * k3 * k2 * k1;
  // Setup Starting Parameters
  double b0  = 0.0;
  double b1  = 0.0;
  double b2  = 0.0;
  double b3  = 0.0;
  double b0z = 1.267;
  for (int ki = 1;ki <= zexp_nterms;ki++){
    b0 += zexp_aterms[ki];
    b1 += ki * zexp_aterms[ki];
    b2 += ki * (ki - 1) * zexp_aterms[ki];
    b3 += ki * (ki - 1) * (ki - 2) * zexp_aterms[ki];
    b0z += zexp_aterms[ki]*pow(z0,ki);
  }
  // A0
  zexp_aterms[0] =
    (- 6.*b0z - b0*(del-6.) + b3*(-z1 + 3.*z2 - 3.*z3 + z4)
     + b2 * (3.*k2*z1 - 3.*(3.*k0+5.)*z2 + 3.*(3.*k0+4.)*z3 - 3.*k1*z4)
     + b1 * (-3.*k3*k2*z1 + 3.*k3*(3.*k0+4.)*z2
	     -3.*k1*(3.*k0+8.)*z3 + 3.*k2*k1*z4) ) / (del);
  // A1
  zexp_aterms[(int)k1] =                        \
    (- (b0-b0z)*k4*k3*k2                                \
     + b3*(1. - 0.5*k4*k3*z2 + k4*k2*z3 - 0.5*k3*k2*z4) \
     + b2*(-3.*k2 + k4*k3*k2*z2                         \
	   - k4*k2*(2.*k0+3.)*z3 + k3*k2*k1*z4)         \
     + b1*(3.*k3*k2 - 0.5*k4*k3*k3*k2*z2                \
	   + k4*k3*k2*k1*z3 - 0.5*k3*k2*k2*k1*z4)       \
     ) / (del) ;
  // A2
  zexp_aterms[(int)k2] =
    ( + 3.*(b0-b0z)*k4*k3*k1                                    \
      + b3*(-3. + 0.5*k4*k3*z1 - (3./2.)*k4*k1*z3 + k3*k1*z4)   \
      + b2*(3.*(3.*k0+5) - k4*k3*k2*z1 + 3*k4*k1*k1*z3          \
	    - k3*k1*(2.*k0+1.)*z4)                              \
      + b1*(-3.*k3*(3.*k0+4.) + 0.5*k4*k3*k3*k2*z1              \
	    -(3./2.)*k4*k3*k1*k0*z3 + k3*k2*k1*k0*z4)           \
      ) / (del);
  // A3
  zexp_aterms[(int)k3] =                                       \
    (- 3.*(b0-b0z)*k4*k2*k1                                    \
     + b3*(3. - k4*k2*z1 + (3./2.)*k4*k1*z2 - 0.5*k2*k1*z4)   \
     + b2*(-3.*(3.*k0+4.) + k4*k2*(2.*k0+3.)*z1               \
	   - 3.*k4*k1*k1*z2 + k2*k1*k0*z4)                    \
     + b1*(3.*k1*(3.*k0+8.) - k4*k3*k2*k1*z1                  \
	   +(3./2.)*k4*k3*k1*k0*z2 - (1./2.)*k2*k1*k1*k0*z4)  \
     ) / (del);
  // A4
  zexp_aterms[(int)k4] =                                      \
    ( + (b0-b0z)*k3*k2*k1                                     \
      + b3*(-1. + (1./2.)*k3*k2*z1 - k3*k1*z2 + 0.5*k2*k1*z3) \
      + b2*(3.*k1 - k3*k2*k1*z1 + k3*k1*(2.*k0+1.)*z2         \
	    - k2*k1*k0*z3)                                    \
      + b1*(-3.*k2*k1 + 0.5*k3*k2*k2*k1*z1                    \
	    -k3*k2*k1*k0*z2 + 0.5*k2*k1*k1*k0*z3)             \
      ) / (del);
  return;
}
void zexp_applyq0limit()
{
  double z = zexp_GetZ(0.0);
  double FA = 0.0;
  for (int i = 1; i <= zexp_nterms; i++)
  {
    FA = FA + pow(z, i)* zexp_aterms[i];
  }
  zexp_aterms[0]= gA - FA;
  return;
}

// Average LQCD axial form factor
double LQCD_FA(const double q2, const double /*bandSwitch*/)
{
    constexpr double tcut = 161604.0;   // (3 m_pi)^2 in MeV^2
    constexpr double t0   = -500000.0;  // MeV^2

    // Ensure physical domain: q² must be < tcut
    if (q2 >= tcut) return 0.0;

    const double sqrt1 = std::sqrt(tcut - q2);
    const double sqrt2 = std::sqrt(tcut - t0);
    const double z = (sqrt1 - sqrt2) / (sqrt1 + sqrt2);

    constexpr double a[] = {
        0.7174202,   // a0
       -1.720897,    // a1
        0.3098271,   // a2
        1.621258,    // a3
       -0.2750699,   // a4
       -1.252979,    // a5
        0.6004408    // a6
    };

    double FA = 0.0, z_power = 1.0;
    for (int k = 0; k <= 6; ++k)
    {
        FA += a[k] * z_power * (-1);
        z_power *= z;
    }
    
//        std::cout << -q2/1e6 << " " << FA << std::endl;

    return FA;
}

// axial form factor determined by MINERvA, Nature 614, 48 (2023) : implementation by A. Ankowski
double MINERvA_parametrization(const double qSq, const double bandWeight)
{
    static const double tCut( 9*std::pow(139.6*MeV, 2) );
    static const double t0( -0.75*GeV2 );
    static const double dummy( sqrt(tCut - t0) );
    const double var( sqrt(tCut - qSq) );
    const double z( (var - dummy)/(var + dummy) );
    ///from the correspondence with Tejin Cai: effectively it gives gA = -1.27212
    const double zExp(   -0.502678 + z*(1.50 + z*(-1.2 + z*(-0.149228 + z*(0.157978 + z*(0.510344 + z*(-0.411635 + z*(0.127889 - z*0.0326702)))))))   );

    if( not bandWeight )
    {
        return zExp;
    }
    static bool firstTime( true );

    static const int parametrizationOrder(8 + 1);
    static const int polynomialOrder( 2*(parametrizationOrder - 1) + 1 );
    static double coefficients[polynomialOrder];

    if ( firstTime )
    {
        for ( int mCnt(0); mCnt < polynomialOrder; ++mCnt)
        {
            coefficients[mCnt] = 0.0;
        }
        const double var0( sqrt(tCut) );
        const double z0( (var0 - dummy)/(var0 + dummy) );///z for Q^2 = 0
        const int paramNum(4);
        ///define the correlation matrix
        const double correlationMatrix[paramNum][paramNum] =
        {
           { 1.000, 0.012,-0.930, 0.520},
           { 0.012, 1.000,-0.320,-0.780},
           {-0.930,-0.320, 1.000,-0.270},
           { 0.520,-0.780,-0.270, 1.000}
        };
        const double parameterUncertainties[paramNum] =
        {
           0.31, 0.70, 1.90, 3.50
        };
        const int funcNum( paramNum + 1 );
        ///f0--f4
        const int fi[funcNum][parametrizationOrder] =
        {
        ///  0      1      2      3      4      5      6      7      8 ///the power of z
           { 1,     0,     0,     0,     0,   -56,   140,  -120,    35 },
           { 0,     1,     0,     0,     0,   -35,    84,   -70,    20 },
           { 0,     0,     1,     0,     0,   -20,    45,   -36,    10 },
           { 0,     0,     0,     1,     0,   -10,    20,   -15,     4 },
           { 0,     0,     0,     0,     1,    -4,     6,    -4,     1 }
        };
        const double f0_for_z0(   fi[0][0] + z0*( fi[0][1] + z0*( fi[0][2] + z0*( fi[0][3] + z0*( fi[0][4] + z0*( fi[0][5] + z0*( fi[0][6] + z0*( fi[0][7] + z0*fi[0][8] ) ) ) ) ) ) )   );
        const double f1_for_z0(   fi[1][0] + z0*( fi[1][1] + z0*( fi[1][2] + z0*( fi[1][3] + z0*( fi[1][4] + z0*( fi[1][5] + z0*( fi[1][6] + z0*( fi[1][7] + z0*fi[1][8] ) ) ) ) ) ) )   );
        const double f2_for_z0(   fi[2][0] + z0*( fi[2][1] + z0*( fi[2][2] + z0*( fi[2][3] + z0*( fi[2][4] + z0*( fi[2][5] + z0*( fi[2][6] + z0*( fi[2][7] + z0*fi[2][8] ) ) ) ) ) ) )   );
        const double f3_for_z0(   fi[3][0] + z0*( fi[3][1] + z0*( fi[3][2] + z0*( fi[3][3] + z0*( fi[3][4] + z0*( fi[3][5] + z0*( fi[3][6] + z0*( fi[3][7] + z0*fi[3][8] ) ) ) ) ) ) )   );
        const double f4_for_z0(   fi[4][0] + z0*( fi[4][1] + z0*( fi[4][2] + z0*( fi[4][3] + z0*( fi[4][4] + z0*( fi[4][5] + z0*( fi[4][6] + z0*( fi[4][7] + z0*fi[4][8] ) ) ) ) ) ) )   );
        ///fi - fi_for_z0/f0_for_z0 (i = 1, 2, 3, 4)
        const double coefficientsMatrix[paramNum][parametrizationOrder] =
        {
        ///                     0      1      2      3      4      5      6      7      8
           { -f1_for_z0/f0_for_z0,     1,     0,     0,     0,   -35,    84,   -70,    20},
           { -f2_for_z0/f0_for_z0,     0,     1,     0,     0,   -20,    45,   -36,    10},
           { -f3_for_z0/f0_for_z0,     0,     0,     1,     0,   -10,    20,   -15,     4},
           { -f4_for_z0/f0_for_z0,     0,     0,     0,     1,    -4,     6,    -4,     1},
        };
        for ( int iCnt(0); iCnt < paramNum; ++iCnt)
        {
            for ( int jCnt(0); jCnt < paramNum; ++jCnt)
            {
                    for ( int ia(0); ia < parametrizationOrder; ++ia)
                    {
                        for ( int jb=0; jb < parametrizationOrder; ++jb)
                            coefficients[ia + jb] += correlationMatrix[iCnt][jCnt]*coefficientsMatrix[iCnt][ia]*parameterUncertainties[iCnt]*coefficientsMatrix[jCnt][jb]*parameterUncertainties[jCnt];
                    }
            }
        }
        for (int i = 0; i < polynomialOrder; ++i)
        {
    //std::cout << "Coefficient " << i << ": " << coefficients[i] << std::endl;
        }
        firstTime = false;
    }

    const double uncSq(   coefficients[0] + z*( coefficients[1] + z*( coefficients[2] + z*( coefficients[3] + z*( coefficients[4] + z*( coefficients[5] + z*( coefficients[6] + z*( coefficients[7] + z*( coefficients[8] + z*( coefficients[9] + z*( coefficients[10] + z*( coefficients[11] + z*( coefficients[12] + z*( coefficients[13] + z*( coefficients[14] + z*( coefficients[15] + z*coefficients[16] ) ) ) ) ) ) ) ) ) ) ) ) ) ) )   );

    return zExp - bandWeight * (sqrt( uncSq ));
}

//double mva_errorBar( 0.0 ); // Upper limit = 1, Lower limit = -1
double MINERvA_FA(const double qSq, const double ma)
{
    double weight = mva_errorBar;
    return MINERvA_parametrization( qSq, weight );
}

// Standard Dipole
double dipole_FA(double q2, double ma)
{
  return gA / pow2(1 - q2 / ma / ma);
}

// axial form factor from neutrino-deuteron scattering data, Meyer et al., Phys. Rev. D 93, 113015 (2016)
double deuterium_FA_parametrization(const double qSq, const double bandSwitch)
{
    static const double tCut( 9*std::pow(140.0*MeV, 2) );
    static const double t0( -0.28*GeV2 );
    static const double dummy( sqrt(tCut - t0) );

    const double var( sqrt(tCut - qSq) );
    const double z( (var - dummy)/(var + dummy) );

    const double zExp(   -0.759 + z*(2.30 + z*(-0.6 + z*(-3.8 + z*(2.3 + z*(2.16 + z*(-0.896 + z*(-1.58 + z*0.823)))))))   );

    if ( not bandSwitch )
        return zExp;

    static bool firstTime( true );

    static const int parametrizationOrder(8 + 1);
    static const int polynomialOrder( 2*(parametrizationOrder - 1) + 1 );

    static double coefficients[polynomialOrder];

    if ( firstTime )
    {
        for ( int mCnt(0); mCnt < polynomialOrder; ++mCnt)
        {
            coefficients[mCnt] = 0.0;
            //std::cout<<mCnt<<" "<<coefficients[mCnt]<<std::endl;
        }

        const double var0( sqrt(tCut) );
        const double z0( (var0 - dummy)/(var0 + dummy) );///z for Q^2 = 0


        const int paramNum(4);

        ///define the correlation matrix
        const double correlationMatrix[paramNum][paramNum] =
        {
           { 1.000, 0.350,-0.678, 0.611},
           { 0.350, 1.000,-0.898, 0.367},
           {-0.678,-0.898, 1.000,-0.685},
           { 0.611, 0.367,-0.685, 1.000}
        };

        const double parameterUncertainties[paramNum] =
        {
           0.13, 1.0, 2.5, 2.7
        };

        const int funcNum( paramNum + 1 );

        ///f0--f4
        const int fi[funcNum][parametrizationOrder] =
        {
        ///  0      1      2      3      4      5      6      7      8 ///the power of z
           { 1,     0,     0,     0,     0,   -56,   140,  -120,    35 },
           { 0,     1,     0,     0,     0,   -35,    84,   -70,    20 },
           { 0,     0,     1,     0,     0,   -20,    45,   -36,    10 },
           { 0,     0,     0,     1,     0,   -10,    20,   -15,     4 },
           { 0,     0,     0,     0,     1,    -4,     6,    -4,     1 }
        };

        const double f0_for_z0(   fi[0][0] + z0*( fi[0][1] + z0*( fi[0][2] + z0*( fi[0][3] + z0*( fi[0][4] + z0*( fi[0][5] + z0*( fi[0][6] + z0*( fi[0][7] + z0*fi[0][8] ) ) ) ) ) ) )   );
        const double f1_for_z0(   fi[1][0] + z0*( fi[1][1] + z0*( fi[1][2] + z0*( fi[1][3] + z0*( fi[1][4] + z0*( fi[1][5] + z0*( fi[1][6] + z0*( fi[1][7] + z0*fi[1][8] ) ) ) ) ) ) )   );
        const double f2_for_z0(   fi[2][0] + z0*( fi[2][1] + z0*( fi[2][2] + z0*( fi[2][3] + z0*( fi[2][4] + z0*( fi[2][5] + z0*( fi[2][6] + z0*( fi[2][7] + z0*fi[2][8] ) ) ) ) ) ) )   );
        const double f3_for_z0(   fi[3][0] + z0*( fi[3][1] + z0*( fi[3][2] + z0*( fi[3][3] + z0*( fi[3][4] + z0*( fi[3][5] + z0*( fi[3][6] + z0*( fi[3][7] + z0*fi[3][8] ) ) ) ) ) ) )   );
        const double f4_for_z0(   fi[4][0] + z0*( fi[4][1] + z0*( fi[4][2] + z0*( fi[4][3] + z0*( fi[4][4] + z0*( fi[4][5] + z0*( fi[4][6] + z0*( fi[4][7] + z0*fi[4][8] ) ) ) ) ) ) )   );

        ///fi - fi_for_z0/f0_for_z0 (i = 1, 2, 3, 4)
        const double coefficientsMatrix[paramNum][parametrizationOrder] =
        {
        ///                     0      1      2      3      4      5      6      7      8
           { -f1_for_z0/f0_for_z0,     1,     0,     0,     0,   -35,    84,   -70,    20},
           { -f2_for_z0/f0_for_z0,     0,     1,     0,     0,   -20,    45,   -36,    10},
           { -f3_for_z0/f0_for_z0,     0,     0,     1,     0,   -10,    20,   -15,     4},
           { -f4_for_z0/f0_for_z0,     0,     0,     0,     1,    -4,     6,    -4,     1},
        };

        for ( int iCnt(0); iCnt < paramNum; ++iCnt)
        {
            for ( int jCnt(0); jCnt < paramNum; ++jCnt)
            {
                    for ( int ia(0); ia < parametrizationOrder; ++ia)
                    {
                        for ( int jb(0); jb < parametrizationOrder; ++jb)
                            coefficients[ia + jb] += correlationMatrix[iCnt][jCnt]*coefficientsMatrix[iCnt][ia]*parameterUncertainties[iCnt]*coefficientsMatrix[jCnt][jb]*parameterUncertainties[jCnt];
                    }
            }
        }

        firstTime = false;
    }

    const double uncSq(   coefficients[0] + z*( coefficients[1] + z*( coefficients[2] + z*( coefficients[3] + z*( coefficients[4] + z*( coefficients[5] + z*( coefficients[6] + z*( coefficients[7] + z*( coefficients[8] + z*( coefficients[9] + z*( coefficients[10] + z*( coefficients[11] + z*( coefficients[12] + z*( coefficients[13] + z*( coefficients[14] + z*( coefficients[15] + z*coefficients[16] ) ) ) ) ) ) ) ) ) ) ) ) ) ) )   );

    return (bandSwitch > 0) ? zExp - sqrt( uncSq ) : zExp + sqrt( uncSq );
}

double deuterium_FA( const double qSq, const double ma )
{
    double weight = deut_errorBar;
    return deuterium_FA_parametrization(qSq, weight);
}

// Calculate the axial form factors
pair<double, double> fap(double q2, int kind)
{
  double ksi = 3.706;
  //C Thorpe added Dec 2018  hyperon mass and kaon mass required for g3 calculation  using Phys Rev D98 (2018) no.3 033005, eq. 48
  double hyp_mass;
  double kmass = PDG::mass_K;
  double Ga, Fpa, Gas, Fpas;
  double Fp = 0, Fa = 0;
  switch (kind)
  {
    case 0:  // cc
      Fa = Axialfromq2(q2, rew.qel_cc_axial_mass.val); //Fa = Axialfromq2(q2, MA_cc);
      Fa *= axialcorr(axialFFset, q2);
      Fp = 2 * Mass2 * Fa / (piMass2 - q2);
      break;
    case 1:  // nc proton
      Fa = 0.5 * Axialfromq2(q2, rew.qel_nc_axial_mass.val); //Fa = 0.5 * Axialfromq2(q2, MA_nc);
      // Fp=2.0*Mass2*Fa/(piMass2 - q2) ;
      break;
    case 2:  // nc neutron
      Fa = -0.5 * Axialfromq2(q2, rew.qel_nc_axial_mass.val); //Fa = -0.5 * Axialfromq2(q2, MA_nc);
      // Fp=2.0*Mass2*Fa/(piMass2 - q2) ;
      break;
    case 3:
    case 6:  // mec cc
      Fa = Axialfromq2(q2, MA_cc_mec);
      Fa *= axialcorr(axialFFset, q2);
      Fp = 2 * Mass2 * Fa / (piMass2 - q2);
      break;
    case 4:
    case 7:  // mec nc proton
      Fa = 0.5 * Axialfromq2(q2, MA_nc_mec);
      // Fp=2.0*Mass2*Fa/(piMass2 - q2) ;
      break;
    case 5:
    case 8:  // mec nc neutron
      Fa = -0.5 * Axialfromq2(q2, MA_nc_mec);
      // Fp=2.0*Mass2*Fa/(piMass2 - q2) ;
      break;
    //hyperon production channels
    //Lambda zero;
    case 12:
      // need to add x=F/(F+D) value to constants file
      // x = 0.36543014996
      //Fa = Axialfromq2(q2, MA_hyp);
      Fa = Axialfromq2(q2,MA_hyp);
      Fa *= (-1)*(1+2*Axial_x)/sqrt(6);
      //SU(3) symmetry breaking
      if(sym_break == true)
      {
        Fa *= 1.072;
      }
      hyp_mass = PDG::mass_Lambda;
      //Fp = Fa*(M+hyp_mass)*(M+hyp_mass)/(2*(kmass*kmass-q2));
      Fp = Fa*(M+hyp_mass)*(M+hyp_mass)/(2*(kmass*kmass-q2));
      break;
    //Sigma zero
    case 13:
      //need to add x=F/(F+D) value to constants file
      // x = 0.36543014996
      Fa = Axialfromq2(q2, MA_hyp);
      //Fa  = 1.267/((1-q2/(MA_hyp*MA_hyp))*(1-q2/(MA_hyp*MA_hyp)));
      Fa *= (1-2*Axial_x)/(pow(2,0.5));
      //SU(3) symmetry breaking
      if(sym_break == true)
      {
        Fa *= 1.051;
      }
      hyp_mass = PDG::mass_Sigma;
      Fp = Fa*(M+hyp_mass)*(M+hyp_mass)/(2*(kmass*kmass-q2));
      break;
    //Sigma minus
    case 14:
      // need to add x=F/(F+D) value to constants file
      // x = 0.36543014996
      Fa = Axialfromq2(q2, MA_hyp);
      //Fa  = 1.267/((1-q2/(MA_hyp*MA_hyp))*(1-q2/(MA_hyp*MA_hyp)));
      Fa *= (1-2*Axial_x);
      //SU(3) symmetry breaking
      if(sym_break == true)
      {
        Fa *= 1.056;
      }
      hyp_mass = PDG::mass_SigmaM;
      Fp = Fa*(M+hyp_mass)*(M+hyp_mass)/(2*(kmass*kmass-q2));
      break;
  }
  if ((kind == 1 or kind == 2) and strange) {
    switch (strange) {
      case 1: {
        Fa -= -0.5 * delta_s / pow2(1 - q2 / MA_nc / MA_nc);
      } break;
      case 2:  // new implementation
      {
        Fa -= -0.5 * delta_s / pow2(1 - q2 / MA_s / MA_s);
      }
      default:
        break;  // no strange content
    }
  }
  return pair<double, double>(Fa, Fp);
}
double axialcorr(int axialFF, double q2)
{
  double min;   // maximal reduction
  double max;   // maximal enhancement
  double szer;  // reduction in Q2
  double dlug;  // enhancement range
  switch (axialFF) {
    case 2:
      min = 0.9;   // maximal reduction
      max = 1.1;   // maximal enhancement
      szer = 0.2;  // reduction in Q2
      dlug = 2.0;  // enhancement range
      break;
    case 3:
      min = 0.8;   // maximal reduction
      max = 1.2;   // maximal enhancement
      szer = 0.2;  // reduction in Q2
      dlug = 2.0;  // enhancement range
    default:
      return 1;
  }
  // parabolic approximation
  double x = -q2 / GeV2;
  if (x < szer)
    return 1.0 + 4.0 / szer / szer * (1.0 - min) * x * (x - szer);
  else if (x < dlug)
    return 1.0 -
           4.0 * (max - 1.0) / pow2(dlug - szer) * (x - szer) * (x - dlug);
  else
    return 1.0;
}
// C Thorpe Added Dec 2018
// Second class current Form Factors in Hyperon production
pair<double,double>g2(double q2,int kind){
  // real and imaginary parts of form factor g2
  // imaginary values corresponding to TRV
  // assume dipole form with same axial mass as g1
  double Rg2 = (-1)*Rg20/(pow(1-q2/(MA_hyp*MA_hyp),2));
  double Ig2 = (-1)*Ig20/(pow(1-q2/(MA_hyp*MA_hyp),2));
  switch (kind)
  {
  //Lambda zero production
  case 12:
    Rg2 *= (-1)*(1+2*Axial_x)/pow(6,0.5);
    Ig2 *= (-1)*(1+2*Axial_x)/pow(6,0.5);
    break;
  //Sigma zero production
  case 13:
    Rg2 *= (1-2*Axial_x)/(pow(2,0.5));
    Ig2 *= (1-2*Axial_x)/(pow(2,0.5));
    break;
  //Sigma minus production
  case 14:
    Rg2 *= (1-2*Axial_x);
    Ig2 *= (1-2*Axial_x);
    break;
  // for ds=0 quasielastic do not include SCC for the time being
  default:
    Rg2 =0;
    Ig2 =0;
    break;
  }
  return pair<double,double>(Rg2,Ig2);
}
/// Form Factor Configuration
void ff_configure(params & p)
{
  // Select vector form factor model based on input parameters
  switch (p.qel_vector_ff_set)
  {
    case 1:
      FFfromq2 = DipoleFF;
      break;
    case 2:
      FFfromq2 = bbba05_FF;
      break;
    case 3:
      FFfromq2 = bba03_FF;
      break;
    case 4:
      FFfromq2 = JLab_FF;
      break;
    case 5:
      FFfromq2 = kg_FF;
      break;
    case 6:
      FFfromq2 = npar_FF;
      break;
    case 7:
      FFfromq2 = bbba07_FF;
      break;
    default:
      throw("bad ffset");
  };
  // Set BBBA07 Vector Pars
  if (p.qel_vector_ff_set == 7) {
    for(int i=0;i<7;i++)
    {
      p_AEp[i] = (&rew.bba07_AEp1)[i].val;
      p_AMp[i] = (&rew.bba07_AMp1)[i].val;
      p_AEn[i] = (&rew.bba07_AEn1)[i].val;
      p_AMn[i] = (&rew.bba07_AMn1)[i].val;
    }
  }
  // AXIAL FF SET
  axialFFset = p.qel_axial_ff_set;
  switch(axialFFset)
  {
  // correct?
  case 1:
  case 2:
  case 3: Axialfromq2 = dipole_FA;  break;
  case 4: Axialfromq2 = bbba07_FA;
    for(int i=0;i<7;i++)
      p_AAx[i] = (&rew.bba07_AAx1)[i].val;
          break;
  case 5: Axialfromq2 = comp2_FA;
    axial_ff_gamma = rew.qel_axial_2comp_gamma.val;
    axial_ff_alpha = rew.qel_axial_2comp_alpha.val;
    break;
  case 6: Axialfromq2 = comp3_FA;
    axial_ff_gamma = rew.qel_axial_2comp_gamma.val;
    axial_ff_alpha = rew.qel_axial_2comp_alpha.val;
    axial_ff_beta = rew.qel_axial_3comp_beta.val;
    axial_ff_theta = rew.qel_axial_3comp_theta.val;
    break;
  case 7:  Axialfromq2 = zexp_FA; // ZEXPANSION
    zexp_nterms = rew.zexp_nterms.val;
    if (zexp_nterms > 10)  // Truncate at 10 terms
      zexp_nterms = 10;
    // Get T values
    zexp_tc = rew.zexp_tc.val;
    zexp_t0 = rew.zexp_t0.val;
    // Get Terms
    for(int i=0;i<10;i++)
       zexp_aterms[i] = (&rew.zexp_a0)[i].val;
    // Set Limits
    zexp_q4limit = bool(rew.zexp_q4limit.val);
    if (zexp_q4limit) zexp_applysumrules();
    else zexp_applyq0limit();
    //    PrintZExpTerms(true);
    //    sleep(5);
    break;
    case 8:
      mva_errorBar = rew.qel_minerva_ff_scale.val;//jancheck
      Axialfromq2 = MINERvA_FA;
      break;
      case 9:
      deut_errorBar = rew.qel_deuterium_ff_scale.val;//jancheck
      Axialfromq2 = deuterium_FA;
      //deutFA_errorBar = rew.deutFA_errorBar.val;
      break;
     case 10:
      Axialfromq2 = LQCD_FA;
      break;
  default:
    throw("bad axial ffset");
    break;
  }
  strange = p.qel_strange;
  strangeEM = p.qel_strangeEM;
  delta_s = rew.delta_s.val;
  MA_cc = rew.qel_cc_axial_mass.val;
  MA_nc = rew.qel_nc_axial_mass.val;
  MA_s = rew.qel_s_axial_mass.val;
  //locate parameter values in params
  //hyperon axial mass
  MA_hyp = p.hyp_axial_mass;
  //SCC Setup
  //values of g2 at Q2 = 0
  Rg20 = p.hyp_g2_Re_part;
  Ig20 = p.hyp_g2_Im_part;
  //symmetry breaking setup
  if(p.hyp_su3_sym_breaking)
  sym_break = true;
}
