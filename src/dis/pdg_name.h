#ifndef _pdg_name_h_
#define _pdg_name_h_

namespace DIS_PDG {

// quarks and diquark pdg from Big Book
const int quark_d = 1;
const int quark_u = 2;
const int quark_s = 3;
const int quark_c = 4;
const int anti_quark_d = -1 * quark_d;
const int anti_quark_u = -1 * quark_u;
const int anti_quark_s = -1 * quark_s;
const int anti_quark_c = -1 * quark_c;
const int diquark_dd_1 = 1103;
const int diquark_ud_0 = 2101;
const int diquark_ud_1 = 2103;
const int diquark_uu_1 = 2203;
const int diquark_sd_0 = 3101;
const int diquark_sd_1 = 3103;
const int diquark_su_0 = 3201;
const int diquark_su_1 = 3203;
const int diquark_ss_1 = 3303;
const int diquark_cd_0 = 4101;
const int diquark_cd_1 = 4103;
const int proton = 2212;
const int neutron = 2112;
const int pizero = 111;
const int piplus = 211;
const int piminus = -211;
const int Kplus = 321;
const int Kzero = 311;
const int Kminus = -321;
const int Dsplus = 431;
const int Dsminus = -431;
const int Dzero = 421;
const int Dzero_bar = -421;
const int Dplus = 411;
const int Dminus = -411;
const int Delta_plus = 2214;
const int Delta_plusplus = 2224;
const int Delta_minus = 1114;

const int Lambda = 3122;
const int Sigma = 3212;
const int SigmaP = 3222;
const int SigmaM = 3212;

const int e = 11;
const int nu_e = 12;
const int mu = 13;
const int nu_mu = 14;
const int tau = 15;
const int nu_tau = 16;

// mass of the particle listed above

const double quark_mass = 0.33;
const double quark_s_mass = 0.5;
const double quark_c_mass = 1.5;
const double diquark_mass = 0.77133;
const double diquark_mass_0 = 0.579;
const double diquark_s_mass = 0.805;
const double diquark_c_mass = 1.969;
const double Kplus_mass = 0.4936;
const double Kzero_mass = 0.498;
const double Dsplus_mass = 1.968;
const double Dzero_mass = 1.8646;
const double Dplus_mass = 1.8694;
const double Lambda_mass = 1.1157;
const double SigmaP_mass = 1.18937;
const double SigmaM_mass = 1.19745;
const double Sigma_mass = 1.19264;
const double mass_e = 0.000510999 * 1000;
const double mass_mu = 0.105658 * 1000;
const double mass_tau = 1.77699 * 1000;
const double m2_mu = mass_mu * mass_mu;

const double mc2 = quark_c_mass * quark_c_mass * 1000 * 1000;

// masa leptonu powstalego
}

using namespace DIS_PDG;

#endif
