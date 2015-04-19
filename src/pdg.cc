#include "pdg.h"
#include <iostream>
#include <map>
#include "jednostki.h"
#include <cstdlib>

using namespace std;

namespace PDG
{
  map < int, double >masstable;

  int init ()
  {
      masstable.insert (pair < int, double >(pdg_B_cP, 6.4));
      masstable.insert (pair < int, double >(pdg_B_s, 5.3696));
      masstable.insert (pair < int, double >(pdg_BstarP, 5.325));
      masstable.insert (pair < int, double >(pdg_Bstar, 5.325));
      masstable.insert (pair < int, double >(pdg_BP, 5.279));
      masstable.insert (pair < int, double >(pdg_B, 5.2794));
      masstable.insert (pair < int, double >(pdg_D_1_2420, 2.4222));
      masstable.insert (pair < int, double >(pdg_D_2_star_2460P, 2.459));
      masstable.insert (pair < int, double >(pdg_D_2_star_2460, 2.4589));
      masstable.insert (pair < int, double >(pdg_D_s_starP, 2.1121));
      masstable.insert (pair < int, double >(pdg_D_sP, 1.9683));
      masstable.insert (pair < int, double >(pdg_D_s1_2536P, 2.5353));
      masstable.insert (pair < int, double >(pdg_D_s2_2573P, 2.5724));
      masstable.insert (pair < int, double >(pdg_D_sJ_2317P, 2.3174));
      masstable.insert (pair < int, double >(pdg_D_sJ_2460P, 2.4593));
      masstable.insert (pair < int, double >(pdg_Dstar_2007, 2.0067));
      masstable.insert (pair < int, double >(pdg_Dstar_2010P, 2.01));
      masstable.insert (pair < int, double >(pdg_DP, 1.8694));
      masstable.insert (pair < int, double >(pdg_D, 1.8646));
      masstable.insert (pair < int, double >(pdg_Delta_1232P, 1.232));
      masstable.insert (pair < int, double >(pdg_Delta_1232PP, 1.232));
      masstable.insert (pair < int, double >(pdg_Delta_1232M, 1.232));
      masstable.insert (pair < int, double >(pdg_Delta_1232, 1.232));
      masstable.insert (pair < int, double >(pdg_Delta_1600P, 1.6));
      masstable.insert (pair < int, double >(pdg_Delta_1600PP, 1.6));
      masstable.insert (pair < int, double >(pdg_Delta_1600M, 1.6));
      masstable.insert (pair < int, double >(pdg_Delta_1600, 1.6));
      masstable.insert (pair < int, double >(pdg_Delta_1620P, 1.62));
      masstable.insert (pair < int, double >(pdg_Delta_1620PP, 1.62));
      masstable.insert (pair < int, double >(pdg_Delta_1620M, 1.62));
      masstable.insert (pair < int, double >(pdg_Delta_1620, 1.62));
      masstable.insert (pair < int, double >(pdg_Delta_1700P, 1.7));
      masstable.insert (pair < int, double >(pdg_Delta_1700PP, 1.7));
      masstable.insert (pair < int, double >(pdg_Delta_1700M, 1.7));
      masstable.insert (pair < int, double >(pdg_Delta_1700, 1.7));
      masstable.insert (pair < int, double >(pdg_Delta_1905P, 1.905));
      masstable.insert (pair < int, double >(pdg_Delta_1905PP, 1.905));
      masstable.insert (pair < int, double >(pdg_Delta_1905M, 1.905));
      masstable.insert (pair < int, double >(pdg_Delta_1905, 1.905));
      masstable.insert (pair < int, double >(pdg_Delta_1910P, 1.91));
      masstable.insert (pair < int, double >(pdg_Delta_1910PP, 1.91));
      masstable.insert (pair < int, double >(pdg_Delta_1910M, 1.91));
      masstable.insert (pair < int, double >(pdg_Delta_1910, 1.91));
      masstable.insert (pair < int, double >(pdg_Delta_1920P, 1.92));
      masstable.insert (pair < int, double >(pdg_Delta_1920PP, 1.92));
      masstable.insert (pair < int, double >(pdg_Delta_1920M, 1.92));
      masstable.insert (pair < int, double >(pdg_Delta_1920, 1.92));
      masstable.insert (pair < int, double >(pdg_Delta_1930P, 1.93));
      masstable.insert (pair < int, double >(pdg_Delta_1930PP, 1.93));
      masstable.insert (pair < int, double >(pdg_Delta_1930M, 1.93));
      masstable.insert (pair < int, double >(pdg_Delta_1930, 1.93));
      masstable.insert (pair < int, double >(pdg_Delta_1950P, 1.95));
      masstable.insert (pair < int, double >(pdg_Delta_1950PP, 1.95));
      masstable.insert (pair < int, double >(pdg_Delta_1950M, 1.95));
      masstable.insert (pair < int, double >(pdg_Delta_1950, 1.95));
      masstable.insert (pair < int, double >(pdg_Jpsi_1S, 3.09692));
      masstable.insert (pair < int, double >(pdg_K_0_star_1430P, 1.412));
      masstable.insert (pair < int, double >(pdg_K_0_star_1430, 1.412));
      masstable.insert (pair < int, double >(pdg_K_1_1270P, 1.272));
      masstable.insert (pair < int, double >(pdg_K_1_1270, 1.272));
      masstable.insert (pair < int, double >(pdg_K_1_1400P, 1.402));
      masstable.insert (pair < int, double >(pdg_K_1_1400, 1.402));
      masstable.insert (pair < int, double >(pdg_K_2_1770P, 1.773));
      masstable.insert (pair < int, double >(pdg_K_2_1770, 1.773));
      masstable.insert (pair < int, double >(pdg_K_2_1820P, 1.816));
      masstable.insert (pair < int, double >(pdg_K_2_1820, 1.816));
      masstable.insert (pair < int, double >(pdg_K_2_star_1430P, 1.4256));
      masstable.insert (pair < int, double >(pdg_K_2_star_1430, 1.4324));
      masstable.insert (pair < int, double >(pdg_K_3_star_1780P, 1.776));
      masstable.insert (pair < int, double >(pdg_K_3_star_1780, 1.776));
      masstable.insert (pair < int, double >(pdg_K_4_star_2045P, 2.045));
      masstable.insert (pair < int, double >(pdg_K_4_star_2045, 2.045));
      masstable.insert (pair < int, double >(pdg_K_L, 0.497648));
      masstable.insert (pair < int, double >(pdg_K_S, 0.497648));
      masstable.insert (pair < int, double >(pdg_Kstar_1410P, 1.414));
      masstable.insert (pair < int, double >(pdg_Kstar_1410, 1.414));
      masstable.insert (pair < int, double >(pdg_Kstar_1680P, 1.717));
      masstable.insert (pair < int, double >(pdg_Kstar_1680, 1.717));
      masstable.insert (pair < int, double >(pdg_Kstar_892P, 0.89166));
      masstable.insert (pair < int, double >(pdg_Kstar_892, 0.8961));
      masstable.insert (pair < int, double >(pdg_KP, 0.493677));
      masstable.insert (pair < int, double >(pdg_K, 0.497648));
      masstable.insert (pair < int, double >(pdg_Lambda_1405, 1.407));
      masstable.insert (pair < int, double >(pdg_Lambda_1520, 1.5195));
      masstable.insert (pair < int, double >(pdg_Lambda_1600, 1.6));
      masstable.insert (pair < int, double >(pdg_Lambda_1670, 1.67));
      masstable.insert (pair < int, double >(pdg_Lambda_1690, 1.69));
      masstable.insert (pair < int, double >(pdg_Lambda_1800, 1.8));
      masstable.insert (pair < int, double >(pdg_Lambda_1810, 1.81));
      masstable.insert (pair < int, double >(pdg_Lambda_1820, 1.82));
      masstable.insert (pair < int, double >(pdg_Lambda_1830, 1.83));
      masstable.insert (pair < int, double >(pdg_Lambda_1890, 1.89));
      masstable.insert (pair < int, double >(pdg_Lambda_2100, 2.1));
      masstable.insert (pair < int, double >(pdg_Lambda_2110, 2.11));
      masstable.insert (pair < int, double >(pdg_Lambda_b, 5.624));
      masstable.insert (pair < int, double >(pdg_Lambda_c_2593P, 2.5939));
      masstable.insert (pair < int, double >(pdg_Lambda_cP, 2.2849));
      masstable.insert (pair < int, double >(pdg_Lambda, 1.11568));
      masstable.insert (pair < int, double >(pdg_N_1440P, 1.44));
      masstable.insert (pair < int, double >(pdg_N_1440, 1.44));
      masstable.insert (pair < int, double >(pdg_N_1520P, 1.52));
      masstable.insert (pair < int, double >(pdg_N_1520, 1.52));
      masstable.insert (pair < int, double >(pdg_N_1535P, 1.535));
      masstable.insert (pair < int, double >(pdg_N_1535, 1.535));
      masstable.insert (pair < int, double >(pdg_N_1650P, 1.65));
      masstable.insert (pair < int, double >(pdg_N_1650, 1.65));
      masstable.insert (pair < int, double >(pdg_N_1675P, 1.675));
      masstable.insert (pair < int, double >(pdg_N_1675, 1.675));
      masstable.insert (pair < int, double >(pdg_N_1680P, 1.68));
      masstable.insert (pair < int, double >(pdg_N_1680, 1.68));
      masstable.insert (pair < int, double >(pdg_N_1700P, 1.7));
      masstable.insert (pair < int, double >(pdg_N_1700, 1.7));
      masstable.insert (pair < int, double >(pdg_N_1710P, 1.71));
      masstable.insert (pair < int, double >(pdg_N_1710, 1.71));
      masstable.insert (pair < int, double >(pdg_N_1720P, 1.72));
      masstable.insert (pair < int, double >(pdg_N_1720, 1.72));
      masstable.insert (pair < int, double >(pdg_N_2190P, 2.19));
      masstable.insert (pair < int, double >(pdg_N_2190, 2.19));
      masstable.insert (pair < int, double >(pdg_Omega_c, 2.6975));
      masstable.insert (pair < int, double >(pdg_OmegaM, 1.67245));
      masstable.insert (pair < int, double >(pdg_Sigma_1385P, 1.3828));
      masstable.insert (pair < int, double >(pdg_Sigma_1385M, 1.3872));
      masstable.insert (pair < int, double >(pdg_Sigma_1385, 1.3837));
      masstable.insert (pair < int, double >(pdg_Sigma_1660P, 1.66));
      masstable.insert (pair < int, double >(pdg_Sigma_1660M, 1.66));
      masstable.insert (pair < int, double >(pdg_Sigma_1660, 1.66));
      masstable.insert (pair < int, double >(pdg_Sigma_1670P, 1.67));
      masstable.insert (pair < int, double >(pdg_Sigma_1670M, 1.67));
      masstable.insert (pair < int, double >(pdg_Sigma_1670, 1.67));
      masstable.insert (pair < int, double >(pdg_Sigma_1750P, 1.75));
      masstable.insert (pair < int, double >(pdg_Sigma_1750M, 1.75));
      masstable.insert (pair < int, double >(pdg_Sigma_1750, 1.75));
      masstable.insert (pair < int, double >(pdg_Sigma_1775P, 1.775));
      masstable.insert (pair < int, double >(pdg_Sigma_1775M, 1.775));
      masstable.insert (pair < int, double >(pdg_Sigma_1775, 1.775));
      masstable.insert (pair < int, double >(pdg_Sigma_1915P, 1.915));
      masstable.insert (pair < int, double >(pdg_Sigma_1915M, 1.915));
      masstable.insert (pair < int, double >(pdg_Sigma_1915, 1.915));
      masstable.insert (pair < int, double >(pdg_Sigma_1940P, 1.94));
      masstable.insert (pair < int, double >(pdg_Sigma_1940M, 1.94));
      masstable.insert (pair < int, double >(pdg_Sigma_1940, 1.94));
      masstable.insert (pair < int, double >(pdg_Sigma_2030P, 2.03));
      masstable.insert (pair < int, double >(pdg_Sigma_2030M, 2.03));
      masstable.insert (pair < int, double >(pdg_Sigma_2030, 2.03));
      masstable.insert (pair < int, double >(pdg_Sigma_c_2455P, 2.4513));
      masstable.insert (pair < int, double >(pdg_Sigma_c_2455PP, 2.4525));
      masstable.insert (pair < int, double >(pdg_Sigma_c_2455, 2.4522));
      masstable.insert (pair < int, double >(pdg_SigmaP, 1.18937));
      masstable.insert (pair < int, double >(pdg_SigmaM, 1.19745));
      masstable.insert (pair < int, double >(pdg_Sigma, 1.19264));
      masstable.insert (pair < int, double >(pdg_Theta_1540P, 1.5392));
      masstable.insert (pair < int, double >(pdg_Upsilon_10860, 10.865));
      masstable.insert (pair < int, double >(pdg_Upsilon_11020, 11.019));
      masstable.insert (pair < int, double >(pdg_Upsilon_1S, 9.4603));
      masstable.insert (pair < int, double >(pdg_Upsilon_2S, 10.0233));
      masstable.insert (pair < int, double >(pdg_Upsilon_3S, 10.3552));
      masstable.insert (pair < int, double >(pdg_Upsilon_4S, 10.58));
      masstable.insert (pair < int, double >(pdg_WP, 80.42));
      masstable.insert (pair < int, double >(pdg_Xi_1530M, 1.535));
      masstable.insert (pair < int, double >(pdg_Xi_1530, 1.5318));
      masstable.insert (pair < int, double >(pdg_Xi_1820M, 1.823));
      masstable.insert (pair < int, double >(pdg_Xi_1820, 1.823));
      masstable.insert (pair < int, double >(pdg_Xi_c_primP, 2.5741));
      masstable.insert (pair < int, double >(pdg_Xi_c_prim, 2.5788));
      masstable.insert (pair < int, double >(pdg_Xi_cP, 2.4663));
      masstable.insert (pair < int, double >(pdg_Xi_c, 2.4718));
      masstable.insert (pair < int, double >(pdg_XiM, 1.32131));
      masstable.insert (pair < int, double >(pdg_Xi, 1.31483));
      masstable.insert (pair < int, double >(pdg_Z, 91.1876));
      masstable.insert (pair < int, double >(pdg_a_0_1450P, 1.474));
      masstable.insert (pair < int, double >(pdg_a_0_1450, 1.474));
      masstable.insert (pair < int, double >(pdg_a_0_980P, 0.9847));
      masstable.insert (pair < int, double >(pdg_a_0_980, 0.9847));
      masstable.insert (pair < int, double >(pdg_a_1_1260P, 1.23));
      masstable.insert (pair < int, double >(pdg_a_1_1260, 1.23));
      masstable.insert (pair < int, double >(pdg_a_2_1320P, 1.3183));
      masstable.insert (pair < int, double >(pdg_a_2_1320, 1.3183));
      masstable.insert (pair < int, double >(pdg_a_4_2040P, 2.01));
      masstable.insert (pair < int, double >(pdg_a_4_2040, 2.01));
      masstable.insert (pair < int, double >(pdg_b_1_1235P, 1.2295));
      masstable.insert (pair < int, double >(pdg_b_1_1235, 1.2295));
      masstable.insert (pair < int, double >(pdg_chi_b0_1P, 9.8599));
      masstable.insert (pair < int, double >(pdg_chi_b0_2P, 10.2321));
      masstable.insert (pair < int, double >(pdg_chi_b1_1P, 9.8927));
      masstable.insert (pair < int, double >(pdg_chi_b1_2P, 10.2552));
      masstable.insert (pair < int, double >(pdg_chi_b2_1P, 9.9126));
      masstable.insert (pair < int, double >(pdg_chi_b2_2P, 10.2685));
      masstable.insert (pair < int, double >(pdg_chi_c0_1P, 3.41519));
      masstable.insert (pair < int, double >(pdg_chi_c1_1P, 3.51059));
      masstable.insert (pair < int, double >(pdg_chi_c2_1P, 3.55626));
      masstable.insert (pair < int, double >(pdg_e, 0.000510999));
      masstable.insert (pair < int, double >(pdg_etaprim_958, 0.95778));
      masstable.insert (pair < int, double >(pdg_eta_1295, 1.294));
      masstable.insert (pair < int, double >(pdg_eta_1405, 1.4103));
      masstable.insert (pair < int, double >(pdg_eta_1475, 1.476));
      masstable.insert (pair < int, double >(pdg_eta_2_1645, 1.617));
      masstable.insert (pair < int, double >(pdg_eta_c_1S, 2.9796));
      masstable.insert (pair < int, double >(pdg_eta, 0.54775));
      masstable.insert (pair < int, double >(pdg_f_0_1500, 1.507));
      masstable.insert (pair < int, double >(pdg_f_0_1710, 1.714));
      masstable.insert (pair < int, double >(pdg_f_0_980, 0.98));
      masstable.insert (pair < int, double >(pdg_f_1_1285, 1.2818));
      masstable.insert (pair < int, double >(pdg_f_1_1420, 1.4263));
      masstable.insert (pair < int, double >(pdg_f_2_prim_1525, 1.525));
      masstable.insert (pair < int, double >(pdg_f_2_1270, 1.2754));
      masstable.insert (pair < int, double >(pdg_f_2_1950, 1.945));
      masstable.insert (pair < int, double >(pdg_f_2_2010, 2.01));
      masstable.insert (pair < int, double >(pdg_f_2_2300, 2.297));
      masstable.insert (pair < int, double >(pdg_f_2_2340, 2.34));
      masstable.insert (pair < int, double >(pdg_f_4_2050, 2.034));
      masstable.insert (pair < int, double >(pdg_gamma, 0));
      masstable.insert (pair < int, double >(pdg_h_1_1170, 1.17));
      masstable.insert (pair < int, double >(pdg_mu, 0.105658));
      masstable.insert (pair < int, double >(pdg_n, 0.939565));
      masstable.insert (pair < int, double >(pdg_neutron, 0.939565));
      masstable.insert (pair < int, double >(pdg_nu_e, 0));
      masstable.insert (pair < int, double >(pdg_nu_mu, 0));
      masstable.insert (pair < int, double >(pdg_nu_tau, 0));
      masstable.insert (pair < int, double >(pdg_omega_1650, 1.67));
      masstable.insert (pair < int, double >(pdg_omega_3_1670, 1.667));
      masstable.insert (pair < int, double >(pdg_omega_782, 0.78259));
      masstable.insert (pair < int, double >(pdg_pP, 0.938272));
      masstable.insert (pair < int, double >(pdg_proton, 0.938272));
      masstable.insert (pair < int, double >(pdg_phi_1020, 1.01946));
      masstable.insert (pair < int, double >(pdg_phi_1680, 1.68));
      masstable.insert (pair < int, double >(pdg_phi_3_1850, 1.854));
      masstable.insert (pair < int, double >(pdg_pi_1_1400P, 1.376));
      masstable.insert (pair < int, double >(pdg_pi_1_1400, 1.376));
      masstable.insert (pair < int, double >(pdg_pi_1_1600P, 1.596));
      masstable.insert (pair < int, double >(pdg_pi_1_1600, 1.596));
      masstable.insert (pair < int, double >(pdg_pi_1300P, 1.3));
      masstable.insert (pair < int, double >(pdg_pi_1300, 1.3));
      masstable.insert (pair < int, double >(pdg_pi_1800P, 1.812));
      masstable.insert (pair < int, double >(pdg_pi_1800, 1.812));
      masstable.insert (pair < int, double >(pdg_pi_2_1670P, 1.6724));
      masstable.insert (pair < int, double >(pdg_pi_2_1670, 1.6724));
      masstable.insert (pair < int, double >(pdg_piP, 0.13957));
      masstable.insert (pair < int, double >(pdg_pi, 0.134977));
      masstable.insert (pair < int, double >(pdg_psi_2S, 3.68609));
      masstable.insert (pair < int, double >(pdg_psi_3770, 3.77));
      masstable.insert (pair < int, double >(pdg_psi_4040, 4.04));
      masstable.insert (pair < int, double >(pdg_psi_4160, 4.159));
      masstable.insert (pair < int, double >(pdg_psi_4415, 4.415));
      masstable.insert (pair < int, double >(pdg_rho_1450P, 1.465));
      masstable.insert (pair < int, double >(pdg_rho_1450, 1.465));
      masstable.insert (pair < int, double >(pdg_rho_1700P, 1.72));
      masstable.insert (pair < int, double >(pdg_rho_1700, 1.72));
      masstable.insert (pair < int, double >(pdg_rho_3_1690P, 1.6888));
      masstable.insert (pair < int, double >(pdg_rho_3_1690, 1.6888));
      masstable.insert (pair < int, double >(pdg_rho_770P, 0.7758));
      masstable.insert (pair < int, double >(pdg_rho_770, 0.7758));
      masstable.insert (pair < int, double >(pdg_t23, 174));
      masstable.insert (pair < int, double >(pdg_tau, 1.77699));

      return 1;
  }

  double mass (int code)
  {
    static int i = init ();
    map < int, double >::iterator x = masstable.find (abs (code));
    if (x != masstable.end ())
      return x->second * GeV;
    else
      return -1;
  };

  int charge (int code1)
  {
    int sum;
    int code = abs (code1);
    if (code < 20)		//  lepton
      {
	sum = -(code & 1);
      }
    else			//bariony i mezony
      {
	int q1 = (code / 1000) % 10;
	int q2 = (code / 100) % 10;
	int q3 = (code / 10) % 10;

	if (q1 != 0)
	  sum = 2 - (q1 & 1) - (q2 & 1) - (q3 & 1);
	else
	  sum = abs ((q3 & 1) - (q2 & 1));
      }


    return code1 > 0 ? sum : -sum;
  }


};
