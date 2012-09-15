#ifndef _simalations_h_
#define _simalations_h_
#include "fsi.h"

const string events5m = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 5000000' ";
const string events1m = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 1000000' ";
const string events100k = "-p 'number_of_test_events = 1000000' -p 'number_of_events = 100000' ";
const string events500k = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 500000' ";
const string events10k = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 10000' ";
const string events0k = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 10' ";
const string events25k = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 25000' ";
const string evtest    = "-p 'number_of_test_events = 5000' -p 'number_of_events = 25000' ";

const string NCdyn      = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 1' ";
const string NCcoh      = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 1' ";
const string NCwocoh    = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string CCdyn      = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 1' -p 'dyn_coh_nc = 0' ";
const string ALLdyn     = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 1' -p 'dyn_coh_nc = 1' ";
const string CCwocoh    = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string QEL        = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string onlydis    = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string wocoh      = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string resdis     = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string CCresdis   = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string onlyres    = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string NCQEL      = "-p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string CCQEL      = "-p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";

const string numu       = "-p 'beam_particle = 14' ";
const string antinumu   = "-p 'beam_particle = -14' ";
const string piP        = "-p 'beam_particle = 211' ";
const string piM        = "-p 'beam_particle = -211' ";
const string pi0        = "-p 'beam_particle = 111' ";
const string proton     = "-p 'beam_particle = 2212' ";

const string pos_ran    = "-p 'beam_placement = 1' ";
const string surface    = "-p 'beam_placement = 2' ";

const string K2Kbeam    = "-p 'beam_type = 0' -p 'beam_energy = 0 4000 0.102189781022 0.953771289538 1.95133819951 2.70559610706 3.92214111922 5.6496350365 7.54744525547 9.22627737226 10.5158150852 11.4647201946 12.1946472019 12.3163017032 11.9270072993 11.1727493917 10.102189781 8.83698296837 7.49878345499 6.20924574209 5.0900243309 4.04379562044 3.21654501217 2.51094890511 1.95133819951 1.48905109489 1.12408759124 0.856447688564 0.661800486618 0.491484184915 0.394160583942 0.296836982968 0.199513381995 0.150851581509 0.126520681265 0.102189781022 0.102189781022 0.102189781022 0.0778588807786 0.0778588807786 0.0778588807786 0.029197080292' ";
const string MBbeam     = "-p 'beam_type = 0' -p 'beam_energy = 0 7200 2.272e-12 8.566e-12 1.112e-11 1.335e-11 1.658e-11 1.82e-11 1.946e-11 2.045e-11 2.161e-11 2.241e-11 2.279e-11 2.292e-11 2.275e-11 2.253e-11 2.214e-11 2.156e-11 2.078e-11 1.992e-11 1.894e-11 1.789e-11 1.677e-11 1.558e-11 1.439e-11 1.318e-11 1.193e-11 1.069e-11 9.503e-12 8.356e-12 7.278e-12 6.292e-12 5.396e-12 4.601e-12 3.902e-12 3.285e-12 2.76e-12 2.312e-12 1.932e-12 1.616e-12 1.355e-12 1.138e-12 9.589e-13 8.15e-13 6.928e-13 5.937e-13 5.147e-13 4.478e-13 3.935e-13 3.5e-13 3.15e-13 2.867e-13 2.615e-13 2.409e-13 2.273e-13 2.11e-13 1.995e-13 1.92e-13 1.815e-13 1.726e-13 1.665e-13 1.601e-13 1.554e-13 1.493e-13 1.442e-13 1.412e-13 1.363e-13 1.323e-13 1.265e-13 1.217e-13 1.183e-13 1.14e-13 1.102e-13 1.06e-13 1.014e-13 9.7e-14 9.34e-14 9.001e-14 8.641e-14 8.19e-14 7.867e-14 7.464e-14 7.146e-14 6.812e-14 6.499e-14 6.185e-14 5.858e-14 5.614e-14 5.32e-14 5.016e-14 4.765e-14 4.561e-14 4.281e-14 4.087e-14 3.841e-14 3.632e-14 3.432e-14 3.263e-14 3.016e-14 2.857e-14 2.689e-14 2.529e-14 2.372e-14 2.227e-14 2.103e-14 1.957e-14 1.834e-14 1.73e-14 1.615e-14 1.513e-14 1.406e-14 1.303e-14 1.214e-14 1.129e-14 1.047e-14 9.569e-15 8.87e-15 8.148e-15 7.429e-15 6.765e-15 6.097e-15 5.492e-15 4.977e-15 4.445e-15 3.967e-15 3.492e-15 3.037e-15 2.595e-15 2.225e-15 1.854e-15 1.537e-15 1.22e-15 9.78e-16 7.842e-16 6.198e-16 4.786e-16 3.334e-16 1.971e-16 9.391e-17 2.738e-17 6.065e-18 4.135e-18 1.933e-18 9.888e-19 4.494e-20 0' ";
const string MBbeamanti = "-p 'beam_type = 0' -p 'beam_energy = 0 7600 2.157e-12 7.84e-12 9.731e-12 1.141e-11 1.319e-11 1.438e-11 1.477e-11 1.479e-11 1.5e-11 1.485e-11 1.447e-11 1.406e-11 1.345e-11 1.287e-11 1.221e-11 1.152e-11 1.075e-11 9.98e-12 9.177e-12 8.411e-12 7.658e-12 6.907e-12 6.18e-12 5.505e-12 4.877e-12 4.269e-12 3.686e-12 3.151e-12 2.678e-12 2.262e-12 1.898e-12 1.58e-12 1.311e-12 1.083e-12 8.917e-13 7.285e-13 5.941e-13 4.834e-13 3.937e-13 3.18e-13 2.577e-13 2.066e-13 1.665e-13 1.346e-13 1.081e-13 8.837e-14 7.136e-14 5.707e-14 4.62e-14 3.778e-14 3.028e-14 2.412e-14 1.977e-14 1.638e-14 1.323e-14 1.038e-14 8.707e-15 6.981e-15 6.078e-15 5.111e-15 3.919e-15 3.328e-15 2.861e-15 2.382e-15 2.295e-15 2.269e-15 1.828e-15 1.613e-15 1.537e-15 1.375e-15 1.247e-15 1.044e-15 9.532e-16 8.148e-16 7.485e-16 7.837e-16 6.202e-16 6.614e-16 5.495e-16 5.132e-16 4.9e-16 4.984e-16 4.479e-16 2.406e-16 2.196e-16 1.973e-16 1.762e-16 1.339e-16 1.256e-16 9.669e-17 7.546e-17 5.98e-17 6.291e-17 4.74e-17 3.581e-17 3.289e-17 2.808e-17 2.706e-17 2.439e-17 1.477e-17 1.24e-17 1.035e-17 9.001e-18 6.87e-18 8.575e-18 4.392e-18 4.736e-18 3.038e-18 2.413e-18 2.072e-18 2.061e-18 1.113e-18 9.256e-19 6.533e-19 8.079e-19 6.662e-19 3.475e-19 3.808e-19 3.495e-19 2.521e-19 1.967e-19 1.264e-19 6.892e-20 7.387e-20 3.467e-20 1.026e-19 1.603e-20 1.62e-20 1.13e-20 1.395e-20 2.114e-21 3.943e-21 1.888e-21 6.712e-21 5.944e-21 0 1.436e-21 0 0 0 0 0 0 1.436e-21 0 0 0 0 0 0 0 0' ";
const string E1         = "-p 'beam_type = 0' -p 'beam_energy = 1000' ";
const string PNSbeam    = "-p 'beam_type = 0' -p 'beam_energy = 150 700' ";
const string PrThe      = "-p 'beam_type = 0' -p 'beam_energy = 1000 5000' ";
const string Pr2        = "-p 'beam_type = 0' -p 'beam_energy = 1000 2300' ";
const string PrTle      = "-p 'beam_type = 0' -p 'beam_energy = 940 1400' ";
const string PiThe      = "-p 'beam_type = 0' -p 'beam_energy = 2500 5000' ";
const string PiTle      = "-p 'beam_type = 0' -p 'beam_energy = 150 700' ";
const string Nomadbeam  = "-p 'beam_type = 0' -p 'beam_energy = 0 240000 124849.136271 378715.650687 920143.201528 1319720.39306 1283611.90082 1086787.35912 846663.511319 659594.625742 513858.297294 400322.166666 311871.654049 256825.743597 205708.866932 155873.068275 131972.039306 108678.735912 92014.3201528 77905.1673922 67814.9274083 59031.5704765 54317.4982495 48612.3916248 43506.5070308 40032.2166666 36835.3720079 33893.8171288 32064.4722602 30333.8622916 27911.4971638 25682.5743597 23631.6461948 21744.4985871 19460.6179421 18410.2735256 16476.5951193 15160.8269833 13950.1318782 12484.9136271 11487.9083956 10570.52081 9460.2718066 8704.80478745 7790.51673922 7168.39101399 6595.94625742 6239.94444169 5741.64245594 5138.58297294 4728.23218804 4473.03616636 3893.69066354 3683.53720079 3296.64580918 3033.38622916 2869.66582223 2714.78186724 2429.64134136 2235.6179624 2114.95535801 1946.06179421 1790.65552969 1602.57837124 1474.60159321 1356.84463096 1283.61190082 1181.10676394 1086.78735912 1000 946.02718066 894.967426547 779.051673922 737.004058669 697.225871758 641.547645827 574.164245594 528.313327144 486.123916248 447.303616636 435.065070308 378.715650687' ";
const string Atmbeam    = "-p 'beam_type = 0' -p 'beam_energy = 400 6200 1.22e+11 1.41e+11 1.16e+11 8.89e+10 6.27e+10 4.00e+10 2.49e+10 1.38e+10 7.40e+09 4.21e+09 2.72e+09 1.93e+09 1.44e+09 1.13e+09 9.34e+08 7.95e+08 7.07e+08 6.50e+08 6.18e+08 5.91e+08 5.47e+08 5.13e+08 4.52e+08 4.33e+08 3.89e+08 3.46e+08 3.02e+08 2.59e+08 2.17e+08' ";
const string pitbeam    = "-p 'beam_type = 0' -p 'beam_energy = 5000' ";
const string t2knu      = "-p 'beam_type = 0' -p 'beam_energy = 50 9450 38571.4 82467.5 134286 208016 291786 389464 496786 601825 747946 945238 1.09036e+06 1.13119e+06 1.07349e+06 915000 676786 455000 293571 194048 141464 110714 91904.8 78285.7 68015.9 60918.4 53954.1 48452.4 43035.7 39365.1 35754 33992.3 32251.3 30510.2 28769.1 27028.1 25287 24085.1 22989.6 21894.1 20798.6 19703 18607.5 17512 16416.5 15321 14266.5 13917.1 13567.6 13218.1 12868.7 12519.2 12169.7 11820.3 11470.8 11121.4 10771.9 10507.9 10260.8 10013.6 9766.52 9519.39 9272.26 9025.14 8778.01 8530.88 8308.78 8206.78 8104.78 8002.77 7900.77 7798.77 7696.77 7594.77 7492.77 7390.76 7288.76 7196.35 7105.94 7015.53 6925.12 6834.71 6744.3 6653.89 6563.47 6473.06 6382.65 6245.64 6098.92 5952.2 5805.48 5658.75 5512.03 5365.31 5218.59 5071.87 4925.15 4787.18 4651.05 4514.91 4378.77 4242.63 4106.49 3970.35 3834.22 3698.08 3561.94 3425.8 3289.66 3153.53 3017.39 2881.25 2767.77 2704.64 2641.51 2578.38 2515.25 2452.12 2389 2325.87 2262.74 2199.61 2136.48 2073.35 2010.22 1947.09 1883.97 1820.84 1757.71 1694.58 1631.45 1568.32 1524.02 1489.62 1455.22 1420.82 1386.42 1352.02 1317.63 1283.23 1248.83 1214.43 1180.03 1145.63 1111.24 1076.84 1042.44 1008.04 973.643 939.245 904.847 870.448 847.568 830.749 813.93 797.111 780.292 763.473 746.654 729.835 713.016 696.197 679.377 662.558 645.739 628.92 612.101 595.282 578.463 561.644 544.825 528.006 514.286 502.196 490.106 478.017 465.927 453.837 441.748 429.658 417.568 405.478 393.389 381.299 369.209 357.12 345.03 332.94 320.85 308.761 296.671' ";
const string to1500     = "-p 'beam_type = 0' -p 'beam_energy = 0 1500 ' ";
const string E40k     = "-p 'beam_type = 0' -p 'beam_energy = 40000 ' ";

const string hydrogen   = "-p 'target_type = 0' -p 'nucleus_p = 1' -p 'nucleus_n = 0' -p 'nucleus_target = 0' -p 'kaskada_on = 0' -p 'pauli_blocking = 0' -p 'sf_method = 0' ";
const string oxygen     = "-p 'target_type = 0' -p 'nucleus_p = 8' -p 'nucleus_n = 8' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string carbon     = "-p 'target_type = 0' -p 'nucleus_p = 6' -p 'nucleus_n = 6' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string calcium    = "-p 'target_type = 0' -p 'nucleus_p = 22' -p 'nucleus_n = 26' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string iron       = "-p 'target_type = 0' -p 'nucleus_p = 26' -p 'nucleus_n = 30' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string lithium    = "-p 'target_type = 0' -p 'nucleus_p = 3' -p 'nucleus_n = 4' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string aluminium  = "-p 'target_type = 0' -p 'nucleus_p = 13' -p 'nucleus_n = 14' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string copper     = "-p 'target_type = 0' -p 'nucleus_p = 29' -p 'nucleus_n = 34' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string neon       = "-p 'target_type = 0' -p 'nucleus_p = 10' -p 'nucleus_n = 10' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string barium     = "-p 'target_type = 0' -p 'nucleus_p = 56' -p 'nucleus_n = 82' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' -p 'sf_method = 0' ";
const string fneutron   = "-p 'target_type = 0' -p 'nucleus_p = 0' -p 'nucleus_n = 1' -p 'kaskada_on = 0' -p 'pauli_blocking = 0' -p 'sf_method = 0' ";
const string nomdet     = "-p 'target_type = 0' -p 'target_content = 1 0 514x 0 0 0' -p 'target_content += 6 6 6430x 34 220 1' -p 'target_content += 7 7 592x 34 220 1' -p 'target_content += 8 8 2213x 34 220 1' -p 'target_content += 13 14 171x 34 220 1' -p 'target_content += 14 14 27x 34 220 1' -p 'target_content += 17 18 30x 34 220 1' -p 'target_content += 18 18 19x 34 220 1' -p 'target_content += 29 34 3x 34 220 1' -p 'kaskada_on = 1' -p 'pauli_blocking = 1' ";

const string fsioff     = "-p 'kaskada_on = 0' ";
const string sf         = "-p 'sf_method = 1' ";//-p 'nucleus_target = 4' ";

std::string fzp (int fz);
std::string xpar (int xs);
void simlog(string command);
void nuint_sim();
void krenz_sim ();
void ft_sim ();
void ks();
void PiT2le (int fz, int xs);
void K2K (int fz, int xs);
void MB (int fz, int xsprob, bool anti);
void simMBback (int fz, int xsprob, bool anti);
void MBCC (int fz, int xsprob);
void simMBCCtotal(int fz, int xs);
void simMBCCratio(int fz, int xs);
void simSBCCtotal(int fz, int xs);
void simNOMADCCtotal(int fz, int xs);
void simsfg();
void simMINOSCCtotal(int fz, int xs);
void PNS (int fz, int xs);
void PrT (int fz, int xs, char *which);
void PiT2 (int fz, int xs);
void PiT3(int fz);
void PiT (int fz, int xs, char *which);
void Nomad (int fz, int xs);
void AtmNu (int fz, int xs);
void Multiplicity (int fz, int xs);
void roman_sim ();
void oset_sim ();
void t2k_anal_sim();
void hayato_sim1();
void hayato_sim2();
void hayato_sim3();
void hayato_sim4();
void dens_test_sim();
void ptsim();
void angle_test();
void kendall_sim(string pdg, string p, string n, string pf, string eb);
void xsec_sim();
void xsec_sim2();
void xsec_sim3();
void xsec_sim4();
void towork_sim();
void ccpip_js_sim();
void hayato_sim0812();
#endif
