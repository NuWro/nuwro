#include "event1.h"
#include "particle.h"
#include "mb_nce_fit.h"

#include <TFile.h>
#include <TTree.h>

#include <string>

using namespace std;

const int nof_ma = 16;
const int ma[nof_ma] = {1100, 1200, 1300, 1400, 1410, 1420, 1430, 1440, 1450, 1460, 1470, 1480, 1490, 1500, 1600, 1700};

const int nof_ds = 2;//3;
const double xds[nof_ds] = {0, 0.5};

const int events = 5000000;

const string parQEL = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 5000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 1' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string parRES = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 5000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 0' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 1' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 1' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string parCC = "-p 'number_of_test_events = 10000000' -p 'number_of_events = 5000000' -p '@data/beam/newMB.txt' -p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string H = "-p '@data/target/H.txt' ";
const string C = "-p '@data/target/C.txt' ";
const string fg = "-p 'sf_method = 0' ";
const string sf = "-p 'sf_method = 1' ";
const string ccma = "-p 'qel_cc_axial_mass = 1350' ";
const string ccmaH = "-p 'qel_cc_axial_mass = 1030' ";

void sim();
void re_sim();
void calc();
void calcH(string in, double *res);
void calcC(string in, double res[5][51]);

void delta_s_sim(double xMa);
void delta_s_calc();
void calcH_he(string in, double *res);
void calcC_he(string in, double res[5][30]);

void chi2();
double calc_chi(string in);

void q2_calc(string rootC);
void q2_calc(string in, double *res);

void sim_cc();
void ratio_calc(string rootCnc);
void ratio_calcH(string in, double *res);
void ratio_calcC(string in, double *res, double *reslike);

void norm(double *tab, double x);
void normHE(double *tab, double x);

void tm_sim(double ma);
void tm_calc(double ma);
double tm_chi(double ma);

void run(string com); //run external program "com"
double crosssection (string filename); //read xsec from *.root.txt
string ncma(double val); //create input string from double for Ma
string wf(double val); //create input string from double for work function
string ds(double val); //create input string from double for delta_s
bool noFile(string filename); //returns true if a file does not exist
