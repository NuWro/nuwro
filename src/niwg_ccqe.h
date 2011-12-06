#ifndef _niwg_ccqe_h_
#define _niwg_ccqe_h_
#include <string>
#include <fstream>
using namespace std;

int ccqe_sim_help(string *tarpar, std::string en);
int ccqe_sim();
int ccqe_calc_help (ofstream *files, string en);
int ccqe_sfg_ratio();
int ccqe_calc();
int ccqe_sfg_plot();
int ccqe_plot();

#endif
