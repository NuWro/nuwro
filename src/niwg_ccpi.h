#ifndef niwg_ccpi_h_
#define niwg_ccpi_h_

#include <string>
#include <fstream>
using namespace std;

int ccpi_sim_help(string *tarpar, string en);
int ccpi_sim();
double pifact(string file, double fac[][2]);
int ccpi_calc_help (ofstream files[][3][2], string en);
int ccpi_calc();
int ccpi_plot();

#endif
