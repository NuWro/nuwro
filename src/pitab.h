#ifndef _pitab_h_
#define _pitab_h_

void maketable(double tab1[][26], double tab2[][26], double result[][26], double a);
void make_app(double Ek, int k, double *res, double pdataHE[][29]);
double *pion_params (double Ek, int xsec, double dens);

#endif
