#ifndef _plots_h_
#define _plots_h_

#include <string>

using namespace std;

void make_data(const double *x, const double *y, const double *xerr, const double *yerr, const int bins, bool rest);
void plotlog(string text);
void findratio(string file, double &ratio, double &ratiofsi);
void findnof(string file, double *back);
int plotK2K(int fz, int xs);
int plotMBangle(int fz, int xs);
int plotMB(int fz, int xs);
int plotMBCC(int fz, int xs);
int plotMBCCtotal(int fz, int xs);
int plotSB(int fz, int xs);
int plotPNS(int fz, int xs);
int plotPrThe(int fz, int xs);
int plotPrTle(int fz, int xs);
int plotPiTle(int fz, int xs);
void plotFZ ();
int plotNomad(int fz, int xs);
int plotPiThe(int fz, int xs);
int plotAtmNu(int fz, int xs);
int plotMBCCrat(int fz, int xs);
int plotnuintPr(int fz, int xs);
void t2k_anal_plot();
#endif
