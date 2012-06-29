#ifndef _fsi_h_
#define _fsi_h_

#include <string>

#include "particle.h"
#include "params.h"
#include "event1.h"

using namespace std;

const int nof_expr = 19;
const int nof_fz = 10;
const int nof_opt = 5;
const int nof_xsec = 2;
const string expr[nof_expr] = {"K2K", "MiniBooNE", "SciBooNE", "Pion-Nucleus Scattering", "Proton Transparency (high energy)", "Proton Transparency (low energy)", "Pion Transparency (high energy)", "Pion Transparency (low energy)", "Nomad", "Atmospheric Neutrinos", "Multiplicity", "MiniBooNE (CC)", "MiniBooNE (CC total)", "SciBooNE (CC total)", "NOMAD (CC total)", "MINOS (CC total)", "MB CCpi+ to CCQE", "NUINT11 proton tr", "MiniBooNE back"};
const string fzname[nof_fz] = {"Without formation zone", "Skat ({/Symbol m}^{2} = 0.08 GeV^{2})", "Stodolsky", "Ranft ({/Symbol t} = 0.342fm)", "Ranft-like ({/Symbol t} = 0.342fm)", "Cosyn", "Delta", "Coherence length", "Formation zone", "pitrans"};
const string fzwork[nof_fz] = {"nofz", "skat8", "stod", "ranft", "rl", "cosyn", "delta", "cohl", "fz", "trans"};
const string options[nof_opt] = {"Simulation", "Calculation", "Plot", "Multiplot", "Vivisection"};
const string xsec[nof_xsec] = {"Metropolis", "Oset"};

extern bool expr_on[nof_expr];
extern bool fz_on[nof_fz];
extern bool options_on[nof_opt];
extern bool xsec_on[nof_xsec];
extern bool allinone_on;
extern string allinone;

extern char day[10];
extern char month[10];
extern char fullm[10];
extern char hour[10];
extern char minute[10];
extern char year[10];
extern string date;

const string dot   = ".";
const string sep   = "_";
const string rot   = ".root ";
const string space = " ";

extern ofstream logfile;

using namespace std;
int getch (void);

void get_date();
void run(string com);
bool noFile(string filename);
string find_last (string name);
double formation_zone (particle &p, params &par, event &e);
double formation_zone (particle &p, params &par);
#endif
