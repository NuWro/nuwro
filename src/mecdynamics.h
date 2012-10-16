#ifndef _mecdynamics_h_
#define _mecdynamics_h_

double mecweight (double cohE, bool cohprocess, int cohA, particle &lepton, particle &nuc1, particle &nuc2, bool fsi, double pot);
inline double los(){ return frandom();}
#endif


