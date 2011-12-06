#ifndef _kaskada7_h_
#define _kaskada7_h_

#include "event1.h"
#include "params.h"
#include "fsi.h"

int kaskadaevent(params &p,event &e, vector<particle> &input, vector<particle>&output,vector<particle>& all, int nod[], vect q, bool qel, vect p0, int place[][20], bool res, int nofpi, double &pabsen, bool &nopp, bool &abs, int &nofi, double *pen, double &absr, double &odl, double *density_hist, double &radius_hist, int &protons_hist, int &neutrons_hist);

inline int kaskadaevent(params &p,event &e)
{
   return  kaskadaevent(p,e,e.out,e.post,e.all,e.nod,e.q(),e.flag.qel, e.in[1].p4(), e.place, e.flag.res, e.nof(211)+e.nof(-211)+e.nof(111), e.pabsen, e.nopp, e.abs, e.nofi, e.pen, e.absr, e.odl, e.density_hist, e.radius_hist, e.protons_hist, e.neutrons_hist);
}

#endif
