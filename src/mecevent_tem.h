#ifndef _mecevent_tem_h_
#define _mecevent_tem_h_

#include <cassert>
#include "particle.h"
#include "jednostki.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "params.h"
#include "beam.h"
#include "pdg.h"
#include "nucleus.h"
#include "event1.h"
#include "qel_sigma.h"

const double M1 = (PDG::mass_proton);// + PDG::mass_neutron) / 2.0;
const double M2 = M1 * M1;
static double Q2;
static double dif;

inline bool blocked (double E, double m)		//determines a lower bound of the neutrino energy
{
	if (2.0 * M1 * E > m * (2.0 * M1 + m)) //(W2-mecm2-mecM2)*(W2-mecm2-mecM2) - 4.0*mecm2*mecM2 > 0
		return false;
	else
		return true;
}

inline double setQ2 (double E, double m2)
{
	double W2 = 2.0 * M1 * E + M2;	//s variable
	double W = sqrt (W2);
	double E_cmf = (W2 - M2) / 2.0 / W;
	double Eprim_cmf = (W2 + m2 - M2) / 2.0 / W;
	double kprim_cmf = sqrt( (W2 - m2 - M2) * (W2 - m2 - M2) - 4.0 * m2 * M2 ) / 2.0 / W;
	
	double Q2min = 2.0 * E_cmf * Eprim_cmf - m2 - 2.0 * E_cmf * kprim_cmf;
	double Q2max = 2.0 * E_cmf * Eprim_cmf - m2 + 2.0 * E_cmf * kprim_cmf;
	
	dif = Q2max - Q2min;
	
	//double Q2mintrue = max(Q2min, 4.0 * 2.0 * M1);
	
	return Q2min + (Q2max - Q2min) * frandom();
}

void tem_kin (double E, particle &meclep, particle *nucleon, nucleus &t, double central, double smearing, 
	      double binding, bool &kinematics, bool czy_pb, int ile_pb);
double mec_do_cc (double &w, double E, particle *p, double m, double ratio, bool nu);
double mec_do_nc (double &w, double E, particle *p, double m, double ratio, bool nu);

#endif
