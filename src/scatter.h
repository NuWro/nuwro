#ifndef _scatter_h_
#define _scatter_h_

#include "particle.h"

/// Scatter p1,p2 into p[0],..,p[n-1] uniformly in phase space
int scatter_n (int n, particle p1, particle p2, 
		    particle p[]);

/// Scatter p1,p2 into p3,p4  isotropically in cms 
/// if f is given returns f(s,q2) for chosen kinematics
double scatter_2 (particle p1,  particle p2, 
	        particle & p3,  particle & p4,
	        double (*f) (double, double) = NULL);

/// Scatter p1,p2 into p3,p4 occording to distribution 
/// given by f=Ax^3+Bx where x=\cos \theta
bool scatterAB (particle p1, particle p2, 
	   particle & p3, particle & p4,
	   double A, double B, double C, double D, double E, double F, double G, double H);

//added by CThorpe
double scatter2_with_BE_SC(particle p1,particle p2, particle &p3, particle &p4, double Y_Eb);

double scatter2_with_BE(particle p1,particle p2, particle &p3, particle &p4, double Y_Eb, vec &cms_dir);

//bool rescale_momenta(vect pcms, vec cms_dir, particle &p3, particle  &p4, double Y_Eb);
bool rescale_momenta(vect pcms, vec cms_dir, particle &p3, particle  &p4);

#endif

