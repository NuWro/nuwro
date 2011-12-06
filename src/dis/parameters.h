#ifndef _parameters_h_
#define _parameters_h_

#include <math.h>
#include "jednostki.h"


const double mpi=138;//masa pionu 
const double m_pi=mpi;//masa pionu 
const double mpi2=mpi*mpi;
const double pi=M_PI; 
const double pi2=pi*pi;
const double MD=1232;
const double MD2=MD*MD;
const double GD0=115;
const double Mv2=710000;
const double Ma=1030;
const double Ma2=Ma*Ma;
const double cos_theta_C = 0.975;
const double cos_2_theta_C = cos_theta_C*cos_theta_C;
const double sin_theta_C = sqrt(1-cos_2_theta_C);
const double sin_2_theta_C = sin_theta_C*sin_theta_C;
const double sin_2_theta_W = 0.23120;

const double stala = G*G;

double kwad (double a);
double max2(double a, double b);
double min2(double a, double b);

double x_d2c(double W, double nu);
double x_s2c(double W, double nu);

#endif
