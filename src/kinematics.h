#ifndef _kinematics_h_
#define _kinematics_h_

double V(double,double); // Momentum dependent potential

double jakob(vect nu,vect N0,vect lepton);

double bodek_kinematics(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double & jakobian);

double CT_bodek_kinematics(particle N0, int Z, int N);

double bodek_binding_energy(particle N0, int Z, int N);

double czarek_kinematics(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double &jakobian);

double czarek_kinematics2(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double &jakobian);

double momentum_dependent_potential_kinematics(particle nu, particle N0, particle &lepton, particle &N1, double &jakobian);

#endif
