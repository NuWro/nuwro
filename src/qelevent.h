#include <cassert>          
#include "params.h"
#include "particle.h"
#include "event1.h"
#include "nucleus.h"

using namespace NUWRO;
//double qelevent(params &p,event &e,nucleus &t,bool);     
double qelevent1(params &p,event &e,nucleus &t,bool);     

//double qelevent2(params &p,event &e,target &t,bool);     

double momentum_dependent_potential_kinematics(particle nu, particle N0, particle &lepton, particle &N1, double & jakobian);

double czarek_kinematics(double Eb,particle nu,particle N0, particle &lepton, particle &N1,double &jakobian);

