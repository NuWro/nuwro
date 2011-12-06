#ifndef _dis_cr_sec_h_
#define _dis_cr_sec_h_


////////////////////////////////decomposition //////////////////////////

double hadr();
double prob_fey();
//cross section for nukleons with bodek corrections
//double cr_sec_dis_temp(double E, double W, double nu, int lepton, int nukleon, bool current)

double cr_sec_dis(double E, double W, double nu, int lepton, int nukleon, bool current);
///////////////////////////////only grv////////////////////////////////
//double cr_sec_dis(double E, double W, double nu, int lepton, int nukleon, bool current)
double cr_sec_dis_grv(double E, double W, double nu, int lepton, int nukleon, bool current);


#endif
