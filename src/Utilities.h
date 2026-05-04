#ifndef UTILITIES_H
#define UTILITIES_H

#include "../data/F1F209/F1F209Wrapper.hh"
#include "event1.h"
#include "nucleus.h"
#include "params.h"
#include "pdg.h"
#include "sfevent.h"

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>

void ComputeQuasiElasticPeak(params &p, nucleus &t, F1F209Wrapper &BostedFit);
std::vector<double> ComputeDerivative(const std::vector<double> &x,
                                      const std::vector<double> &fx);
double PeakExtractor(const std::vector<double> &x,
                     const std::vector<double> &fx,
                     const std::vector<double> &f_prime);

extern double w_peak;
extern double quasielastic_peak_scaling_factor;

#endif
