#ifndef EMPERICAL_FITS_H
#define EMPERICAL_FITS_H

class params;

/// Ensemble mean over `NumberOfNetworks` ONNX models (see EmpericalFits.cc).
/// Returns 1.0 when `NeuralNetwork_guided` is false, nucleus is unsupported, or on load/inference failure.
double emperical_fit_ensemble_mean(const params &p, int mass_number_A,
                                   const double inputs5[5]);

/// Drop cached ensemble (e.g. after changing model directory at runtime).
void emperical_fit_reset_cache();

#endif
