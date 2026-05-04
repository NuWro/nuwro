#include "DNN_event.h"
#include "NeuralNetworks/EmpericalFits.h"
#include "Utilities.h"

void DNN_event(params &p, event &e, nucleus &t, F1F209Wrapper &BostedFit,
               bool nc) {

  event *qelevent = new event();
  event *mecevent = new event();
  event *sppevent = new event();

  // The sfevent returns non-zero weighted only in the selected slice.
  double qel_weight = sfevent(p, *qelevent, t);
  if (!qel_weight) {
    e.weight = 0;
    return;
  }

  sppevent->in[0] = qelevent->in[0];
  sppevent->in[1] =
      qelevent
          ->in[1]; // Momentum of the initial nucleon from spectral function ?!
  sppevent->out[0] = qelevent->out[0];
  double spp_weight = 0;
  //============================================================================//
  double Z = t.Z();                         // Atomic number
  double A = t.A();                         // Mass number
  double Ein = qelevent->in[0].E() / GeV;   // Incoming electron energy in GeV
  double Eout = qelevent->out[0].E() / GeV; // Outgoing electron energy in GeV
  double cos_theta = p.el_costh_lab;        // Cosine angle
  double theta = std::acos(p.el_costh_lab) * 180. / Pi; // Angle in degree
  double w = Ein - Eout; // Energy transfer In GeV
  double Q2 =
      -qelevent->q2() / GeV / GeV; // Square of energy transfer in GeV*GeV
  //============================================================================//

  //=================== P. Bosted & M. Christy fit
  //=============================//
  //=================== F1F209.f
  //===============================================//

  /*
   *
   *      Differential cross section d2sigma / domega dOmega
   *      in units of mb / GeV sr
   *
   *
   */
  double xsec_quasielastic =
      BostedFit.GetXS_QE(Z, A, Ein, Eout, theta); // In mb / GeV sr
  double xsec_mec =
      BostedFit.GetXS_MEC(Z, A, Ein, Eout, theta); // In mb / GeV sr
  xsec_mec = (xsec_mec >= 0) * xsec_mec;           // Asserting xsec_mec > 0
  double xsec_quasielastic_like = xsec_quasielastic + xsec_mec; // In mb/ GeV sr

  //============================================================================//

  std::vector<double> intermediate_weight(3);
  if (w < w_peak) { /*
                       Before quasielastic peak
                    */
    intermediate_weight[0] = xsec_quasielastic_like;
    intermediate_weight[1] = spp_weight;
    intermediate_weight[2] = 0.0;
  } else { /*
              After quasielastic peak
             */
    intermediate_weight[0] = quasielastic_peak_scaling_factor * qel_weight;
    intermediate_weight[1] = spp_weight;
    intermediate_weight[2] = xsec_quasielastic_like - intermediate_weight[0];
  }

  double total_weight = 0.0;
  for (size_t i = 0; i < intermediate_weight.size(); i++) {
    total_weight += intermediate_weight[i];
  }

  std::vector<double> fractions(3);
  for (size_t i = 0; i < fractions.size(); i++) {
    fractions[i] = intermediate_weight[i] / total_weight;
  }

  // ===============  Inference from Neural Networks ======== //
  const double input[5] = {Ein / 20.0, w / 20.0, theta / 180.0, cos_theta,
                           Q2 / 100.0};
  double inclusive_weight_DNN = emperical_fit_ensemble_mean(p, A, input);
  auto DNN_scaling_factor = [](double energy, double angle) {
    double cos_theta_by2 = std::cos(angle * Pi / 360.0);
    double sin_theta_by2 = std::sin(angle * Pi / 360.0);
    double s = std::pow(10, 9) * cos_theta_by2 * cos_theta_by2;
    s = s / (137.0 * 137.0 * energy * cos_theta_by2 * 4.0 * energy * energy *
             std::pow(sin_theta_by2, 4));
    return s;
  };
  inclusive_weight_DNN *=
      DNN_scaling_factor(Ein, theta); // Inclusive cross section
                                      // in nb / GeV sr
  // ========================================================= //

  std::vector<double> new_weight(3);
  for (size_t i = 0; i < fractions.size(); i++) {
    new_weight[i] = fractions[i] * inclusive_weight_DNN;
  }

  switch (e.dyn) {
  case 20:
    for (size_t i = 0; i < qelevent->in.size(); i++) {
      e.in[i] = qelevent->in[i];
    }
    for (size_t i = 0; i < qelevent->out.size(); i++) {
      e.out[i] = qelevent->out[i];
    }
    e.weight = new_weight[0];
    break;
  case 21:
    for (size_t i = 0; i < sppevent->in.size(); i++) {
      e.in[i] = sppevent->in[i];
    }
    for (size_t i = 0; i < sppevent->out.size(); i++) {
      e.out[i] = sppevent->out[i];
    }
    e.weight = new_weight[1];
    break;
  case 22:
    for (size_t i = 0; i < qelevent->in.size(); i++) {
      e.in[i] = qelevent->in[i];
    }
    Generate_outgoing_nucleons(p, e, t);
    e.weight = new_weight[2];
    break;
  }

  return;
}

void Generate_outgoing_nucleons(params &p, event &e, nucleus &t) {}
