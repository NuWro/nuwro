#include "Utilities.h"

double w_peak = 0;
double quasielastic_peak_scaling_factor = 1.0;

void ComputeQuasiElasticPeak(params &p, nucleus &t, F1F209Wrapper &BostedFit) {

  double Z = t.Z();
  double A = t.A();
  double Ein = std::stod(p.beam_energy) / GeV;
  double cos_theta = p.el_costh_lab;
  double theta = std::acos(cos_theta) * 180.0 / Pi;
  static const int Npoints = 500;
  static const int Ml = PDG::mass(p.beam_particle) / GeV;
  std::vector<double> x;
  std::vector<double> fx;

  for (int index = 0; index < Npoints; index++) {
    double width = (Ein - Ml) / Npoints;
    double lowerEdge = width * index;
    double upperEdge = lowerEdge + width;
    double binCenter = 0.5 * (lowerEdge + upperEdge);
    double Eout = binCenter;
    double w = Ein - Eout;

    double xsec_quasielastic =
        BostedFit.GetXS_QE(Z, A, Ein, Eout, theta) * 1e3; // In nb/GeV sr
    double xsec_mec =
        BostedFit.GetXS_MEC(Z, A, Ein, Eout, theta) * 1e3; // In nb/GeV sr

    x.push_back(w);
    fx.push_back(xsec_quasielastic + xsec_mec);
  }

  std::vector<double> f_prime = ComputeDerivative(x, fx);
  w_peak = PeakExtractor(x, fx, f_prime); // <- In GeV
  w_peak = w_peak * GeV;                  // <- In MeV

  // Cross section at peak
  auto itr = std::find(x.begin(), x.end(), w_peak);
  int index = itr - x.begin();
  double Bosted_quasielasticlike_weight = fx.at(index); // In nb/GeV sr
  Bosted_quasielasticlike_weight *= 2 * Pi * Ein;       // In nb
  Bosted_quasielasticlike_weight *=
      1.0; // NOTE: CHeck this later !! for now I multiplied 1

  // Compute scaling factor at peak
  double w = 0;
  double sfevent_weight = 0;
  do {
    event *sampleevent = new event();
    sfevent_weight = sfevent(p, *sampleevent, t);
    w = sampleevent->in[0].E() - sampleevent->out[0].E();
    delete sampleevent;
  } while (std::fabs(w - w_peak) > 5 * MeV);

  quasielastic_peak_scaling_factor =
      Bosted_quasielasticlike_weight / sfevent_weight;

  return;
}

std::vector<double> ComputeDerivative(const std::vector<double> &x,
                                      const std::vector<double> &fx) {
  if (x.size() != fx.size() || x.size() < 3) {
    std::cerr << "size x: " << x.size() << "\n";
    std::cerr << "size fx: " << fx.size() << "\n";
    std::cerr << "Vector must be of same size and size < 3\n";
    std::exit(EXIT_FAILURE);
  }

  size_t n = x.size();
  std::vector<double> f_prime(n, 0.0f);

  /**/
  // 1. Central difference for all internal points
  for (size_t i = 1; i < n - 1; ++i) {
    f_prime[i] = (fx[i + 1] - fx[i - 1]) / (x[i + 1] - x[i - 1]);
  }

  // 2. Forward difference for the first element
  f_prime[0] = (fx[1] - fx[0]) / (x[1] - x[0]);

  // 3. Backward difference for the last element
  f_prime[n - 1] = (fx[n - 1] - fx[n - 2]) / (x[n - 1] - x[n - 2]);
  /**/

  /*
  float h = x[1] - x[0];
  // 1. 6th Order Richardson's Extrapolation for all internal points
  for (size_t i=4; i < n-4; ++i) {
    f_prime[i] = (fx[i+4] - 8.*fx[i+2] + 8.*fx[i-2] - fx[i-4]
    -32*fx[i+2] + 256.fx[i+1] - 256.*fx[i-1] + 32.*fx[i-2]) / (360.0 * h);
  }

  // 2. 4th Order Richardson's Extrapolation for intermediate points
  for (size_t i=2; i<4; ++i) {
    f_prime[i] = (-fx[i+2] + 8.*fx[i+1] - 8.*fx[i-1] + fx[i-2]) / (12.*h);
  }
  for (size_t i=n-4; i<n-2; ++i) {
    f_prime[i] = (-fx[i+2] + 8.*fx[i+1] - 8.*fx[i-1] + fx[i-2]) / (12.*h);
  }

  // 3. Central difference for 2nd and second to last point
  f_prime[1] = (fx[2] - fx[0]) / (2.*h);
  f_prime[n-2] = (fx[n-1] - fx[n-3]) / (2.*h);

  // 4. Forward difference for first point
  f_prime[0] = (fx[1] - fx[0]) / h;

  // 5. Backward difference for the last point
  f_prime[n-1] = (fx[n-1] - fx[n-2]) / h;
  */
  return f_prime;
}

double PeakExtractor(const std::vector<double> &x,
                     const std::vector<double> &fx,
                     const std::vector<double> &f_prime) {
  // 1. Safety check
  if (x.size() != fx.size() || x.size() != f_prime.size()) {
    std::cerr << x.size() << " " << fx.size() << " " << f_prime.size() << "\n";
    std::exit(EXIT_FAILURE);
  }

  // 2. Zero-crossing detection
  for (size_t i = 0; i < f_prime.size() - 1; ++i) {
    // A peak occurs when the slope changes from positive to negative (or zero).
    if (f_prime[i] > 0.0f && f_prime[i + 1] <= 0.0f) {

      // Optional: We can do a quick sanity check to ensure it's higher than the
      // starting baseline to avoid triggering on tiny ripples in the noise.
      if (fx[i] > fx[0]) {
        // If f_prime[i+1] is exactly 0, the peak is at i+1. Otherwise, it's
        // roughly at i. For simplicity with discrete dense data, returning x[i]
        // is usually highly accurate.
        return static_cast<double>(x[i]);
      }
    }
  }

  std::cerr << "No peak found !!\n";
  std::exit(1);
}
