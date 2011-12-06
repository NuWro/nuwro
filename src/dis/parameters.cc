#include <math.h>
#include "jednostki.h"
#include "masses.h"
#include "pdg_name.h"

double kwad (double a)
{
  return a * a;
}

double max2 (double a, double b)
{
  return (a > b ? a : b);
}

double min2 (double a, double b)
{
  return (a < b ? a : b);
}

double x_d2c (double W, double nu)
{
  double Q2 = M2 + 2 * M12 * nu - W * W;
  double Deltadc = kwad (M12 + Dplus_mass) - M2;
  return Q2 / (Q2 + Deltadc * GeV2);
}

double x_s2c (double W, double nu)
{
  double Q2 = M2 + 2 * M12 * nu - W * W;
  double Deltasc = kwad (Dplus_mass + Lambda_mass) - M2;
  return Q2 / (Q2 + Deltasc * GeV2);
}
