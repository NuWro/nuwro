#include "shell_sampler.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "jednostki.h"


shell_sampler::shell_sampler(const std::string& filename)
{
  load_data(filename);

  if (data.size() < 2) {
    throw std::runtime_error("not enough data points in shell file");
  }

  // ensure sorted (important for correctness)
  std::sort(data.begin(), data.end(),
            [](const point& a, const point& b) {
              return a.r < b.r;
            });

  build_cdf();
}

double shell_sampler::r()
{
  double u = frandom();

  auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
  if (it == cdf.end())
    return data.back().r;

  size_t i = std::distance(cdf.begin(), it);
  if (i == 0) return data[0].r;

  return 0.5 * (data[i-1].r + data[i].r);
}

double shell_sampler::dens(double r) const
{
  if (r <= data.front().r) return data.front().rho / (norm * 4 * Pi);
  if (r >= data.back().r)  return 0.0;

  auto it = std::lower_bound(data.begin(), data.end(), r,
                 [](const point& p, double val) { return p.r < val; });

  size_t i = std::distance(data.begin(), it);
  double r0 = data[i-1].r,  r1 = data[i].r;
  double f0 = data[i-1].rho, f1 = data[i].rho;
  double t = (r - r0) / (r1 - r0);
  return ((1-t)*f0 + t*f1) / (norm * 4 * Pi);
}

void shell_sampler::load_data(const std::string& filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("unable to open shell data file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);

    point p;
    if (!(iss >> p.r >> p.rho)) {
        throw std::runtime_error("invalid line format in shell data file: " + line);
    }
    p.r = p.r / 197.3;  // convert fm → 1/MeV at load time

    data.push_back(p);
  }
  file.close();
}

void shell_sampler::build_cdf()
{
  size_t n = data.size();
  cdf.resize(n);

  std::vector<double> weight(n);

  double sum = 0.0;

  for (size_t i = 0; i < n; ++i) {
    double r_left, r_right;

    if (i == 0) {
      r_left = std::max(0.0, data[i].r - 0.5 * (data[i+1].r - data[i].r));
    } else {
      r_left  = 0.5 * (data[i].r + data[i-1].r);
    }

    if (i == n - 1) {
      r_right = data[i].r + 0.5 * (data[i].r - data[i-1].r);
    } else {
      r_right = 0.5 * (data[i].r + data[i+1].r);
    }

    double dr = r_right - r_left;

    double w = data[i].rho * data[i].r * data[i].r * dr;

    sum += w;
    weight[i] = sum;
  }

  if (sum <= 0.0) {
    throw std::runtime_error("invalid zero total weight in CDF");
  }

  for (size_t i = 0; i < n; ++i) {
    cdf[i] = weight[i] / sum;
  }

  norm = sum;
}