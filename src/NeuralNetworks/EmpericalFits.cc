#include "EmpericalFits.h"
#include "NN.h"

#include "../dirs.h"
#include "../params.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <glob.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

std::string trim_slash(std::string s) {
  while (!s.empty() && (s.back() == '/' || s.back() == '\\'))
    s.pop_back();
  return s;
}

std::string default_model_root(const params &p) {
  if (!p.NeuralNetwork_model_dir.empty())
    return trim_slash(p.NeuralNetwork_model_dir);
  return trim_slash(get_data_dir() + "NeuralNetworks/EmpericalFits");
}

std::string nuclear_subdir(int A) {
  if (A == 12)
    return "C12";
  if (A == 56)
    return "Fe56";
  return "";
}

std::string file_basename(const std::string &path) {
  const auto p = path.find_last_of("/\\");
  if (p == std::string::npos)
    return path;
  return path.substr(p + 1);
}

/// Ensemble index from `model_electron_A_12_7.onnx`; exclude `_ver_` names.
int ensemble_index_from_basename(const std::string &base) {
  if (base.find("_ver_") != std::string::npos)
    return -1;
  const auto dot = base.rfind('.');
  const auto usc = base.rfind('_');
  if (dot == std::string::npos || usc == std::string::npos || usc >= dot)
    return -1;
  const std::string num = base.substr(usc + 1, dot - usc - 1);
  char *end = nullptr;
  const long v = std::strtol(num.c_str(), &end, 10);
  if (!end || end != num.c_str() + num.size() || v < 0)
    return -1;
  return static_cast<int>(v);
}

std::vector<std::string> discover_sorted_onnx(const std::string &abs_dir,
                                              int max_networks) {
  const std::string pattern = abs_dir + "/model_electron_*_*.onnx";
  glob_t g;
  memset(&g, 0, sizeof(g));
  if (glob(pattern.c_str(), 0, nullptr, &g) != 0) {
    globfree(&g);
    return {};
  }
  std::vector<std::pair<int, std::string>> items;
  for (size_t i = 0; i < g.gl_pathc; ++i) {
    const std::string path = g.gl_pathv[i];
    const int idx = ensemble_index_from_basename(file_basename(path));
    if (idx < 0)
      continue;
    items.push_back({idx, path});
  }
  globfree(&g);
  std::sort(items.begin(), items.end(),
            [](const std::pair<int, std::string> &a,
               const std::pair<int, std::string> &b) { return a.first < b.first; });
  std::vector<std::string> out;
  int last = -1;
  for (const auto &pr : items) {
    if (pr.first == last)
      continue;
    last = pr.first;
    out.push_back(pr.second);
    if (static_cast<int>(out.size()) >= max_networks)
      break;
  }
  return out;
}

class EmpiricalEnsemble {
public:
  explicit EmpiricalEnsemble(const params &p, const std::string &subdir) {
    const std::string root = default_model_root(p);
    const std::string dir = root + "/" + subdir;
    const std::vector<std::string> paths =
        discover_sorted_onnx(dir, p.NumberOfNetworks);
    if (paths.empty()) {
      std::cerr << "[EmpericalFits] no ONNX models under " << dir << std::endl;
      return;
    }
    runtime_ = std::make_shared<OrtRuntime>("NuWroEmpericalFits");
    nets_.reserve(paths.size());
    for (const std::string &path : paths) {
      try {
        nets_.push_back(
            std::make_unique<NeuralNetwork>(runtime_, path, p, true));
      } catch (const std::exception &ex) {
        std::cerr << "[EmpericalFits] skip " << path << ": " << ex.what()
                  << std::endl;
      }
    }
    std::cerr << "[EmpericalFits] loaded " << nets_.size() << " / "
              << paths.size() << " sessions from " << dir << std::endl;
  }

  bool ok() const { return !nets_.empty(); }

  double predict_mean(const params &p, const double x[5]) const {
    if (nets_.empty())
      return 1.0;
    const int cap = p.NumberOfCores > 0 ? p.NumberOfCores : 1;
    const int workers = std::max(1, std::min(cap, static_cast<int>(nets_.size())));
    std::vector<double> partial(nets_.size());
    std::atomic<size_t> next{0};
    auto worker = [&]() {
      for (;;) {
        const size_t i = next.fetch_add(1);
        if (i >= nets_.size())
          return;
        partial[i] =
            nets_[i]->run_scalar_io(x, "input", "Identity:0");
      }
    };
    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(workers));
    for (int t = 0; t < workers; ++t)
      threads.emplace_back(worker);
    for (auto &th : threads)
      th.join();
    double sum = 0;
    for (double v : partial)
      sum += v;
    return sum / static_cast<double>(partial.size());
  }

private:
  std::shared_ptr<OrtRuntime> runtime_;
  std::vector<std::unique_ptr<NeuralNetwork>> nets_;
};

std::mutex g_cache_mu;
std::unique_ptr<EmpiricalEnsemble> g_ensemble;
int g_ensemble_A = -1;
std::string g_ensemble_subdir;
std::string g_ensemble_cache_key;

static std::string ensemble_cache_key(const params &p, const std::string &sub) {
  return default_model_root(p) + "/" + sub + "/" +
         std::to_string(p.NumberOfNetworks);
}

} // namespace

void emperical_fit_reset_cache() {
  std::lock_guard<std::mutex> lock(g_cache_mu);
  g_ensemble.reset();
  g_ensemble_A = -1;
  g_ensemble_subdir.clear();
  g_ensemble_cache_key.clear();
}

double emperical_fit_ensemble_mean(const params &p, int mass_number_A,
                                   const double inputs5[5]) {
  if (!p.NeuralNetwork_guided)
    return 1.0;
  const std::string sub = nuclear_subdir(mass_number_A);
  if (sub.empty())
    return 1.0;
  try {
    std::lock_guard<std::mutex> lock(g_cache_mu);
    const std::string want_key = ensemble_cache_key(p, sub);
    if (!g_ensemble || g_ensemble_A != mass_number_A ||
        g_ensemble_subdir != sub || g_ensemble_cache_key != want_key) {
      g_ensemble = std::make_unique<EmpiricalEnsemble>(p, sub);
      g_ensemble_A = mass_number_A;
      g_ensemble_subdir = sub;
      g_ensemble_cache_key = want_key;
    }
    if (!g_ensemble->ok())
      return 1.0;
    const double m = g_ensemble->predict_mean(p, inputs5);
    if (!std::isfinite(m) || m <= 0.0)
      return 1.0;
    return m;
  } catch (const std::exception &ex) {
    std::cerr << "[EmpericalFits] " << ex.what() << std::endl;
    return 1.0;
  }
}
