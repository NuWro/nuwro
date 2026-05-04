#ifndef NN_H
#define NN_H

#include "../params.h"
#include "../params_all.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <onnxruntime/onnxruntime_c_api.h>

/// Shared ONNX Runtime environment (one OrtEnv per process/ensemble).
class OrtRuntime {
public:
  explicit OrtRuntime(const std::string &instance_name,
                      OrtLoggingLevel log_level = ORT_LOGGING_LEVEL_WARNING);
  ~OrtRuntime();

  OrtRuntime(const OrtRuntime &) = delete;
  OrtRuntime &operator=(const OrtRuntime &) = delete;

  const OrtApi *api() const { return api_; }
  OrtEnv *env() const { return env_; }

private:
  const OrtApi *api_;
  OrtEnv *env_;
};

/// One loaded ONNX model: session + I/O metadata + optional scalar Run.
class NeuralNetwork {
public:
  /// Single model: creates its own OrtRuntime.
  NeuralNetwork(const std::string &model_path, const params &par,
                const std::string &instance_name = "NuWroNN",
                OrtLoggingLevel log_level = ORT_LOGGING_LEVEL_WARNING);

  /// Ensemble member: shares OrtRuntime; use ensemble_member=true for intra-op
  /// threads = 1 when outer parallelism runs over models.
  NeuralNetwork(std::shared_ptr<OrtRuntime> runtime,
                const std::string &model_path, const params &par,
                bool ensemble_member);

  ~NeuralNetwork();

  NeuralNetwork(const NeuralNetwork &) = delete;
  NeuralNetwork &operator=(const NeuralNetwork &) = delete;

  OrtSession *session() const { return session_; }
  OrtMemoryInfo *memory_info() const { return memory_info_; }

  const std::vector<std::string> &input_names() const { return input_names_; }
  const std::vector<std::string> &output_names() const { return output_names_; }

  bool using_gpu() const { return using_gpu_; }
  const std::string &selected_backend() const { return selected_backend_; }
  const std::string &backend_reason() const { return backend_reason_; }

  /// Run with tensor names "input" and "Identity:0" (or overrides). Uses model
  /// element type (double or float). Returns first scalar output element.
  double run_scalar_io(const double input5[5],
                       const char *input_tensor_name = "input",
                       const char *output_tensor_name = "Identity:0");

private:
  const OrtApi *api_;
  std::shared_ptr<OrtRuntime> runtime_;
  OrtSession *session_;
  OrtSessionOptions *session_options_;
  OrtMemoryInfo *memory_info_;

  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;

  bool using_gpu_;
  std::string selected_backend_;
  std::string backend_reason_;

  std::vector<int64_t> input_shape_;
  ONNXTensorElementDataType input_element_type_{
      ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED};
  bool input_shape_cached_{false};

  void init_session_options(const params &par, bool ensemble_member);
  void configure_execution_provider(const params &par);
  void create_session(const std::string &model_path);
  void cache_io_names();
  void cache_input_shape_and_type();
  void release_session_resources();
};

#endif
