#include "NN.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#if defined(_WIN32)
#include <codecvt>
#include <locale>
#endif
#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
#include <onnxruntime/coreml_provider_factory.h>
#endif

namespace {
void ThrowOnOrtError(const OrtApi *api, OrtStatus *status) {
  if (status == nullptr)
    return;
  const char *msg = api->GetErrorMessage(status);
  std::string err = msg ? msg : "unknown ONNX Runtime error";
  api->ReleaseStatus(status);
  throw std::runtime_error(err);
}

#if defined(_WIN32)
std::wstring ToWide(const std::string &s) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> conv;
  return conv.from_bytes(s);
}
#endif

static size_t product_shape(const std::vector<int64_t> &shape) {
  size_t n = 1;
  for (int64_t d : shape) {
    if (d <= 0)
      return 0;
    n *= static_cast<size_t>(d);
  }
  return n;
}

static double read_scalar_from_tensor(const OrtApi *api, OrtValue *v) {
  OrtTensorTypeAndShapeInfo *info = nullptr;
  ThrowOnOrtError(api, api->GetTensorTypeAndShape(v, &info));
  enum ONNXTensorElementDataType etype = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;
  ThrowOnOrtError(api, api->GetTensorElementType(info, &etype));
  size_t dim_count = 0;
  ThrowOnOrtError(api, api->GetDimensionsCount(info, &dim_count));
  if (dim_count > 0) {
    std::vector<int64_t> dims(dim_count);
    ThrowOnOrtError(api, api->GetDimensions(info, dims.data(), dim_count));
    size_t prod = product_shape(dims);
    if (prod != 1) {
      api->ReleaseTensorTypeAndShapeInfo(info);
      throw std::runtime_error("ONNX output: expected scalar or [1] tensor");
    }
  }
  void *base = nullptr;
  ThrowOnOrtError(api, api->GetTensorMutableData(v, &base));
  double out = 0;
  if (etype == ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE) {
    out = *static_cast<double *>(base);
  } else if (etype == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT) {
    out = static_cast<double>(*static_cast<float *>(base));
  } else {
    api->ReleaseTensorTypeAndShapeInfo(info);
    throw std::runtime_error("ONNX output: expected float32 or float64 tensor");
  }
  api->ReleaseTensorTypeAndShapeInfo(info);
  return out;
}
} // namespace

OrtRuntime::OrtRuntime(const std::string &instance_name, OrtLoggingLevel log_level)
    : api_(OrtGetApiBase()->GetApi(ORT_API_VERSION)), env_(nullptr) {
  ThrowOnOrtError(api_,
                  api_->CreateEnv(log_level, instance_name.c_str(), &env_));
}

OrtRuntime::~OrtRuntime() {
  if (env_) {
    api_->ReleaseEnv(env_);
    env_ = nullptr;
  }
}

NeuralNetwork::NeuralNetwork(const std::string &model_path, const params &par,
                             const std::string &instance_name,
                             OrtLoggingLevel log_level)
    : NeuralNetwork(std::make_shared<OrtRuntime>(instance_name, log_level),
                    model_path, par, false) {}

NeuralNetwork::NeuralNetwork(std::shared_ptr<OrtRuntime> runtime,
                             const std::string &model_path, const params &par,
                             bool ensemble_member)
    : api_(runtime->api()), runtime_(std::move(runtime)), session_(nullptr),
      session_options_(nullptr), memory_info_(nullptr), using_gpu_(false),
      selected_backend_("CPU"), backend_reason_("") {
  try {
    init_session_options(par, ensemble_member);
    configure_execution_provider(par);
    create_session(model_path);
    cache_io_names();
    ThrowOnOrtError(api_, api_->CreateCpuMemoryInfo(OrtArenaAllocator,
                                                    OrtMemTypeDefault,
                                                    &memory_info_));
    cache_input_shape_and_type();

    std::cerr << "[NeuralNetwork] backend=" << selected_backend_
              << " (HardwareSupport=" << par.HardwareSupport
              << ", ensemble=" << (ensemble_member ? "1" : "0")
              << ", NumberOfCores=" << par.NumberOfCores << ")";
    if (!backend_reason_.empty())
      std::cerr << " reason=\"" << backend_reason_ << "\"";
    std::cerr << std::endl;
  } catch (...) {
    release_session_resources();
    throw;
  }
}

NeuralNetwork::~NeuralNetwork() { release_session_resources(); }

void NeuralNetwork::init_session_options(const params &par,
                                         bool ensemble_member) {
  ThrowOnOrtError(api_, api_->CreateSessionOptions(&session_options_));
  ThrowOnOrtError(api_, api_->SetSessionGraphOptimizationLevel(session_options_,
                                                               ORT_ENABLE_ALL));
  int intra = 1;
  if (!ensemble_member && par.NumberOfCores > 0)
    intra = par.NumberOfCores;
  ThrowOnOrtError(
      api_, api_->SetIntraOpNumThreads(session_options_, intra));
}

void NeuralNetwork::configure_execution_provider(const params &par) {
  if (par.HardwareSupport == 0) {
    selected_backend_ = "CPU";
    backend_reason_ = "HardwareSupport=0";
    return;
  }

#if defined(_WIN32)
  OrtStatus *cuda_status =
      OrtSessionOptionsAppendExecutionProvider_CUDA(session_options_, 0);
  if (cuda_status == nullptr) {
    using_gpu_ = true;
    selected_backend_ = "CUDA";
    backend_reason_ = "Windows GPU: CUDA execution provider initialized";
    return;
  }
  backend_reason_ =
      std::string("CUDA unavailable: ") +
      (api_->GetErrorMessage(cuda_status) ? api_->GetErrorMessage(cuda_status)
                                          : "unknown");
  api_->ReleaseStatus(cuda_status);
#elif defined(__linux__)
  OrtStatus *cuda_status =
      OrtSessionOptionsAppendExecutionProvider_CUDA(session_options_, 0);
  if (cuda_status == nullptr) {
    using_gpu_ = true;
    selected_backend_ = "CUDA";
    backend_reason_ = "Linux GPU: CUDA execution provider initialized";
    return;
  }
  backend_reason_ =
      std::string("CUDA unavailable: ") +
      (api_->GetErrorMessage(cuda_status) ? api_->GetErrorMessage(cuda_status)
                                          : "unknown");
  api_->ReleaseStatus(cuda_status);
#elif defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
  OrtStatus *coreml_status =
      OrtSessionOptionsAppendExecutionProvider_CoreML(session_options_, 0);
  if (coreml_status == nullptr) {
    using_gpu_ = true;
    selected_backend_ = "CoreML";
    backend_reason_ =
        "macOS Apple Silicon GPU: CoreML execution provider initialized";
    return;
  }
  backend_reason_ = std::string("CoreML unavailable: ") +
                    (api_->GetErrorMessage(coreml_status)
                         ? api_->GetErrorMessage(coreml_status)
                         : "unknown");
  api_->ReleaseStatus(coreml_status);
#else
  backend_reason_ = "No GPU provider configured for this platform";
#endif

  using_gpu_ = false;
  selected_backend_ = "CPU";
  if (backend_reason_.empty())
    backend_reason_ = "GPU requested but provider not available";
}

void NeuralNetwork::create_session(const std::string &model_path) {
#if defined(_WIN32)
  const std::wstring wide_path = ToWide(model_path);
  ThrowOnOrtError(api_, api_->CreateSession(runtime_->env(), wide_path.c_str(),
                                            session_options_, &session_));
#else
  ThrowOnOrtError(api_, api_->CreateSession(runtime_->env(), model_path.c_str(),
                                            session_options_, &session_));
#endif
}

void NeuralNetwork::cache_io_names() {
  OrtAllocator *allocator = nullptr;
  ThrowOnOrtError(api_, api_->GetAllocatorWithDefaultOptions(&allocator));

  size_t input_count = 0;
  ThrowOnOrtError(api_, api_->SessionGetInputCount(session_, &input_count));
  for (size_t i = 0; i < input_count; ++i) {
    char *name = nullptr;
    ThrowOnOrtError(api_,
                    api_->SessionGetInputName(session_, i, allocator, &name));
    input_names_.push_back(name ? name : "");
    allocator->Free(allocator, name);
  }

  size_t output_count = 0;
  ThrowOnOrtError(api_, api_->SessionGetOutputCount(session_, &output_count));
  for (size_t i = 0; i < output_count; ++i) {
    char *name = nullptr;
    ThrowOnOrtError(api_,
                    api_->SessionGetOutputName(session_, i, allocator, &name));
    output_names_.push_back(name ? name : "");
    allocator->Free(allocator, name);
  }
}

void NeuralNetwork::cache_input_shape_and_type() {
  OrtTypeInfo *typeinfo = nullptr;
  ThrowOnOrtError(api_,
                  api_->SessionGetInputTypeInfo(session_, 0, &typeinfo));
  const OrtTensorTypeAndShapeInfo *tensor_info = nullptr;
  ThrowOnOrtError(api_, api_->CastTypeInfoToTensorInfo(typeinfo, &tensor_info));
  if (!tensor_info) {
    api_->ReleaseTypeInfo(typeinfo);
    throw std::runtime_error("ONNX model: input 0 is not a tensor");
  }

  ThrowOnOrtError(api_,
                    api_->GetTensorElementType(tensor_info, &input_element_type_));
  size_t dim_count = 0;
  ThrowOnOrtError(api_, api_->GetDimensionsCount(tensor_info, &dim_count));
  input_shape_.resize(dim_count);
  if (dim_count > 0) {
    ThrowOnOrtError(api_, api_->GetDimensions(tensor_info, input_shape_.data(),
                                              dim_count));
  }
  api_->ReleaseTypeInfo(typeinfo);
  size_t nel = product_shape(input_shape_);
  if (nel != 5)
    throw std::runtime_error("ONNX model: expected exactly 5 input elements");
  input_shape_cached_ = true;
}

void NeuralNetwork::release_session_resources() {
  if (memory_info_) {
    api_->ReleaseMemoryInfo(memory_info_);
    memory_info_ = nullptr;
  }
  if (session_) {
    api_->ReleaseSession(session_);
    session_ = nullptr;
  }
  if (session_options_) {
    api_->ReleaseSessionOptions(session_options_);
    session_options_ = nullptr;
  }
}

double NeuralNetwork::run_scalar_io(const double input5[5],
                                     const char *input_tensor_name,
                                     const char *output_tensor_name) {
  if (!input_shape_cached_)
    cache_input_shape_and_type();

  const size_t nelem = product_shape(input_shape_);
  if (nelem != 5)
    throw std::runtime_error("run_scalar_io: input size mismatch");

  std::vector<float> float_buf(5);
  std::vector<double> double_buf(5);
  void *data_ptr = nullptr;
  size_t bytes = 0;
  ONNXTensorElementDataType tensor_type = input_element_type_;

  if (tensor_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE) {
    std::memcpy(double_buf.data(), input5, 5 * sizeof(double));
    data_ptr = double_buf.data();
    bytes = 5 * sizeof(double);
  } else if (tensor_type == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT) {
    for (int i = 0; i < 5; ++i)
      float_buf[i] = static_cast<float>(input5[i]);
    data_ptr = float_buf.data();
    bytes = 5 * sizeof(float);
  } else {
    throw std::runtime_error(
        "ONNX model: input tensor must be float32 or float64");
  }

  OrtValue *input_tensor = nullptr;
  ThrowOnOrtError(
      api_,
      api_->CreateTensorWithDataAsOrtValue(
          memory_info_, data_ptr, bytes, input_shape_.data(),
          input_shape_.size(), tensor_type, &input_tensor));

  const char *in_names[1] = {input_tensor_name};
  const char *out_names[1] = {output_tensor_name};
  OrtValue *output_tensor = nullptr;

  const OrtValue *const_input_vals[1] = {input_tensor};
  ThrowOnOrtError(api_, api_->Run(session_, nullptr, in_names, const_input_vals, 1,
                                  out_names, 1, &output_tensor));

  api_->ReleaseValue(input_tensor);
  double y = read_scalar_from_tensor(api_, output_tensor);
  api_->ReleaseValue(output_tensor);
  return y;
}
