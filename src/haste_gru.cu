#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

#include "haste.hpp"

#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>

template<typename T>
struct device_ptr {
  static constexpr size_t ElemSize = sizeof(typename T::Scalar);

  static device_ptr<T> NewByteSized(size_t bytes) {
    return device_ptr<T>((bytes + ElemSize - 1) / ElemSize);
  }

  explicit device_ptr(size_t size_)
      : data(nullptr), size(size_) {
    void* tmp;
    cudaMalloc(&tmp, size * ElemSize);
    data = static_cast<typename T::Scalar*>(tmp);
  }

  explicit device_ptr(const T& elem)
      : data(nullptr), size(elem.size()) {
    void* tmp;
    cudaMalloc(&tmp, size * ElemSize);
    data = static_cast<typename T::Scalar*>(tmp);
    ToDevice(elem);
  }

  device_ptr(device_ptr<T>&& other) : data(other.data), size(other.size) {
    other.data = nullptr;
    other.size = 0;
  }

  device_ptr& operator=(const device_ptr<T>&& other) {
    if (&other != this) {
      data = other.data;
      size = other.size;
      other.data = nullptr;
      other.size = 0;
    }
    return *this;
  }

  device_ptr(const device_ptr<T>& other) = delete;
  device_ptr& operator=(const device_ptr<T>& other) = delete;

  void ToDevice(const T& src) {
    assert(size == src.size());
    cudaMemcpy(data, src.data(), src.size() * ElemSize, cudaMemcpyHostToDevice);
  }

  void ToHost(T& target) const {
    assert(size == target.size());
    cudaMemcpy(target.data(), data, target.size() * ElemSize, cudaMemcpyDeviceToHost);
  }

  size_t Size() const {
    return size;
  }

  void zero() {
    cudaMemset(data, 0, size * ElemSize);
  }

  ~device_ptr() {
    cudaFree(data);
  }

  typename T::Scalar* data;
  size_t size;
};

using haste::v0::gru::ForwardPass;
using haste::v0::gru::BackwardPass;
using std::string;

using Tensor1 = Eigen::Tensor<float, 1>;
using Tensor2 = Eigen::Tensor<float, 2>;
using Tensor3 = Eigen::Tensor<float, 3>;

constexpr int BATCH_SIZE = 64;
constexpr int SEQUENCE_LEN = 1000;
constexpr int HIDDEN_DIMS = 512;
constexpr int INPUT_DIMS = 512;

static cublasHandle_t g_blas_handle;

class ScopeTimer {
  public:
    ScopeTimer(const string& msg) : msg_(msg) {
      cudaEventCreate(&start_);
      cudaEventCreate(&stop_);
      cudaDeviceSynchronize();
      cudaEventRecord(start_);
    }

    ~ScopeTimer() {
      float elapsed_ms;
      cudaEventRecord(stop_);
      cudaEventSynchronize(stop_);
      cudaEventElapsedTime(&elapsed_ms, start_, stop_);
      printf("%s %fms\n", msg_.c_str(), elapsed_ms);
      cudaEventDestroy(start_);
      cudaEventDestroy(stop_);
    }

  private:
    string msg_;
    cudaEvent_t start_, stop_;
};

void GruInference(
    const Tensor2& W,
    const Tensor2& R,
    const Tensor1& bx,
    const Tensor1& br,
    const Tensor3& x) {
  const int time_steps = x.dimension(2);
  const int batch_size = x.dimension(1);
  const int input_size = x.dimension(0);
  const int hidden_size = R.dimension(1);

  // Copy weights over to GPU.
  device_ptr<Tensor2> W_dev(W);
  device_ptr<Tensor2> R_dev(R);
  device_ptr<Tensor1> bx_dev(bx);
  device_ptr<Tensor1> br_dev(br);
  device_ptr<Tensor3> x_dev(x);

  device_ptr<Tensor2> h_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> tmp_Wx_dev(time_steps * batch_size * hidden_size * 3);
  device_ptr<Tensor2> tmp_Rh_dev(batch_size * hidden_size * 3);

  h_dev.zero();

  ScopeTimer t("Inference:");

  ForwardPass<float> forward = ForwardPass<float>(
      false,  // training
      batch_size,
      input_size,
      hidden_size,
      g_blas_handle);

  forward.Run(
      time_steps,
      W_dev.data,
      R_dev.data,
      bx_dev.data,
      br_dev.data,
      x_dev.data,
      h_dev.data,
      nullptr,
      tmp_Wx_dev.data,
      tmp_Rh_dev.data,
      0.0f,
      nullptr);
}

void GruTrain(
    const Tensor2& W,
    const Tensor2& R,
    const Tensor1& bx,
    const Tensor1& br,
    const Tensor3& x,
    const Tensor3& dh_new) {
  const int time_steps = x.dimension(2);
  const int batch_size = x.dimension(1);
  const int input_size = x.dimension(0);
  const int hidden_size = R.dimension(1);

  // Copy weights over to GPU.
  device_ptr<Tensor2> W_dev(W);
  device_ptr<Tensor2> R_dev(R);
  device_ptr<Tensor1> bx_dev(bx);
  device_ptr<Tensor1> br_dev(br);
  device_ptr<Tensor3> x_dev(x);
  device_ptr<Tensor3> dh_new_dev(dh_new);

  device_ptr<Tensor2> h_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> tmp_Wx_dev(time_steps * batch_size * hidden_size * 3);
  device_ptr<Tensor2> tmp_Rh_dev(batch_size * hidden_size * 3);
  device_ptr<Tensor3> v_dev(time_steps * batch_size * hidden_size * 4);

  h_dev.zero();

  {
    ScopeTimer t("Train forward:");
    ForwardPass<float> forward = ForwardPass<float>(
        true,  // training
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    forward.Run(
        time_steps,
        W_dev.data,
        R_dev.data,
        bx_dev.data,
        br_dev.data,
        x_dev.data,
        h_dev.data,
        v_dev.data,
        tmp_Wx_dev.data,
        tmp_Rh_dev.data,
        0.0f,
        nullptr);
  }

  device_ptr<Tensor3> dx_dev(time_steps * batch_size * input_size);
  device_ptr<Tensor2> dW_dev(input_size * hidden_size * 3);
  device_ptr<Tensor2> dR_dev(hidden_size * hidden_size * 3);
  device_ptr<Tensor1> dbx_dev(hidden_size * 3);
  device_ptr<Tensor1> dbr_dev(hidden_size * 3);
  device_ptr<Tensor2> dh_dev(batch_size * hidden_size);
  device_ptr<Tensor3> dp_dev(time_steps * batch_size * hidden_size * 3);
  device_ptr<Tensor3> dq_dev(time_steps * batch_size * hidden_size * 3);

  {
    ScopeTimer t("Train backward:");
    BackwardPass<float> backward(
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    backward.Run(
        time_steps,
        W_dev.data,
        R_dev.data,
        bx_dev.data,
        br_dev.data,
        x_dev.data,
        h_dev.data,
        v_dev.data,
        dh_new_dev.data,
        dx_dev.data,
        dW_dev.data,
        dR_dev.data,
        dbx_dev.data,
        dbr_dev.data,
        dh_dev.data,
        dp_dev.data,
        dq_dev.data,
        nullptr);
  }
}

int main() {
  srand(time(0));

  cublasCreate(&g_blas_handle);

  // Weights.
  Tensor2 W(HIDDEN_DIMS * 3, INPUT_DIMS);
  Tensor2 R(HIDDEN_DIMS * 3, HIDDEN_DIMS);
  Tensor1 bx(HIDDEN_DIMS * 3);
  Tensor1 br(HIDDEN_DIMS * 3);

  // Input.
  Tensor3 x(INPUT_DIMS, BATCH_SIZE, SEQUENCE_LEN);

  // Gradients from upstream layers.
  Tensor3 dh(HIDDEN_DIMS, BATCH_SIZE, SEQUENCE_LEN + 1);

  W.setRandom();
  R.setRandom();
  bx.setRandom();
  br.setRandom();
  x.setRandom();
  dh.setRandom();

  GruInference(W, R, bx, br, x);
  GruTrain(W, R, bx, br, x, dh);

  cublasDestroy(g_blas_handle);

  return 0;
}