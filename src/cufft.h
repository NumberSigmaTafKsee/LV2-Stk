
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <complex>
#include <ccomplex>

#include "cuda_runtime.h"
#include "cublas_v2.h"
#include <cufftw.h>
#include "viper.hpp"
// 1D = audio/signals
// 2D = images matrix
// 3D = vision 

namespace Viper::cuFFT
{
    template<typename T>
    struct Window 
    {   
        using Vector = Viper::Vector<T>;

        Vector window;

        Window(size_t i) { window.resize(i); }
        virtual ~Window() = default;

        T& operator[](size_t i) { return window[i]; }

        Vector operator * (const Vector& v) { return window * v; }
    };
    template<typename T>
    struct ComplexWindow 
    {   
        using Vector = Viper::ComplexVector<T>;
        Vector window;

        Window(size_t i) { window.resize(i); }
        virtual ~Window() = default;

        T& operator[](size_t i) { return window[i]; }

        Vector operator * (const Vector& v) { return window * v; }
    };

    template<typename T>
    struct Rectangle: public Window<T>
    {
        Rectangle(size_t i) : Window<T>(i) { 
            fill(this->window,1.0f);
            } 
    };
    template<typename T>
    struct Hamming: public Window<T>
    {
        Hamming(size_t n) : Window<T>(n) {            
            for(size_t i = 0; i < this->window.size(); i++)
            {
                this->window[i] = 0.54 - (0.46 * std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
    };
    template<typename T>
    struct Hanning: public Window<T>
    {
        Hanning(size_t n) : Window<T>(n) {            
            for(size_t i = 0; i < this->window.size(); i++)
            {
                this->window[i] = 0.5*(1 - std::cos(2*M_PI*(double)i/(double)n));
            }        
        }
    };
    template<typename T>
    struct Blackman: public Window<T>
    {
        Blackman(size_t n) : Window<T>(n)    
        {            
            for(size_t i = 0; i < this->window.size(); i++)                    
                this->window[i] = 0.42 - (0.5* std::cos(2*M_PI*i/(n)) + (0.08*std::cos(4*M_PI*i/n)));        
        }
    };
    template<typename T>
    struct BlackmanHarris: public Window<T>
    {
        BlackmanHarris(size_t n) : Window<T>(n)    
        {            
            for(size_t i = 0; i < this->window.size(); i++)            
            {   
                double ci = (double) i / (double) n;
                this->window[i] = 0.35875 
                        - 0.48829*std::cos(2*M_PI*(ci))
                        + 0.14128*std::cos(4.0*M_PI*(ci)) 
                        - 0.01168*std::cos(6.0*M_PI*(ci));
            }
        }
    };
    template<typename T>
    struct Gaussian: public Window<T>
    {
        Gaussian(size_t i) : Window<T>(i)
        {
            T a,b,c=0.5;
            for(size_t n = 0; n < this->window.size(); n++)
            {
                a = ((double)n - c*(this->window.size()-1)/(std::sqrt(c)*this->window.size()-1));
                b = -c * std::sqrt(a);
                this->window(n) = std::exp(b);
            }
        }
    };
    template<typename T>
    struct Welch: public Window<T>
    {
        Welch(size_t n) : Window<T>(n)
        {
            for(size_t i = 0; i < this->window.size(); i++)
                this->window[i] = 1.0 - std::sqrt((2.0*(double)i-(double)this->window.size()-1)/((double)this->window.size()));        
        }
    };
    template<typename T>
    struct Parzen: public Window<T>
    {

        Parzen(size_t n) : Window<T>(n)
        {
            for(size_t i = 0; i < this->window.size(); i++)
                this->window[i] = 1.0 - std::abs((2.0*(double)i-this->window.size()-1)/(this->window.size()));        
        }    
    };
    template<typename T>
    struct Tukey: public Window<T>
    {
        Tukey(size_t num_samples, T alpha) : Window<T>(num_samples)
        {            
            T value = (-1*(num_samples/2)) + 1;
            double n2 = (double)num_samples / 2.0;
            for(size_t i = 0; i < this->window.size(); i++)
            {    
                if(value >= 0 && value <= (alpha * (n2))) 
                    this->window[i] = 1.0; 
                else if(value <= 0 && (value >= (-1*alpha*(n2)))) 
                    this->vector.vector[i] = 1.0;
                else 
                    this->vector.vector[i] = 0.5 * (1 + std::cos(M_PI *(((2.0*value)/(alpha*(double)num_samples))-1)))        ;
                value = value + 1;
            }     
        }
    };

    fftw_complex* fftw_alloc_complex(size_t n) {
        fftw_complex * p = nullptr;
        cudaMalloc(&p,n * sizeof(fftw_complex));
        return p;
    }
    fftwf_complex* fftwf_alloc_complex(size_t n) {
        fftwf_complex * p = nullptr;
        cudaMalloc(&p,n * sizeof(fftwf_complex));
        return p;
    }
    double* fftw_alloc_real(size_t n) {
        double * p = nullptr;
        cudaMalloc(&p,n * sizeof(double));
        return p;
    }
    float* fftwf_alloc_real(size_t n) {
        float * p = nullptr;
        cudaMalloc(&p,n * sizeof(float));
        return p;
    }
    ////////////////////////////////////////////////////////////////
    // FFTW Complex 2 Complex
    ////////////////////////////////////////////////////////////////
    struct C2CD
    {
        fftw_complex * in;    
        fftw_complex * out;
        fftw_complex * host;
        size_t size;
        fftw_plan p;

        enum Direction {
            BACKWARD= FFTW_BACKWARD,
            FORWARD = FFTW_FORWARD,
        };

        C2CD(size_t n, Direction dir = FORWARD) {
            in = fftw_alloc_complex(n);
            out= fftw_alloc_complex(n);        
            host= (fftw_complex*)malloc(n * sizeof(fftw_complex));
            size = n;
            p = fftw_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CD() {
            fftw_destroy_plan(p);
            cudaFree(in);
            cudaFree(out);    
            free(host);
        }

        void download_host() {
            cudaMemcpy(host,out,size*sizeof(fftw_complex),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host,size*sizeof(fftw_complex),cudaMemcpyHostToDevice);
        }
        fftw_complex& operator[](size_t index) {
            return host[index];
        }
        
        void set_input(Viper::ComplexVector<double> & input) {
            for(size_t i = 0; i < size; i++) {
                host[i][0] = input[i].real();
                host[i][1] = input[i].imag();
            }
            upload_device();
        }
                
        Viper::ComplexVector<double> get_output() {
            Viper::ComplexVector<double> r(size);            
            for(size_t i = 0; i < size; i++ )
            {
                r[i].real(host[i][0]);
                r[i].imag(host[i][1]);
            }
            return r;
        }
        void get_output(complex<double> * output)
        {            
            download_host();
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host[i][0]);
                output[i].imag(host[i][1]);
            }
        }        
        void get_output(Viper::ComplexVector<double>&  output)
        {
            download_host();
            if(output.size() != size) output.resize(size);            
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host[i][0]);
                output[i].imag(host[i][1]);
            }
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                host[i][0] /= (double)size;    
                host[i][1] /= (double)size;
            }            
        }
        void Execute() {
            fftwf_execute(p);            
        }
        Viper::ComplexVector<double> execute(Viper::ComplexVector<double> & in) {
            if(dir == FORWARD) return Forward(in);
            return Backward(in);
        }        

        Viper::ComplexVector<double> Forward(Viper::ComplexVector<double> & in, bool normalize=true)
        {
            set_input(in);
            fftw_execute(p);
            if(normalize) normalize();
            return get_output();
        }
        Viper::ComplexVector<double> Backward(Viper::ComplexVector<double> & in)
        {
            set_input(in);
            Execute();
            return get_output();
        }

    };

    
    struct C2CF
    {
        fftwf_complex * in;    
        fftwf_complex * out;
        fftwf_complex * host;
        size_t size;
        fftwf_plan p;

        enum Direction {
            BACKWARD=FFTW_BACKWARD,
            FORWARD=FFTW_FORWARD,
        };

        C2CF(size_t n, Direction dir = FORWARD) {
            in = fftwf_alloc_complex(n);
            out= fftwf_alloc_complex(n);        
            host =(fftwf_complex*)malloc(n*sizeof(fftwf_complex));
            size = n;
            p = fftwf_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
        }
        ~C2CF() {
            fftwf_destroy_plan(p);
            cudaFree(in);
            cudaFree(out);    
            free(host);
        }
        void download_host() {
            cudaMemcpy(host,out,size*sizeof(fftwf_complex),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host,size*sizeof(fftwf_complex),cudaMemcpyHostToDevice);
        }
        fftwf_complex& operator[](size_t index) {            
            return host[index];
        }
        void set_input(Viper::ComplexVector<float>> & input) {
            for(size_t i = 0; i < size; i++) {
                host[i][0] = input[i].real();
                host[i][1] = input[i].imag();
            }
            upload_device();
        }
        void set_input(complex<float> * buffer) {
            cudaMemcpy(in,buffer,size*sizeof(complex<float>),cudaMemcpyHostToDevice);
        }
        Viper::ComplexVector<float>> get_output() {
            Viper::ComplexVector<float>> r(size);
            for(size_t i = 0; i < size; i++ )
            {
                r[i].real(host[i][0]);
                r[i].imag(host[i][1]);
            }
            return r;
        }
        void get_output(complex<float> * output)
        {
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host[i][0]);
                output[i].imag(host[i][1]);
            }
        }
        void get_output(Viper::ComplexVector<float>>& output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host[i][0]);
                output[i].imag(host[i][1]);
            }
        }        
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                host[i][0] /= (float)size;    
                host[i][1] /= (float)size;
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }
        Viper::ComplexVector<float> execute(Viper::ComplexVector<float> & in) {
            if(dir == FORWARD) return Forward(in);
            return Backward(in);
        }        

        Viper::ComplexVector<float> Forward(Viper::ComplexVector<float> & in, bool normalize=true)
        {
            set_input(in);
            fftw_execute(p);
            if(normalize) normalize();
            return get_output();
        }
        Viper::ComplexVector<float> Backward(Viper::ComplexVector<float> & in)
        {
            set_input(in);
            Execute();
            return get_output();
        }
    };

    
    ////////////////////////////////////////////////////////////////
    // FFTW Complex 2 Real
    ////////////////////////////////////////////////////////////////
    struct C2RD
    {
        fftw_complex * in;    
        double * out;
        fftw_complex * host_c;
        double *host_d;
        size_t size;
        fftw_plan p;

        C2RD() {
            in = NULL;
            out = NULL;
            host_c = NULL;
            host_d = NULL;
            size = 0;
        }
        C2RD(size_t n) {
            init(n);
        }
        ~C2RD() {
            fftw_destroy_plan(p);
            if(in) cudaFree(in);
            if(out) cudaFree(out);    
            if(host_c) free(host_c);
            if(host_d) free(host_d);
        }
        void download_host() {
            cudaMemcpy(host_d,out,size*sizeof(double),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host_c,size*sizeof(fftw_complex),cudaMemcpyHostToDevice);
        }        
        void init(size_t n) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) cudaFree(in);
            if(out!= NULL) cudaFree(out);
            if(host_c != nullptr) free(host_c);
            if(host_d != nullptr) free(host_d);
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_d = nullptr;
            size = n;
            in = fftw_alloc_complex(n);
            out= fftw_alloc_real(n);       
            host_c= (fftw_complex*)malloc(n*sizeof(fftw_complex));
            host_d= (double*)malloc(n*sizeof(double));
            p = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(Viper::ComplexVector<double>> & input) {
            for(size_t i = 0; i < size; i++) {
                host_c[i][0] = input[i].real();
                host_c[i][1] = input[i].imag();
            }
            upload_device();
        }
        void set_input(complex<double> * buffer) {
            memcpy(in,buffer,size*sizeof(complex<double>));
            upload_device();
        }
        Viper::Vector<double> get_output() {
            Viper::Vector<double> r(size);
            memcpy(r.data(),host_d, size * sizeof(double));
            return r;
        }
        void get_output(double * output)
        {
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = host_d[i];                
            }
        }
        void get_output( Viper::Vector<double> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = host_d[i];                
            }
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                host_d[i] /= (double)size;                    
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }        
        Viper::Vector<double> execute(Viper::ComplexVector<double> & in) {
            fftw_execute(p);
            return get_output();
        }                    
    };

    struct C2RF
    {
        fftwf_complex * in;  
        fftwf_complex * host_c;
        float * host_f;  
        float * out;
        size_t size;
        fftwf_plan p;

        C2RF() {
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = 0;
        }
        C2RF(size_t n) {
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = 0;
            init(n);
        }
        ~C2RF() {
            fftwf_destroy_plan(p);
            if(in) cudaFree(in);
            if(out) cudaFree(out);    
            if(host_c) free(host_c);
            if(host_f) free(host_f);
        }
        void download_host() {
            cudaMemcpy(host_f,out,size*sizeof(float),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host_c,size*sizeof(fftwf_complex),cudaMemcpyHostToDevice);
        }        
        void init(size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) cudaFree(in);
            if(out != NULL) cudaFree(out);
            if(host_c != nullptr) free(host_c);
            if(host_f != nullptr) free(host_f);
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = n;
            in = fftwf_alloc_complex(n);
            out= fftwf_alloc_real(n);                    
            host_c = (fftwf_complex*)malloc(n*sizeof(fftwf_complex));
            host_f = (float*)malloc(n*sizeof(float));
            p = fftwf_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(Viper::ComplexVector<float>> & input) {
            for(size_t i = 0; i < size; i++) {
                host_c[i][0] = input[i].real();
                host_c[i][1] = input[i].imag();
            }
            upload_device();
        }
        void set_input(complex<float> * buffer) {
            memcpy(host_c,buffer,size*sizeof(fftwf_complex));
            upload_device();
        }
        Viper::Vector<float> get_output() {
            Viper::Vector<float> r(size);
            memcpy(r.data(),host_f, size*sizeof(float));
            return r;
        }
        void get_output(float * output)
        {
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = host_f[i];
            }
        }
        void get_output( Viper::Vector<float> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = host_f[i];
            }
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) 
                host_f[i] /= (float)size;                
        }
        void Execute() {
            fftwf_execute(p);            
        }        
        Viper::Vector<double> execute(Viper::ComplexVector<double> & in) {
            fftw_execute(p);
            return get_output();
        }                    
    };


    ////////////////////////////////////////////////////////////////
    // FFTW Real 2 Complex
    ////////////////////////////////////////////////////////////////
    struct R2CD
    {
        double       * in;    
        fftw_complex * out;
        double       * host_d;
        fftw_complex * host_c;
        size_t size;
        fftw_plan p;

        R2CD() {
            in = NULL;
            out = NULL;
            host_d = nullptr;
            host_c = nullptr;
            size= 0;
        }
        R2CD(size_t n) {
            in = NULL;
            out = NULL;
            host_d = nullptr;
            host_c = nullptr;
            size= 0;
            init(n);            
        }
        ~R2CD() {
            fftw_destroy_plan(p);
            if(in) cudaFree(in);
            if(out) cudaFree(out);    
            if(host_c) free(host_c);
            if(host_d) free(host_d);
        }        
        void init(size_t n) {
            if(in != NULL) fftw_destroy_plan(p);
            if(in != NULL) cudaFree(in);
            if(out != NULL) cudaFree(out);
            if(host_c != nullptr) free(host_c);
            if(host_d != nullptr) free(host_d);
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_d = nullptr;
            size = n;
            in = fftw_alloc_real(n);
            out= fftw_alloc_complex(n);                                
            host_c = (fftw_complex*) malloc(n*sizeof(fftw_complex));
            host_d = (double*) malloc(n*sizeof(double));
            p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
        }
        void download_host() {
            cudaMemcpy(host_c,out,size*sizeof(fftw_complex),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host_d,size*sizeof(double),cudaMemcpyHostToDevice);
        }        
        void set_input(Viper::Vector<double> & input) {
            memcpy(host_d,input.data(),size*sizeof(double));
            upload_device();
        }
        void set_input(double * buffer) {
            memcpy(host_d,buffer,size*sizeof(double));
            upload_device();
        }
        Viper::ComplexVector<double>> get_output() {
            Viper::ComplexVector<double>> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i].real(host_c[i][0]);
                r[i].imag(host_c[i][1]);
            }
            return r;
        }
        void get_output(complex<double> * output)
        {
            for(size_t i = 0; i < size; i++)
            {
                output[i].real(host_c[i][0]);
                output[i].imag(host_c[i][1]);
            }
        }
        void get_output(Viper::ComplexVector<double>> & output) {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++)
            {
                output[i].real(host_c[i][0]);
                output[i].imag(host_c[i][1]);
            }            
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                host_c[i][0] /= (double)size;    
                host_c[i][1] /= (double)size;
            }
        }
        void Execute() {
            fftw_execute(p);            
        }
        Viper::ComplexVector<double> execute(Viper::Vector<double> & in, bool normal=true) {
            fftw_execute(p);
            if(normal) normalize();
            return get_output();
        }                    
    };

    struct R2CF
    {
        float * in;    
        fftwf_complex * out;
        float * host_f;
        fftwf_complex * host_c;
        size_t size;
        fftwf_plan p;

        R2CF() {
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = 0;
        }
        R2CF(size_t n) {
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = 0;
            init(n);            
        }
        ~R2CF() {
            fftwf_destroy_plan(p);
            if(in) cudaFree(in);
            if(out) cudaFree(out);    
            if(host_c) free(host_c);
            if(host_f) free(host_f);
            
        }
        void download_host() {
            cudaMemcpy(host_c,out,size*sizeof(fftwf_complex),cudaMemcpyDeviceToHost);
        }
        void upload_device() {
            cudaMemcpy(in,host_f,size*sizeof(float),cudaMemcpyHostToDevice);
        }        
        void init(size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) cudaFree(in);
            if(out != NULL) cudaFree(out);
            if(host_c != nullptr) free(host_c);
            if(host_f != nullptr) free(host_f);
            in = NULL;
            out = NULL;
            host_c = nullptr;
            host_f = nullptr;
            size = n;
            in = fftwf_alloc_real(n);
            out= fftwf_alloc_complex(n);                    
            host_c = (fftwf_complex*)malloc(n*sizeof(fftwf_complex));
            host_f = (float*)malloc(n*sizeof(float));
            p = fftwf_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(Viper::Vector<float> & input) {
            memcpy(host_f,input.data(),size*sizeof(float));
            upload_device();
        }
        void set_input(float * buffer) {
            memcpy(host_f,buffer,size*sizeof(float));
            upload_device();
        }
        Viper::ComplexVector<float>> get_output() {
            Viper::ComplexVector<float>> r(size);
            for(size_t i = 0; i < size; i++) {
                r[i].real(host_c[i][0]);
                r[i].imag(host_c[i][1]);
            }                
            return r;
        }    
        void get_output(complex<float> * output)
        {
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host_c[i][0]);
                output[i].imag(host_c[i][1]);
            }
        }
        void get_output( Viper::ComplexVector<float>> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i].real(host_c[i][0]);
                output[i].imag(host_c[i][1]);
            }
        }
        void normalize() {
            for(size_t i = 0; i < size; i++) {
                host_c[i][0] /= (float)size;    
                host_c[i][1] /= (float)size;
            }
        }
        void Execute() {
            fftwf_execute(p);            
        }
        Viper::ComplexVector<float> execute(Viper::Vector<float> & in, bool normal=true) {
            fftw_execute(p);
            if(normal) normalize();
            return get_output();
        }                    
    };

    Viper::Vector<float> conv(Viper::Vector<float> x, Viper::Vector<float> y)
    {
        int size = x.size() + y.size() - 1;
        if(size % 2 != 0) size =  pow(2,log2(size)+1):size;
        R2CF forward(size);
        C2RF backward(size);
        int size1 = x.size % 2 != 0? pow(2,log2(x.size())+1):x.size();
        int size2 = y.size % 2 != 0? pow(2,log2(y.size())+1):y.size();
        Viper::Vector<float> t1(log2(size1));
        Viper::Vector<float> t2(log2(size2));
        forward.set_input(t1);
        t1 = forward.execute(t1);
        forward.set_input(t2);
        t2 = forward.execute(t2);
        Viper::Vector<float> c = t1*t2;
        backward.set_input(c);
        c = backward.get_output();
        return c;
    }
    void blockconv(Viper::Vector<float> h, Viper::Vector<float> x, Viper::Vector<float>& y, Viper::Vector<float> & ytemp)    
    {
        int i;
        int M = h.size();
        int L = x.size();
        y = conv(h,x);      
        for (i=0; i<M; i++) {
            y[i] += ytemp[i]; 
            ytemp[i] = y[i+L];
        }        
    }

    Viper::Vector<float> deconv(Viper::Vector<float> x, Viper::Vector<float> y)
    {
        int size = x.size() + y.size() - 1;
        if(size % 2 != 0) size =  pow(2,log2(size)+1):size;
        R2CF forward(size);
        C2RF backward(size);
        int size1 = x.size % 2 != 0? pow(2,log2(x.size())+1):x.size();
        int size2 = y.size % 2 != 0? pow(2,log2(y.size())+1):y.size();
        Viper::Vector<float> t1(log2(size1));
        Viper::Vector<float> t2(log2(size2));
        forward.set_input(t1);
        t1 = forward.execute(t1);
        forward.set_input(t2);
        t2 = forward.execute(t2);
        Viper::Vector<float> c = t1/t2;
        backward.set_input(c);
        c = backward.get_output();
        return c;
    }
    /*
    Viper::Vector<float> xcorr(Viper::Vector<float> x, Viper::Vector<float> y)
    {
        int size = x.size() + y.size() - 1;
        if(size % 2 != 0) size =  pow(2,log2(size)+1):size;
        R2CF forward(size);
        C2RF backward(size);
        int size1 = x.size % 2 != 0? pow(2,log2(x.size())+1):x.size();
        int size2 = y.size % 2 != 0? pow(2,log2(y.size())+1):y.size();
        Viper::Vector<float> t1(log2(size1));
        Viper::Vector<float> t2(log2(size2));
        forward.set_input(t1);
        t1 = forward.execute(t1);
        forward.set_input(t2);
        t2 = forward.execute(t2);
        Viper::Vector<float> c = conj(t1)*t2;
        backward.set_input(c);
        c = backward.get_output();
        return c;
    }
    */
}
