#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <valarray>

typedef float SampleType; 

#include "TinyEigen.hpp"

template<typename T>
using EigenArray   = Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>;
template<typename T>
using EigenArray2D = Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

namespace Casino::TinyEigen
{        
    template<typename T>
    struct SampleVector : public EigenArray<T>
    {        
        
        size_t channels;

        SampleVector() {
            channels = 1;
        }
        SampleVector(size_t size, size_t channels = 1) {
            resize(size * channels);
            this->channels = channels;
        }
        SampleVector(const std::vector<T> & s, size_t chans = 1) {
            resize(s.size());
            memcpy(data(),s.data(),s.size()*sizeof(T));
            channels = chans;
        }
        SampleVector(const T * ptr, size_t n, size_t chans = 1) {
            resize(n*chans);
            memcpy(data(),ptr,n*chans*sizeof(T));
            channels = chans;
        }
        SampleVector(const SampleVector<T> & s) {
            (*this) = s;
            channels = s.channels;
        }    
        SampleVector(const EigenArray<T> & v, size_t ch = 1) {
            (*this) = v;
            channels = ch;
        }
        SampleVector(size_t channels, const EigenArray<T>& a) {
            (*this) = a;
            this->channels = channels;
        }

        SampleVector<T> clamp(T min = (T)0.0, T max = (T)1.0) {
            SampleVector<T> r(*this,channels);
            r = SampleVector<T>(r.cwiseMin(max).cwiseMax(min));
            return r;
        }
        
        T&   get_stride(size_t ch, size_t pos) { return (*this)[pos*channels + ch]; }
        void set_stride(size_t ch, size_t pos, const T val) { (*this)[pos*channels + ch] = val; } 
        
        void swap_stereo_channels() {
            assert(channels == 2);
            for(size_t i = 0; i < size(); i+= channels) {
                T temp = (*this)[i];
                (*this)[i] = (*this)[i+1];
                (*this)[i+1] = temp;
            }
        }

        void set_data(size_t n, size_t channels, T * samples)
        {
            resize(n);
            this->channels = channels;
            memcpy(data(),samples,n*sizeof(T));
        }
        void copy_data(T * samples) {
            memcpy(samples,data(),size()*sizeof(T));
        }
        void set_channel(size_t ch, T * samples) {
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = samples[x++];
        }

        size_t size() const { return EigenArray<T>::size(); }
        size_t num_channels() const { return channels; } 
        void zero() { EigenArray<T>::setZero(); }
        void ones() { EigenArray<T>::setOnes(); }
        void random() { EigenArray<T>::setRandom(); }
        void fill(T s) { EigenArray<T>::fill(s); }
        size_t samples_per_channel() const { return size() / num_channels(); }
        T sum() { return EigenArray<T>::sum(); }
        T min() { return EigenArray<T>::minCoeff(); }
        T max() { return EigenArray<T>::maxCoeff(); }
        size_t min_index() { size_t i; EigenArray<T>::minCoeff(&i); return i; }
        size_t max_index() { size_t i; EigenArray<T>::maxCoeff(&i); return i; }
        T* data() { return EigenArray<T>::data(); }
        
        void normalize() { EigenArray<T>::matrix().normalize(); }
        SampleVector<T> normalized() { return SampleVector<T>(EigenArray<T>::matrix().normalized().array()); }

        SampleVector<T>& operator = (const SampleVector<T> & v) {
            (*this) = v;
            channels  = v.channels;
            return *this;
        }    
        T& operator()(size_t i, size_t ch) {
            EigenArray<T> & r = *this;
            return r[i*channels + ch];
        }
        T operator()(size_t i, size_t ch) const {
            const EigenArray<T> & r = *this;
            return r[i*channels + ch];
        }
        T& operator[](size_t i) {
            EigenArray<T> & r = *this;
            return r[i];
        }
        T operator[](size_t i) const {
            const EigenArray<T> & r = *this;
            return r[i];
        }
         
        SampleVector<T> get_channel(size_t channel) {            
            SampleVector<T> r(samples_per_channel());
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];
            return r;
        }
        void get_channel(size_t channel, T * r) {                        
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];            
        }
        SampleVector<T> get_channel(size_t channel) const {            
            SampleVector<T> r(samples_per_channel());
            size_t x=0;
            for(size_t i = channel; i < size(); i+=channels)
                r[x++] = (*this)[i];
            return r;
        }
        void set_channel(const SampleVector<T> & v, size_t ch) {            
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = v[x++];
        }
        void set_channel(const T* v, size_t ch) {            
            size_t x = 0;
            for(size_t i = ch; i < size(); i+=channels) (*this)[i] = v[x++];
        }
        void make_stereo(const T * in) {
            set_channel(in,0);
            set_channel(in,1);
        }
        void make_stereo(const T * left, const T * right) {
            set_channel(left,0);
            set_channel(right,1);
        }

        void set_channels(size_t c) {
            channels = c;
            resize(size()*c);
        }
        void resize(size_t n) {
            EigenArray<T>::resize(n);
        }
        void resize(size_t s, size_t c) {
            channels = c;
            EigenArray<T>::resize(s * c);
        }

        SampleVector<T> get_channel_count(size_t ch, size_t pos, size_t n) {
            SampleVector<T> r(n);
            memcpy(r.data(), get_channel(ch).data()+pos, n*sizeof(T));
            return r;
        }

        bool operator == (const SampleVector<T> & s) {
            return s.channels == channels && size() == s.size();
        }

        void copy(SampleVector<T> & in, size_t block) {
            resize(in.size());
            memcpy(data(),in.data(),block*sizeof(T));
        }
        void copy(SampleVector<T> & in, size_t pos, size_t n) {
            resize((pos+n));
            memcpy(data(),in.data()+pos,n*sizeof(T));
        }
        void copy_from(const T* in, size_t block) {            
            memcpy(data(),in,block*sizeof(T));
        }
        void copy_to(T* in, size_t block) {            
            memcpy(in,data(),block*sizeof(T));
        }

        SampleVector<T> slice(size_t i, size_t n) {
            SampleVector<T> r(n);
            memcpy(r.data(), data()+i, n * sizeof(T));
            return r;
        }
        
        SampleVector<T> pan(SampleType pan) {
            SampleType pan_mapped = ((pan+1.0)/2.0) * (M_PI/2.0);
            SampleVector<T> r(size(),2);
            for(size_t i = 0; i < size(); i+=channels)
            {
                r[i] = (*this)[i] * std::sin(pan_mapped);
                if(num_channels() == 2)
                    r[i+1] = (*this)[i+1] * std::cos(pan_mapped);
                else
                    r[i+1] = (*this)[i] * std::cos(pan_mapped);
            }
            return r;        
        }
        SampleVector<T> stride_slice(const SampleVector<T> & in, size_t stride, size_t pos, size_t n)
        {
            SampleVector<T> r(n/stride);
            size_t x = pos;
            for(size_t i = 0; i < n; i++) {
                r[i] = (*this)[x];
                x += stride;
            } 
            return r;
        }
        SampleVector<T> operator + (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r += s;
            return r;
        }
        SampleVector<T> operator - (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r -= s;
            return r;
        }
        SampleVector<T> operator * (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r *= s;
            return r;
        }
        SampleVector<T> operator / (const SampleVector<T> &s) {
            assert(*this == s);
            SampleVector<T> r(*this);
            r /= s;
            return r;
        }
        
        SampleVector<T> operator + (const T s) {    
            SampleVector<T> r(*this);
            r += s;
            return r;
        }
        SampleVector<T> operator - (const T s) {        
            SampleVector<T> r(*this);
            r -= s;
            return r;
        }
        SampleVector<T> operator / (const T s) {        
            SampleVector<T> r(*this);
            r /= s;
            return r;
        }
        SampleVector<T> operator * (const T s) {        
            SampleVector<T> r(*this);
            r *= s;
            return r;
        }
        

        SampleVector<T>& operator += (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) += s;
            return *this;
        }
        SampleVector<T>& operator -= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) -= s;
            return *this;
        }
        SampleVector<T>& operator *= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) *= s;
            return *this;
        }
        SampleVector<T>& operator /= (const SampleVector<T>& s) {        
            assert(*this == s);
            (*this) /= s;
            return *this;
        }
        
        SampleVector<T> operator += (const T s) {
            (*this) += s;
            return *this;
        }
        SampleVector<T> operator -= (const T s) {
            (*this) -= s;
            return *this;
        }
        SampleVector<T> operator /= (const T s) {
            (*this) /= s;
            return *this;
        }
        SampleVector<T> operator *= (const T s) {
            (*this) *= s;
            return *this;
        }
        
        
        SampleVector<T> abs() { 
            SampleVector<T> r(*this);
            r = this->abs();
            return r;
        }
        SampleVector<T> abs2() { 
            SampleVector<T> r(*this);
            r = this->abs2();
            return r;
        }
        SampleVector<T> inverse() { 
            SampleVector<T> r(*this);
            r = this->inverse();
            return r;
        }
        SampleVector<T> exp() { 
            SampleVector<T> r(*this);
            r = this->exp();
            return r;
        }
        SampleVector<T> log() { 
            SampleVector<T> r(*this);
            r = this->log();
            return r;
        }
        SampleVector<T> log1p() { 
            SampleVector<T> r(*this);
            r = this->log1p();
            return r;
        }
        SampleVector<T> log10() { 
            SampleVector<T> r(*this);
            r = this->log10();
            return r;
        }
        SampleVector<T> pow(const SampleVector<T> & s) { 
            SampleVector<T> r(*this);
            r = this->pow(s);
            return r;
        }
        SampleVector<T> pow(const T s) { 
            SampleVector<T> r(*this);
            r = this->pow(s);        
            return r;        
        }
        SampleVector<T> sqrt() {
            SampleVector<T> r(*this);
            r = this->sqrt();
            return r;
        }
        SampleVector<T> rsqrt() {
            SampleVector<T> r(*this);
            r = this->rsqrt();
            return r;
        }
        SampleVector<T> square() {
            SampleVector<T> r(*this);
            r = this->square();
            return r;
        }
        SampleVector<T> sin() {
            SampleVector<T> r(*this);
            r = this->sin();
            return r;
        }
        SampleVector<T> cos() {
            SampleVector<T> r(*this);
            r = this->cos();
            return r;
        }
        SampleVector<T> tan() {
            SampleVector<T> r(*this);
            r = this->tan();
            return r;
        }
        SampleVector<T> asin() {
            SampleVector<T> r(*this);
            r = this->asin();
            return r;
        }
        SampleVector<T> acos() {
            SampleVector<T> r(*this);
            r = this->acos();
            return r;
        }
        SampleVector<T> atan() {
            SampleVector<T> r(*this);
            r = this->atan();
            return r;
        }        
        SampleVector<T> sinh() {
            SampleVector<T> r(*this);
            r = this->sinh();
            return r;
        }
        SampleVector<T> cosh() {
            SampleVector<T> r(*this);
            r = this->cosh();
            return r;
        }
        SampleVector<T> tanh() {
            SampleVector<T> r(*this);
            r = this->tanh();
            return r;
        }
        SampleVector<T> ceil() { 
            SampleVector<T> r(*this);
            r = this->ceil();
            return r;        
        }
        SampleVector<T> floor() { 
            SampleVector<T> r(*this);
            r = this->floor();
            return r;
        }
        SampleVector<T> round() { 
            SampleVector<T> r(*this);
            r = this->round();
            return r;
        }

        void println() {
            std::cout << "Size=" << size() << std::endl;
            std::cout << "Channel=" << num_channels() << std::endl;
            std::cout << "Samples per channel=" << samples_per_channel() << std::endl;
            for(size_t i = 0; i < size(); i++)
                std:: cout << (*this)[i] << ",";
            std::cout << std::endl;
        }

    };

    template<typename T>
    struct SampleMatrix : public EigenArray2D<T>
    {        
        SampleMatrix() = default;
        
        SampleMatrix(size_t chan, size_t samps) {
            resize(chan,samps);
        }
        SampleMatrix(const SampleVector<T> & v) {
            resize(v.num_channels(), v.size()/v.num_channels());
            for(size_t i = 0; i < v.num_channels(); i++)             
                this->row(i) = v.get_channel(i);
        }
        SampleMatrix(const EigenArray2D<T> & c) {
            (*this) = c;
        }
        SampleMatrix(const SampleMatrix<T> & m) {
            (*this) = m.channels;
        }

        T * data() { return EigenArray2D<T>::data(); }
        
        void deinterleave(const SampleVector<T> & s) {
            
            esize(s.num_channels(),s.samples_per_channel());            
            
            for(size_t i = 0; i < s.num_channels(); i++) {
                this->row(i) = s.get_channel(i);
            }
            
        }    

        SampleVector<T> interleave() { 
            SampleVector<T> r;
            r.resize(cols() * cols());
            for(size_t i = 0; i < rows(); i++)            
            {
                SampleVector<T> tmp(1,(*this)(i));            
                r.set_channel(tmp,i);
            }
            return r;
        }

        size_t num_channels() const { return EigenArray2D<T>::rows(); }
        size_t rows() const { return EigenArray2D<T>::rows(); }
        size_t cols() const { return EigenArray2D<T>::cols(); }
        size_t size() const { return rows()*cols(); }
        void   resize(size_t r, size_t c) { EigenArray2D<T>::resize(r,c); }
        
        SampleVector<T>& row(size_t i) {
            return EigenArray2D<T>::row(i);
        }
        SampleVector<T>& col(size_t i) {
            return EigenArray2D<T>::col(i);
        }
        SampleVector<T> operator[](size_t ch) { 
            SampleVector<T> tmp(1,row(ch));
            return tmp;
        }
        SampleVector<T> get_channel(size_t c) { 
            SampleVector<T> tmp(1,row(c));
            return tmp;
        }
        void set_channel(size_t c, const SampleVector<T> & s) 
        { 
            row(c) = s;
        }

        SampleMatrix<T>& operator = (const SampleMatrix<T> & v) {        
            (*this)  = v;
            return *this;
        }
        T& operator()(size_t c, size_t n) {
            return (*this)(c,n);
        }
        bool operator == (const SampleMatrix<T> & b) {
            return rows() == b.rows() && cols() == b.cols();
        }
        SampleMatrix<T> operator + (const SampleMatrix<T> & m) {        
            SampleMatrix<T> r(*this);
            r = r + m;
            return r;
        }
        SampleMatrix<T> operator - (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r = r - m;
            return r;
        }
        SampleMatrix<T> operator * (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r = (r.matrix() * m.matrix()).array();
            return r;
        }
        SampleMatrix<T> operator / (const SampleMatrix<T> & m) {
            assert(*this == m);
            SampleMatrix<T> r(*this);
            r =  r / m;
            return r;
        }
        
        SampleMatrix<T> operator + (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r + s;        
            return r;
        }
        SampleMatrix<T> operator - (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r - s;
            return r;
        }
        SampleMatrix<T> operator * (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r * s;
            return r;
        }
        SampleMatrix<T> operator / (const T s) {            
            SampleMatrix<T> r(*this);        
            r = r / s;
            return r;
        }
        

        SampleMatrix<T> abs() { 
            SampleMatrix<T> r(*this);        
            r = abs();
            return r;
        }
        SampleMatrix<T> abs2() { 
            SampleMatrix<T> r(*this);        
            r = abs2();
            return r;
        }
        SampleMatrix<T> inverse() { 
            SampleMatrix<T> r(*this);        
            r = inverse();
            return r;
        }
        SampleMatrix<T> exp() { 
            SampleMatrix<T> r(*this);        
            r = exp();
            return r;
        }
        SampleMatrix<T> log() { 
            SampleMatrix<T> r(*this);        
            r = log();
            return r;
        }
        SampleMatrix<T> log1p() { 
            SampleMatrix<T> r(*this);        
            r = log1p();
            return r;
        }
        SampleMatrix<T> log10() { 
            SampleMatrix<T> r(*this);        
            r = log10();
            return r;
        }
        SampleMatrix<T> pow(const T s) { 
            SampleMatrix<T> r(*this);        
            r = pow(s);
            return r;
        }
        SampleMatrix<T> sqrt() {
            SampleMatrix<T> r(*this);        
            r = sqrt();
            return r;
        }
        SampleMatrix<T> rsqrt() {
            SampleMatrix<T> r(*this);        
            r = rsqrt();
            return r;
        }
        SampleMatrix<T> square() {
            SampleMatrix<T> r(*this);        
            r = square();
            return r;
        }
        SampleMatrix<T> sin() {
            SampleMatrix<T> r(*this);        
            r = sin();
            return r;
        }
        SampleMatrix<T> cos() {
            SampleMatrix<T> r(*this);        
            r = cos();
            return r;
        }
        SampleMatrix<T> tan() {
            SampleMatrix<T> r(*this);        
            r = tan();
            return r;
        }
        SampleMatrix<T> asin() {
            SampleMatrix<T> r(*this);        
            r = asin();
            return r;
        }
        SampleMatrix<T> acos() {
            SampleMatrix<T> r(*this);
            r = acos();
            return r;
        }
        SampleMatrix<T> atan() {
            SampleMatrix<T> r(*this);
            r = atan();
            return r;
        }
        SampleMatrix<T> sinh() {
            SampleMatrix<T> r(*this);        
            r = sinh();
            return r;
        }
        SampleMatrix<T> cosh() {
            SampleMatrix<T> r(*this);        
            r = cosh();
            return r;
        }
        SampleMatrix<T> tanh() {
            SampleMatrix<T> r(*this);        
            r = tanh();
            return r;
        }
        SampleMatrix<T> ceil() {
            SampleMatrix<T> r(*this);        
            r = ceil();
            return r;
        }
        SampleMatrix<T> floor() {
            SampleMatrix<T> r(*this);        
            r = floor();
            return r;
        }
        SampleMatrix<T> round() {
            SampleMatrix<T> r(*this);        
            r = round();
            return r;
        }
    };

    template<typename T>  SampleVector<T> abs( SampleVector<T> & m) { return m.abs(); }
    template<typename T>  SampleVector<T> abs2( SampleVector<T> & m) { return m.abs2(); }
    template<typename T>  SampleVector<T> inverse( SampleVector<T> & m) { return m.inverse(); }
    template<typename T>  SampleVector<T> exp( SampleVector<T> & m) { return m.exp(); }
    template<typename T>  SampleVector<T> log( SampleVector<T> & m) { return m.log(); }
    template<typename T>  SampleVector<T> log1p( SampleVector<T> & m) { return m.log1p(); }
    template<typename T>  SampleVector<T> log10( SampleVector<T> & m) { return m.log10(); }
    template<typename T>  SampleVector<T> pow( SampleVector<T> & m,  SampleVector<T> & p) { return m.pow(p); }
    template<typename T>  SampleVector<T> pow( SampleVector<T> & m,  T p) { return m.pow(p); }
    template<typename T>  SampleVector<T> sqrt( SampleVector<T> & m) { return m.sqrt(); }
    template<typename T>  SampleVector<T> rsqrt( SampleVector<T> & m) { return m.rsqrt(); }
    template<typename T>  SampleVector<T> square( SampleVector<T> & m) { return m.square(); }
    template<typename T>  SampleVector<T> sin( SampleVector<T> & m) { return m.sin(); }
    template<typename T>  SampleVector<T> cos( SampleVector<T> & m) { return m.cos(); }
    template<typename T>  SampleVector<T> tan( SampleVector<T> & m) { return m.tan(); }
    template<typename T>  SampleVector<T> asin( SampleVector<T> & m) { return m.asin(); }
    template<typename T>  SampleVector<T> acos( SampleVector<T> & m) { return m.acos(); }
    template<typename T>  SampleVector<T> atan( SampleVector<T> & m) { return m.atan(); }
    template<typename T>  SampleVector<T> sinh( SampleVector<T> & m) { return m.sinh(); }
    template<typename T>  SampleVector<T> cosh( SampleVector<T> & m) { return m.cosh(); }
    template<typename T>  SampleVector<T> tanh( SampleVector<T> & m) { return m.tanh(); }
    template<typename T>  SampleVector<T> ceil( SampleVector<T> & m) { return m.ceil(); }
    template<typename T>  SampleVector<T> floor( SampleVector<T> & m) { return m.floor(); }
    template<typename T>  SampleVector<T> round( SampleVector<T> & m) { return m.round(); }

    template<typename T>  SampleMatrix<T> abs( SampleMatrix<T> & m) { return m.abs(); } 
    template<typename T>  SampleMatrix<T> abs2( SampleMatrix<T> & m) { return m.abs2(); }
    template<typename T>  SampleMatrix<T> inverse( SampleMatrix<T> & m) { return m.inverse(); }
    template<typename T>  SampleMatrix<T> exp( SampleMatrix<T> & m) { return m.exp(); }
    template<typename T>  SampleMatrix<T> log( SampleMatrix<T> & m) { return m.log(); }
    template<typename T>  SampleMatrix<T> log1p( SampleMatrix<T> & m) { return m.log1p(); }
    template<typename T>  SampleMatrix<T> log10( SampleMatrix<T> & m) { return m.log10(); }
    template<typename T>  SampleMatrix<T> pow( SampleMatrix<T> & m,  T &p) { return m.pow(p); }
    template<typename T>  SampleMatrix<T> sqrt( SampleMatrix<T> & m) { return m.sqrt(); }
    template<typename T>  SampleMatrix<T> rsqrt( SampleMatrix<T> & m) { return m.rsqrt(); }
    template<typename T>  SampleMatrix<T> square( SampleMatrix<T> & m) { return m.square(); }
    template<typename T>  SampleMatrix<T> sin( SampleMatrix<T> & m) { return m.sin(); }
    template<typename T>  SampleMatrix<T> cos( SampleMatrix<T> & m) { return m.cos(); }
    template<typename T>  SampleMatrix<T> tan( SampleMatrix<T> & m) { return m.tan(); }
    template<typename T>  SampleMatrix<T> asin( SampleMatrix<T> & m) { return m.asin(); }
    template<typename T>  SampleMatrix<T> acos( SampleMatrix<T> & m) { return m.acos(); }
    template<typename T>  SampleMatrix<T> atan( SampleMatrix<T> & m) { return m.atan(); }
    template<typename T>  SampleMatrix<T> sinh( SampleMatrix<T> & m) { return m.sinh(); }
    template<typename T>  SampleMatrix<T> cosh( SampleMatrix<T> & m) { return m.cosh(); }
    template<typename T>  SampleMatrix<T> tanh( SampleMatrix<T> & m) { return m.tanh(); }
    template<typename T>  SampleMatrix<T> ceil( SampleMatrix<T> & m) { return m.ceil(); }
    template<typename T>  SampleMatrix<T> floor( SampleMatrix<T> & m) { return m.floor(); }
    template<typename T>  SampleMatrix<T> round( SampleMatrix<T> & m) { return m.round(); }
    
    template<typename T> 
    SampleMatrix<T> deinterleave(const SampleVector<T> & v) {
        SampleMatrix<T> m;
        m.deinterleave(v);
        return m;
    }

    template<typename T> 
    SampleVector<T> interleave(SampleMatrix<T> & m) {
        return m.interleave();
    }

    
    template<typename T>
    struct Window : public SampleVector<T>
    {       
        Window(size_t i) : SampleVector<T>(i) {}
        virtual ~Window() = default;     
    };

    template<typename T>
    struct Rectangle: public Window<T>
    {
        Rectangle(size_t i) : Window<T>(i) { this->fill(1); } 
    };

    template<typename T>
    struct Hamming: public Window<T>
    {
        Hamming(size_t i) : Window<T>(i) {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)(i) = 0.54 - (0.46 * std::cos(2*M_PI*i/n));
            }        
        }
    };

    template<typename T>
    struct Hanning: public Window<T>
    {
        Hanning(size_t i) : Window<T>(i) {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)
            {
                (*this)(i) = 0.5*(1 - std::cos(2*M_PI*i/n));
            }        
        }
    };

    template<typename T>
    struct Blackman: public Window<T>
    {
        Blackman(size_t i) : Window<T>(i)    
        {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)                    
                (*this)(i) = 0.42 - (0.5* std::cos(2*M_PI*i/(n)) + (0.08*std::cos(4*M_PI*i/n)));        
        }
    };

    template<typename T>
    struct BlackmanHarris: public Window<T>
    {
        BlackmanHarris(size_t i) : Window<T>(i)    
        {
            T n = (*this).size()-1;
            for(size_t i = 0; i < this->size(); i++)            
                (*this)(i) = 0.35875 
                        - 0.48829*std::cos(2*M_PI*(i/n))
                        + 0.14128*std::cos(4.0*M_PI*(i/n)) 
                        - 0.01168*std::cos(6.0*M_PI*(i/n));
        }
    };

    template<typename T>
    struct Gaussian: public Window<T>
    {
        Gaussian(size_t i) : Window<T>(i)
        {
            T a,b,c=0.5;
            for(size_t n = 0; n < this->size(); n++)
            {
                a = (n - c*((*this).size()-1)/(std::sqrt(c)*(*this).size()-1));
                b = -c * std::sqrt(a);
                (*this)(n) = std::exp(b);
            }
        }
    };
    template<typename T>
    struct Welch: public Window<T>
    {
        Welch(size_t i) : Window<T>(i)
        {
            for(size_t i = 0; i < (*this).size(); i++)
                (*this)(i) = 1.0 - std::sqrt((2*i-(*this).size())/((*this).size()+1));        
        }
    };
    template<typename T>
    struct Parzen: public Window<T>
    {

        Parzen(size_t i) : Window<T>(i)
        {
            for(size_t i = 0; i < (*this).size(); i++)
                (*this)(i) = 1.0 - std::abs((2.0f*i-(*this).size())/((*this).size()+1.0f));        
        }    
    };
    template<typename T>
    struct Tukey: public Window<T>
    {
        Tukey(size_t i, T alpha) : Window<T>(i)
        {
            size_t num_samples = (*this).size();
            T value = (-1*(num_samples/2)) + 1;
            for(size_t i = 0; i < num_samples; i++)
            {    
                if(value >= 0 && value <= (alpha * (num_samples/2))) 
                    (*this)[i] = 1.0; 
                else if(value <= 0 && (value >= (-1*alpha*(num_samples/2)))) 
                    (*this)[i] = 1.0;
                else 
                    (*this)[i] = 0.5 * (1 + std::cos(M_PI *(((2*value)/(alpha*num_samples))-1)))        ;
                value = value + 1;
            }     
        }
    };

    template<typename T>
    void swap(SampleVector<T> & vector,size_t i, size_t j)
    {
        T x = vector(i);
        vector(i) = vector(j);
        vector(j) = x;
    }        

    template<typename T>
    void shift(SampleVector<T> & vector) 
    {
        size_t half = vector.size() / 2; 
        size_t start= half; 
        if(2*half < vector.size()) start++;
        for(size_t i = 0; i < half; i++)
        {
            swap(vector,i, i+start);
            
        }
        if(start != half)
        {
            for(size_t i = 0; i < half; i++)            
            {
                swap(vector,i+start-1,i+start);
            }
        }
    }

    template<typename T>    
    void ishift(SampleVector<T> & vector) 
    {
        size_t half = vector.size() / 2; 
        size_t start= half; 
        if(2*half < vector.size()) start++;
        for(size_t i = 0; i < half; i++)
        {
            swap(vector,i,i+start);
        }
        if(start != half)
        {
            for(size_t i = 0; i < half; i++)            
            {
                swap(vector,half,i);
            }
        }
    }

    template<typename T>    
    T quadratic_peak_pos(SampleVector<T> & vector,size_t pos)
    {
        T s0,s1,s2;
        size_t x0,x2;        
        if(pos == 0 || pos == vector.size()-1) return pos; 
        x0 = (pos < 1)? pos:pos-1;
        x2 = (pos+1 < vector.size()) ? pos+1:pos;
        if(x0 == pos) return vector(pos) <= vector(x2)? pos:x2;
        if(x2 == pos) return vector(pos) <= vector(x0)? pos:x0;
        s0 = vector(x0);
        s1 = vector(pos);
        s2 = vector(x2);
        return pos + 0.5 *(s0-s2) / (s0 - 2*s1 + s2);
    }

    template<typename T>    
    T quadratic_peak_mag(SampleVector<T> & vector,size_t pos)
    {
        T x0,x1,x2;
        size_t index = (size_t)(pos - 0.5) + 1;
        if(pos >= vector.size() || pos < (T)0) return (T)0;
        if((T)index == pos) return vector(index);
        x0 = vector(index-1);
        x1 = vector(index);
        x2 = vector(index+1);
        return x1 - 0.25 * (x0 - x2) * (pos - index);
    }

    template<typename T>    
    T median(SampleVector<T> & vector)
    {        
        size_t n = vector.size();
        size_t high,low;
        size_t median;
        size_t middle,ll,hh;

        low = 0;
        high = n-1;
        median = (low+high)/2;

        for(;;)        
        {
            if(high <= low)
                return vector(median);

            if(high == low+1)
            {
                if(vector(low) > vector(high))
                    swap(vector,low,high);
                return vector(median);
            }
            middle = (low+high)/2;
            if(vector(middle) > vector(high)) swap(vector,middle,high);
            if(vector(low) > vector(high)) swap(vector,low,high);
            if(vector(middle) > vector(low)) swap(vector,middle,low);
            swap(vector,middle,low+1);

            ll=low+1;
            hh=high;
            for(;;)
            {
                do ll++; while( vector(low) > vector(ll));
                do hh--; while( vector(hh) > vector(low));
                if(hh < ll) break;
                swap(vector,ll,hh);
            }
            swap(vector,low,hh);
            if(hh <= median) low = ll;
            if(hh >= median) high = hh-1;
        }
    }

    template<typename T>    
    T moving_threshold(SampleVector<T> & input, size_t post, size_t pre, size_t pos)
    {
        size_t length = input.size();
        size_t win_length = post+pre+1;
        SampleVector<T> vector = input;
        if(pos < post+1)
        {
            for(size_t k = 0; k < post +1 - pos; k++)
                vector(k) = 0;
            for(size_t k = post+1-pos; k < win_length; k++)
                vector(k) = vector(k+pos-post);            
        }
        else if(pos + pre < length)
        {
            for(size_t k = 0; k < win_length; k++)
                vector(k) = vector(k+pos-post);
        }
        else 
        {
            for(size_t k = 0; k < length - pos + post; k++)
                vector(k) = vector(k+pos-post);
            for(size_t k = length - pos + post; k < win_length; k++)
                vector(k) = 0;
        }        
        return median(vector);
    }

    template<typename T> 
    T zero_crossing_rate(SampleVector<T> & vector)
    {
        T zcr = 0;
        for(size_t i = 1; i < vector.size(); i++)
        {
            bool current = vector(i) > 0;
            bool prev    = vector(i-1) > 0;
            if(current != prev) zcr++;
        }    
        return zcr;   
    }
    
    template<typename T> 
    void autocorr(SampleVector<T> & vector,SampleVector<T> & output)
    {
        T tmp;
        output.resize(vector.size());
        for(size_t i = 0; i < vector.size(); i++)
        {
            tmp = (T)0;
            for(size_t j = 0; j < vector.size(); j++)
                tmp += vector(j-i) * vector(j);                
            
            output.vector(i) = tmp / (T)(vector.size()-1);
        }
    }

    
    template<typename T>
    void push(SampleVector<T> & vector,const T& new_elem)
    {
        for(size_t i = 0; i < vector.size()-1; i++) vector(i) = vector(i+1);
        vector(vector.size()-1) = new_elem;
    }

    template<typename T>
    void clamp(SampleVector<T> & vector,T absmax) { 
        for(size_t i = 0; i < vector.size(); i++)
        {
            if( vector(i) > absmax) vector(i) = absmax;
            else if(vector(i) < -absmax) vector(i) = -absmax;
        }
    }

    template<typename T> 
    void  normalize(SampleVector<T> & vector) { 
            Eigen::Matrix<T,1,Eigen::Dynamic> r = vector.matrix();
            r.normalize();
            vector = r.array();         
    }


    template<typename T> 
    bool peakpick(SampleVector<T> & vector,size_t pos)
    {
        bool tmp = false;
        tmp = (vector(pos) > vector(pos-1)) && (vector(pos) > vector(pos+1)) && vector(pos) > 0;
        return tmp;
    }
    
    template<typename T> 
    T RMS(SampleVector<T> & vector)
    {
        T sum = vector.pow(2).sum();        
        return std::sqrt(sum/vector.size());
    }

    template<typename T> 
    T peak_energy(SampleVector<T> & vector)
    {
        T peak = -10000;
        T as;         
        for(size_t i = 0; i < vector.size(); i++)
        {
            T as = std::fabs(vector(i));
            if(as > peak) peak = as;
        }
        return peak;
    }

    template<typename T> 
    T min(SampleVector<T> & vector) { 
        T min = 1e120;
        for(size_t i = 0; i < vector.size(); i++)
            if(vector(i) < min) min = vector(i);
        return min;            
        //return vector(vector.minCoeff()); 
    }

    template<typename T> 
    T max(SampleVector<T> & vector) { 
        T max = -1e120;
        for(size_t i = 0; i < vector.size(); i++)
            if(vector(i) > max) max = vector(i);
        return max;
        
    }

    template<typename T> 
    SampleVector<T>& weighted_copy(SampleVector<T> & vector,const SampleVector<T> & weight, SampleVector<T> & out)
    {
        out.vector = vector * weight.vector;
        return out;
    }    

    template<typename T> 
    T level_lin(SampleVector<T> & vector) 
    {
        T energy = vector.sqrt().sum();
        return energy/vector.size();
    }

    template<typename T> 
    T local_hfc(SampleVector<T> & vector)
    {
        T hfc = 0;
        for(size_t j = 0; j < vector.size(); j++)
            hfc += (j+1) * vector(j);
        return hfc;
    }

    template<typename T> 
    void min_removal(SampleVector<T> & vector)
    {
        T m = min(vector);
        vector -= m;
    }

    template<typename T> 
    T alpha_norm(SampleVector<T> & vector,SampleVector<T> & out, T alpha)
    {
        T tmp = vector.abs().pow(alpha).sum();
        return std::pow(tmp/vector.size(),1.0/alpha);
    }

    template<typename T> 
    void alpha_normalize(SampleVector<T> & vector,SampleVector<T> & out, T alpha)
    {
        T tmp = alpha_norm(out,alpha);
        out.vector = out.vector / tmp;         
    }    

    template<typename T> 
    T unwrap2pi(T phase)
    {
        return phase + 2*M_PI*(1.0 + std::floor(-(phase+M_PI)/(2*M_PI)));
    }

    template<typename T> 
    T quadfrac(T s0, T s1, T s2, T pf)
    {
        T tmp = s0 + (pf/2)*(pf*(s0-2.0*s1+s2)-3.0*s0+4*s1-s2);
        return tmp;
    }

    template<typename T> 
    T freqtomidi(T freq)
    {
        T midi;
        if(freq < 2.0 || freq > 100000) return 0;        
        midi = 12*std::log(midi/6.875) / 0.6931471805599453 -3;
        return midi;
    }

    template<typename T> 
    T miditofreq(T midi)
    {
        if(midi > 140) return 0;
        T freq = std::exp(((midi + 3) / 12)*0.6931471805599453)*6.875;
        return freq;
    }

    template<typename T> 
    T bintofreq(T bin, T sample_rate, T fft_size)
    {
        T freq = sample_rate / fft_size;
        return freq * std::max(bin,(T)0);
    }


    template<typename T> 
    T bintomidi(T bin, T sample_rate, T fft_size)
    {
        T midi = bintofreq(bin,sample_rate,fft_size);
        return freqtomidi(midi);
    }

    template<typename T> 
    T freqtobin(T freq, T sample_rate, T fft_size)
    {
        T bin = fft_size / sample_rate;
        return std::max(freq,(T)0)*bin;
    }

    template<typename T> 
    T miditobin(T midi, T sample_rate, T fft_size)
    {
        T freq = miditofeq(midi);
        return freqtobin(freq,sample_rate, fft_size);
    }
    template<typename T> 
    bool is_power_of_two(uint64_t a)
    {
        if((a & (a-1))==0)
            return true;
        return false;
    }
    template<typename T> 
    uint64_t next_power_of_two(uint64_t a)
    {
        uint64_t i = 1;
        while(i < a) i <<= 1;
        return i;
    }

    template<typename T> 
    T hztomel(T freq)
    {
        T lin_space = (T)3/(T)200;
        T split_hz = (T)1000;
        T split_mel = split_hz * lin_space;
        T log_space = (T)27 / std::log(6.4);
        assert(freq >= (T)0);
        if(freq < split_hz)        
            return freq * lin_space;        
        else
            return split_mel + log_space * std::log(freq/split_hz);
    }
    template<typename T> 
    T hztomel_htk(T freq)
    {        
        T split_hz = (T)700;        
        T log_space = (T)1127;
        assert(freq >= (T)0);
        return log_space * std::log(1 + freq/split_hz);
    }
        
    template<typename T> 
    T meltohz(T mel)
    {
        T lin_space = (T)200/(T)3;
        T split_hz  = (T)1000;
        T split_mel = split_hz / lin_space;
        T log_spacing = std::pow(6.4,1.0/27.0);
        assert(mel >= 0);
        if( mel < split_mel) return lin_space * mel;
        return split_hz * std::pow(log_spacing, mel-split_mel);
    }

    template<typename T> 
    T meltohz_htk(T mel)
    {
        T split_hz = (T)700;
        T log_space = (T)1/(T)1127;
        assert(mel >= 0);
        return split_hz * (std::exp(mel *log_space) -(T)1);
    }

    template<typename T> 
    T db_spl(SampleVector<T> & vector) { return 10 * std::log10(level_lin(vector)); }

    template<typename T> 
    T level_detect(SampleVector<T> & vector,T threshold) { 
        T db = db_spl(vector);
        if(db < threshold) return 1;
        return db;
    }

    template<typename T> 
    size_t size(SampleVector<T> & v) { return v.size(); }
    
    template<typename T> 
    SampleVector<T> random(size_t n, size_t channels=1) { 
        SampleVector<T> r(n,channels);
        r.random();
        return r;
    }

    template<typename T> 
    void fill(SampleVector<T> & v, T value) { v.fill(value); }

    template<typename T> 
    void resize(SampleVector<T> & v, size_t n, size_t channels=1) { v.resize(n,channels); }

    template<typename T> 
    bool    silence_detect(SampleVector<T> & vector,T threshold) { return db_spl(vector) < threshold; }

    template<typename T>
    void set_left_channel(const std::vector<T> & left, SampleVector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_left_channel(const SampleVector<T> & left, SampleVector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_left_channel(const SampleVector<T> & left, std::vector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const std::vector<T> & right, SampleVector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_right_channel(const SampleVector<T> & right, SampleVector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_right_channel(const SampleVector<T> & right, std::vector<T> & out) {
        size_t x = 0;
        for(size_t i = 0; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const std::vector<T> & in, SampleVector<T> & out) {
        size_t x = 0;
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    SampleVector<T> interleave(size_t n, size_t channels, const std::vector<SampleVector<T>> & in) {
        std::vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    std::vector<SampleVector<T>> deinterleave(size_t n, size_t channels, const SampleVector<T> & in) {
        std::vector<SampleVector<T>> r(n);
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }
};