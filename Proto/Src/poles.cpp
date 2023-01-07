#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "PolynomialRoots.hpp"
#include "LaguerrePolynomialRoots.hpp"
#include "SndFile.hpp"
#include "Undenormal.hpp"


typedef double DspFloatType;


int TestLaguerre()
{
    deque<complex<DspFloatType>> P = {1,2,3}, R;
  LaguerreMethod L(P);
  R = L.solve_roots();

  // Display the equation to solve
  for (int i = P.size()-1; i >=1 ; i--) cout << P[i] <<"*x^" << i << " + ";
  cout << P[0] << "= 0"<< endl;

  //Display roots
  cout << "--------- ROOTS ---------" << endl;
  for (int i = 0; i < R.size(); i++) cout << R[i] << endl;
  
  //In order to gather the quality of the roots found, the lines below evaluate the polynom at every
  //root, R_i, and takes the maximum and the mean deviation to zero 
  cout << endl << endl;
  complex<DspFloatType> P_x(0.,0.);
  DspFloatType mean_err = 0, max_err = -1e12, abs_err;
  cout << "--------- Root analysis ---------" << endl;
  for (int j = 0; j < R.size(); j++)
  { 
    P_x = complex<DspFloatType>(0,0);
    for (int i = 0; i < P.size(); i++) P_x += P[i]*pow(R[j], i);
    abs_err = abs(P_x);
    if (abs(P_x) > max_err) max_err = abs_err;
    mean_err += abs_err;
    cout << "Error{Root[" << j << "]}= " << abs_err << endl;
  }
  cout << endl;
  cout << "Mean error = sum(|P(R_i)|, i={1, N}/N = " << mean_err/R.size() << endl;
  cout << "Maximum error = max(|P(R_i)|), i={1, N}/N = " << max_err << endl;
  return 0;
}


// -(p1 + p2)s
DspFloatType complex_pole1(std::complex<DspFloatType> p1, std::complex<DspFloatType> p2)
{
    return abs(p1 + p2);
}

// (p1*p2)
DspFloatType complex_pole2(std::complex<DspFloatType> p1, std::complex<DspFloatType> p2)
{
    return abs(p1 * p2);
}

std::complex<DspFloatType> butterworthpole(DspFloatType k, DspFloatType n)
{
    DspFloatType p = M_PI * ((2 * k + n - 1) / (2 * n));
    return std::complex<DspFloatType>(-std::cos(p), -std::sin(p));
}
std::complex<DspFloatType> ButterworthPoles(DspFloatType K, DspFloatType N)
{
    DspFloatType theta = ((2*K+N-1)/(2*N))*M_PI;
    DspFloatType sk = cos(theta);
    DspFloatType omegak = -sin(theta);
    std::complex<DspFloatType> p(sk,omegak);
    return p;
}
std::complex<DspFloatType> ChebyshevH0(DspFloatType N, DspFloatType r=1.0)
{        
    DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*1-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> P(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
    for(size_t K=2; K <= N; K++)
    {
        e     = sqrt(pow(10.0,r/10.0)-1.0);
        theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
        phi   = (1.0/N)*asinh(1.0/e);
        std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
        P *= -p;        
    }
    if(fmod(N,2) == 0) return P/sqrt(1 + e*e);
    return P;
}

// Chebyshev2 = 1/pole
std::complex<DspFloatType> ChebyshevPole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return p;
}


// Chebyshev2 = 1/pole
std::complex<DspFloatType> Chebyshev2Zeros(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType uk    = ((2*K-1)/N)*(M_PI/2.0);
    std::complex<DspFloatType> p(0,cos(uk));
    return 1.0/p;
}

// Chebyshev2 = 1/pole
std::complex<DspFloatType> Chebyshev2Pole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
{      
    DspFloatType e     = 1.0/sqrt(pow(10.0,r/10.0)-1.0);
    DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
    DspFloatType phi   = (1.0/N)*asinh(1.0/e);
    std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
    return 1.0/p;
}

std::vector<std::complex<DspFloatType>> ChebyshevPoles(DspFloatType N, DspFloatType r)
{
    std::vector<std::complex<DspFloatType>> out(N);
    for(size_t K = 1; K <= N; K++)
    {
        out.push_back(ChebyshevPole(K,N,r));
    }
    return out;
}

// Hn(s) = H0/d * PROD( (s^2 + ai) / (s^2 + bi*s + ci))
// d = s+p0 n=odd
// d = 1    n=even


DspFloatType factorial(DspFloatType x)
{
    if (x == 0)
        return 1;
    return x * factorial(x - 1);
}
DspFloatType binomial(DspFloatType n, DspFloatType k)
{
    if(n == 0 || k == 0) return 1;
    return factorial(n) / (factorial(k) * factorial(fabs(n - k)));
}

DspFloatType bk(DspFloatType n, DspFloatType k) {
    DspFloatType num = factorial(2*n-k);
    DspFloatType den = pow(2.0,n-k)*factorial(k)*factorial(n-k);
    return num/den;
}

// calculate bessel filter coefficients qn(s) order n 
// it will return the polynomial coefficients starting at b0
// b0,b1*s,b2*s^2 + ... bn*s^n
// these use the polynomial root solver to find the poles (complex)
std::vector<DspFloatType> qn(size_t n) 
{
    std::vector<DspFloatType> r(n+1);
    r[0] = bk(n,0);
    for(size_t k=1; k <= n; k++) {
        DspFloatType b = bk(n,k);
        r[k] = b;
    }
    return r;
}



// this will calculate from 2 to 8
// for higher order need to calculate the bessel functions
// and find the roots of the polynomial for the poles.
// it seems bessel filter does not to preserve group delay in digital
// but the tests for these were ok
// dont know why exactly they say that it doesn't have much problems for audio

#define MAXORDER 9
void bessel_coefficients(int order, std::vector<DspFloatType> & coeff, bool D=true)
{
    int i,N,index,indexM1,indexM2;
    DspFloatType B[3][MAXORDER];
    DspFloatType A,renorm[MAXORDER];
    renorm[2] = 0.72675;
    renorm[3] = 0.57145;
    renorm[4] = 0.46946;
    renorm[5] = 0.41322;
    renorm[6] = 0.37038;
    renorm[7] = 0.33898;
    renorm[8] = 0.31546;
    index = 1;
    indexM1 = 0;
    indexM2 = 2;

    memset(B,0,3*MAXORDER*sizeof(DspFloatType));
    B[0][0] = 1.0;
    B[1][0] = 1.0;
    B[1][1] = 1.0;
    coeff.resize(order+1);
    for(N=2; N <= order; N++)
    {
        index = (index+1)%3;
        indexM1= (indexM1 + 1)%3;
        indexM2= (indexM2 + 1)%3;
        for(i=0; i < N; i++) B[index][i] = (2*N-1) * B[indexM1][i];
        for(i=2; i <= N; i++) B[index][i] = B[index][i] + B[indexM2][i-2];
        if(D) 
            for(i=0; i <= order; i++) coeff[i] = B[index][i];
        else 
            for(i=0; i <= order; i++) coeff[i] = B[index][i] * pow(A,order-i);
    }
}

// precomputed poles for 2 to 8 order bessel filter
std::complex<DspFloatType> bessel_poles_2[] = { std::complex<DspFloatType>(-1.5,0.8660), 
                                          std::complex<DspFloatType>(-1.5,-0.8660)};
std::complex<DspFloatType> bessel_poles_3[] = { std::complex<DspFloatType>(-2.3222,0),
                                          std::complex<DspFloatType>(-1.8390,1.7543),
                                          std::complex<DspFloatType>(-1.8390,-1.7543)};                                          
std::complex<DspFloatType> bessel_poles_4[] = { std::complex<DspFloatType>(-2.1039,2.6575), 
                                          std::complex<DspFloatType>(-2.1039,-2.6575), 
                                          std::complex<DspFloatType>(-2.8961,0.8672),                                          
                                          std::complex<DspFloatType>(-2.8961,-0.8672)};
std::complex<DspFloatType> bessel_poles_5[] = { std::complex<DspFloatType>(-3.6467,0),
                                          std::complex<DspFloatType>(-2.3247,3.5710),
                                          std::complex<DspFloatType>(-2.3247,-3.5710),
                                          std::complex<DspFloatType>(-3.3520,1.7427),
                                          std::complex<DspFloatType>(-3.3520,-1.7427)};
std::complex<DspFloatType> bessel_poles_6[] = { std::complex<DspFloatType>(-2.5158,4.4927),
                                          std::complex<DspFloatType>(-2.5158,-4.4927),
                                          std::complex<DspFloatType>(-3.7357,2.6263),
                                          std::complex<DspFloatType>(-3.7357,-2.6263),
                                          std::complex<DspFloatType>(-4.2484,0.8675),
                                          std::complex<DspFloatType>(-4.2484,-0.8675)};
std::complex<DspFloatType> bessel_poles_7[] = { std::complex<DspFloatType>(-4.9716,0),
                                          std::complex<DspFloatType>(-2.6857,5.4206), 
                                          std::complex<DspFloatType>(-2.6857,5.4206),
                                          std::complex<DspFloatType>(-4.0701,3.5173),
                                          std::complex<DspFloatType>(-4.0701,-3.5173),
                                          std::complex<DspFloatType>(-4.7584,1.7393),
                                          std::complex<DspFloatType>(-4.7584,-1.7393)};
std::complex<DspFloatType> bessel_poles_8[] = { std::complex<DspFloatType>(-5.2049,2.6162),
                                          std::complex<DspFloatType>(-5.2049,2.6162),
                                          std::complex<DspFloatType>(-4.3683,4.4146),
                                          std::complex<DspFloatType>(-4.3683,-4.4146),
                                          std::complex<DspFloatType>(-2.3888,6.3540),
                                          std::complex<DspFloatType>(-2.3888,-6.3540),
                                          std::complex<DspFloatType>(-5.5878,0.8676),
                                          std::complex<DspFloatType>(-5.5878,-0.8676) };




struct FilterCoefficients
{
    DspFloatType a[2];
    DspFloatType b[3];
};

struct BiquadSection
{
    DspFloatType z[3];
    DspFloatType p[3];

    BiquadSection()
    {
        memset(z, 0, sizeof(z));
        memset(p, 0, sizeof(p));
    }
    BiquadSection(const FilterCoefficients &c)
    {
        z[0] = c.b[0];
        z[1] = c.b[1];
        z[2] = c.b[2];
        p[0] = c.a[0];
        p[1] = c.a[1];
    }
    BiquadSection(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
    {
        z[0] = z1;
        z[1] = z2;
        z[2] = z3;
        p[0] = p1;
        p[1] = p2;
    }
    BiquadSection(const BiquadSection &b)
    {
        //memcpy(z, b.z, sizeof(z));
        //memcpy(p, b.p, sizeof(p));
        z[0] = b.z[0];
        z[1] = b.z[1];
        z[2] = b.z[2];
        p[0] = b.p[0];
        p[1] = b.p[1];
        p[2] = b.p[2];
    }
    void setCoefficients(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
    {
        z[0] = z1;
        z[1] = z2;
        z[2] = z3;
        p[0] = p1;
        p[1] = p2;
    }
    void setCoefficients(DspFloatType n[3], DspFloatType d[2])
    {
        //memcpy(z, n, sizeof(z));
        //memcpy(p, d, sizeof(p));
        z[0] = n[0];
        z[1] = n[1];
        z[2] = n[2];
        p[0] = d[0];
        p[1] = d[1];
    }
    void setCoefficients(const FilterCoefficients &c)
    {
        z[0] = c.b[0];
        z[1] = c.b[1];
        z[2] = c.b[2];
        p[0] = c.a[0];
        p[1] = c.a[1];
    }
    BiquadSection &operator=(const BiquadSection &b)
    {
        //memcpy(z, b.z, sizeof(z));
        //memcpy(p, b.p, sizeof(p));
        z[0] = b.z[0];
        z[1] = b.z[1];
        z[2] = b.z[2];
        p[0] = b.p[0];
        p[1] = b.p[1];
        p[2] = b.p[2];
        return *this;
    }

    void print()
    {
        std::cout << z[0] << "," << z[1] << "," << z[2] << std::endl;
        std::cout << p[0] << "," << p[1] << "," << p[2] << std::endl;
        
    }
};


using BiquadSOS = std::vector<BiquadSection>;

struct BiquadTransposedTypeII 
{
    BiquadSection biquad;
    DspFloatType x, y, d1, d2;

    BiquadTransposedTypeII() 
    {
        x = y = 0;
        d1 = d2 = 0;
    }
    BiquadTransposedTypeII(const BiquadSection &b) : biquad(b)
    {
        x = y = 0;
        d1 = d2 = 0;
    }
    BiquadTransposedTypeII &operator=(const BiquadTransposedTypeII &b)
    {
        biquad = b.biquad;
        x = b.x;
        y = b.y;
        d1 = b.d1;
        d2 = b.d2;
        return *this;
    }
    void setCoefficients(const BiquadSection &b)
    {
        biquad = b;
    }        
    void setBiquad(const BiquadSection &b)
    {
        biquad = b;
    }

    // transposed is just flip - to +
    DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
    {
        Undenormal denormal;
        x = I;
        y = biquad.z[0] * x + d1;
        d1 = biquad.z[1] * x - biquad.p[0] * y + d2;
        d2 = biquad.z[2] * x - biquad.p[1] * y;
        return A * y;
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
    }
};

struct FilterBase
{
    
    std::vector<BiquadTransposedTypeII> biquads;
    size_t order;
    DspFloatType fc,sr,R,q,bw,g,ripple,rolloff,stop,pass;
    bool init = false;

    FilterBase() 
    {
    }
    FilterBase(const BiquadSOS &s)
    {
        setCoefficients(s);
    }
    void setCoefficients(const BiquadSOS &s)
    {        
        if(s.size() == 0) {
            init = false;
            return;
        }
        biquads.resize(s.size());
        for (size_t i = 0; i < s.size(); i++)
        {
            biquads[i].setCoefficients(s[i]);
        }
        init = true;
    }
    DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
    {
        if(!init) return 0;
        DspFloatType o = biquads[0].Tick(I, A, X, Y);
        for (size_t i = 1; i < biquads.size(); i++)
            o = biquads[i].Tick(o, A, X, Y);
        return A * o;
    }

    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A) {
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i]);
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A, DspFloatType * X) {
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i],X[i]);
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A, DspFloatType * X, DspFloatType * Y) {
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i],X[i],Y[i]);
    }
};

void prewarp(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs)
{
    DspFloatType wp, pi;

    pi = 4.0 * std::atan(1.0);
    wp = 2.0 * fs * std::tan(pi * fc / fs);

    *a2 = (*a2) / (wp * wp);
    *a1 = (*a1) / wp;
}
void prewarpR(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType R)
{
    DspFloatType wp, pi;

    pi = 4.0 * std::atan(1.0);
    wp = 2.0 * fs * std::tan(pi * fc / fs);

    *a2 = R * R * (*a2) / (wp * wp);
    *a1 = R * (*a1) / wp;
}
void prewarpQ(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType Q)
{
    DspFloatType wp, pi;

    pi = 4.0 * std::atan(1.0);
    wp = 2.0 * fs * std::tan(pi * fc / fs);

    *a2 = (*a2) / (Q * Q * wp * wp);
    *a1 = (*a1) / (Q * wp);
}
void prewarpRQ(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType R, DspFloatType Q)
{
    DspFloatType wp, pi;

    pi = 4.0 * std::atan(1.0);
    wp = 2.0 * fs * std::tan(pi * fc / fs);

    *a2 = R * R * (*a2) / (Q * Q * wp * wp);
    *a1 = R * (*a1) / (Q * wp);
}

void inversebilinear(
    DspFloatType z[3], DspFloatType p[3],
    DspFloatType k,    /* overall gain factor */
    DspFloatType fs,   /* sampling rate */
    DspFloatType *coef /* pointer to 4 iir coefficients */
)
{
    DspFloatType ad, bd;
    DspFloatType b0 = k;
    DspFloatType b1 = coef[0];
    DspFloatType b2 = coef[1];
    DspFloatType a0 = 1;
    DspFloatType a1 = coef[2];
    DspFloatType a2 = coef[3];

    ad = 1 / 4. * a2 / (fs * fs) + a1 / (2 * fs) + a0;
    bd = 1 / 4. * b2 / (fs * fs) + b1 / (2 * fs) + b0;

    z[0] = k * bd / ad;
    z[1] = b1 * bd;
    z[2] = b2 * bd;
    p[0] = 1;
    p[1] = a1 * ad;
    p[2] = a2 * ad;
}

void bilinear(
    DspFloatType a0, DspFloatType a1, DspFloatType a2, /* numerator coefficients */
    DspFloatType b0, DspFloatType b1, DspFloatType b2, /* denominator coefficients */
    DspFloatType *k,                       /* overall gain factor */
    DspFloatType fs,                       /* sampling rate */
    DspFloatType *coef                     /* pointer to 4 iir coefficients */
)
{
    DspFloatType ad, bd;

    
    /* alpha (Numerator in s-domain) */
    ad = 4. * a2 * fs * fs + 2. * a1 * fs + a0;
    /* beta (Denominator in s-domain) */
    bd = 4. * b2 * fs * fs + 2. * b1 * fs + b0;

    /* update gain constant for this section */
    *k *= ad / bd;
   
    
    /* Nominator */
    *coef++ = (2. * a0 - 8. * a2 * fs * fs) / ad;         /* alpha1 */
    *coef++ = (4. * a2 * fs * fs - 2. * a1 * fs + a0) / ad; /* alpha2 */

    *coef++ = (2. * b0 - 8. * b2 * fs * fs) / bd;           /* beta1 */
    *coef++ = (4. * b2 * fs * fs - 2. * b1 * fs + b0) / bd; /* beta2 */

    
}

// frequency transform
// lp(s/wc) => ( (s/wc) - p1) * (s/wc) - p2) => (s/wc)^2 - (s/wc)*p1 - (s/wc)*p2 + p1*p2
// hp(wc/s) => ( (wc/s) - p1) * (wc/s) - p2) => (wc/s)^2 - (wc/s)*p1 - (wc/s) *p2 + p1*p2
//          => (wc^2/s^2 - (wc/s)*p1 - (wc/s)*p2 +p1p2)
//          => s^2 / (wc^2 - wc*s*p1 - wc*s*p2 + p1p2*s^2)


void lp2lp(std::vector<std::complex<DspFloatType>> & poles,
            DspFloatType wc, DspFloatType & gain) {
        
    gain *= pow(wc,poles.size());
    for(size_t i = 0; i < poles.size(); i++)
        poles[i] = wc * poles[i];
        
}
void lp2hp(std::vector<std::complex<DspFloatType>> & zeros,
            std::vector<std::complex<DspFloatType>> & poles,
            DspFloatType wc, DspFloatType & gain)
{
    std::complex<DspFloatType> prodz(1.0,0.0),prodp(1.0,0.0);
    for(size_t i = 0; i < zeros.size(); i++) prodz *= -zeros[i];
    for(size_t i = 0; i < poles.size(); i++) prodp *= -poles[i];
    gain *= prodz.real() / prodp.real();
    for(size_t i = 0; i < poles.size(); i++)
        if(abs(poles[i])) poles[i] = std::complex<DspFloatType>(wc) / poles[i];
    zeros.resize(poles.size());
    for(size_t i = 0; i < zeros.size(); i++)
        zeros[i] = std::complex<DspFloatType>(0.0);
}               
void lp2bp(std::vector<std::complex<DspFloatType>> & zeros,
            std::vector<std::complex<DspFloatType>> & poles,
            DspFloatType wu, DspFloatType wl, DspFloatType & gain)
{
    DspFloatType wc = sqrt(wu*wl);
    DspFloatType bw = wu-wl;
    gain      *= pow(bw,poles.size()-zeros.size());
    std::vector<std::complex<DspFloatType>> temp;
    for(size_t i = 0; i < poles.size(); i++) 
    {
        if(abs(poles[i])) {
            std::complex<DspFloatType> first = DspFloatType(0.5) * poles[i] * bw;
            std::complex<DspFloatType> second= DspFloatType(0.5) * sqrt(bw*bw) * (poles[i]*poles[i]-DspFloatType(4.0)*wc*wc);
            temp.push_back(first + second);
        }
    }
    for(size_t i = 0; i < poles.size(); i++) {
        if(abs(poles[i])) {
            std::complex<DspFloatType> first = DspFloatType(0.5) * poles[i] * bw;
            std::complex<DspFloatType> second= DspFloatType(0.5) * sqrt(bw*bw) * (poles[i]*poles[i]-DspFloatType(4.0)*wc*wc);
            temp.push_back(first - second);
        }
    }
    zeros.resize(poles.size());
    for(size_t i = 0; i < zeros.size(); i++) {
        zeros[i] = std::complex<DspFloatType>(0);
    }
    size_t index = 0;
    poles.resize(temp.size());
    for(auto i = temp.begin(); i != temp.end(); i++) {
        poles[index] = *i;
        index++;
    }        
}       
void lp2bs(std::vector<std::complex<DspFloatType>> & zeros,
            std::vector<std::complex<DspFloatType>> & poles,
            DspFloatType wu, DspFloatType wl, DspFloatType & gain)
{ 
    DspFloatType bw = wu-wl;
    DspFloatType Wc = sqrt(wu*wl);
    std::complex<DspFloatType> prodz(1.0,0.0);
    std::complex<DspFloatType> prodp(1.0,0.0);
    for(size_t i = 0; i < zeros.size(); i++)
        prodz *= -zeros[i];
    for(size_t i = 0; i < poles.size(); i++)
        prodp *= -poles[i];
    gain *= prodz.real() / prodp.real();
    std::vector<std::complex<DspFloatType>> ztmp;
    for(size_t i = 0; i < zeros.size(); i++) {
        ztmp.push_back(std::complex<DspFloatType>(0.0,Wc));
        ztmp.push_back(std::complex<DspFloatType>(0.0,-Wc));            
    }
    std::vector<std::complex<DspFloatType>> ptmp;
    for(size_t i = 0; i < poles.size(); i++) {
        if(abs(poles[i])) {
            std::complex<DspFloatType> term1 = DspFloatType(0.5) * bw / poles[i];
            std::complex<DspFloatType> term2 = DspFloatType(0.5) * sqrt((bw*bw) / (poles[i]*poles[i]) - (DspFloatType(4)*Wc*Wc));
            ptmp.push_back(term1+term2);
        }
    }
    size_t index = 0;
    for(auto i = ztmp.begin(); i != ztmp.end(); i++) {
        zeros[index++] = *i;
    }
    index = 0;
    for(auto i = ptmp.begin(); i != ptmp.end(); i++) {
        poles[index++] = *i;
    }
} 


// convert analog setion to biquad type I
BiquadSection AnalogBiquadSection(const BiquadSection &section, DspFloatType fc, DspFloatType fs)
{
    BiquadSection ns = section;
    prewarp(&ns.z[0], &ns.z[1], &ns.z[2], fc, fs);
    prewarp(&ns.p[0], &ns.p[1], &ns.p[2], fc, fs);
    
    //std::vector<DspFloatType> coeffs(4);
    DspFloatType coeffs[4] = {0,0,0,0};
    DspFloatType k =1;
    
    bilinear(ns.z[0], ns.z[1], ns.z[2], ns.p[0], ns.p[1], ns.p[2], &k, fs, coeffs);
    ns.z[0] = k;
    ns.z[1] = k*coeffs[0];
    ns.z[2] = k*coeffs[1];
    ns.p[0] = coeffs[2];
    ns.p[1] = coeffs[3];
    ns.p[2] = 0;    
    return ns;
}
// H(s) => Bilinear/Z => H(z)
// convert analog sos to biquad cascade type I
BiquadSOS AnalogBiquadCascade(const BiquadSOS &sos, DspFloatType fc, DspFloatType fs)
{
    BiquadSOS nsos = sos;
    for (size_t i = 0; i < sos.size(); i++)
    {
        BiquadSection b = AnalogBiquadSection(sos[i], fc, fs);
        nsos[i] = b;
    }
    return nsos;
}

struct Biquad
{
    DspFloatType z[3];
    DspFloatType p[3];
    DspFloatType x[2];
    DspFloatType y[2];

    Biquad() {
        x[0] = x[1] = 0;
        y[0] = y[1] = 0;
    }
    void setCoefficients(DspFloatType _z[3], DspFloatType _p[3]) {
        memcpy(z,_z,3*sizeof(DspFloatType));
        memcpy(p,_p,3*sizeof(DspFloatType));
    }
    DspFloatType Tick(DspFloatType I) {
        DspFloatType out = z[0]*I + z[1]*x[0] + z[2] * z[1];
        out = out - p[0]*y[0] - p[1]*y[1];
        x[1] = x[0];
        x[0] = I;
        y[1] = y[0];
        y[1] = out;
        return out;
    }
};

BiquadSOS butter(int order, double fc, double sr, double Q=1.0)
{
    BiquadSOS sos;    
    size_t n = 1;
    if(order %2 != 0) {
        BiquadSection c;
        std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
        DspFloatType x1 = abs(p1);
        DspFloatType x2 = 0;
                
        // (s-p1)        
        c.z[0] = 1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1.0;;        
        c.p[1] = 1/x1;
        c.p[2] = 0.0;
        sos.push_back(c);        
        n++;
    }
        
    for(size_t i = n; i < order; i += 2)
    {
        std::complex<DspFloatType> p1  = ButterworthPoles(n,order);
        std::complex<DspFloatType> p2  = ButterworthPoles(n+1,order);
        
        DspFloatType x1 = abs(p1*p2);
        DspFloatType x2 = abs(-p1-p2);
        
        // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
        BiquadSection c;
        c.z[0] = 1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1;        
        c.p[1] = Q*x2/x1;
        c.p[2] = 1/x1;    
        sos.push_back(c);
    }
    return AnalogBiquadCascade(sos,fc,sr);
}

BiquadSection butter2(double Q=1.0)
{    
    std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
    std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
  
    DspFloatType x1 = abs(p1*p2);
    DspFloatType x2 = abs(-p1-p2);
    std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
    // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
    BiquadSection c;
    c.z[0] = 1.0/x1;
    c.z[1] = 0.0;
    c.z[2] = 0.0;
    c.p[0] = 1;    
    c.p[1] = Q*x2/x1;
    c.p[2] = 1.0/x1;        

    return c;
}

void test_butterworth()
{
    BiquadSection c = butter2();
    SndFileReaderDouble r("baby_elephant.wav");
    std::vector<DspFloatType> v(r.size()),q(r.size());
    r.read(v.size(),v.data());      
    auto x = AnalogBiquadSection(c,r.samplerate()/2-100,r.samplerate());    
    x.print();
    BiquadTransposedTypeII filter(x);
    
    filter.ProcessBlock(v.size(),v.data(),q.data());
    SndFileWriterDouble w("test.wav",0x10006,r.channels(),r.samplerate());
    w << q;
}


void Besselchubs()
{
    std::complex<DspFloatType> p1  = bessel_poles_2[0];
    std::complex<DspFloatType> p2  = bessel_poles_2[1];
  
    DspFloatType x1 = abs(p1*p2);
    DspFloatType x2 = abs(-p1-p2);
    std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
    // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
    BiquadSection c;
    c.z[0] = 1.0/x1;
    c.z[1] = 0.0;
    c.z[2] = 0.0;
    c.p[0] = 1;
    c.p[1] = (1/200.0)*x2/x1;
    c.p[2] = 1/x1;    
    
    SndFileReaderDouble r("baby_elephant.wav");
    std::vector<DspFloatType> v(r.size()),q(r.size());
    r.read(v.size(),v.data());      
    auto x = AnalogBiquadSection(c,400,r.samplerate());    
    x.print();
    BiquadTransposedTypeII filter(x);
    
    filter.ProcessBlock(v.size(),v.data(),q.data());
    SndFileWriterDouble w("test.wav",0x10006,r.channels(),r.samplerate());
    w << q;
}

void Cheby1chubs()
{
    std::complex<DspFloatType> h0  = ChebyshevH0(2);
    std::complex<DspFloatType> p1  = ChebyshevPole(1,2);
    std::complex<DspFloatType> p2  = ChebyshevPole(2,2);
  
    DspFloatType x1 = abs(p1*p2);
    DspFloatType x2 = abs(-p1-p2);
    std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
    // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
    BiquadSection c;

    c.z[0] = abs(h0)/x1;
    c.z[1] = 0.0;
    c.z[2] = 0.0;
    c.p[0] = 1;
    c.p[1] = (1/200.0)*x2/x1;
    c.p[2] = 1/x1;    
    
    SndFileReaderDouble r("baby_elephant.wav");
    std::vector<DspFloatType> v(r.size()),q(r.size());
    r.read(v.size(),v.data());      
    auto x = AnalogBiquadSection(c,400,r.samplerate());    
    x.print();
    BiquadTransposedTypeII filter(x);
    
    filter.ProcessBlock(v.size(),v.data(),q.data());
    SndFileWriterDouble w("test.wav",0x10006,r.channels(),r.samplerate());
    w << q;
}

int Cheby2Chubs()
{
    std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,2,1.0);
    std::complex<DspFloatType> H1  = Chebyshev2Zeros(2,2,1.0);
    std::complex<DspFloatType> p1  = Chebyshev2Pole(1,2,1);
    std::complex<DspFloatType> p2  = Chebyshev2Pole(2,2,1);
  
    DspFloatType x1 = abs(p1*p2);
    DspFloatType x2 = abs(-p1-p2);
    DspFloatType z1 = abs(H0*H1);
    DspFloatType z2 = abs(-H0-H1);

    std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
    // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
    BiquadSection c;

    c.z[0] = z1/x1;
    c.z[1] = z2/x1;
    c.z[2] = 0;
    c.p[0] = 1;
    // radius is the same thing but goes from 0..1
    // 0 = most resonant
    // 1 = least resonant
    c.p[1] = (1/200.0)*x2/x1;
    c.p[2] = 1/x1;    
    
    SndFileReaderDouble r("baby_elephant.wav");
    std::vector<DspFloatType> v(r.size()),q(r.size());
    r.read(v.size(),v.data());      
    auto x = AnalogBiquadSection(c,400,r.samplerate());    
    x.print();
    BiquadTransposedTypeII filter(x);
    
    filter.ProcessBlock(v.size(),v.data(),q.data());
    SndFileWriterDouble w("test.wav",0x10006,r.channels(),r.samplerate());
    w << q;
}

// (s-p1)(s-p2) = s^2 -(p1+p2)s + p1*p2

// this doesn't work correctly use CauerFilter.hpp
// but elliptical filter isn't very good for audio anyways it's too difficult to calculate the cutoff
// it has the steepest cutoff but it has too much ripples to be useful for anything else
void elliptical_order_estimate(DspFloatType omegaPass, DspFloatType omegaStop, DspFloatType maxPassLoss, DspFloatType minStopLoss,
                                int & order, DspFloatType& actualMinStopLoss)
{
    DspFloatType k,u,q,dd,kk,lambda,w,mu,om;
    DspFloatType sum,term,denom,numer,sigma,v;
    int i,m,r;

    k = omegaPass/omegaStop;
    kk = sqrt(sqrt(1.0-k*k));
    u = 0.5 * (1.0-kk)/(1.0+kk);
    q = 150.0 * pow(u,13.0);
    q = q + 15.0*pow(u,9.0);
    q = q + 2.0*pow(u,5.0);
    q = q + u;
    dd = pow(10.0,minStopLoss/10.0)-1.0;
    dd = dd/(pow(10.0,maxPassLoss/10.0)-1.0);
    order = ceil(log10(16.0*dd)/log10(1.0/q));
    numer = pow(10.0,(maxPassLoss/10.0))-1.0;
    actualMinStopLoss = 10.0*log10(numer/(16.0*pow(q,order))+1.0);    
}                                


void elliptical_coeffs(DspFloatType omegaPass, DspFloatType omegaStop, DspFloatType maxPassLoss, int order,
                    std::vector<DspFloatType> &aa,
                    std::vector<DspFloatType> &bb,
                    std::vector<DspFloatType> &cc,
                    int& numSecs,
                    DspFloatType& hZero,
                    DspFloatType& pZero)
{
    DspFloatType k,kk,u,q,vv,ww,mu,xx,yy,sum,term,denom,numer;
    int i,m,r;

    k = omegaPass/omegaStop;
    kk= sqrt(sqrt(1.0 - k*k));
    u = 0.5*(1.0-kk)/(1.0+kk);

    q = 150.0 * pow(u,13.0);
    q = q + 15.0* pow(u,9.0);
    q = q + 2.0 * pow(u,5.0);
    q = q + u;

    numer = pow(10.0,maxPassLoss/20.0) + 1.0;
    vv = log(numer / (pow(10.0,maxPassLoss/20.0)-1.0))/(2.0*order);


    sum = 0.0;
    for(m = 0; m < 5; m++) {
        term = pow(-1.0,m);
        term = term * pow(q,m*(m+1));
        term = term * sinh((2*m+1)*vv);
        sum  = sum + term;
    }
    numer = 2.0 * sum * sqrt(sqrt(q));
    sum = 0.0;
    for(m = 1; m < 5; m++) {
        term = pow(-1.0,m);
        term = term * pow(q,m*m);
        term = term * cosh(2.0*m*vv);
        sum = sum + term;
    }
    denom = 1.0 + 2.0*sum;
    pZero = fabs(numer/denom);
    
    ww = 1.0 + k * pZero*pZero;
    ww = sqrt(ww * (1.0 + pZero*pZero/k));

    r = (order - (order %2))/2;
    numSecs = r;

    aa.resize(r+1);
    bb.resize(r+1);
    cc.resize(r+1);

    for(i=1; i <= r; i++) {
        if(order % 2) 
            mu = i;
        else
            mu = i-0.5;

        sum = 0.0;
        for(m=0; m < 5; m++) {
            term = pow(-1.0,m);
            term = term * pow(q,m*(m+1));
            term = term * sin((2*m+1)*M_PI*mu/order);
            sum  += term;
        }
        numer = 2.0 * sum * sqrt(sqrt(q));

        sum = 0.0;
        for(m=1; m < 5; m++) {
            term = pow(-1.0,m);
            term = term * pow(q,m*m);
            term = term * cos(2.0*M_PI*m*mu/order);
            sum += term;
        }
        denom = 1.0 + 2.0*sum;
        xx = numer/denom;
        yy = 1.0 - k*xx*xx;
        yy = sqrt(yy * (1.0 - (xx*xx/k)));
        aa[i] = 1.0/(xx*xx);
        denom = 1.0 + pow(pZero*xx,2.0);
        bb[i] = 2.0 * pZero * yy/denom;
        denom = pow(denom,2.0);
        numer = pow(pZero*yy,2.0) + pow(xx*ww,2.0);
        cc[i] = numer/denom;
    }
    term = 1.0;
    for(i=1; i <= r; i++) {            
        term = term * cc[i]/aa[i];
    }
    if(order % 2) {
        term = term * pZero;
    }
    else {
        term = term * pow(10.0,maxPassLoss/(-20.0));
    }
    hZero = term;
}

void cauerRescale(int order, std::vector<DspFloatType> & aa,
                            std::vector<DspFloatType> & bb,
                            std::vector<DspFloatType> & cc,
                            DspFloatType *hZero,
                            DspFloatType *pZero,
                            DspFloatType alpha)
{
    DspFloatType alphaSqrd;
    int r,i;

    alphaSqrd = alpha*alpha;
    if(order % 2 ) {
        r = (order-1)/2;
        *hZero = *hZero * alpha;
       *pZero = *pZero * alpha;
    }
    else {
        r = order/2;
    }
    
    for(i=1; i <= r; i++) {
        aa[i] = aa[i] * alphaSqrd;
        cc[i] = cc[i] * alphaSqrd;
        bb[i] = bb[i] * alpha;
    }
}                 

#include "CauerFilter.hpp"

int main()
{
    Filters::Elliptical::CauerFilter filter(0.1,51.556,30,5200);
    filter.calculate_coefficients();
    std::cout << filter.getHo() << "\n";
    for(size_t i = 0; i < filter.Ai.size(); i++) std::cout << filter.Ai[i] << "\n";
    for(size_t i = 0; i < filter.Bi.size(); i++) std::cout << filter.Bi[i] << "\n";
    for(size_t i = 0; i < filter.Ci.size(); i++) std::cout << filter.Ci[i] << "\n";

    
    double hZero = filter.getHo();
    BiquadSOS sos;
    for(size_t i = 0; i < filter.r; i++)
    {
        BiquadSection c;
        c.z[0] = hZero*filter.Ai[i];
        c.z[1] = 0;
        c.z[2] = hZero;    
        c.p[0] = filter.Ci[i];
        c.p[1] = filter.Bi[i];
        c.p[2] = 1;       
        c.z[0] /= c.p[0];    
        c.z[2] /= c.p[0];    
        c.p[1] /= c.p[0];
        c.p[2] /= c.p[0];
        c.p[0] /= c.p[0];
        sos.push_back(c);
    }
    SndFileReaderDouble r("baby_elephant.wav");
    std::vector<DspFloatType> v(r.size()),q(r.size());
    r.read(v.size(),v.data());      
    auto x = AnalogBiquadCascade(sos,8000,r.samplerate());    
    
    FilterBase filt(x);
    
    filt.ProcessBlock(v.size(),v.data(),q.data());
    SndFileWriterDouble w("test.wav",0x10006,r.channels(),r.samplerate());
    w << q;
}

void tests()
{
    int order,numSections;
    std::vector<DspFloatType> a,b,c;
    DspFloatType pZero,hZero;
    DspFloatType omega_p = 3000.0;
    DspFloatType omega_s = 3200.0;

    
    elliptical_coeffs(omega_p, omega_s, 0.1, 9, a,b,c,numSections,hZero,pZero);
    std::cout << numSections << std::endl;
    
    for(size_t i = 1; i < 5; i++)
        std::cout << a[i] << "," << b[i] << "," << c[i] << std::endl;

    DspFloatType h0;
    if(order % 2 != 0) {        
        
        // d = s + p0
        // h0/(s + p0)
        
        // b0=h0/p0
        // b1=0
        // b2=0
        // a1=1
        // a2=/p0
        // a3=0
        
    }
    else 
    {
        h0 = hZero;
        // b0 = h0
        // b1 = 0;
        // b2 = 0;        
    }
    
    
    std::complex<DspFloatType> p1  = ButterworthPoles(1,6);
    std::complex<DspFloatType> p2  = ButterworthPoles(2,6);
    std::complex<DspFloatType> p3  = ButterworthPoles(3,6);
    std::complex<DspFloatType> p4  = ButterworthPoles(4,6);
    std::complex<DspFloatType> p5  = ButterworthPoles(5,6);
    std::complex<DspFloatType> p6  = ButterworthPoles(6,6);
    
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << p3 << std::endl;
    std::cout << p4 << std::endl;
    std::cout << p5 << std::endl;
    std::cout << p6 << std::endl;
    
    std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,2,1.0);
    std::complex<DspFloatType> H1  = Chebyshev2Zeros(2,2,1.0);
    p1  = Chebyshev2Pole(1,2,1);
    p2  = Chebyshev2Pole(2,2,1);

    std::cout << H0 << std::endl;
    std::cout << H1 << std::endl;
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    
    std::vector<DspFloatType> bessel;
    bessel_coefficients(3,bessel);
    for(size_t i = 0; i < bessel.size(); i++)
        std::cout << bessel[i] << ",";
    std::cout << std::endl;

    bessel = qn(7);
    for(size_t i = 0; i < bessel.size(); i++)
        std::cout << bessel[i] << ",";

}