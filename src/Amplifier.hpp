#pragma once

#include "Amplifiers.hpp"

// Amplifiers.hpp
// ClipFunctions.hpp
// DistortionFunctions.hpp
// WaveShapers.hpp
struct Amplifier
{
    // ToneStack
    double G;
    double bias;
    double prv;
    double min = -1;
    double max = 1;
    double postgain = 1;
    double pregain = 1;
    

    enum {
        TANH,
        TANH_NORMAL,
        POSITIVE,
        NEGATIVE,
        SIGMOID,
        SIGMOID_MINUS,
        BPSIGMOID,
        FULLRECTIFY,
        HALFRECTIFY,
        ASIGMOID,
        ASIGMOID2,
        DISTORTION,
        CUBIC,
        ASIN,
        ACOS,
        ATAN,
        ASINH,
        ACOSH,
        ATANH,
        EXP,
        DC,
        BIPOLAR,
        QUADRATIC,
        QUADRATIC2,
        QUADRATIC3,
        PARAMETRIC,
        ARCTAN,
        SOFT,
        ERF,
        SIGMOID,
        HARDCLIP,
        HYPERTAN,
        DIODE,
        FUZZ,
        PIECEWISE,
        TUBE,
        ARRAYA,
        GALLO,
        DOUBLESOFT,
        CRUSH,
        TUBOID,
        YEH,        
    };
    int distortion = TANH;

    Amplifier(double Gain = 1.0f, double b = 0.0f)
    {
        G = Gain;
        bias = b;        
    }
    
    
    double Integrator(double in) {
        double r = in + prv;
        prv = in;
        return r;
    }
    double Differencer(double in) {
        double r = in - prv;
        prv = in;
        return r;
    }
    void SetBias(double b) {
        bias = b;
    }

    double Tick(double I, double A = 1, double X = 1, double Y = 1, double B=0)
    {
        double r = 0;
        double K = G;
        extern float ::pregain = pregain;
        extern float ::postgain= postgain;
        
        switch(distortion)   
        {
            case TANH: r = tanh(I,K); break;
            case TANH_NORMAL: r = tanh_normal(I,K*A); break;
        }
        r = clamp(A*r,min*X,max*Y);
        return r;
    }
};