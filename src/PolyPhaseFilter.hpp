#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef enum {
    FIR_LOWPASS,
    FIR_HIGHPASS,
    FIR_BANDPASS,
    FIR_BANDSTOP
} FIR_Response;

typedef struct {
    float* pfb;
    int num_phases;
    int taps_per_phase;
    int interpolation;
    int decimation;
    int phase_index;
    int input_deficit;
    float* history;
    int history_length;
} FIR_Filter;

float* src_generate_fir_coeffs(int num_taps, float cutoff);
FIR_Filter* src_generate_fir_filter(float* coefficients, int num_taps,
        int interpolation, int decimation);
FIR_Filter* src_generate(int interpolation, int decimation);
int src_filt(FIR_Filter* filter, float* samples, int count, float* p_output);

struct PolyPhaseFilter : public FilterProcessor
{
    std::vector<float> taps;
    FIR_Filter * filter;
    PolyPhaseFilter(int num_taps, float cutoff, int interpolation, int decimation) : FilterProcessor()
    {
        taps.resize)num_taps;
        float * coeffs = src_generate_fir_coeffs(num_taps,cutoff);
        memcpy(taps.size(),coeffs,num_taps * sizeof(float));
        filter = src_generate_fir_filter(coeffs,num_taps,interpolation,decimation);        
    }
    PolyPhaseFilter(int interpolation, int decimation) : FilterProcessor()
    {        
        filter = src_generate(interpolation,decimation);        
    }
    ~PolyPhaseFilter() {
        if(filter) {
            free(filter->history);
            free(filter);
        }
    }
    void ProcessBlock(size_t n, float * input, float * output)
    {        
        src_filt(filter,input,n,output);
    }
};


#define sinc(x)  (x == 0.0) ? (1.0) : (sin(M_PI * x)/(M_PI * x))
#define sincf(x) (x == 0.0) ? (1.0) : (sinf(M_PI * x)/(M_PI * x))

// n/d rounding up
int iceil(int n, int d) {
    int rounding = (n % d != 0) ? 1 : 0;
    return n / d + rounding;
}

// Dot product of two vectors
float dot(float* a, int a_length, float* history, float* b, int b_last_index) {
    float dotprod = 0.0f;
    int i = 0;
    for (; i < (a_length - b_last_index - 1); i++) {
        dotprod += a[i] * history[b_last_index + i];
    }
    for (; i < a_length; i++) {
        dotprod += a[i] * b[b_last_index - a_length + 1 + i];
    }
    return dotprod;
}

// Shift b into a
void src_shiftin(float* a, int a_length, float* b, int b_length) {
    if (b_length > a_length) 
        memcpy(a, &b[b_length - a_length], a_length * sizeof(float));
    else {
        for (int i = 0; i < (a_length - b_length); i++) {
            a[i] = a[i + b_length];
        }
        for (int i = 0; i < b_length; i++) {
            a[i + a_length - b_length] = b[i];
        }
    }
}

// Hamming Window
float* src_hamming(int num_taps) {
    float* window = malloc(num_taps * sizeof(float));
    for (int i = 0; i < num_taps; i++) {
        float alpha = 0.54;
        float beta = 0.46;
        window[i] = alpha - 
                beta * cos(2.0f * M_PI * (float)i / (float)(num_taps - 1));
    }
    return window;
}

// Hann Window
float* src_hann(int num_taps) {
    float* window = malloc(num_taps * sizeof(float));
    for (int i = 0; i < num_taps; i++) {
        float alpha = 0.5;
        float beta = 0.5;
        window[i] = alpha - 
                beta * cos(2.0f * M_PI * (float)i / (float)(num_taps - 1));
    }
    return window;
}

// Blackman Window
float* src_blackman(int num_taps) {
    float* window = malloc(num_taps * sizeof(float));
    for (int i = 0; i < num_taps; i++) {
        float alpha = 0.42;
        float beta = 0.5;
        window[i] = alpha - 
                beta * cos(2.0f * M_PI * (float)i / (float)(num_taps - 1)) +
                (beta - alpha) * 
                cos(4.0f * M_PI * (float)i / (float)(num_taps - 1));
    }
    return window;
}

// Generate an FIR prototype
// For LPF and HPF, feed cutoff_low with the cutoff frequency
float* src_fir_prototype(int num_taps, float cutoff_low, float cutoff_high, 
        FIR_Response response) {
    float* proto = malloc(num_taps * sizeof(float));
    // some alias to help make the code cleaner
    float f = cutoff_low;
    float f1 = cutoff_low;
    float f2 = cutoff_high;
    int m = num_taps - 1;
    switch(response) {
        case FIR_LOWPASS:
            for (int i = 0; i < num_taps; i++) 
                proto[i] = 2.f*f*sincf(2.f*f*(i-m/2.f));
            break;
        case FIR_HIGHPASS:
            for (int i = 0; i < num_taps; i++)
                proto[i] = sincf(i-m/2.f)-2.f*f*sincf(2.f*f*(i-m/2.f));
            break;
        case FIR_BANDPASS:
            for (int i = 0; i < num_taps; i++)
                proto[i] = 2.f*(f1*sincf(2.f*f1*(i-m/2.f)) - 
                        f2*sincf(2.f*f2*(i-m/2.f)));
            break;
        case FIR_BANDSTOP:
            for (int i = 0; i < num_taps; i++)
                proto[i] = 2.f*(f2*sincf(2.f*f2*(i-m/2.f)) - 
                        f1*sincf(2*f1*(i-m/2.f)));
            break;
        default:
            free(proto);
            proto = NULL;
            fprintf(stderr, "Invalid response type!\n");
            break;
    }
    return proto;
}

// Design an FIR filter with windowing
float* src_generate_fir_coeffs(int num_taps, float cutoff) {
    float* proto = src_fir_prototype(num_taps, cutoff, 0, FIR_LOWPASS);
    float* window = src_hann(num_taps);
    for (int i = 0; i < num_taps; i++) {
        proto[i] *= window[i];
    }
    free(window);

    return proto;
}

// Generate Polyphase filter bank for interpolator and rational
void src_taps_to_pfb(float* coefficients, int num_taps, int interpolation, 
        float** p_pfb, int* p_num_phases,
        int* p_taps_per_phase) {
    int num_phases = interpolation;
    int taps_per_phase = iceil(num_taps, num_phases); // iceil(a/b)
    int pfb_size = taps_per_phase * num_phases;
    float* pfb = malloc(pfb_size * sizeof(float));
    int c_index = 0;

    for (int phase = 0; phase < num_phases; phase++)
        for (int tap = 0; tap < taps_per_phase; tap++) {
            pfb[phase * taps_per_phase + taps_per_phase - 1 - tap] = 
                    coefficients[tap * num_phases + phase];
        }

    *p_pfb = pfb;
    *p_num_phases = num_phases;
    *p_taps_per_phase = taps_per_phase;
}

// Generate a FIR filter with given coefficients
// It is will make a local copy of coefficients
FIR_Filter* src_generate_fir_filter(float* coefficients, int num_taps,
        int interpolation, int decimation) {
    FIR_Filter* filter = malloc(sizeof(FIR_Filter));

    src_taps_to_pfb(coefficients, num_taps, interpolation, 
            &filter->pfb,
            &filter->num_phases,
            &filter->taps_per_phase);
    filter->interpolation = interpolation;
    filter->decimation = decimation;
    filter->phase_index = 0;
    filter->input_deficit = 0;

    filter->history_length = filter->taps_per_phase - 1;
    filter->history = malloc(filter->history_length * sizeof(float));
    memset(filter->history, 0, filter->history_length * sizeof(float));

    return filter;
}

// Generate a filter
FIR_Filter* src_generate(int interpolation, int decimation) {
    int num_taps = 24;
    float cutoff_freq = 0.5 / num_taps;
    float* coefficients = src_generate_fir_coeffs(num_taps, cutoff_freq);
    FIR_Filter* filter = src_generate_fir_filter(coefficients, num_taps, 
            interpolation, decimation);
    free(coefficients);
    return filter;
}

// Filt incoming sample block
int src_filt(FIR_Filter* filter, float* samples, int count, 
        float* output) {
    int phase_index_step = filter->decimation % filter->interpolation;
    int phase = filter->phase_index;
    int i = filter->input_deficit;
    int input_step;
    int j = 0;
    
    while (i < count) {
        output[j++] = dot(
                &filter->pfb[phase * filter->taps_per_phase],
                filter->taps_per_phase, filter->history, samples, i);
        i += (phase + filter->decimation) / filter->interpolation;
        phase = (phase + phase_index_step) % filter->interpolation;
    }

    filter->input_deficit = i - count;
    filter->phase_index = phase;
    src_shiftin(filter->history, filter->history_length, samples, count);

    return j;
}

