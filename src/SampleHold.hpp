#pragma once

struct SampleHold
{
    double rate,sr;
    double phase;
    double latch;

    SampleHold(float F, float Fs=44100.0f)
    {
        sr   = Fs;
        rate = F / Fs;
        phase = 0;
        latch = 0;
    }
    void setRate(float frequency) {
        rate = frequency/sr;
    }
    float Tick(float in)
    {
        if(phase >= 1.0) latch = in;
        phase = (phase + rate,1);
        return latch;
    }
};
