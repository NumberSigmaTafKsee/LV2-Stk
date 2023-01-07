#pragma once

struct Stereofyier
{
    float pan = 0.5;
    SinewaveGenerator lfo;
    float WideCoeff = 1;
    float rotAmt = 0.5;
    float lfoFreq = 0.01;

    Stereofyier() : lfo(lfoFreq,sampleRate)
    {
        
    }
    void setFrequency(float f) {
        lfo.setFrequency(f);
    }
    void setPan(float p) {
        pan = p;
    }
    float panLeft(float x, float lfoTick) {

        return x + rotAmt *x * lfoTick * sin((1-pan)*M_PI/2);        
    }
    float panRight(float x, float lfoTick) {
        return x + rotAmt * x * lfoTick * cos(pan*M_PI/2);
    }
    void ProcessBuffer(size_t n, float ** inputs, float ** outputs) 
    {            
        for(size_t i = 0; i < n; i++)
        {
            float inL,inR;
            float lfoTick = lfo.Tick();
            inL = inputs[0][i];
            inR = inputs[1][i];
        
            float MonoSign   = (inL + inR)/2.0;
            float DeltaLeft  = inL - MonoSign;
            DeltaLeft *= WideCoeff;
            float DeltaRight = inR - MonoSign;
            DeltaRight *= WideCoeff;

            inL += DeltaLeft;
            inR += DeltaRight;

            // calculate scale coefficient
            float coef_S = WideCoeff*0.5;

            // then do this per sample
            float m = MonoSign;
            float s = (inL - inR )*coef_S;

            float out_left  = m - s;
            float out_right = m + s;

            outputs[0][i] = panRight(out_left,lfoTick);
            outputs[1][i] = panLeft(out_right,lfoTick);  
        }
    }

    void InplaceProcess(size_t numSamples, float ** buffer)
    {
        ProcessBuffer(numSamples,buffer,buffer);
    }
};

