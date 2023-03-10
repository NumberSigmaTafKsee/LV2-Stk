//
//  ADRS.h
//
//  Created by Nigel Redmon on 12/18/12.
//  EarLevel Engineering: earlevel.com
//  Copyright 2012 Nigel Redmon
//
//  For a complete explanation of the ADSR envelope generator and code,
//  read the series of articles by the author, starting here:
//  http://www.earlevel.com/main/2013/06/01/envelope-generators/
//
//  License:
//
//  This source code is provided as is, without warranty.
//  You may copy and distribute verbatim copies of this document.
//  You may modify and use this source code to create binary code for your own purposes, free or commercial.
//

#pragma once

// experiment with filter
// not needed

#include "OnePole.hpp"
namespace SoundWave {
class ADSR {
public:
	
    ADSR(float sample_rate=44100.0f, float cutoff=1000.0f)
    {
        reset();
        setAttackRate(0);
        setDecayRate(0);
        setReleaseRate(0);
        setSustainLevel(1.0);
        setTargetRatioA(0.3);
        setTargetRatioDR(0.0001);
        filter.setFc(cutoff);
        sr = sample_rate;
    }

    ~ADSR(void) {
    }


    ADSR(float a, float d, float s, float r, float sample_rate=44100.0f, float cutoff=1000.0f) {
        filter.setFc(cutoff);
        sr = sample_rate;
        reset();
        setAttackRate(a * sample_rate);
        setDecayRate(d * sample_rate);
        setSustainLevel(s);
        setReleaseRate(r * sample_rate);                
        setTargetRatioA(0.3);
        setTargetRatioDR(0.0001);
    }
	
	float process(void);
    float getOutput(void);
    int getState(void);
	void gate(int on);
    
    void setAllTimes(float a, float d, float s, float r) {
        setAttackTime(a);
        setDecayTime(d);
        setSustainLevel(s);
        setReleaseTime(r);
    }

    void setAttackTime(float rate) { attackRate=rate/sr;}
    void setDecayTime(float rate) { decayRate=rate/sr; }
    void setReleaseTime(float rate) { releaseRate=rate/sr;}

    void setAttackRate(float rate);
    void setDecayRate(float rate);
    void setReleaseRate(float rate) ;
	void setSustainLevel(float level);

    void setTargetRatioA(float targetRatio);
    void setTargetRatioDR(float targetRatio);

    void reset(void);

    float Tick(float in) {
        return process() * in;
    }
    float Tick() {
        return process();
    }
    void Process(size_t n, float * input, float * output) {
        for(size_t i = 0; i < n; i++)
        {
            output[i] = process()*input[i];
        }
    }
    void Process(float * samples, size_t n) {
        for(size_t i = 0; i < n ; i++)
            samples[i] = process();
    }

    enum envState {
        env_idle = 0,
        env_attack,
        env_decay,
        env_sustain,
        env_release
    };

protected:
    OnePole filter;
	int state;
	float output;
	float attackRate;
	float decayRate;
	float releaseRate;
	float attackCoef;
	float decayCoef;
	float releaseCoef;
	float sustainLevel;
    float targetRatioA;
    float targetRatioDR;
    float attackBase;
    float decayBase;
    float releaseBase;
    float sr;

    float calcCoef(float rate, float targetRatio);
};

inline float ADSR::process() {
	switch (state) {
        case env_idle:
            break;
        case env_attack:
            output = attackBase + output * attackCoef;
            if (output >= 1.0) {
                output = 1.0;
                state = env_decay;
            }
            break;
        case env_decay:
            output = decayBase + output * decayCoef;
            if (output <= sustainLevel) {
                output = sustainLevel;
                state = env_sustain;
            }
            break;
        case env_sustain:
            break;
        case env_release:
            output = releaseBase + output * releaseCoef;
            if (output <= 0.0) {
                output = 0.0;
                state = env_idle;
            }
	}
	return filter.process(output);
}

inline void ADSR::gate(int gate) {
	if (gate)
		state = env_attack;
	else if (state != env_idle)
        state = env_release;
}

inline int ADSR::getState() {
    return state;
}

inline void ADSR::reset() {
    state = env_idle;
    output = 0.0;
}

inline float ADSR::getOutput() {
	return output;
}

inline void ADSR::setAttackRate(float rate) {
    attackRate = rate;
    attackCoef = calcCoef(rate, targetRatioA);
    attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
}

inline void ADSR::setDecayRate(float rate) {
    decayRate = rate;
    decayCoef = calcCoef(rate, targetRatioDR);
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
}

inline void ADSR::setReleaseRate(float rate) {
    releaseRate = rate;
    releaseCoef = calcCoef(rate, targetRatioDR);
    releaseBase = -targetRatioDR * (1.0 - releaseCoef);
}

inline float ADSR::calcCoef(float rate, float targetRatio) {
    return (rate <= 0) ? 0.0 : exp(-log((1.0 + targetRatio) / targetRatio) / rate);
}

inline void ADSR::setSustainLevel(float level) {
    sustainLevel = level;
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
}

inline void ADSR::setTargetRatioA(float targetRatio) {
    if (targetRatio < 0.000000001)
        targetRatio = 0.000000001;  // -180 dB
    targetRatioA = targetRatio;
    attackCoef = calcCoef(attackRate, targetRatioA);
    attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
}

inline void ADSR::setTargetRatioDR(float targetRatio) {
    if (targetRatio < 0.000000001)
        targetRatio = 0.000000001;  // -180 dB
    targetRatioDR = targetRatio;
    decayCoef = calcCoef(decayRate, targetRatioDR);
    releaseCoef = calcCoef(releaseRate, targetRatioDR);
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    releaseBase = -targetRatioDR * (1.0 - releaseCoef);
}
}