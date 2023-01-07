#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <list>
#include <pthread.h>
#include <unistd.h>
#include <glob.h>


#include "AudioMidi/audiosystem.h"
#include "Undenormal.hpp"
#include "Std/StdObject.h"
#include "Std/StdRandom.h"
#include "samples/sample.hpp"
#include "samples/sample_dsp.hpp"
#include "SndFile.hpp"

#include "DSP/Decimators.hpp"
#include "DSP/SimpleResampler.hpp"


#include "Envelopes/Modulation.hpp"

#include "Oscillators/FunctionGenerator.hpp"
#include "Oscillators/FunctionLFO.hpp"

#include "Envelopes/ADSR.hpp"

#include "Filters/Filters.h"
#include "Filters/Filters.hpp"
#include "Filters/Biquads.hpp"
//#include "Filters/BiquadFilters.hpp"
//#include "Filters/MoogFilters.hpp"
//#include "Filters/StateVariableFilters.hpp"


#include "Distortion/Amplifier.hpp"
#include "Distortion/Amplifiers.hpp"
#include "Distortion/ClipFunctions.hpp"
#include "Distortion/DistortionFunctions.hpp"
#include "Distortion/Diode.hpp"
#include "Distortion/SimpleDynamics.hpp"


#include "Delay/BasicDelay.hpp"
#include "Delay/DelayLines.hpp"
#include "Delay/Delays.hpp"


/*
#include "Delay/AllPass.hpp"
#include "Delay/BasicDelayLine.hpp"
#include "Delay/BasicDelayLineStereo.hpp"
#include "Delay/CircularDelay.hpp"
#include "Delay/Comb.hpp"
#include "Delay/CombFilter.hpp"
#include "Delay/CrossDelayLine.hpp"
#include "Delay/DCBlock.h"
#include "Delay/DelayLine.hpp"
#include "Delay/ERTapDelayLine.hpp"
#include "Delay/FractionalDelayBuffer.hpp"
#include "Delay/LowPass.hpp"
#include "Delay/Moorer.hpp"
#include "Delay/MoorerStereo.hpp"
#include "Delay/Multiverb.hpp"
#include "Delay/PPDelayLine.hpp"
#include "Delay/reverb.hpp"
#include "Delay/Schroeder.hpp"
#include "Delay/SchroederImproved.hpp"
#include "Delay/SchroederStereo.hpp"
#include "Delay/SyncedTapDelayLine.hpp"
*/

#include "Plugins/LV2Plugin.hpp"

#include "FX/AudioEffects.hpp"
#include "FX/MVerb.h"
#include "FX/AudioDSP_Delay.hpp"
#include "FX/AudioDSP_VibraFlange.hpp"
#include "FX/AudioDSP_Chorus.hpp"
#include "FX/AudioDSP_Phaser.hpp"
#include "FX/FXChorusDelay.hpp"
#include "FX/FXYKChorus.hpp"
#include "FX/FXCE2Chorus.hpp"

//#include "SstFilters.hpp"
//#include "SstWaveshaper.hpp"

#include "MusicFunctions.hpp"
#include "Threads.hpp"



using namespace std;

Std::RandomMersenne noise;
int    BufferSize=256;
double sampleRate = 44100.0f;
double sampleRate2 = 44100.0f/2.0f;
double inverseSampleRate = 1.0/sampleRate;

#include "Oscillators/FunctionGenerator.hpp"
#include "Oscillators/WaveGenerators.hpp"
#include "FX/Stereofiyer.hpp"

vector<std::vector<uint8_t>> patterns;
vector<uint8_t>  global1;
vector<uint8_t>  patch;

// true if want to dump sysex to edit buffer
// not all parameters will respond in the keyboard like LFOs
bool send_sysex = true;
bool auto_send  = true;

// smooth sysex
// morph from one patch to another over time 
// filter parameter changes 
// filter cc/nrpn 
// CSmoothFilter
// BiquadFilter 

uint8_t* PackBits(uint8_t * data)
{
    uint8_t * out = (uint8_t*)malloc(1024);
    uint8_t * p = out;
    uint8_t a7,b7,c7,e7,f7,g7;
    uint8_t d1,d2,d3,d4,d5,d6,d7,d8;

    for(size_t i=0; i < 252; i+=7)
    {	

		a7 = (data[i+0] & 0x80) >> 7;
		b7 = (data[i+1] & 0x80) >> 6;
		c7 = (data[i+2] & 0x80) >> 5;
		d7 = (data[i+3] & 0x80) >> 4;
        e7 = (data[i+4] & 0x80) >> 3;
		f7 = (data[i+5] & 0x80) >> 2;
		g7 = (data[i+6] & 0x80) >> 1;

		d1 = a7 | b7 | c7 | d7 | e7 | f7 | g7;
		d2 = data[i+0] & 0x7f;
		d3 = data[i+1] & 0x7f;
		d4 = data[i+2] & 0x7f;
		d5 = data[i+3] & 0x7f;
		d6 = data[i+4] & 0x7f;
		d7 = data[i+5] & 0x7f;
		d8 = data[i+6] & 0x7f;

        *p++ = d1;
        *p++ = d2;
        *p++ = d3;
        *p++ = d4;
        *p++ = d5;
        *p++ = d6;
        *p++ = d7;        
        *p++ = d8;
    }
	int i = 252;
	a7 = (data[i+0] & 0x80) >> 7;
	b7 = (data[i+1] & 0x80) >> 6;
	c7 = (data[i+2] & 0x80) >> 5;
	d7 = (data[i+3] & 0x80) >> 4;

	d1 = a7 | b7 | c7 | d7;
	d2 = data[i+0] & 0x7f;
	d3 = data[i+1] & 0x7f;
	d4 = data[i+2] & 0x7f;
	d5 = data[i+3] & 0x7f;

    *p++ = d1;
    *p++ = d2;
    *p++ = d3;
    *p++ = d4;
    *p++ = d5;
	
	return out;
}

uint8_t* PackTetraBits(uint8_t * data)
{
    uint8_t * out = (uint8_t*)malloc(1024);
    uint8_t * p = out;
    uint8_t a7,b7,c7,e7,f7,g7;
    uint8_t d1,d2,d3,d4,d5,d6,d7,d8;

    for(size_t i=0; i < 378; i+=7)
    {	

		a7 = (data[i+0] & 0x80) >> 7;
		b7 = (data[i+1] & 0x80) >> 6;
		c7 = (data[i+2] & 0x80) >> 5;
		d7 = (data[i+3] & 0x80) >> 4;
        e7 = (data[i+4] & 0x80) >> 3;
		f7 = (data[i+5] & 0x80) >> 2;
		g7 = (data[i+6] & 0x80) >> 1;

		d1 = a7 | b7 | c7 | d7 | e7 | f7 | g7;
		d2 = data[i+0] & 0x7f;
		d3 = data[i+1] & 0x7f;
		d4 = data[i+2] & 0x7f;
		d5 = data[i+3] & 0x7f;
		d6 = data[i+4] & 0x7f;
		d7 = data[i+5] & 0x7f;
		d8 = data[i+6] & 0x7f;

        *p++ = d1;
        *p++ = d2;
        *p++ = d3;
        *p++ = d4;
        *p++ = d5;
        *p++ = d6;
        *p++ = d7;        
        *p++ = d8;
    }
	int i = 378;
	a7 = (data[i+0] & 0x80) >> 7;
	b7 = (data[i+1] & 0x80) >> 6;
	c7 = (data[i+2] & 0x80) >> 5;
	d7 = (data[i+3] & 0x80) >> 4;
    e7 = (data[i+2] & 0x80) >> 3;
	f7 = (data[i+3] & 0x80) >> 2;

	d1 = a7 | b7 | c7 | d7;
	d2 = data[i+0] & 0x7f;
	d3 = data[i+1] & 0x7f;
	d4 = data[i+2] & 0x7f;
	d5 = data[i+3] & 0x7f;
    d6 = data[i+4] & 0x7f;
    d7 = data[i+5] & 0x7f;

    *p++ = d1;
    *p++ = d2;
    *p++ = d3;
    *p++ = d4;
    *p++ = d5;
    *p++ = d6;
    *p++ = d7;
	
	return out;
}

uint8_t* UnpackBits(uint8_t * data)
{	
    uint8_t *buf = (uint8_t*)malloc(1024);
    uint8_t a7,b7,c7,e7,f7,g7;
    uint8_t d1,d2,d3,d4,d5,d6,d7,d8;
	uint8_t o1,o2,o3,o4,o5,o6,o7,o8;
    uint8_t *p = buf;

    for(size_t i = 0; i < 288; i+= 8)
    {
	
		g7 = (data[i+0] & 0x40) << 1;
		f7 = (data[i+0]  & 0x20) << 2;
		e7 = (data[i+0] & 0x10) << 3;
		d7 = (data[i+0] & 0x8) << 4;
		c7 = (data[i+0] & 0x4) << 5;
		b7 = (data[i+0] & 0x2) << 6;
		a7 = (data[i+0] & 0x01) << 7;

		o1 = data[i+1] | a7;
		o2 = data[i+2] | b7;
		o3 = data[i+3] | c7;
		o4 = data[i+4] | d7;
		o5 = data[i+5] | e7;
		o6 = data[i+6] | f7;
		o7 = data[i+7] | g7;

        *p++ = o1;        
        *p++ = o2;        
        *p++ = o3;        
        *p++ = o4;        
        *p++ = o5;        
        *p++ = o6;        
        *p++ = o7;
    }

	size_t i = 288;
	e7 = (data[i+0] & 0x10) << 1;
	d7 = (data[i+0] & 0x8) << 2;
	c7 = (data[i+0] & 0x4) << 3;
	b7 = (data[i+0] & 0x2) << 4;
	a7 = (data[i+0] & 0x01) << 5;

	o1 = data[i+1] | a7;
	o2 = data[i+2] | b7;
	o3 = data[i+3] | c7;
	o4 = data[i+4] | d7;
	o5 = data[i+5] | e7;
    
    *p++ = o1;    
    *p++ = o2;    
    *p++ = o3;    
    *p++ = o4;    
    *p++ = o5;
	
	return buf;
}

#define ITERATE(index,start,end) for(size_t index = start; index < end; index += 1)
#define STEP(index,start,end,step) for(size_t index = start; index < end; index += step)


void print_pattern(uint8_t * p)
{
    for(size_t i = 0; i < 293; i++) cout << (int)p[i] << ",";
    cout << endl;
}


Random& GetRandom()
{
    static Random r;
    return r;
}

void sendNRPN(int nrpn, int value)
{
    int nrpn_mb = (nrpn >> 7) & 0x7F;
    int nrpn_lb = (nrpn & 0x7F);
    int val_mb  = (value >> 7) & 0x7F;
    int val_lb  = value & 0x7F;
    
    SendMidiMessage(0xB0,99,nrpn_mb);    
    SendMidiMessage(0xB0,98,nrpn_lb);    
    SendMidiMessage(0xB0,6,val_mb);
    SendMidiMessage(0xB0,38,val_lb);
    
} 


void LoadPatterns()
{
    FILE * f;
    glob_t g;        
    glob("Mopho/*.syx",GLOB_TILDE,NULL,&g);
    printf("%d\n",g.gl_pathc);
    for(size_t i = 0; i < g.gl_pathc; i++)
    {
        string path = string(g.gl_pathv[i]);            
        uint8_t buffer[300];            
        memset(buffer,0x0,sizeof(buffer));
        f = fopen(path.c_str(),"rb");                    
        fseek(f,0,SEEK_END);
        int fsize = ftell(f);            
        fseek(f,0,SEEK_SET);            
        size_t r = fread(buffer,1,fsize-7,f);            
        fclose(f);
        uint8_t* pattern = UnpackBits(&buffer[6]);            
        std::vector<uint8_t> p(293);
        memcpy(p.data(),pattern,293);
        patterns.push_back(p);
        free(pattern);
    }
    printf("Patterns=%d\n",patterns.size());
}

int random_parameter(int n)
{
    Random r;
    switch(n)
    {
    case 0: return r.randint(0,120);
    case 1: return r.randint(0,100);
    case 2: return r.randint(0,103);
    case 3: return r.randint(0,127);
    case 4: return r.randint(0,1);
    case 5: return r.randint(0,127);    
    case 6: return r.randint(0,120);
    case 7: return r.randint(0,100);
    case 8: return r.randint(0,103);
    case 9: return r.randint(0,127);
    case 10: return r.randint(0,1);
    case 11: return r.randint(0,127);    
    case 12: return r.randint(0,1);
    case 13: return r.randint(0,3);
    case 14: return r.randint(0,5);
    case 15: return r.randint(0,12);
    case 16: return r.randint(0,5);
    case 17: return r.randint(0,127);
    case 18: return r.randint(0,127);
    case 19: return r.randint(0,127);
    
    case 20: return r.randint(0,164);
    case 21: return r.randint(0,127);
    case 22: return r.randint(0,127);
    case 23: return r.randint(0,127);
    case 24: return r.randint(0,1);
    case 25: return r.randint(0,254);
    case 26: return r.randint(0,127);
    case 27: return r.randint(0,127);
    case 28: return r.randint(0,127);
    case 29: return r.randint(0,127);
    case 30: return r.randint(0,127);
    case 31: return r.randint(0,127);

    case 32: return r.randint(0,127);
    case 33: return r.randint(0,127);
    case 34: return r.randint(0,127);
    case 35: return r.randint(0,127);
    case 36: return r.randint(0,127);
    case 37: return r.randint(0,127);
    case 38: return r.randint(0,127);
    case 39: return r.randint(0,127);
    case 40: return r.randint(0,127);
    case 41: return r.randint(0,127);

    case 42: return r.randint(0,166);
    case 43: return r.randint(0,4);    
    case 44: return r.randint(0,127);    
    case 45: return r.randint(0,47);    
    case 46: return r.randint(0,1);

    case 47: return r.randint(0,166);
    case 48: return r.randint(0,4);    
    case 49: return r.randint(0,127);    
    case 50: return r.randint(0,47);    
    case 51: return r.randint(0,1);

    case 52: return r.randint(0,166);
    case 53: return r.randint(0,4);    
    case 54: return r.randint(0,127);    
    case 55: return r.randint(0,47);    
    case 56: return r.randint(0,1);

    case 57: return r.randint(0,166);
    case 58: return r.randint(0,4);    
    case 59: return r.randint(0,127);    
    case 60: return r.randint(0,47);    
    case 61: return r.randint(0,1);

    case 62: return r.randint(0,47);
    case 63: return r.randint(0,254);
    case 64: return r.randint(0,127);
    case 65: return r.randint(0,127);
    case 66: return r.randint(0,127);
    case 67: return r.randint(0,127);
    case 68: return r.randint(0,127);
    case 69: return r.randint(0,127);
    case 70: return r.randint(0,1);

    case 71: return r.randint(0,22);
    case 72: return r.randint(0,254);
    case 73: return r.randint(0,47);

    case 74: return r.randint(0,22);
    case 75: return r.randint(0,254);
    case 76: return r.randint(0,47);

    case 77: return r.randint(0,22);
    case 78: return r.randint(0,254);
    case 79: return r.randint(0,47);

    case 80: return r.randint(0,22);
    case 81: return r.randint(0,254);
    case 82: return r.randint(0,47);

    case 83: return r.randint(0,254);
    case 84: return r.randint(0,47);

    case 85: return r.randint(0,254);
    case 86: return r.randint(0,47);

    case 87: return r.randint(0,254);
    case 88: return r.randint(0,47);

    case 89: return r.randint(0,254);
    case 90: return r.randint(0,47);

    case 91: return r.randint(0,254);
    case 92: return r.randint(0,47);

    case 101: return r.randint(30,250);
    case 102: return r.randint(0,12);

    case 103: return r.randint(0,14);
    case 104: return r.randint(0,1);

    case 105: return r.randint(0,5);    
    case 106: return r.randint(0,1);    
    case 107: return r.randint(0,48);    
    case 108: return r.randint(0,48);    
    case 109: return r.randint(0,48);    
    case 110: return r.randint(0,48); 

    
    case 120: 
    case 121:
    case 122:
    case 123:
    case 124:
    case 125:
    case 126:
    case 127:
    case 128:
    case 129:
    case 130: 
    case 131:
    case 132:
    case 133:
    case 134:
    case 135:
    case 136:
    case 137:
    case 138:
    case 139:
    case 140: 
    case 141:
    case 142:
    case 143:
    case 144:
    case 145:
    case 146:
    case 147:
    case 148:
    case 149:
    case 150: 
    case 151:
    case 152:
    case 153:
    case 154:
    case 155:
    case 156:
    case 157:
    case 158:
    case 159:
    case 160: 
    case 161:
    case 162:
    case 163:
    case 164:
    case 165:
    case 166:
    case 167:
    case 168:
    case 169:
    case 170: 
    case 171:
    case 172:
    case 173:
    case 174:
    case 175:
    case 176:
    case 177:
    case 178:
    case 179:
    case 180: 
    case 181:
    case 182:
    case 183:
        return r.randint(0,127);    
    }
    return 0;
}
void Osc1(uint8_t * p)
{
    Random r;
    p[0] = r.randint(0,120);
    p[1] = r.randint(0,100);
    p[2] = r.randint(0,103);
    p[3] = r.randint(0,127);
    p[4] = r.randint(0,1);
    p[5] = r.randint(0,127);    
}
void InterpOsc1(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[0] = setter(0);
    p[1] = setter(1);
    p[2] = setter(2);
    p[3] = setter(3);
    p[4] = setter(4);
    p[5] = setter(5);
}
void Osc2(uint8_t * p)
{
    Random r;
    p[6] = r.randint(0,120);
    p[7] = r.randint(0,100);
    p[8] = r.randint(0,103);
    p[9] = r.randint(0,127);
    p[10] = r.randint(0,1);
    p[11] = r.randint(0,127);    
}
void InterpOsc2(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{    
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[6] = setter(6);
    p[7] = setter(7);
    p[8] = setter(8);
    p[9] = setter(9);
    p[10] = setter(10);
    p[11] = setter(11);
}
void OscMisc(uint8_t * p)
{
    Random r;
    p[12] = r.randint(0,1);
    p[13] = r.randint(0,3);
    p[14] = r.randint(0,5);
    p[15] = r.randint(0,12);
    p[16] = r.randint(0,5);
    p[17] = r.randint(0,127);
    p[18] = r.randint(0,127);
    p[19] = r.randint(0,127);
}
void InterpOscMisc(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[12] = setter(12);
    p[13] = setter(13);
    p[14] = setter(14);
    p[15] = setter(15);
    p[16] = setter(16);
    p[17] = setter(17);
    p[18] = setter(18);
    p[19] = setter(19);
}
void Filter(uint8_t * p)
{
    Random r;
    p[20] = r.randint(0,164);
    p[21] = r.randint(0,127);
    p[22] = r.randint(0,127);
    p[23] = r.randint(0,127);
    p[24] = r.randint(0,1);
    p[25] = r.randint(0,254);
    p[26] = r.randint(0,127);
    p[27] = r.randint(0,127);
    p[28] = r.randint(0,127);
    p[29] = r.randint(0,127);
    p[30] = r.randint(0,127);
    p[31] = r.randint(0,127);
}
void InterpFilter(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[20] = setter(20);
    p[21] = setter(21);
    p[22] = setter(22);
    p[23] = setter(23);
    p[24] = setter(24);
    p[25] = setter(25);
    p[26] = setter(26);
    p[27] = setter(27);
    p[28] = setter(28);
    p[29] = setter(29);
    p[30] = setter(30);
    p[31] = setter(31);
}
void VCA(uint8_t * p)
{
    Random r;
    p[32] = r.randint(0,127);
    p[33] = r.randint(0,127);
    p[34] = r.randint(0,127);
    p[35] = r.randint(0,127);
    p[36] = r.randint(0,127);
    p[37] = r.randint(0,127);
    p[38] = r.randint(0,127);
    p[39] = r.randint(0,127);
    p[40] = r.randint(0,127);
    p[41] = 127;
}
void InterpVCA(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[32] = setter(32);
    p[33] = setter(33);
    p[34] = setter(34);
    p[35] = setter(35);
    p[36] = setter(36);
    p[37] = setter(37);
    p[38] = setter(38);
    p[39] = setter(39);
    p[40] = setter(40);
    p[41] = 127;
}
void LFO(uint8_t * p)
{
    Random r;
    
    //1
    p[42] = r.randint(0,166);
    p[43] = r.randint(0,4);    
    p[44] = r.randint(0,127);    
    p[45] = r.randint(0,47);    
    p[46] = r.randint(0,1);

    //2
    p[47] = r.randint(0,166);
    p[48] = r.randint(0,4);
    p[49] = r.randint(0,127);
    p[50] = r.randint(0,47);
    p[51] = r.randint(0,1);
    

    //3
    p[52] = r.randint(0,166);
    p[53] = r.randint(0,4);
    p[54] = r.randint(0,127);
    p[55] = r.randint(0,47);
    p[56] = r.randint(0,1);

    
    //4
    p[57] = r.randint(0,166);
    p[58] = r.randint(0,4);
    p[59] = r.randint(0,127);
    p[60] = r.randint(0,47);
    p[61] = r.randint(0,1);
    
}
void InterpLFO(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };

    //1
    p[42] = setter(42);
    p[43] = setter(43);
    p[44] = setter(44);
    p[45] = setter(45);
    p[46] = setter(46);
    
    //2
    p[47] = setter(47);
    p[48] = setter(48);
    p[49] = setter(49);
    p[50] = setter(50);
    p[51] = setter(51);
    
    //3
    p[52] = setter(52);
    p[53] = setter(53);
    p[54] = setter(54);
    p[55] = setter(55);
    p[56] = setter(56);

    //4
    p[57] = setter(57);
    p[58] = setter(58);
    p[59] = setter(59);
    p[60] = setter(60);
    p[61] = setter(61);   
}
void Envelope3(uint8_t * p)
{
    Random r;
    p[62] = r.randint(0,47);
    p[63] = r.randint(0,254);
    p[64] = r.randint(0,127);
    p[65] = r.randint(0,127);
    p[66] = r.randint(0,127);
    p[67] = r.randint(0,127);
    p[68] = r.randint(0,127);
    p[69] = r.randint(0,127);
    p[70] = r.randint(0,1);
}
void InterpEnvelope3(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };
    p[62] = setter(62);
    p[63] = setter(63);
    p[64] = setter(64);
    p[65] = setter(65);
    p[66] = setter(66);
    p[67] = setter(67);
    p[68] = setter(68);
    p[69] = setter(69);
    p[70] = setter(70);
}
void Mod(uint8_t * p)
{
    Random r;
    // 1
    p[71] = r.randint(0,22);
    p[72] = r.randint(0,254);
    p[73] = r.randint(0,47);

    // 2
    p[74] = r.randint(0,22);
    p[75] = r.randint(0,254);
    p[76] = r.randint(0,47);

    // 3
    p[77] = r.randint(0,22);
    p[78] = r.randint(0,254);
    p[79] = r.randint(0,47);

    // 4
    p[80] = r.randint(0,22);
    p[81] = r.randint(0,254);
    p[82] = r.randint(0,47);

    // Mod Wheel
    p[83] = r.randint(0,254);
    p[84] = r.randint(0,47);
    
    // pressure
    p[85] = r.randint(0,254);
    p[86] = r.randint(0,47);

    // breath
    p[87] = r.randint(0,254);
    p[88] = r.randint(0,47);

    // velocity
    p[89] = r.randint(0,254);
    p[90] = r.randint(0,47);

    // foot control
    p[91] = r.randint(0,254);
    p[92] = r.randint(0,47);
    
}
void InterpMod(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    // 1
    p[71] = setter(71);
    p[72] = setter(72);
    p[73] = setter(73);

    // 2
    p[74] = setter(74);
    p[75] = setter(75);
    p[76] = setter(76);

    // 3
    p[77] = setter(77);
    p[78] = setter(78);
    p[79] = setter(79);

    // 4
    p[80] = setter(80);
    p[81] = setter(81);
    p[82] = setter(82);

    // Mod Wheel
    p[83] = setter(83);
    p[84] = setter(84);
    
    // pressure
    p[85] = setter(85);
    p[86] = setter(86);

    // breath
    p[87] = setter(87);
    p[88] = setter(88);

    // velocity
    p[89] = setter(89);
    p[90] = setter(90);

    // foot control
    p[91] = setter(91);
    p[92] = setter(92);
    
}
void Clock(uint8_t * p)
{
    Random rand;
    p[101] = rand.randint(30,250);
    p[102] = rand.randint(0,12);
}
void InterpClock(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    // 1
    p[101] = setter(101);
    p[102] = setter(102);
}
void Arpeggiator(uint8_t * p)
{
    Random rand;
    p[103] = rand.randint(0,14);
    p[104] = rand.randint(0,1);
    
}
void InterpArpeggiator(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    p[103] = setter(103);
    p[104] = setter(104);
    
}
void Sequencer(uint8_t * p)
{
    Random rand;
    p[105] = rand.randint(0,5);    
    p[106] = rand.randint(0,1);    
    p[107] = rand.randint(0,48);    
    p[108] = rand.randint(0,48);    
    p[109] = rand.randint(0,48);    
    p[110] = rand.randint(0,48);    
}
void InterpSequencer(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    p[105] = setter(105);
    p[106] = setter(106);
    p[107] = setter(107);
    p[108] = setter(108);
    p[109] = setter(109);
    p[110] = setter(110);
}
void Seq1(uint8_t * p)
{
    Random rand;
    for(size_t i = 0; i < 16; i++)
        p[120 + i] = rand.randint(0,127);    
}
void InterpSeq1(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    for(size_t i = 0; i < 16; i++)
        p[120 + i] = setter(120+i);
}
void Seq2(uint8_t * p)
{
    Random rand;
    for(size_t i = 0; i < 16; i++)
        p[136+ i] = rand.randint(0,126);    
}
void InterpSeq2(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    for(size_t i = 0; i < 16; i++)
        p[136+ i] = setter(136+i);
}
void Seq3(uint8_t * p)
{
    Random rand;
    for(size_t i = 0; i < 16; i++)
        p[152 + i] = rand.randint(0,126);    
}
void InterpSeq3(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    for(size_t i = 0; i < 16; i++)
        p[152 + i] = setter(152);
}
void Seq4(uint8_t * p)
{
    Random rand;
    for(size_t i = 0; i < 16; i++)
        p[168 + i] = rand.randint(0,126);    
}
void InterpSeq4(uint8_t * p, float f, uint8_t * p1, uint8_t * p2)
{
    auto setter = [f,p1,p2](int n) { return p1[n] + round(f*(p2[n]-p1[n])); };    
    for(size_t i = 0; i < 16; i++)
        p[168 + i] = setter(168+i);
}
void write(uint8_t * sysex)
{    
    FILE * file = fopen("output.syx","wb");
    fwrite(sysex,1,261,file);
    fclose(file);        
}
std::vector<uint8_t> doRandom()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Osc1(pattern2);
    Osc2(pattern2);
    OscMisc(pattern2);
    Filter(pattern2);
    VCA(pattern2);
    LFO(pattern2);
    Envelope3(pattern2);
    Mod(pattern2);    
    Clock(pattern2);
    Arpeggiator(pattern2);
    Sequencer(pattern2);
    Seq1(pattern2);
    Seq2(pattern2);
    Seq3(pattern2);
    Seq4(pattern2);
    return r;
}
std::vector<uint8_t> doInterp(float f = 0.5,int patchno=-1)
{
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpOsc1(pattern2,f,pattern1,pattern2);
    InterpOsc2(pattern2,f,pattern1,pattern2);
    InterpOscMisc(pattern2,f,pattern1,pattern2);
    InterpFilter(pattern2,f,pattern1,pattern2);
    InterpVCA(pattern2,f,pattern1,pattern2);
    InterpLFO(pattern2,f,pattern1,pattern2);
    InterpEnvelope3(pattern2,f,pattern1,pattern2);
    InterpMod(pattern2,f,pattern1,pattern2);
    InterpClock(pattern2,f,pattern1,pattern2);
    InterpArpeggiator(pattern2,f,pattern1,pattern2);
    InterpSequencer(pattern2,f,pattern1,pattern2);
    InterpSeq1(pattern2,f,pattern1,pattern2);
    InterpSeq2(pattern2,f,pattern1,pattern2);
    InterpSeq3(pattern2,f,pattern1,pattern2);
    InterpSeq4(pattern2,f,pattern1,pattern2);
    return r;
}
std::vector<uint8_t> doRandomOsc1()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Osc1(pattern2);    
    return r;
}
std::vector<uint8_t> doSwap(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    for(size_t i = 0; i < r.size(); i++)
        if(rand.frand() < f) r[i] = pattern1[i];
    return r;
}
std::vector<uint8_t> doRandomInjection(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    for(size_t i = 0; i < r.size(); i++)
        if(rand.frand() < f) r[i] = random_parameter(i);
    return r;
}
std::vector<uint8_t> doInterpProb(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    for(size_t i = 0; i < r.size(); i++)
        if(rand.frand() < f) r[i] = r[i] + round(f*(pattern1[i] - r[i]));
    return r;
}
std::vector<uint8_t> doInterp1(int index, float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    r[index] = r[index] + round(f*(pattern1[index] - r[index]));
    return r;
}
std::vector<uint8_t> doInterpBlock(int index1, int index2, float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    for(size_t i = index1; i < index2; i++)
        r[i] = r[i] + round(f*(pattern1[i] - r[i]));
    return r;
}
std::vector<uint8_t> doInterpOsc1(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpOsc1(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomOsc2()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Osc2(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpOsc2(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpOsc2(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomOscMisc()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    OscMisc(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpOscMisc(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpOscMisc(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomFilter()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Filter(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpFilter(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpFilter(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomVCA()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    VCA(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpVCA(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpVCA(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomLFOs()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    LFO(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpLFOs(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpLFO(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomEnv3()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Envelope3(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpEnv3(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpEnvelope3(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomMod()
{
    std::vector<uint8_t> r(256);    
    r = patch;
    uint8_t * pattern2  = r.data();        
    Mod(pattern2);        
    return r;
}
std::vector<uint8_t> doInterpMod(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpMod(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomClock()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Clock(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpClock(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpClock(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomArpeggiator()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Arpeggiator(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpArpeggiator(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpArpeggiator(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomSequencer()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Sequencer(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpSequencer(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpSequencer(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomSeq1()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Seq1(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpSeq1(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpSeq1(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomSeq2()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Seq2(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpSeq2(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpSeq2(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomSeq3()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Seq3(pattern2);    
    return r;
}
std::vector<uint8_t> doInterpSeq3(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpSeq3(pattern2,f,pattern1,pattern2);    
    return r;
}
std::vector<uint8_t> doRandomSeq4()
{
    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    Seq4(pattern2);    
    return r;
}

std::vector<uint8_t> doInterpSeq4(float f = 0.5,int patchno=-1)
{    
    Random rand;
    uint8_t* pattern1 = patterns[rand.randint(0,patterns.size()-1)].data();
    if(patchno != -1) pattern1 = patterns[patchno].data();

    std::vector<uint8_t> r(256);
    r = patch;
    uint8_t * pattern2  = r.data();        
    InterpSeq4(pattern2,f,pattern1,pattern2);    
    return r;
}

void get_patch(int n)
{    
    patch = patterns[n];
}
std::vector<uint8_t> get_random_patch()
{
    Random rand;    
    int n = rand.randint(0,patterns.size()-1);
    return patterns[n];
}
void random_patch()
{
    Random rand;    
    patch = patterns[rand.randint(0,patterns.size()-1)];
}


std::vector<uint8_t> mopho_patch(float f = 0.5)
{
    
    Random r;
    std::vector<uint8_t> pattern1 = patterns[r.randint(0,patterns.size()-1)];
    std::vector<uint8_t> pattern2 = patterns[r.randint(0,patterns.size()-1)];
    
    
    std::vector<uint8_t> pattern(293);
    for(size_t i = 0; i < 293; i++)
    {
        float x1 = float(pattern1[i]);
        float x2 = float(pattern2[i]);
        pattern[i] = round( x1 + f*(x2-x1));
    }
 
    //pattern[24] = 1;
    pattern[32] = 0;   
    pattern[40] = 127;
    pattern[41] = 127;
    
    /*
    std::vector<uint8_t> sysex(300);
    uint8_t * p = sysex.data();
    *p++ = 0xF0;
    *p++ = 0x1;
    *p++ = 0x25;
    *p++ = 0x3;
    
    uint8_t* packed = PackBits(pattern.data());
    for(size_t i = 0; i < 293; i++) {
        *p++ = packed[i];
    }
    *p = 0xF7;
    free(packed);
    printf("%d\n",sysex.size());
    
    SendSysex((char*)sysex.data()); 
    write(sysex.data());
    */
    return pattern;
}

// it's better to send NRPN
// not all the parameters in the sysex can be modified in the edit buffer like LFOs
// it will just ignore it on the mopho keyboard.
void sendNRPNPatch()
{
    uint8_t * p = patch.data();
    
    for(size_t  x= 0; x < 1; x++)
    {
    int i = 0;
    if(x == 1) i = 200;
    
    sendNRPN(0+i,p[0]);
    sendNRPN(1+i,p[1]);
    sendNRPN(2+i,p[2]);
    sendNRPN(3+i,p[3]);
    sendNRPN(4+i,p[4]);
    sendNRPN(114+i,p[5]);

    sendNRPN(5+i,p[6]);
    sendNRPN(6+i,p[7]);
    sendNRPN(7+i,p[8]);
    sendNRPN(8+i,p[9]);
    sendNRPN(9+i,p[10]);
    sendNRPN(115+i,p[11]);
    sendNRPN(10+i,p[12]);
    sendNRPN(11+i,p[13]);
    sendNRPN(12+i,p[14]);
    sendNRPN(93+i,p[15]);
    sendNRPN(13+i,p[16]);
    sendNRPN(14+i,p[17]);

    sendNRPN(116+i,p[18]);
    sendNRPN(110+i,p[19]);

    sendNRPN(15+i,p[20]);
    sendNRPN(16+i,p[21]);
    sendNRPN(17+i,p[22]);
    sendNRPN(18+i,p[23]);
    sendNRPN(19+i,p[24]);
    sendNRPN(20+i,p[25]);
    sendNRPN(21+i,p[26]);
    sendNRPN(22+i,p[27]);
    sendNRPN(23+i,p[28]);
    sendNRPN(24+i,p[29]);
    sendNRPN(25+i,p[30]);
    sendNRPN(26+i,p[31]);

    sendNRPN(27+i,p[32]);
    sendNRPN(30+i,p[33]);
    sendNRPN(31+i,p[34]);    
    sendNRPN(32+i,p[35]);    
    sendNRPN(33+i,p[36]);
    sendNRPN(34+i,p[37]);
    sendNRPN(35+i,p[38]);
    sendNRPN(36+i,p[39]);

    sendNRPN(29+i,127);
    

    // fuck sakes
    sendNRPN(37+i,p[42]);
    sendNRPN(38+i,p[43]);
    sendNRPN(39+i,p[44]);
    sendNRPN(40+i,p[45]);
    sendNRPN(41+i,p[46]);
    
    sendNRPN(42+i,p[47]);
    sendNRPN(43+i,p[48]);
    sendNRPN(44+i,p[49]);
    sendNRPN(45+i,p[50]);
    sendNRPN(46+i,p[51]);

    sendNRPN(47+i,p[52]);
    sendNRPN(48+i,p[53]);
    sendNRPN(49+i,p[54]);
    sendNRPN(50+i,p[55]);
    sendNRPN(51+i,p[56]);

    sendNRPN(52+i,p[57]);
    sendNRPN(53+i,p[58]);
    sendNRPN(54+i,p[59]);
    sendNRPN(55+i,p[60]);
    sendNRPN(56+i,p[61]);

    sendNRPN(57+i,p[62]);
    sendNRPN(58+i,p[63]);
    sendNRPN(59+i,p[64]);
    sendNRPN(60+i,p[65]);
    sendNRPN(61+i,p[66]);
    sendNRPN(62+i,p[67]);
    sendNRPN(63+i,p[68]);
    sendNRPN(64+i,p[69]);
    
    sendNRPN(98+i,p[70]);

    sendNRPN(65+i,p[71]);
    sendNRPN(66+i,p[72]);
    sendNRPN(67+i,p[73]);
    sendNRPN(68+i,p[74]);
    sendNRPN(69+i,p[75]);
    sendNRPN(70+i,p[76]);
    sendNRPN(71+i,p[77]);
    sendNRPN(72+i,p[78]);
    sendNRPN(73+i,p[79]);
    sendNRPN(74+i,p[80]);
    sendNRPN(75+i,p[81]);
    sendNRPN(76+i,p[82]);

    sendNRPN(81+i,p[83]);
    sendNRPN(82+i,p[84]);
    sendNRPN(83+i,p[85]);
    sendNRPN(84+i,p[86]);
    sendNRPN(85+i,p[87]);
    sendNRPN(86+i,p[88]);
    sendNRPN(87+i,p[89]);
    sendNRPN(88+i,p[90]);
    sendNRPN(89+i,p[91]);
    sendNRPN(90+i,p[92]);

    
    sendNRPN(91+i,p[101]);
    sendNRPN(92+i,p[102]);

    sendNRPN(97+i,p[103]);
    sendNRPN(100+i,p[104]);
    sendNRPN(94+i,p[105]);
    sendNRPN(101+i,p[106]);
    sendNRPN(77+i,p[107]);
    sendNRPN(78+i,p[108]);
    sendNRPN(79+i,p[109]);
    sendNRPN(80+i,p[110]);
    
    for(size_t j = 120; j < 184; j++)
        sendNRPN(i+j,p[j]);
    
       
    }
}


void SendPatch(uint8_t * ptr)
{
    std::vector<uint8_t> sysex(600);
    uint8_t * p = sysex.data();
    ptr[40]=127;
    *p++ = 0xF0;
    *p++ = 0x1;
    *p++ = 0x25;
    *p++ = 0x3;
    uint8_t* packed = PackBits(ptr);
    for(size_t i = 0; i < 293; i++) {
        *p++ = packed[i];
    }
    *p = 0xF7;
    
    free(packed);    
    SendSysex((char*)sysex.data()); 
    Pa_Sleep(5);

    std::vector<uint8_t> tetrap(500);
    // for tetra 
    memcpy(tetrap.data(),ptr,200);
    for(size_t i = 200; i < 383; i++)
        tetrap[i] = ptr[i-200];
    
    
    
    packed = PackTetraBits(tetrap.data());
    p = sysex.data();
    *p++ = 0xF0;
    *p++ = 0x1;
    *p++ = 0x26;
    *p++ = 0x3;
    for(size_t i = 0; i < 383; i++) {
        *p++ = packed[i];
    }
    *p = 0xF7;
    free(packed);    
    SendSysex((char*)sysex.data()); 
    Pa_Sleep(5);
    
}


//////////////////////////////////////////////////////////////////
// Audio Stream
//////////////////////////////////////////////////////////////////


float ** temp1;
float ** temp2;

bool flangerOn=false;
LV2Plugins * plugin;
LV2Plugin * lv2flangerL,* lv2flangerR;

JoonasFX::DelayEffect jdelayfx(sampleRate);
JoonasFX::VibraFlange jflanger;
JoonasFX::Chorus jchorus;
JoonasFX::Phaser jphaser;
MVerb<float> mverb;
FX::Chorus chorus;

/*

RingModulator ring;
Tremolo trem;
DistortionEffect dist(DistortionEffect::_dcSigmoid,12);
Moog::MoogFilter1 filter(sampleRate,100,0.5);
AnalogSVF svf(sampleRate,1000,0.5);
AutoWah awah;
StereoCompressor compressor;
DelayEffect delay;
Flanger flanger;
Phaser phaser;
PingPongDelay pingpong;
PVPitchShifter pitch;
Vibrato vibrato;
ChorusDelay cdelay;
EarlyReflectionReverb rev1(sampleRate);
HallReverb rev2(sampleRate);
PlateReverb rev3(sampleRate);
EarlyDelayReverb rev4(sampleRate);
FX::ChorusEngine ykchorus(sampleRate);
FX::Ce2Chorus ce2L(sampleRate),ce2R(sampleRate);
Stereofyier stereofy;

*/
// all filters mostly work
/*
Filters::FilterBase::FilterType type = Filters::FilterBase::FilterType::LOWPASS ;
Filters::BiquadFilter     filter(type,1000,sampleRate,0.5);
Filters::Butterworth24db  butter(sampleRate,800.0,0.5);
Filters::TwoZeroFilter    twozero(type,sampleRate,1000.0,0.5);
Filters::TwoPoleFilter    twopole(type,sampleRate,1000.0,0.5);
Filters::OnePoleOneZeroFilter polezero(type,sampleRate,1000.0,0.5);
Filters::OnePoleFilter    onepole(type,sampleRate,1000.0,0.5);
Filters::OneZeroFilter    onezero(type,sampleRate,1000.0,0.5);
*/

/*
void process_lv2(size_t framesPerBuffer,float ** temp1, float ** temp2)
{
    if(flangerOn)
    {    
        lv2flangerL->Run(framesPerBuffer,temp1[0],temp1[0]);
        lv2flangerR->Run(framesPerBuffer,temp1[1],temp1[1]);
    }
}
*/
void process_effects(size_t framesPerBuffer,float ** temp1, float ** temp2)
{
    
    jchorus.InplaceProcess(framesPerBuffer,temp1);
    jflanger.InplaceProcess(framesPerBuffer,temp1);
    jphaser.InplaceProcess(framesPerBuffer,temp1);
    jdelayfx.InplaceProcess(framesPerBuffer,temp1);

    //pitch.ProcessBlock(framesPerBuffer,temp1,temp2);
    //cdelay.ProcessBlock(framesPerBuffer,temp1,temp2);
    //vflanger.ProcessBlock(framesPerBuffer,temp1,temp2);
    //jphaser.ProcessBlock(framesPerBuffer,temp1,temp2);       
    //phaser.InplaceProcess(framesPerBuffer,temp2);
    //vibrato.InplaceProcess(framesPerBuffer,temp2);
    //pingpong.InplaceProcess(framesPerBuffer,temp2);
    //ring.ProcessBlock(framesPerBuffer,temp2,temp2);   
    //chorus.ProcessBlock(framesPerBuffer,temp2,temp2);   
    //trem.ProcessBlock(framesPerBuffer,temp2,temp2);   
    //awah.ProcessBlock(framesPerBuffer,temp1,temp2);   
    //flanger.InplaceProcess(framesPerBuffer,temp2);    
    //ykchorus.ProcessBlock(framesPerBuffer,temp1,temp2);
    //compressor.InplaceProcess(framesPerBuffer,temp2);       
    
    //rev2.run(framesPerBuffer,temp2,temp1);        
    //convolver.InplaceProcess(framesPerBuffer,temp2);
    //stereofy.InplaceProcess(framesPerBuffer,temp2);    

    mverb.process(temp1,temp2,framesPerBuffer);
    //limiter.Process(temp1,temp2,framesPerBuffer);
}
int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    
    float * output = (float*)outputBuffer;
    float * input  = (float*)inputBuffer;
    float * p = output;
    float gain = pow(10,-1);
    for(size_t i = 0; i < framesPerBuffer; i++)
    {        
        float in   =  gain*input[i*2];
        temp1[0][i] = in;
        temp1[1][i] = in;
    }
    
    process_effects(framesPerBuffer,temp1,temp2);

    for(size_t i = 0; i < framesPerBuffer; i++)
    {        
        *p++ = temp2[0][i];
        *p++ = temp2[1][i];
    }
    return 0;
}            


#include "template_midi"

void print_vector(std::vector<uint8_t> & p)
{
    for(size_t i = 0; i < p.size(); i++)
        printf("%d,",p[i]);
    printf("\n");
}
void sysex_callback(MidiMsg * m)
{
    extern uint8_t sysex_buffer[8192];   
        
    if(m->msg == 240)
    {
        printf("sysex\n");
        std::vector<uint8_t> buffer(293);
        for(size_t i = 4; i < 293; i++)
        {
            buffer[i-4] = sysex_buffer[i];
            
        }        
        uint8_t * unpacked = UnpackBits(buffer.data());
        patch.resize(256);
        memcpy(patch.data(),unpacked,256);
        print_vector(patch);
        free(unpacked);        
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lua REPL
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "LuaJIT.hpp"
LuaJIT * lua;


int test(lua_State *L)
{
    printf("test\n");
    return 0;
}

bool is_cmd(const char * cmd, const char * key)
{
    return !strcmp(cmd,key);
}

void strupr(char * s)
{
    for(size_t i = 0; i < strlen(s); i++)
        s[i] = toupper(s[i]);
}
int CreatePatch(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    patch = mopho_patch(f);
    return 0;
}
int Crossover(lua_State * L)
{
    Random rand;
    std::vector<uint8_t> r = get_random_patch();
    size_t n = lua_tonumber(L,-1);
    for(size_t i = 0; i < n; i++)
    {
        size_t x1 = rand.randint(0,256/2);
        size_t x2 = rand.randint(x1+1,256);
        memcpy(patch.data()+x1, r.data()+x1, x2-x1);        
    }
    return 0;
}
int Mutate(lua_State * L)
{
    Random r;    
    size_t n = lua_tonumber(L,-2);
    float  prob = lua_tonumber(L,-1);
    for(size_t i = 0; i < n; i++)
    {
        size_t x1 = r.randint(0,256/2);
        size_t x2 = r.randint(x1+1,256);
        for(size_t j = x1; j < x2; j++)
            if(r.frand() < prob) patch[j] = random_parameter(j);
    }
    return 0;
}
int RandomPatch(lua_State *L)
{
    std::vector<uint8_t> r = doRandom();
    patch = r;
    return 0;
}
int RandomOsc1(lua_State *L)
{
    std::vector<uint8_t> r = doRandomOsc1();
    patch = r;
    return 0;
}
int RandomOsc2(lua_State *L)
{
    std::vector<uint8_t> r = doRandomOsc2();
    patch = r;
    return 0;
}
int RandomOscMisc(lua_State *L)
{
    std::vector<uint8_t> r = doRandomOscMisc();
    patch = r;
    return 0;
}
int RandomFilter(lua_State *L)
{
    std::vector<uint8_t> r = doRandomFilter();
    patch = r;
    return 0;
}
int RandomVCA(lua_State *L)
{
    std::vector<uint8_t> r = doRandomVCA();
    patch = r;
    return 0;
}
int RandomEnv3(lua_State *L)
{
    std::vector<uint8_t> r = doRandomEnv3();
    patch = r;
    return 0;
}
int RandomLFOs(lua_State *L)
{
    std::vector<uint8_t> r = doRandomLFOs();
    patch = r;
    return 0;
}
int RandomMod(lua_State *L)
{
    std::vector<uint8_t> r = doRandomMod();
    patch = r;
    return 0;
}
int RandomClock(lua_State *L)
{
    std::vector<uint8_t> r = doRandomClock();
    patch = r;
    return 0;
}
int RandomArpeggiator(lua_State *L)
{
    std::vector<uint8_t> r = doRandomArpeggiator();
    patch = r;
    return 0;
}
int RandomSequencer(lua_State *L)
{
    std::vector<uint8_t> r = doRandomSequencer();
    patch = r;
    return 0;
}
int RandomSeq1(lua_State *L)
{
    std::vector<uint8_t> r = doRandomSeq1();
    patch = r;
    return 0;
}
int RandomSeq2(lua_State *L)
{
    std::vector<uint8_t> r = doRandomSeq2();
    patch = r;
    return 0;
}
int RandomSeq3(lua_State *L)
{
    std::vector<uint8_t> r = doRandomSeq3();
    patch = r;
    return 0;
}
int RandomSeq4(lua_State *L)
{
    std::vector<uint8_t> r = doRandomSeq4();
    patch = r;
    return 0;
}
int InterpPatch(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterp(f);
    patch = r;
    return 0;
}
int SwapPatch(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doSwap(f);
    patch = r;
    return 0;
}
int InjectPatch(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doRandomInjection(f);
    patch = r;
    return 0;
}
int InterpProb(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpProb(f);
    patch = r;
    return 0;
}
int Interp1(lua_State *L)
{
    int i = lua_tonumber(L,-2);
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterp1(i,f);
    patch = r;
    return 0;
}
int InterpBlock(lua_State *L)
{
    int i1 = lua_tonumber(L,-3);
    int i2 = lua_tonumber(L,-2);
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpBlock(i1,i2,f);
    patch = r;
    return 0;
}
int InterpOsc1(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpOsc1(f);
    patch = r;
    return 0;
}
int InterpOsc2(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpOsc2(f);
    patch = r;
    return 0;
}
int InterpOscMisc(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpOscMisc(f);
    patch = r;
    return 0;
}
int InterpFilter(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpFilter(f);
    patch = r;
    return 0;
}
int InterpVCA(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpVCA(f);
    patch = r;
    return 0;
}
int InterpLFOs(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpLFOs(f);
    patch = r;
    return 0;
}
int InterpEnv3(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpEnv3(f);
    patch = r;
    return 0;
}
int InterpMod(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpMod(f);
    patch = r;
    return 0;
}
int InterpClock(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpClock(f);
    patch = r;
    return 0;
}
int InterpArpeggiator(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpArpeggiator(f);
    patch = r;
    return 0;
}
int InterpSequencer(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpSequencer(f);
    patch = r;
    return 0;
}
int InterpSeq1(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpSeq1(f);
    patch = r;
    return 0;
}
int InterpSeq2(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpSeq2(f);
    patch = r;
    return 0;
}
int InterpSeq3(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpSeq3(f);
    patch = r;
    return 0;
}
int InterpSeq4(lua_State *L)
{
    float f = lua_tonumber(L,-1);
    std::vector<uint8_t> r = doInterpSeq4(f);
    patch = r;
    return 0;
}
int getPatch(lua_State * L)
{
    uint8_t sysex[] = {0xF0,0x01,0x27,0x6,0xF7};
    //uint8_t sysex[] = {0xF0,0x7E,0x01,0x6,0x1,0xF7};
    SendSysex((char*)&sysex[0]);  
    Pa_Sleep(10);    
    return 0;
}
int sendPatch(lua_State * L)
{    
    if(send_sysex == true)    
        SendPatch(patch.data());            
    else
        sendNRPNPatch();
    Pa_Sleep(10);    
    return 0;
}
int patchPoke(lua_State * L)
{
    uint32_t i = lua_tonumber(L,-2);
    uint8_t v = lua_tonumber(L,-1);
    printf("%d %d\n",i,v);
    patch[i] = v;    
    return 0;
}
int patchPeek(lua_State * L)
{
    uint32_t i = lua_tonumber(L,-1);    
    printf("patch[%d]=%d",i,patch[i]);    
    return 0;
}
int LoadPatch(lua_State * L)
{    
    int i = lua_tonumber(L,-1);
    get_patch(i);
    SendPatch(patch.data());
    return 0;
}
int LoadRandomPatch(lua_State * L)
{
    uint32_t i = lua_tonumber(L,-1);     
    random_patch();
    SendPatch(patch.data());    
    return 0;
    
}
int SavePatch(lua_State * L)
{
    const char * filename = lua_tostring(L,-1);
    return 0;
    
}

int setValueCmd(lua_State *L)
{
    char cmd[128];
    const char * cmdS = lua_tostring(L,-2);
    strcpy(cmd,cmdS);
    strupr((char*)cmd);
    
                          
    float val = lua_tonumber(L,-1);
    if(is_cmd(cmd,"LIST"))
    {
        std::cout << "DELAY,BPM,INVERT,FEEDBACK,WAVESHAPE\n";
    }
    if(is_cmd(cmd,"DELAY")) 
    {
        if(val < 0.0f) val = 0.0f;
        if(val > 1.0f) val = 1.0f;
        lv2flangerL->connections[0][0] = val;
        lv2flangerR->connections[0][0] = val;
    }
    if(is_cmd(cmd,"BPM")) {
        if(val < 5) val = 5.0f;
        if(val > 300.0f) val = 300.0f;
        lv2flangerL->connections[1][0] = val;
        lv2flangerR->connections[1][0] = val;
    }
    if(is_cmd(cmd,"DEPTH")) {
        if(val < 0) val = 0;
        if(val > 1) val = 1;
        lv2flangerL->connections[2][0] = val;
        lv2flangerR->connections[2][0] = val;
    }
    if(is_cmd(cmd,"INVERT"))
    {
        if(val < 0.5) val = 0;
        if(val > 0.5) val = 1.0;
        lv2flangerL->connections[3][0] = val;
        lv2flangerR->connections[3][0] = val;
    }
    if(is_cmd(cmd,"FEEDBACK"))
    {
        if(val < -0.995) val = -0.995;
        if(val > 0.995) val = 0.995;
        lv2flangerL->connections[4][0] = val;
        lv2flangerR->connections[4][0] = val;
    }
    if(is_cmd(cmd,"WAVESHAPE")) 
    {
        if(val < 0) val = 0;
        if(val > 1) val = 1;
        lv2flangerL->connections[5][0] = val;
        lv2flangerR->connections[5][0] = val;
    }    
    return 0;
}

int setFlangerBPM(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < 5) val = 5.0f;
    if(val > 300.0f) val = 300.0f;
    lv2flangerL->connections[1][0] = val;
    lv2flangerR->connections[1][0] = val;
    return 0;
}    
int setFlangerDelay(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < 0.0f) val = 0.0f;
    if(val > 1.0f) val = 1.0f;
    lv2flangerL->connections[0][0] = val;
    lv2flangerR->connections[0][0] = val;
    return 0;
}    
int setFlangerDepth(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < 0) val = 0;
    if(val > 1) val = 1;
    lv2flangerL->connections[2][0] = val;
    lv2flangerR->connections[2][0] = val;    
    return 0;
}    
int setFlangerInvert(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < 0.5) val = 0;
    if(val > 0.5) val = 1.0;
    lv2flangerL->connections[3][0] = val;
    lv2flangerR->connections[3][0] = val;
    return 0;
}    
int setFlangerFeedback(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < -0.995) val = -0.995;
    if(val > 0.995) val = 0.995;
    lv2flangerL->connections[4][0] = val;
    lv2flangerR->connections[4][0] = val;
    return 0;
}    
int setFlangerWaveshape(lua_State *L)
{
    float val = lua_tonumber(L,-1);
    if(val < 0) val = 0;
    if(val > 1) val = 1;
    lv2flangerL->connections[5][0] = val;
    lv2flangerR->connections[5][0] = val;
    return 0;
}    
int setFlangerOn(lua_State * L)
{
    flangerOn = !flangerOn;
    return 0;
}
int saveFlanger(lua_State * L)
{
    const char * file = lua_tostring(L,-1);
    FILE * f = fopen(file,"wb");
    for(size_t i = 0; i < 6; i++)
    {
        float val = lv2flangerL->connections[i][0];
        fwrite(&val,sizeof(val),1,f);
    }
    fclose(f);
    return 0;
}
int loadFlanger(lua_State * L )
{
    const char * file = lua_tostring(L,-1);
    FILE * f = fopen(file,"rb");
    float vals[6];
    fread(vals,sizeof(float),6,f);
    fclose(f);
    for(size_t i = 0; i < 6; i++)
    {
        lv2flangerL->connections[i][0] = vals[i];
        lv2flangerR->connections[i][0] = vals[i];        
    }
    return 0;
}
int randomFlanger(lua_State * L)
{
    lv2flangerL->connections[0][0] = noise.rand();
    lv2flangerR->connections[0][0] = noise.rand();

    lv2flangerL->connections[1][0] = 5 + 295*noise.rand();
    lv2flangerR->connections[1][0] = 5 + 295*noise.rand();

    lv2flangerL->connections[2][0] = noise.rand();
    lv2flangerR->connections[2][0] = noise.rand();

    lv2flangerL->connections[3][0] = noise.rand();
    lv2flangerR->connections[3][0] = noise.rand();

    lv2flangerL->connections[4][0] = noise.rand();
    lv2flangerR->connections[4][0] = noise.rand();

    lv2flangerL->connections[5][0] = noise.rand();
    lv2flangerR->connections[5][0] = noise.rand();

    return 0;

}
int sendnrpn(lua_State * L)
{
    int a = lua_tonumber(L,-2);
    int b = lua_tonumber(L,-1);
    sendNRPN(a,b);    
    return 0;
}
int nrpnsend(lua_State * L)
{
    /*
    for(size_t i = 42; i < 101; i++)
        patch[i] = mod[i];
    SendPatch(patch.data());            
    */
    Pa_Sleep(5);            
    sendNRPNPatch();    
    Pa_Sleep(5);
    return 0;
}

int jchorus_power(lua_State * L)
{
    bool val = lua_toboolean(L,-1);
    jchorus.power(val);
    return 0;
}
int jchorus_randomize(lua_State * L)
{
    jchorus.randomize();
    return 0;
}
int jchorus_freq(lua_State * L)
{
    float x = lua_tonumber(L,-1);
    jchorus.freqs[0] = x;
    jchorus.freqs[1] = x;
    jchorus.freqs[2] = x;
    return 0;
}
int jchorus_depth(lua_State * L)
{
    float x = lua_tonumber(L,-1);
    jchorus.depths[0] = x;
    jchorus.depths[1] = x;
    jchorus.depths[2] = x;
    return 0;
}
int jchorus_wet(lua_State * L)
{
    float x = lua_tonumber(L,-1);
    for(size_t i = 0; i < 3; i++) jchorus.wets[i] = x;
    return 0;
}
int jchorus_feedback(lua_State * L)
{
    float x = lua_tonumber(L,-1);
    for(size_t i = 0; i < 3; i++) jchorus.feedbacks[i] = x;
    return 0;
}
int jphaser_power(lua_State * L)
{
    bool on = lua_toboolean(L,-1);
    jphaser.power(on);
    return 0;
}
int jphaser_rand(lua_State * L)
{    
    jphaser.randomize();
    return 0;
}
int jflanger_power(lua_State * L)
{
    bool on = lua_toboolean(L,-1);
    jflanger.power(on);
    return 0;
}
int jflanger_rand(lua_State * L)
{    
    jflanger.randomize();
    return 0;
}
int mverb_power(lua_State * L)
{
    bool on = lua_toboolean(L,-1);
    mverb.power(on);
    return 0;
}
int mverb_rand(lua_State * L)
{
    mverb.randomize();
    return 0;
}
int jdelay_power(lua_State * L)
{
    bool on = lua_toboolean(L,-1);
    jdelayfx.power(on);
    return 0;
}
int jdelay_rand(lua_State * L)
{
    jdelayfx.randomize();
    return 0;
}
void connectLua()
{
    lua = new LuaJIT("main.lua");

    /*
    lua->CreateCFunction("flanger",setValueCmd);    
    lua->CreateCFunction("fbpm", setFlangerBPM);
    lua->CreateCFunction("fdelay", setFlangerDelay);
    lua->CreateCFunction("fdepth", setFlangerDepth);
    lua->CreateCFunction("finvert", setFlangerInvert);
    lua->CreateCFunction("ffeedback", setFlangerFeedback);
    lua->CreateCFunction("fwaveshape", setFlangerWaveshape);
    lua->CreateCFunction("fpower", setFlangerOn);
    lua->CreateCFunction("fsave", saveFlanger);
    lua->CreateCFunction("fload", loadFlanger);
    lua->CreateCFunction("frandom", randomFlanger);
    */

    lua->CreateCFunction("jchorus_power",jchorus_power);
    lua->CreateCFunction("jchorus_rand",jchorus_randomize);
    lua->CreateCFunction("jchorus_freq",jchorus_freq);    
    lua->CreateCFunction("jchorus_depth",jchorus_depth);
    lua->CreateCFunction("jchorus_wet",jchorus_wet);
    lua->CreateCFunction("jchorus_feedback",jchorus_feedback);

    lua->CreateCFunction("jdelay_power",jdelay_power);
    lua->CreateCFunction("jdelay_rand",jdelay_rand);

    lua->CreateCFunction("jphaser_power",jphaser_power);
    lua->CreateCFunction("jphaser_rand",jphaser_rand);

    lua->CreateCFunction("jflanger_power",jflanger_power);
    lua->CreateCFunction("jflanger_rand",jflanger_rand);

    lua->CreateCFunction("mverb_power",mverb_power);
    lua->CreateCFunction("mverb_rand",mverb_rand);

    lua->CreateCFunction("load",LoadPatch);    
    lua->CreateCFunction("save",SavePatch);    

    lua->CreateCFunction("create",CreatePatch);    
    lua->CreateCFunction("random",RandomPatch);    
    lua->CreateCFunction("swap",SwapPatch);
    lua->CreateCFunction("inject",InjectPatch);
    lua->CreateCFunction("cross",Crossover);
    lua->CreateCFunction("mutate",Mutate);

    lua->CreateCFunction("interp",InterpPatch);
    lua->CreateCFunction("interprob",InterpProb);
    lua->CreateCFunction("interp1",Interp1);
    lua->CreateCFunction("interpblock",InterpBlock);

    lua->CreateCFunction("get",getPatch);
    lua->CreateCFunction("send",sendPatch);

    lua->CreateCFunction("random_osc1",RandomOsc1);
    lua->CreateCFunction("random_osc2",RandomOsc2);
    lua->CreateCFunction("random_oscmisc",RandomOscMisc);
    lua->CreateCFunction("random_vcf",RandomFilter);
    lua->CreateCFunction("random_vca",RandomVCA);
    lua->CreateCFunction("random_lfos",RandomLFOs);
    lua->CreateCFunction("random_env3",RandomEnv3);
    lua->CreateCFunction("random_mod",RandomMod);
    lua->CreateCFunction("random_clock",RandomClock);
    lua->CreateCFunction("random_arp",RandomArpeggiator);
    lua->CreateCFunction("random_sequencer",RandomSequencer);
    lua->CreateCFunction("random_seq1",RandomSeq1);
    lua->CreateCFunction("random_seq2",RandomSeq2);
    lua->CreateCFunction("random_seq3",RandomSeq3);
    lua->CreateCFunction("random_seq4",RandomSeq4);


    lua->CreateCFunction("interp_osc1",InterpOsc1);
    lua->CreateCFunction("interp_osc2",InterpOsc2);
    lua->CreateCFunction("interp_oscmisc",InterpOscMisc);
    lua->CreateCFunction("interp_vcf",InterpFilter);
    lua->CreateCFunction("interp_vca",InterpVCA);
    lua->CreateCFunction("interp_lfos",InterpLFOs);
    lua->CreateCFunction("interp_env3",InterpEnv3);
    lua->CreateCFunction("interp_mod",InterpMod);
    lua->CreateCFunction("interp_clock",InterpClock);
    lua->CreateCFunction("interp_arp",InterpArpeggiator);
    lua->CreateCFunction("interp_sequencer",InterpSequencer);
    lua->CreateCFunction("interp_seq1",InterpSeq1);
    lua->CreateCFunction("interp_seq2",InterpSeq2);
    lua->CreateCFunction("interp_seq3",InterpSeq3);
    lua->CreateCFunction("interp_seq4",InterpSeq4);

    lua->CreateCFunction("nrpn",sendnrpn);
    lua->CreateCFunction("nrpn_send",nrpnsend);
    lua->CreateCFunction("poke",patchPoke);
    lua->CreateCFunction("peek",patchPeek);
    
}



void repl() {
    
    
    if(send_sysex) {
        uint8_t sysex[] = {0xF0,0x01,0x27,0x6,0xF7};    
        SendSysex((char*)&sysex[0]); 
    }
    
    std::string cmd;
    std::cin >> cmd;
    lua->DoCmd(cmd);
    
    Pa_Sleep(5);

    if(patch.size() > 0)
        if(send_sysex == true) 
        {
            //for(size_t i = 42; i < 101; i++)
            //    patch[i] = mod[i];
            //SendPatch(patch.data());            
            Pa_Sleep(5);              
            sendNRPNPatch();                        
            Pa_Sleep(5);
        }        
}

int main()
{    
    LoadPatterns();
    patch.resize(300);    
    //mod.resize(300);
    Init();
    
    connectLua();
    noise.seed_engine();     

    //
    //jdelayfx.power(false);
    jchorus.power(false);
    jflanger.power(false);
    jphaser.power(false);

    //sst::waveshapers::initializeWaveshaperRegister(sst::waveshapers::WaveshaperType::wst_fuzz, waveState.R);
    
    temp1 = (float**)calloc(2,sizeof(float*));
    temp1[0] = (float*)calloc(1024,sizeof(float));
    temp1[1] = (float*)calloc(1024,sizeof(float));

    temp2 = (float**)calloc(2,sizeof(float*));
    temp2[0] = (float*)calloc(1024,sizeof(float));
    temp2[1] = (float*)calloc(1024,sizeof(float));
            
    mverb.setSampleRate(sampleRate);
    mverb.setParameter(MVerb<float>::MIX,0.5);
    mverb.setParameter(MVerb<float>::SIZE,0.7);
    mverb.setParameter(MVerb<float>::PREDELAY,0.01);

    int num_midi = GetNumMidiDevices();
    ITERATE(i,0,num_midi)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    int num_audio = GetNumAudioDevices();
    int pulse = 8;

    ITERATE(i, 0, num_audio)    
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }
    
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    set_system_exclusive_func(sysex_callback);

    //plugin = new LV2Plugins;  
    /*
    lv2flangerL =plugin->LoadPlugin("http://polyeffects.com/lv2/flanger");
    lv2flangerL->connections[0][0] = 0.05;
    lv2flangerL->connections[1][0] = 5;
    lv2flangerL->connections[2][0] = 0.9;
    lv2flangerL->connections[3][0] = 1;
    lv2flangerL->connections[4][0] = 0.9;
    lv2flangerL->connections[5][0] = 0;   

    lv2flangerR =plugin->LoadPlugin("http://polyeffects.com/lv2/flanger");
    lv2flangerR->connections[0][0] = 0.05;
    lv2flangerR->connections[1][0] = 5;
    lv2flangerR->connections[2][0] = 0.9;
    lv2flangerR->connections[3][0] = 1;
    lv2flangerR->connections[4][0] = 0.1;
    lv2flangerR->connections[5][0] = 0;                            
    */  
    InitMidiDevice(1,3,2);
    InitAudioDevice(pulse,pulse,2,sampleRate,BufferSize);
    //Thread audio(run,NULL);
    //runfltk();
    //SendPatch(0);
    RunAudio();
    StopAudio();
}

