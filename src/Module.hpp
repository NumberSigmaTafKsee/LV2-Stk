#pragma once

struct Port
{
    float * port;
    size_t  n;
};

struct Connect
{
    Port * from;
    Port * to;
};

struct ControlMap
{
    std::map<std::string,float> map;

    ControlMap() {
        map["input"] = 0;
        map["output"] = 0;
    }
    void addControl(const std::string & name, float value=0.0f) {
        map[name] = value;
    }
    float getControl(const std::string & name)
    {
        if(map.find(name) == map.end()) return 0;
        return map[name];
    }
    void setValue(const std::string & name, float value) {
        if(map.find(name) == map.end()) return;
        map[name] = value;
    }
};

struct ControlVector
{
    std::vector<std::vector<float>> data;
    int channels;

    ControlVector(size_t n, size_t block) {
        data.resize(n);
        for(size_t i = 0; i < n; i++)
            data[i].resize(block);
    }

    float getValue(size_t channel, size_t i) {
        return data[channel][i];
    }
    void setValue(size_t channel, size_t i, float v) {
        data[channel][i] = v;
    }
};

struct Module
{
    ControlMap controls;
    
    virtual float Tick()   = 0;

    virtual void noteOn(float note, float velocity) {

    }
    virtual void noteOff(float note = 0) {

    }    
};

struct Connection
{
    Module * from;
    std::string from_port;
    Module * to;
    std::string to_port;
};

struct Rack
{
    std::vector<Module*> modules;
    std::vector<Connection*> connections;

    void Run()
    {
        for(auto i = modules.begin(); i != modules.end(); i++)
            *i->Tick();
        for(auto i = connections.begin(); i != connections.end(); i++)
        {
            *i->to->controls.setValue(*i->to_port, *i->from->controls.getValue(*i->from_port));
        }
    }
};

struct minBLEPModule : public Module
{
    minBLEP osc;
    osc_t * slave = NULL;

    minBLEPModule() : osc(sampleRate)
    {
        controls.map["Amp"] = 0;
        controls.map["Fm"]  = 0;
        controls.map["Am"]  = 0;        
        controls.map["Freq"] = 440.0f;
        controls.map["Sync"] = 0;
        controls.map["Waveform"] = OT_SAW;
    }
    ~minBLEPModule() {

    }
    void setWaveform( float wave )
    {
        float waveform = OT_SAW;
        if(wave == 0.0f) waveform = OT_SAW;
        if(wave == 1.0f) waveform = OT_SQUARE;
        if(wave == 2.0f) waveform = OT_TRIANGLE;
        if(wave == 3.0f) waveform = OT_RSAW;
        controls.map["Waveform"] = waveform;
    }
    void setFreq(float freq) {
        controls.map["Freq"] = freq;
    }
    void slaveOsc(minBLEPModule m) {
        slave = (m.osc.lpO);
    }
    float Tick()
    {
        float I = controls.map["In"];        
        float A = controls.map["Amp"];
        float Fm= controls.map["Fm"];
        float Pm= controls.map["Pm"];
        bool  sync = (bool)controls.map["Sync"];
        float Freq = controls.map["Freq"];
        int   Waveform = (int)controls.map["Waveform"];
        if(sync) {
            osc.setSlave(slave);
        }
        else {
            osc.setSlave(NULL);
        }
        osc.setFrequency(Freq);
        float out = osc.Tick(I,A,Fm,Pm);
        controls.map["Out"] = out;
        return out;
    }
};
