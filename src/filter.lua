require('soundwave')
require('audiosystem')
require('Filter')

osc = soundwave.BandlimitedOscillator(44100,soundwave.SAWTOOTH)
osc:setFrequency(440)
filt = Filter.DinkyFilter()
filt:setCutoff(0.25)
filt:setResonance(.09)

function new_buffer(p)
    local b = {}
    b.buffer = p

    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function sound(input,output,frames)
    local out = new_buffer(output)
    for i=0,frames-1 do 
        local sample = osc:tick()
        sample = filt:process(sample)
        out[i*2] = sample
        out[i*2+1] = sample        
    end 
end

function cb_run()
	print('hi')
end

audiosystem.set_audio_func(sound)
audiosystem.Pa_Initialize()
device = 14
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    local dev = audiosystem.GetAudioDeviceName(i)
    if(dev == 'pulse') then device = i end     
end

audiosystem.InitAudioDevice(device,-1,2,44100,256)
audiosystem.set_callback_func(cb_run)
audiosystem.RunAudio()
audiosystem.Sleep(5000)
audiosystem.StopAudio()
