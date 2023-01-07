--[[
1 Pole/1 zero
         b0           b1          a1
LP      K/(K+1)     K/(K+1)     (K-1)/(K+1)   
HP      1/(K+1)     -1/(K+1)    (K-1)/(K+1)
AP     (K-1)/(K+1)     1        (K-1)/(K+1)

wc = 2 * PI * fc/fs
K  = tan(wc/2)    
a  = 1+K

a0 = 1
a1 = -[(1-K)/a]

LP:
b0 = b1 = K/a

HP:
b0 = 1/a
b1 = -1/a

2nd order 
RBJ lol http://www.apogeebio.com/ddx/PDFs/AN-06.pdf
A = sqrt(10(G/20))
fc = 2*PI*fc/fs
wS = sin(fc)
wC = cos(fc)
a = wS/(2*Q)
B = sqrt(A)/Q (shelving only)

LP:
b0 = (1-Wc)/2
b1 = 1-Wc
b2 = (1-Wc)/2
a0 = 1 + a
a1 = -2*wc
a2 = 1 - a

HP:
b0 = (1+wC)/2
b1 = -(1+wC)
b2 = (1+wC)/2
a0 = 1 + a
a1 = -2*wC
a2 = 1-a

BP:
Peak gain = Q
b0 = wS/2 = Q*a
b1 = 0
b2 = -wS/2 = -Q*a
a0 = 1+a
a1 = -2*wC
a2 = 1-a

Constant 0-db peak BP
b0 = a
b1 = 0
b2 = -a
a0 = 1+a
a1 = -2*wC
a2 = 1-a

Notch:
b0=1
b1=-2*wC
b2 = 1
a0=1+a
a1=-2*wC
a2=1-a

All pass
b0 = 1-a
b1 = -2*wC
b2 = 1+a
a0 = 1+a
a1 = -2*wC
a2 = 1-a

Peak EQ
b0 = 1 + (a*A)
b1 = -2*wC
b2 = 1-(a*A)
a0 = 1 - (a/A)
a1 = -2*wC
a2 = 1-(a/A)

Low Shelf
b0 = A * {(A+1)-[(A-1)*wC]+(b*wS)}
b1 = 2*A*{(A-1)-[(A+1)*wC]}
b2 = A * {(A+1)-[(A-1)*wC]-(b*wS)}
a0 = {(A+1)+[(A-1)*wC]-(b*wS)}
a1 = -2*{(A-1)+[(A+1)*wC]
a2 = {(A+1)+[(A-1)*wC]-(b*wS)}

High SHelf
b0 = A * {(A+1)+[(A-1)*wC] + (b*wS)}
b1 = -2*A*{(A-1)+[(A-1)*wC]}
b2 = A * {(A+1)+[(A-1)*wC]-(b*wS)}
a0 = {(A+1)-[(A-1)*wC]}
a1 = 2*{(A-1)-[(A+1)*wC]}
a2 = {(A+1)-[(A-1)*wc]-(b*ws)}
]]

require('audiosystem')
require('soundwave')
require('biquad')

filt = { 
    x0 = 0;
    x1 = 0;
    y0 = 0;
    y1 = 0;
}

function filter1z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 
    filt.x0 = x     
    return y 
end 

function filter2z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 + B[3]*filt.x1
    filt.x0 = filt.x1     
    filt.x1 = x
    return y 
end 

function filter1p(B,A,x)
    local y = B[1]*x - A[1]*filt.y0
    filt.y0 = y    
    return y 
end 

function filter2p(B,A,x)
    local y = B[1]*x - A[1]*filt.y0 - A[2]*filt.y1
    filt.y0 = filt.y1
    filt.y1 = y
    return y 
end 

function filter1p1z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 - A[1]*filt.y0     
    filt.x0 = x 
    filt.y0 = y  
    return y 
end 

function filter1p2z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 -B[3]*filt.x1 - A[1]*filt.y0 
    filt.x0 = file.x1
    file.x1 = x 
    filt.y0 = y
    return y 
end 

function filter2p1z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 - A[1]*filt.y0 - A[2]*filt.y1
    filt.x0 = file.x1
    file.x1 = x 
    filt.y0 = filt.y1
    filt.y1 = y
    return y 
end 


function filter2p2z(B,A,x)
    local y = B[1]*x + B[2]*filt.x0 + B[3]*filt.x1 - A[1]*filt.y0 - A[2]*filt.y1
    filt.x0 = filt.x1
    filt.x1 = x
    filt.y0 = filt.y1
    filt.y1 = y
    return y 
end 

function one_pole_lp(fc)
    -- 1p lp/hp
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {-((1-K)/alpha)}
    local B = {K/alpha,K/alpha}
    return A,B
end

function one_pole_hp(fc)
    -- 1p lp/hp
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {-((1-K)/alpha)}    
    local B = {1/alpha,-1/alpha}
    return A,B
end

function one_zero_lp(fc)
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {1,-((1-K)/alpha)}
    local B= {K/alpha,0}
    return A,B
end

function one_pole1z_lp(fc)
    -- 1p lp/hp
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {1,-((1-K)/alpha)}
    local B = {K/alpha,K/alpha}
    return A,B
end

function one_pole1z_hp(fc)
    -- 1p lp/hp
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {1,-((1-K)/alpha)}    
    local B = {1/alpha,-1/alpha}
    return A,B
end


-- 1p/1z
function one_pole_one_zero_lp(fc)
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {-((1-K)/alpha)}
    local B = {K/alpha,K/alpha}
    return A,B
end


function one_pole_one_zero_hp(fc)
    local wc = 2*math.pi*fc/44100
    local K = math.tan(wc/2)
    local alpha = 1 + K
    local A= {1,-((1-K)/alpha)}
    local B = {1/alpha,-1/alpha}
    return A,B
end

function one_pole_reiss_lp(fc)
    wc = 2*math.pi*40/44100
    K = math.tan(wc/2)
    alpha = 1 + K
    A = {(K-1)/(K+1)}
    B = {K/(K-1),K/(K-1)}
    return A,B
end 

function resonant_2pole_2zero_lp(cutoff,resonanceDB)
    local c, csq, resonance, q, a0, a1, a2, b1, b2;
    c = 1.0 / (math.tan(math.pi * (cutoff / 44100)));
    csq = c * c;
    resonance = math.pow(10.0, -(resonancedB * 0.5));
    q = math.sqrt(2.0) * resonance;
    a0 = 1.0 / (1.0 + (q * c) + (csq));
    a1 = 2.0 * a0;
    a2 = a0;
    b1 = (2.0 * a0) * (1.0 - csq);
    b2 = a0 * (1.0 - (q * c) + csq);
    local A = {b1,b2}
    local B = {a0,a1,a2}
    return A,B
end

-- LP
--//recursion:   tmp = (1-p)*in + p*tmp with output = tmp
--//coefficient: p = (2-cos(x)) - sqrt((2-cos(x))^2 - 1) with x = 2*pi*cutoff/samplerate
--//coeficient approximation: p = (1 - 2*cutoff/samplerate)^2

--p = (2-math.cos(wc)) - math.sqrt((2-math.cos(wc))^2 - 1)
--A = {p}
--B = {1-p,0}

-- HP
-- recursion: tmp = (p-1)*in - p*tmp with output = tmp
-- coefficient: p = (2+cos(x)) - sqrt((2+cos(x))^2 - 1) with x = 2*pi*cutoff/samplerate
-- coeficient approximation: p = (2*cutoff/samplerate)^2
--p = (2+math.cos(wc)) - math.sqrt((2-math.cos(wc))^2 - 1)
--A = {-p}
--B = {p-1,0}

--[[
-- 1p/1z low
w = 2.0 * 44100
fCut = 440.0*2*math.pi
Norm = 1.0/(fCut+1)
b1   = (w-fCut)*Norm
a0   = (w-fCut)*Norm
a1   = a0

A = {b1}
B = {a0,a1}
]]

--[[
-- 1p/1z high
w = 2.0*44100
fCut = 440.0 * 2.0 * math.pi
Norm = 1.0 / (fCut+w)
a0 = w * Norm
a1 = -a0
b1 = (w - fCut) * Norm
]]

--[[
Q = 0.5
cosW = math.cos(wc)
sinW = math.sin(wc)
alpha = sinW/(2*Q)
b0 = (1 - cosW)/2;
b1 = 1 - cosW;
b2 = (1 - cosW)/2;
a0 = 1 + alpha;
a1 = -2*cosW;
a2 = 1 - alpha;
biquad = biquad.BiQuad(b0/a0,b1/a0,b2/a0,a1/a0,a2/a0)

]]

--A,B = one_pole_lp(440)
--A,B = one_zero_lp(440)
--A,B = one_pole1z_lp(440)
--A,B = one_zero_lp(40)
--A,B = resonant_2pole_2zero_lp(440)
K = math.tan(math.pi* 220 / 44100)

A = {1}
B = {K/(K+1),K/(K+1)}

sine = soundwave.BandlimitedOscillator(44100,soundwave.BandlimitedOscillator.SQUARE)

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function noise(input,output,frames)        
    local buffer = new_buffer(output)
    local mono   = soundwave.float_vector(frames)        
    sine:Process(frames,input,mono:data())   
    for i=1,mono:size() do
        --local y = filt:step(mono[i])
        local y = filter1z(B,A,mono[i])
        buffer[i-1] = y
    end
    
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end
function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0            
    freq = f
    sine:setFrequency(f)
    filt:SetCutoff(f + f/2)
    adsr:gate(1)    
end
function note_off(c,n,v)    
    adsr:gate(0)
end

-- app callback, midi handling and logic
-- isAudioRunning shuld be atomic
-- either block audio or wait until finished
-- this is run every 10ms, or can be changed in portaudio.i
function callback()
    print('hi')
end 

audiosystem.Init()
audiosystem.Pm_Initialize()

audiosystem.set_note_on_func(note_on)
audiosystem.set_note_off_func(note_off)

for i=0,audiosystem.GetNumMidiDevices()-1 do 
    print(i,audiosystem.GetMidiDeviceName(i))
end

audiosystem.set_audio_func(noise)
--audiosystem.set_callback_func(callback)
device=14
audiosystem.Pa_Initialize()
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    if( audiosystem.GetAudioDeviceName(i) == 'pulse') then        
        device = i 
        goto done
    end    
end
::done::
audiosystem.InitMidiDevice(1,3,3)
audiosystem.InitAudioDevice(device,-1,1,44100,256)
audiosystem.RunAudio()
audiosystem.StopAudio()
