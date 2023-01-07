require('luapa')
require('sineosc')
require('Cranium')
require('LuaPolyBLEP')
require('luamoog')

NUM=1
sin = LuaPolyBLEP.PolyBLEP(44100,LuaPolyBLEP.PolyBLEP.SAWTOOTH)
network = Network_New(NUM, 1,{{size=16, activation=cranium.active_tanh}}, 1, cranium.active_tanh)

input={}
output = {}
for i=0,255 do
    data = {} 
    data[1] = i/255
    table.insert(input,data)    
    data = {} 
    data[1] = math.sin(2*math.pi*i/255)
    table.insert(output,data)
end 
    
example = DataSet_New(256, 1, input)
classes = DataSet_New(256, 1, output)
 

local params = cranium.ParameterSet()    
params.network = network.network
params.data     = example.ds
params.classes  = classes.ds
params.lossFunction = 1
params.batchSize    = 16
params.learningRate = 1e-1
params.searchTime   = 1000
params.regularizationStrength = 1e-6
params.momentumFactor = 0.9
params.maxIters = 1000
params.shuffle  = 0
params.verbose  = 1

f = 1
sin:setFrequency(1)
sinosc = SineOsc(44100,256)

filt = luamoog.StilsonMoog(44100)
cut = 440.0
filt:SetCutoff(cut)
filt:SetResonance(0.995)

function train()
    cranium.optimize(params)                                                                                                      
end 


train() 

data = cranium.float_buffer(1)
m = cranium.createMatrix(1,1,data)

--[[
for i=0,255 do 
    cranium.float_set(data,0,i/255)       
    cranium.forwardPass(network.network,m)
    local o = network:getOutput()
    local x = o:get_index(0)
    print(x,math.sin(2*math.pi*i/255))
end
os.exit()
]]

f =440
sin:setFrequency(f)
phase = 0
freq = 440/44100
inc = 0.002
q =  0
xn = 0 
yn = 0



function noise(input,output,frames)                
    for i = 0,frames-1 do                          
        cranium.forwardPass(network.network,m)
        local o = network:getOutput()
        local x = o:get_index(0)
        local r = 2*i            
                
        x = math.tanh(x/4)
                
        phase = phase + freq
        if( phase >= 1.0) then phase = 0.0 end 
        phase = 0.5*(1+sinosc:Generate(220))
        cranium.float_set(data,0,phase)
                        
        --[[
        phase = phase + freq
        if( phase >= 1.0) then phase=1.0 freq = -freq end 
        if( phase <= 0.0) then phase=0.0 freq = -freq end 
        q = 1 - math.abs(phase % 2) -1
        table.remove(data)        
        table.insert(data,1,q)
        ]]
        luapa.float_set(output,r,x)
        luapa.float_set(output,r+1,x)                                                   
    end     
end 

function sys()
    print(x)
end 

luapa.set_audio_func(noise)    
luapa.InitAudioDevice(10,-1,2,44100,2*8192)
luapa.RunAudio()
luapa.StopAudio()



