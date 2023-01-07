require('aquila')

size=64
sampleFreq = 2000
f1 = 125
f2 = 700
sine1 = aquila.SineGenerator(sampleFreq)
sine1:setAmplitude(32):setFrequency(f1):setPhase(0.75):generate(size)
sine2 = aquila.SineGenerator(sampleFreq)
sine2:setAmplitude(8):setFrequency(f2):setPhase(0.75):generate(size)
x = sine1+sine2
plot = aquila.TextPlot("Input Signal")
-- unfortunately Aquila returns std::shared_ptr which is not simple to wrap
fft = aquila.FFT(size)
spectrum = fft:fft(x:toArray())
plot:setTitle("Spectrum")
plot:plotSpectrum(spectrum)