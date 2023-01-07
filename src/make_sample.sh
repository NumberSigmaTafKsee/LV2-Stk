swig -lua -c++ -Iinclude/samples samples.i
gcc -std=c++17 -Iinclude -Iinclude/samples -O2 -fPIC -shared -o samples.so samples_wrap.cxx libaudiofft.a -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
