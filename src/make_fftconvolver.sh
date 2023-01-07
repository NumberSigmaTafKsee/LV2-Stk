swig -lua -c++ -Iinclude src/fftconvolver.i
gcc -DAUDIOFFT_FFTW3 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fftconvolver.so src/fftconvolver_wrap.cxx libfftconvolver.a -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
