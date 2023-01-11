swig -lua -c++ btrack.i
gcc -O2 -DUSE_KISS_FFT -I/usr/local/include/kissfft -fPIC -march=native -mavx2 -shared -o btrack.so btrack_wrap.cxx BTrack.cpp OnsetDetectionFunction.cpp libbtrack.a -lstdc++ -lm -lfftw3 -lluajit -lkissfft-float -lsndfile -lsamplerate

