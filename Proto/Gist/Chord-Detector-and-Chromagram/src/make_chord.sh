swig -lua -c++ chord.i
gcc  -DUSE_KISS_FFT -std=c++17 -I/usr/local/include/kissfft -O2 -fPIC -march=native -mavx2 -shared -o chord.so chord_wrap.cxx ChordDetector.cpp Chromagram.cpp -lstdc++ -lm -lluajit -lkissfft-float
