swig -IGist/src -lua -c++ gist.i
gcc -std=c++17 -DUSE_FFTW -IGist/src -O2 -fPIC -march=native -mavx2 -shared -o gist.so gist_wrap.cxx libGist.a -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
