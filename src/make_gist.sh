swig -lua -c++ -Iinclude/Gist gist.i
gcc -I. -Iinclude/Gist -O2 -fPIC -march=native -mavx2 -shared -o gist.so gist_wrap.cxx lib/libGist.a -lstdc++ -lm -lluajit -L. -lfftw3
