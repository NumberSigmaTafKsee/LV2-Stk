swig -lua -c++ eig.i
gcc -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o eig.so  eig_wrap.cxx -lstdc++ -lm -lluajit
