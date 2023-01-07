swig -lua -c++ -Iinclude/uKfr -I/usr/local/include/kfr kfr.i
gcc -std=c++17 -fmax-errors=1  -Iinclude/uKfr -O2 -fPIC -march=native -mavx2 -shared -o kfr.so kfr_wrap.cxx -lstdc++ -lm -lluajit -lkfr_dft -lkfr_io 
