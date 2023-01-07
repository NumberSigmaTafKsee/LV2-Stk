swig -Ikfrcore/include/Kfr -lua -c++ kfr_core.i
gcc -std=c++17 -Ikfrcore/include/Kfr -O2 -march=native -mavx2 -fPIC -shared -o kfr.so kfr_core_wrap.cxx -lstdc++ -lm -lluajit -lkfr_dft -lkfr_io
