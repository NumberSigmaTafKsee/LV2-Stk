swig -lua -c++ mini.i
gcc -fmax-errors=1 -Iinclude/MiniDNN -fPIC -O2 -march=native -mavx2 -shared -o minidnn.so mini_wrap.cxx -lstdc++ -lm -lluajit
