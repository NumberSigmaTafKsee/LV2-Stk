swig -lua -c++ -Iinclude -I/usr/local/cude/include -Iinclude/cuMatDNN cumini.i
gcc -fmax-errors=1 -I/usr/local/cuda/include -Iinclude -Iinclude/cuMatDNN -fPIC -O2 -march=native -mavx2 -shared -o minidnn.so mini_wrap.cxx -lstdc++ -lm -lluajit
