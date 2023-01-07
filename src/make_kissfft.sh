swig -lua -c++ -Iinclude src/kiss_complex.i
swig -lua -c++ -Iinclude src/kiss_real.i
gcc -Iinclude -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared -o kiss_complex.so src/kiss_complex_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
gcc -Iinclude -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared -o kiss_real.so src/kiss_real_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
