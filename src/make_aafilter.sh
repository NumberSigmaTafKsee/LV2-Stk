swig -lua -c++ aafilter.i
gcc -Iinclude -Isrc -O2 -fPIC -march=native -mavx2 -shared -o aafilter.so  aafilter_wrap.cxx -lstdc++ -lm -lluajit
