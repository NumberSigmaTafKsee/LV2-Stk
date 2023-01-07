gcc -std=c++17 -DUSE_FFTW -O2 -fPIC -march=native -mavx2 -c *.cpp
ar -rcv -o libGist.a *.o
cp *.a ../..

