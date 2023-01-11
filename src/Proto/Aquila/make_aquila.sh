clear
echo 'compiling aquila'
swig -lua -c++ -Iaquila aquila.i
gcc -fmax-errors=1 -Iaquila -O2 -fPIC -shared -march=native -o aquila.so aquila_wrap.cxx fft4g.c -lstdc++ -L. -lAquila
