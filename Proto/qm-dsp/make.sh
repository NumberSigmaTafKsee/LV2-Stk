swig -lua -c++ -Iinclude qmdsp.i
gcc -I. -fPIC -O2 -march=native -shared -oqmdsp.so qmdsp_wrap.cxx -lstdc++ -lluajit-5.1 -L. -lqm-dsp -llapack -lblas
