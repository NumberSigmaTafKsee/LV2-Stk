%module chord
%{
#include "ChordDetector.h"
#include "Chromagram.h"
%}
%include "std_vector.i"
%include "ChordDetector.h"
%include "Chromagram.h"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
