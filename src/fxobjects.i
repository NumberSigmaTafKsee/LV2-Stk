%module fxobjects
%{
typedef float DspFloatType;
#include "FXObjects.hpp"
#include "FXOscillators.hpp"
#include "FXDelays.hpp"
#include "FXFilters.hpp"
#include "FXFirFilters.hpp"
#include "FXDistortion.hpp"
#include "FXDynamics.hpp"
#include "FXFFTSpectrum.hpp"
#include "FXPhaseVocoders.hpp"
#include "FXReverbs.hpp"
#include "FXSampleRate.hpp"
#include "FXWaveDigitalFilters.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%ignore makeWindow;
%include "FXObjects.hpp"
%include "FXOscillators.hpp"
%include "FXDelays.hpp"
%include "FXFilters.hpp"
%include "FXFirFilters.hpp"
%include "FXDistortion.hpp"
%include "FXDynamics.hpp"
%include "FXFFTSpectrum.hpp"
%include "FXObjects.hpp"
%include "FXPhaseVocoders.hpp"
%include "FXReverbs.hpp"
%include "FXSampleRate.hpp"
%include "FXWaveDigitalFilters.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;
