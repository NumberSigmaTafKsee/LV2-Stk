%module btrack
%{
#include "BTrack.h"
#include "OnsetDetectionFunction.h"
%}

%include "std_vector.i"

%include "BTrack.h"
#include "OnsetDetectionFunction.h"

%template (float_vector) std::vector<float>;
%template (double_vector) std::vector<double>;

%inline %{
std::vector<double> convert(std::vector<float> x)
{
	std::vector<double> v(x.size());
	for(size_t i = 0; i < x.size(); i++) v[i] = x[i];
	return v;
}
std::vector<float> convert(std::vector<double> x)
{
	std::vector<float> v(x.size());
	for(size_t i = 0; i < x.size(); i++) v[i] = x[i];
	return v;
}
%}

