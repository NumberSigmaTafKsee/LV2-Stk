%module gist 
%{
#define USE_FFTW
#include <Gist.h>    
%}
#define USE_FFTW
%include <CoreFrequencyDomainFeatures.h>
%include <CoreTimeDomainFeatures.h>
%include <MFCC.h>
%include <OnsetDetectionFunction.h>
%include <WindowFunctions.h>
%include <Yin.h>
%include <Gist.h>


%template (GistFloat) Gist<float>;
%template (GistDouble) Gist<double>;

