%module aquila
%{ 
#include "aquila.h"
#include <complex>
#include <vector>
using namespace Aquila;
%}

%include "std_math.i"
%include "std_vector.i"
%include "std_string.i"

// Get the STL typemaps
%include "stl.i"

// Handle standard exceptions
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
}

%inline %{    
    namespace std {
    typedef unsigned long int size_t;
    }
    typedef Aquila::SampleType SampleType;
    typedef double FrequencyType;
    typedef std::complex<double> ComplexType;
    typedef std::vector<ComplexType> SpectrumType;

    struct FFT
    {
        std::shared_ptr<Fft> _fft;

        FFT(size_t n) {
            _fft = Aquila::FftFactory::getFft(n);
        }
        Aquila::SpectrumType fft(Aquila::SampleType * x) {
            return _fft->fft(x);
        }
        void ifft(Aquila::SpectrumType spectrum, double * x) {
            _fft->ifft(spectrum,x);
        }
    };
%}

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(complex_type) std::complex<double>;
%template(spectrum_type) std::vector<std::complex<double>>;

%include "global.h"
%include "functions.h"


namespace Aquila 
{
    class SignalSource
    {
    public:
        
        SignalSource();
        SignalSource(FrequencyType sampleFrequency);        
        SignalSource(SampleType* data, std::size_t dataLength, FrequencyType sampleFrequency = 0);
        SignalSource(const std::vector<SampleType>& data, FrequencyType sampleFrequency = 0);
        virtual ~SignalSource();
        
        virtual FrequencyType getSampleFrequency() const;
        virtual void setSampleFrequency(FrequencyType sampleFrequency);
        virtual unsigned short getBitsPerSample() const;
        virtual std::size_t getSamplesCount() const;
        virtual SampleType sample(std::size_t position) const;
        virtual const SampleType* toArray() const;
        size_t length() const;
        
        
        %extend {
            SampleType __getitem(size_t i) { return $self->sample(i); }
            //void __setitem(size_t i, SampleType val) { (*$self)[i] = val; }
            SignalSource operator+(SampleType x)
            {
                return *$self + x;
            }
            SignalSource operator+(const SignalSource& rhs)
            {
                return *$self + rhs;
            }        
            
            

            SignalSource operator*(SampleType x)
            {
                return *$self * x;
            }
            SignalSource operator*(const SignalSource& rhs)
            {
                return *$self * rhs;
            }
            
        }
    protected:        
        std::vector<SampleType> m_data;
        FrequencyType m_sampleFrequency;
    };

    class Frame
    {
    public:
        Frame(const SignalSource& source, unsigned int indexBegin,
                unsigned int indexEnd);
        Frame(const Frame& other);        
        Frame& operator=(const Frame& other);

        virtual std::size_t getSamplesCount() const;
        virtual unsigned short getBitsPerSample() const;
        virtual SampleType sample(std::size_t position) const;
        virtual const SampleType* toArray() const;

        %extend {
            SampleType __getitem__(size_t i) { return $self->sample(i); }
            //void __setitem__(size_t i, SampleType val) { (*$self)[i] = val; }
        }

    private:
        const SignalSource* m_source;
        unsigned int m_begin, m_end;

    };

    class FramesCollection
    {
    public:

        FramesCollection();
        FramesCollection(const SignalSource& source,
                        unsigned int samplesPerFrame,
                        unsigned int samplesPerOverlap = 0);
        ~FramesCollection();

        static FramesCollection createFromDuration(const SignalSource& source,
                                                double frameDuration,
                                                double overlap = 0.0);

        void divideFrames(const SignalSource& source,
                        unsigned int samplesPerFrame,
                        unsigned int samplesPerOverlap = 0);
        void clear();

        
        size_t count() const;
        unsigned int getSamplesPerFrame() const;
        Frame frame(std::size_t index) const;
        
        %extend {
            Frame __getitem__(size_t i) { return $self->frame(i); }
            //void __setitem__(size_t i, Frame & f) { ($self)[i] = f; }
        }
    private:
        /**
         * Frames container.
         */
        Container m_frames;

        /**
         * Number of samples in each frame.
         */
        unsigned int m_samplesPerFrame;
    };
}


%include "source/PlainTextFile.h"
%include "source/RawPcmFile.h"
%include "source/WaveFile.h"
%include "source/WaveFileHandler.h"
%include "source/generator/Generator.h"
%include "source/generator/SineGenerator.h"
%include "source/generator/SquareGenerator.h"
%include "source/generator/TriangleGenerator.h"
%include "source/generator/PinkNoiseGenerator.h"
%include "source/generator/WhiteNoiseGenerator.h"
%include "source/window/BarlettWindow.h"
%include "source/window/BlackmanWindow.h"
%include "source/window/FlattopWindow.h"
%include "source/window/GaussianWindow.h"
%include "source/window/HammingWindow.h"
%include "source/window/HannWindow.h"
%include "source/window/RectangularWindow.h"
%include "transform/Fft.h"
%include "transform/Dft.h"
%include "transform/AquilaFft.h"
%include "transform/OouraFft.h"
%include "transform/FftFactory.h"
%include "transform/Dct.h"
%include "transform/Mfcc.h"
%include "transform/Spectrogram.h"
%include "filter/MelFilter.h"
%include "filter/MelFilterBank.h"
%include "ml/DtwPoint.h"
%include "ml/Dtw.h"
%include "ml.h"
%include "tools/TextPlot.h"
