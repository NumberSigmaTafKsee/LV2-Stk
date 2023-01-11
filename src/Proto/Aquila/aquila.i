%module aquila
%{ 
#include "aquila.h"
using namespace Aquila;
%}

#define AQUILA_EXPORT
%inline %{
    typedef Aquila::SampleType SampleType;
%}

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
        

            SignalSource operator+(const SignalSource& lhs, SampleType x)
            {
                return lhs + x;
            }
            SignalSource operator+(SampleType x, const SignalSource& rhs)
            {
                return x + rhs;
            }
            SignalSource operator+(const SignalSource& lhs, const SignalSource& rhs)
            {
                return lhs + rhs;
            }

            SignalSource operator*(const SignalSource& lhs, SampleType x)
            {
                return lhs * x;
            }
            SignalSource operator*(SampleType x, const SignalSource& rhs)
            {
                return x * rhs;
            }
            SignalSource operator*(const SignalSource& lhs, const SignalSource& rhs)
            {
                return lhs*rhs;
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
            SampleType __getitem(size_t i) { return $self->sample(i); }
            //void __setitem(size_t i, SampleType val) { SampleType * m_data = toArray(); m_data[i] = val;}
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
            Frame __getitem(size_t i) { return $self->frame(i); }
            //void __setitem(size_t i, Frame & f) { m_frames[i] = f; }
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
