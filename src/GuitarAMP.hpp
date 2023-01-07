#include <ATK/EQ/ButterworthFilter.h>
#include <ATK/EQ/IIRFilter.h>
#include <ATK/EQ/ToneStackFilter.h>
#include <ATK/Preamplifier/Triode3Filter.h>
#include <ATK/Special/ConvolutionFilter.h>
#include <ATK/Tools/DecimationFilter.h>
#include <ATK/Tools/OversamplingFilter.h>
#include <ATK/Tools/SumFilter.h>
#include <ATK/Tools/VolumeFilter.h>
#include <ATK/Core/TypedBaseFilter.h>

template <typename DataType>
class ATKInputFilter : public ATK::TypedBaseFilter<DataType> {
  public:
    using ATK::TypedBaseFilter<DataType>::outputs;

    explicit ATKInputFilter(int channels) : ATK::TypedBaseFilter<DataType>(0, channels), mChannels(channels) {}

    void set_inputs(double** inputs, int size) {
        mInputs = inputs;
        mSize = size;
    }

  protected:
    double** mInputs = nullptr;
    int mSize = 0;
    int mChannels = 1;

    virtual void process_impl(int64_t size) const {
        for (int c = 0; c < mChannels; ++c) {
            for (int64_t i = 0; i < size; ++i) {
                outputs[c][i] = mInputs[c][i];
            }
        }
    }
};

template class ATKInputFilter<float>;
template class ATKInputFilter<double>;




template <typename DataType>
class ATKOutputFilter : public ATK::TypedBaseFilter<DataType> {
  public:
    using ATK::TypedBaseFilter<DataType>::converted_inputs;

    explicit ATKOutputFilter(int channels) : ATK::TypedBaseFilter<DataType>(channels, 0), mChannels(channels) {}

    void set_outputs(double** outputs, int size) {
        mOutputs = outputs;
        mSize = size;
    }

  protected:
    double** mOutputs = nullptr;
    int mSize = 0;
    int mChannels = 1;

    virtual void process_impl(int64_t size) const {
        for (int c = 0; c < mChannels; ++c) {
            for (int64_t i = 0; i < size; ++i) {
                mOutputs[c][i] = converted_inputs[c][i];
            }
        }
    }
};

template class ATKOutputFilter<float>;
template class ATKOutputFilter<double>;


using InType = ATKInputFilter<double>;
using OutType = ATKOutputFilter<double>;
using BaseFilterType = ATK::TypedBaseFilter<double>;
using SumType = ATK::SumFilter<double>;
using VolumeType = ATK::VolumeFilter<double>;
using OversamplingType_2 = ATK::OversamplingFilter<double, ATK::Oversampling6points5order_2<double>>;
using OversamplingType_4 = ATK::OversamplingFilter<double, ATK::Oversampling6points5order_4<double>>;
using OversamplingType_8 = ATK::OversamplingFilter<double, ATK::Oversampling6points5order_8<double>>;
using DownsampingType = ATK::DecimationFilter<double>;
using LowPassType = ATK::IIRFilter<ATK::ButterworthLowPassCoefficients<double>>;
using TubeFilterType = ATK::Triode3Filter<double>;
using ToneStackType = ATK::IIRFilter<ATK::ToneStackCoefficients<double>>;
using CabType = ATK::ConvolutionFilter<double>;
using PeakFilterType = PeakMeterFilter<double>;

int mOversamplingFactor = 4;
ATK::Triode3Type mTubeType = ATK::Triode3Type::ECC83;

std::unique_ptr<InType> mInputs;
std::unique_ptr<OutType> mOutputs;
std::unique_ptr<SumType> mSum;
std::unique_ptr<VolumeType> mGain;
std::unique_ptr<VolumeType> mVolume;
std::unique_ptr<BaseFilterType> mOversample;
std::unique_ptr<DownsampingType> mDownsample;
std::unique_ptr<LowPassType> mLowPass;
std::unique_ptr<TubeFilterType> mTube[NUM_TUBES];
std::unique_ptr<ToneStackType> mToneStack;
std::unique_ptr<CabType> mCabinet;
std::unique_ptr<PeakFilterType> mPeakOut;



void ProcessDoubleReplacing( int nFrames, double** inputs, double** outputs) {
    // Mutex is already locked for us.
    mInputs->set_inputs(inputs, nFrames);
    mOutputs->set_outputs(outputs, nFrames);
    mOutputs->process(nFrames);
}


void Setup() {
    int channels = 2;
    int rate = sampleRate;
    int rateOver = rate * mOversamplingFactor;

    // if (IR_SAMPLERATE != rate) {
    //    mResampler.reset(new r8b::CDSPResampler16IR(IR_SAMPLERATE, rate, RESAMPLE_BLOCK_SIZE));
    //}

    mInputs.reset(new InType(channels));
    mInputs->set_output_sampling_rate(rate);

    if (!mBypass) {
        switch (mOversamplingFactor) {
            case 2:
                mOversample.reset(new OversamplingType_2(1));
                break;
            case 4:
                mOversample.reset(new OversamplingType_4(1));
                break;
            case 8:
                mOversample.reset(new OversamplingType_8(1));
                break;
        }
        mOversample->set_input_sampling_rate(rate);
        mOversample->set_output_sampling_rate(rateOver);

        if (channels > 1) {
            mSum.reset(new SumType());
            mSum->set_input_sampling_rate(rate);
            mSum->set_input_port(0, mInputs.get(), 0);
            mSum->set_input_port(1, mInputs.get(), 1);
            mOversample->set_input_port(0, mSum.get(), 0);
        } else {
            mOversample->set_input_port(0, mInputs.get(), 0);
        }

        mLowPass.reset(new LowPassType(1));
        mLowPass->set_cut_frequency(rate);
        mLowPass->set_order(5);
        mLowPass->set_input_sampling_rate(rateOver);

        mDownsample.reset(new DownsampingType(1));
        mDownsample->set_input_sampling_rate(rateOver);
        mDownsample->set_output_sampling_rate(rate);
        mDownsample->set_input_port(0, mLowPass.get(), 0);

        mToneStack.reset(new ToneStackType(ToneStackType::buildJCM800Stack()));
        mToneStack->set_input_sampling_rate(rate);
        mToneStack->set_low(GetParam(mLowID)->Value());
        mToneStack->set_middle(GetParam(mMidID)->Value());
        mToneStack->set_high(GetParam(mHighID)->Value());
        mToneStack->set_input_port(0, mDownsample.get(), 0);

        mVolume.reset(new VolumeType(1));
        mVolume->set_input_sampling_rate(rate);
        mVolume->set_volume(GetParam(mVolID)->Value());
    }

    mPeakOut.reset(new PeakFilterType(1, mMeterOut));
    mPeakOut->set_input_sampling_rate(rate);
    if (!mBypass) {
        mPeakOut->set_input_port(0, mVolume.get(), 0);
    } else {
        mPeakOut->set_input_port(0, mInputs.get(), 0);
    }

    mOutputs.reset(new OutType(channels));
    mOutputs->set_input_sampling_rate(rate);
    mOutputs->set_input_port(0, mPeakOut.get(), 0);
    if (channels > 1) {
        mOutputs->set_input_port(1, mPeakOut.get(), 0);
    }

    if (!mBypass) {
        SetupTubes();

        // Add Gain control between the first and second tube stages
        mGain.reset(new VolumeType(1));
        mGain->set_input_sampling_rate(rateOver);
        mGain->set_volume(GetParam(mGainID)->Value());
        mGain->set_input_port(0, mTube[0].get(), 0);
        mTube[1]->set_input_port(0, mGain.get(), 0);

        SetupCabinet();
    }
}

void SetupTubes(bool useTubeParameters) {
    int rate = sampleRate;
    int rateOver = rate * mOversamplingFactor;
    for (int i = 0; i < NUM_TUBES; ++i) {
        mTube[i].reset(new TubeFilterType(mTubeType));        
        mTube[i]->set_input_sampling_rate(rateOver);
        if (i > 1) {
            // We connect the gain stage between tube 1 and 2, so we skip this here
            mTube[i]->connect_stage(mTube[i - 1].get());
        }
    }
    mTube[0]->set_input_port(0, mOversample.get(), 0);
    mLowPass->set_input_port(0, mTube[NUM_TUBES - 1].get(), 0);
}

void SetupCabinet() {
    int rate = GetSampleRate();
    mCabinet.reset(new CabType());
    mCabinet->set_split_size(CONVO_SPLIT_SIZE);
    mCabinet->set_input_sampling_rate(rate);
    Impulse& impulse(mCab < gImpulses.size() ? gImpulses[mCab] : *mCustIR);
    if (impulse.SampleRate != rate) {
        std::vector<double> ir;
        Resample(impulse.Frames, ir, impulse.SampleRate, rate);
        mCabinet->set_impulse(std::move(ir));
    } else {
        std::vector<double> ir(impulse.Frames);
        mCabinet->set_impulse(std::move(ir));
    }
    mCabinet->set_input_port(0, mToneStack.get(), 0);
    mVolume->set_input_port(0, mCabinet.get(), 0);
}