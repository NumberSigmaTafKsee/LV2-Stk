#pragma once
#include <cstdint>
#include <cmath>
#include <cassert>
#include <random>
#include <cstdlib>

/** PIs
*/
#define PI_F 3.1415927410125732421875f
#define TWOPI_F (2.0f * PI_F)
#define HALFPI_F (PI_F * 0.5f)
#define DSY_MIN(in, mn) (in < mn ? in : mn)
#define DSY_MAX(in, mx) (in > mx ? in : mx)
#define DSY_CLAMP(in, mn, mx) (DSY_MIN(DSY_MAX(in, mn), mx))
#define DSY_COUNTOF(_arr) (sizeof(_arr) / sizeof(_arr[0]))
#define ROOT2 (1.4142135623730950488)

#define DSY_REVERBSC_MAX_SIZE 98936
#define REVSC_OK 0
#define REVSC_NOT_OK 1

#define DEFAULT_SRATE 48000.0
#define MIN_SRATE 5000.0
#define MAX_SRATE 1000000.0
#define MAX_PITCHMOD 20.0
#define DELAYPOS_SHIFT 28
#define DELAYPOS_SCALE 0x10000000
#define DELAYPOS_MASK 0x0FFFFFFF


#define EXPF expf
    
#ifndef FT_MAXLEN
#define FT_MAXLEN 0x1000000L
#endif

#ifndef FT_PHMASK
#define FT_PHMASK 0x0FFFFFFL
#endif

#define DSY_NLFILT_MAX_DELAY 1024
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifndef max
#define max(a, b) ((a < b) ? b : a)
#endif

#ifndef min
#define min(a, b) ((a < b) ? a : b)
#endif

#define SLOPE(out, in, positive, negative)                \
    {                                                     \
        float error = (in)-out;                           \
        out += (error > 0 ? positive : negative) * error; \
    }


/** Shift can be 30-100 ms lets just start with 50 for now.
0.050 * SR = 2400 samples (at 48kHz)
*/
#define SHIFT_BUFFER_SIZE 16384
//#define SHIFT_BUFFER_SIZE 4800
//#define SHIFT_BUFFER_SIZE 8192
//#define SHIFT_BUFFER_SIZE 1024

#define WUTR_SOUND_DECAY 0.95f
#define WUTR_SYSTEM_DECAY 0.996f
#define WUTR_GAIN 1.0f
#define WUTR_NUM_SOURCES 10.0f
#define WUTR_CENTER_FREQ0 450.0f
#define WUTR_CENTER_FREQ1 600.0f
#define WUTR_CENTER_FREQ2 750.0f
#define WUTR_RESON 0.9985f
#define WUTR_FREQ_SWEEP 1.0001f
#define MAX_SHAKE 2000.0f

#define PLUKMIN 64
#define INTERPFACTOR 256.0f
#define INTERPFACTOR_I 255


namespace Daisy
{
    //////////////////////////////////////////////////////////////////////////////////////
    // DSP
    //////////////////////////////////////////////////////////////////////////////////////

    //Avoids division for random floats. e.g. rand() * kRandFrac
    static constexpr float kRandFrac = 1.f / (float)RAND_MAX;

    //Convert from semitones to other units. e.g. 2 ^ (kOneTwelfth * x)
    static constexpr float kOneTwelfth = 1.f / 12.f;

    /** efficient floating point min/max
    c/o stephen mccaul
    */
    inline float fmax(float a, float b)
    {
        float r;
    #ifdef __arm__
        asm("vmaxnm.f32 %[d], %[n], %[m]" : [d] "=t"(r) : [n] "t"(a), [m] "t"(b) :);
    #else
        r = (a > b) ? a : b;
    #endif // __arm__
        return r;
    }

    inline float fmin(float a, float b)
    {
        float r;
    #ifdef __arm__
        asm("vminnm.f32 %[d], %[n], %[m]" : [d] "=t"(r) : [n] "t"(a), [m] "t"(b) :);
    #else
        r = (a < b) ? a : b;
    #endif // __arm__
        return r;
    }

    /** quick fp clamp
    */
    inline float fclamp(float in, float min, float max)
    {
        return fmin(fmax(in, min), max);
    }

    /** From Musicdsp.org "Fast power and root estimates for 32bit floats)
    Original code by Stefan Stenzel
    These are approximations
    */
    inline float fastpower(float f, int n)
    {
        long *lp, l;
        lp = (long *)(&f);
        l  = *lp;
        l -= 0x3F800000;
        l <<= (n - 1);
        l += 0x3F800000;
        *lp = l;
        return f;
    }

    inline float fastroot(float f, int n)
    {
        long *lp, l;
        lp = (long *)(&f);
        l  = *lp;
        l -= 0x3F800000;
        l >>= (n = 1);
        l += 0x3F800000;
        *lp = l;
        return f;
    }

    /** From http://openaudio.blogspot.com/2017/02/faster-log10-and-pow.html
    No approximation, pow10f(x) gives a 90% speed increase over powf(10.f, x)
    */
    inline float pow10f(float f)
    {
        return expf(2.302585092994046f * f);
    }

    /* Original code for fastlog2f by Dr. Paul Beckmann from the ARM community forum, adapted from the CMSIS-DSP library
    About 25% performance increase over std::log10f
    */
    inline float fastlog2f(float f)
    {
        float frac;
        int   exp;
        frac = frexpf(fabsf(f), &exp);
        f    = 1.23149591368684f;
        f *= frac;
        f += -4.11852516267426f;
        f *= frac;
        f += 6.02197014179219f;
        f *= frac;
        f += -3.13396450166353f;
        f += exp;
        return (f);
    }

    inline float fastlog10f(float f)
    {
        return fastlog2f(f) * 0.3010299956639812f;
    }

    /** Midi to frequency helper
    */
    inline float mtof(float m)
    {
        return powf(2, (m - 69.0f) / 12.0f) * 440.0f;
    }


    /** one pole lpf
    out is passed by reference, and must be retained between
    calls to properly filter the signal
    coeff can be calculated:
    coeff = 1.0 / (time * sample_rate) ; where time is in seconds
    */
    inline void fonepole(float &out, float in, float coeff)
    {
        out += coeff * (in - out);
    }

    /** Curves to use with the fmap function */
    enum class Mapping
    {
        LINEAR,
        EXP,
        LOG,
    };

    /** Maps a float between a specified range, using a specified curve. 
     * 
     *  \param in a value between 0 to 1 that will be mapped to the new range.
     *  \param min the new minimum value
     *  \param max the new maxmimum value
     *  \param curve a Mapping Value to adjust the response curve of the transformation
     *               defaults to Linear. @see Mapping
     * 
     *  When using the log curve min and max, must be greater than zero.
     * 
     *  \retval returns the transformed float within the new range.
    */
    inline float
    fmap(float in, float min, float max, Mapping curve = Mapping::LINEAR)
    {
        switch(curve)
        {
            case Mapping::EXP:
                return fclamp(min + (in * in) * (max - min), min, max);
            case Mapping::LOG:
            {
                const float a = 1.f / log10f(max / min);
                return fclamp(min * powf(10, in / a), min, max);
            }
            case Mapping::LINEAR:
            default: return fclamp(min + in * (max - min), min, max);
        }
    }

    /** Simple 3-point median filter
    c/o stephen mccaul
    */
    template <typename T>
    T median(T a, T b, T c)
    {
        return (b < a) ? (b < c) ? (c < a) ? c : a : b
                    : (a < c) ? (c < b) ? c : b : a;
    }

    /** Ported from pichenettes/eurorack/plaits/dsp/oscillator/oscillator.h
    */
    inline float ThisBlepSample(float t)
    {
        return 0.5f * t * t;
    }

    /** Ported from pichenettes/eurorack/plaits/dsp/oscillator/oscillator.h
    */
    inline float NextBlepSample(float t)
    {
        t = 1.0f - t;
        return -0.5f * t * t;
    }

    /** Ported from pichenettes/eurorack/plaits/dsp/oscillator/oscillator.h
    */
    inline float NextIntegratedBlepSample(float t)
    {
        const float t1 = 0.5f * t;
        const float t2 = t1 * t1;
        const float t4 = t2 * t2;
        return 0.1875f - t1 + 1.5f * t2 - t4;
    }

    /** Ported from pichenettes/eurorack/plaits/dsp/oscillator/oscillator.h
    */
    inline float ThisIntegratedBlepSample(float t)
    {
        return NextIntegratedBlepSample(1.0f - t);
    }

    /** Soft Limiting function ported extracted from pichenettes/stmlib */
    inline float SoftLimit(float x)
    {
        return x * (27.f + x * x) / (27.f + 9.f * x * x);
    }

    /** Soft Clipping function extracted from pichenettes/stmlib */
    inline float SoftClip(float x)
    {
        if(x < -3.0f)
            return -1.0f;
        else if(x > 3.0f)
            return 1.0f;
        else
            return SoftLimit(x);
    }

    /** Quick check for Invalid float values (NaN, Inf, out of range) 
     ** \param x value passed by reference, replaced by y if invalid. 
    ** \param y value to replace x if invalidity is found. 
    ** 
    ** When DEBUG is true in the build, this will halt 
    ** execution for tracing the reason for the invalidity. */
    inline void TestFloat(float &x, float y = 0.f)
    {
        if(!std::isnormal(x) && x != 0)
        {
    #if defined(__arm__) && defined(DEBUG)
            asm("bkpt 255");
    #else
            x = y;
    #endif
        }
    }

    /** Based on soft saturate from:
    [musicdsp.org](musicdsp.org/en/latest/Effects/42-soft-saturation.html)
    Bram de Jong (2002-01-17)
    This still needs to be tested/fixed. Definitely does some weird stuff
    described as:
    x < a:
        f(x) = x
    x > a:
        f(x) = a + (x-a)/(1+((x-a)/(1-a))^2)
    x > 1:
        f(x) = (a + 1)/2
    */
    inline float soft_saturate(float in, float thresh)
    {
        bool  flip;
        float val, out;
        //val = fabsf(in);
        out  = 0.f;
        flip = in < 0.0f;
        val  = flip ? -in : in;
        if(val < thresh)
        {
            out = in;
        }
        else if(val > 1.0f)
        {
            out = (thresh + 1.0f) / 2.0f;
            if(flip)
                out *= -1.0f;
        }
        else if(val > thresh)
        {
            float temp;
            temp = (val - thresh) / (1 - thresh);
            out  = thresh + (val - thresh) / (1.0f + (temp * temp));
            if(flip)
                out *= -1.0f;
        }
        return out;
        //    return val < thresh
        //               ? val
        //               : val > 1.0f
        //                     ? (thresh + 1.0f) / 2.0f
        //                     : thresh
        //                           + (val - thresh)
        //                                 / (1.0f
        //                                    + (((val - thresh) / (1.0f - thresh))
        //                                       * ((val - thresh) / (1.0f - thresh))));
    }
    constexpr bool is_power2(uint32_t x)
    {
        return ((x - 1) & x) == 0;
    }


    
    // This causes with infinity with certain curves,
    // which then causes NaN erros...
    //#define EXPF expf_fast

    // To resolve annoying bugs when using this you can:
    // if (val != val)
    //     val = 0.0f; // This will un-NaN the value.

    // Fast Exp approximation
    // 8x multiply version
    //inline float expf_fast(float x)
    //{
    //   x = 1.0f + x / 256.0f;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   x *= x;
    //   return x;
    //}

    // 10x multiply version
    inline float expf_fast(float x)
    {
        x = 1.0f + x / 1024.0f;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        return x;
    }

    /** Prior to C++14 constexpr functions were required to be a single return statement.
     *  So this clause guards against that behavior to allow the library, and this function
     *  to continue to work with C++11.
     *  The function itself is not currently (12 May 2021) used within the library itself.
     */
    #if __cplusplus <= 201103L
    inline uint32_t get_next_power2(uint32_t x)
    #else
    constexpr uint32_t get_next_power2(uint32_t x)
    #endif
    {
        x--;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        x++;

        assert(is_power2(x));
        return x;
    }



    //////////////////////////////////////////////////////////////////////////////////////
    /** Simple Delay line.
    November 2019

    Converted to Template December 2019

    declaration example: (1 second of floats)

    DelayLine<float, SAMPLE_RATE> del;

    By: shensley
    */
    //////////////////////////////////////////////////////////////////////////////////////
  
    template <typename T, size_t max_size>
    class DelayLine
    {
    public:
        DelayLine() {}
        ~DelayLine() {}
        /** initializes the delay line by clearing the values within, and setting delay to 1 sample.
        */
        void Init() { Reset(); }
        /** clears buffer, sets write ptr to 0, and delay to 1 sample.
        */
        void Reset()
        {
            for(size_t i = 0; i < max_size; i++)
            {
                line_[i] = T(0);
            }
            write_ptr_ = 0;
            delay_     = 1;
        }

        /** sets the delay time in samples
            If a float is passed in, a fractional component will be calculated for interpolating the delay line.
        */
        inline void SetDelay(size_t delay)
        {
            frac_  = 0.0f;
            delay_ = delay < max_size ? delay : max_size - 1;
        }

        /** sets the delay time in samples
            If a float is passed in, a fractional component will be calculated for interpolating the delay line.
        */
        inline void SetDelay(float delay)
        {
            int32_t int_delay = static_cast<int32_t>(delay);
            frac_             = delay - static_cast<float>(int_delay);
            delay_ = static_cast<size_t>(int_delay) < max_size ? int_delay
                                                            : max_size - 1;
        }

        /** writes the sample of type T to the delay line, and advances the write ptr
        */
        inline void Write(const T sample)
        {
            line_[write_ptr_] = sample;
            write_ptr_        = (write_ptr_ - 1 + max_size) % max_size;
        }

        /** returns the next sample of type T in the delay line, interpolated if necessary.
        */
        inline const T Read() const
        {
            T a = line_[(write_ptr_ + delay_) % max_size];
            T b = line_[(write_ptr_ + delay_ + 1) % max_size];
            return a + (b - a) * frac_;
        }

        /** Read from a set location */
        inline const T Read(float delay) const
        {
            int32_t delay_integral   = static_cast<int32_t>(delay);
            float   delay_fractional = delay - static_cast<float>(delay_integral);
            const T a = line_[(write_ptr_ + delay_integral) % max_size];
            const T b = line_[(write_ptr_ + delay_integral + 1) % max_size];
            return a + (b - a) * delay_fractional;
        }

        inline const T ReadHermite(float delay) const
        {
            int32_t delay_integral   = static_cast<int32_t>(delay);
            float   delay_fractional = delay - static_cast<float>(delay_integral);

            int32_t     t     = (write_ptr_ + delay_integral + max_size);
            const T     xm1   = line_[(t - 1) % max_size];
            const T     x0    = line_[(t) % max_size];
            const T     x1    = line_[(t + 1) % max_size];
            const T     x2    = line_[(t + 2) % max_size];
            const float c     = (x1 - xm1) * 0.5f;
            const float v     = x0 - x1;
            const float w     = c + v;
            const float a     = w + v + (x2 - x0) * 0.5f;
            const float b_neg = w + a;
            const float f     = delay_fractional;
            return (((a * f) - b_neg) * f + c) * f + x0;
        }

        inline const T Allpass(const T sample, size_t delay, const T coefficient)
        {
            T read  = line_[(write_ptr_ + delay) % max_size];
            T write = sample + coefficient * read;
            Write(write);
            return -write * coefficient + read;
        }

    private:
        float  frac_;
        size_t write_ptr_;
        size_t delay_;
        T      line_[max_size];
    };

    //////////////////////////////////////////////////////////////////////////////////////
    // DC Block
    //////////////////////////////////////////////////////////////////////////////////////
    /** Removes DC component of a signal
        */
    class DcBlock
    {
    public:
        DcBlock(){};
        ~DcBlock(){};

        /** Initializes DcBlock module
        */
        void Init(float sample_rate);

        /** performs DcBlock Process 
        */
        float Process(float in);

    private:
        float input_, output_, gain_;
    };
    void DcBlock::Init(float sample_rate)
    {
        output_ = 0.0;
        input_  = 0.0;
        gain_   = 0.99;
    }

    float DcBlock::Process(float in)
    {
        float out;
        out     = in - input_ + (gain_ * output_);
        output_ = out;
        input_  = in;
        return out;
    }


    //////////////////////////////////////////////////////////////////////////////////////
    // Jitter CSound Opcode 
    //////////////////////////////////////////////////////////////////////////////////////

    class Jitter
    {
    public:
        Jitter() {}
        ~Jitter() {}

        /** Initializes Jitter module
        \param sample_rate Audio engine sample rate
        */
        void Init(float sample_rate);

        /** Get next floating point jitter sample */
        float Process();

        /** Set the minimum speed of the jitter engine.
        \param cps_min Number of new jitters per second
        */
        void SetCpsMin(float cps_min);

        /** Set the maximum speed of the jitter engine.
        \param cps_max Maximum number of jitters per second.
        */
        void SetCpsMax(float cps_max);

        /** Set the amplitude of the jitter. Jitters fall from -amp to +amp
        \param amp Jitter amplitude
        */
        void SetAmp(float amp);

    private:
        float   amp_, cps_min_, cps_max_, cps_, sample_rate_;
        int32_t phs_;
        bool    init_flag_;
        float   num1_, num2_, dfd_max_;
        float   randGab();
        float   biRandGab();
        void    Reset();
    };

    

    float Jitter::randGab()
    {
        return (float)((rand() >> 1) & 0x7fffffff)
            * (4.656612875245796924105750827168e-10);
    }

    float Jitter::biRandGab()
    {
        return (float)(rand() & 0x7fffffff)
            * (4.656612875245796924105750827168e-10);
    }

    void Jitter::SetAmp(float amp)
    {
        amp_ = amp;
        Reset();
    }

    void Jitter::SetCpsMin(float cps_min)
    {
        cps_min_ = cps_min;
        Reset();
    }

    void Jitter::SetCpsMax(float cps_max)
    {
        cps_max_ = cps_max;
        Reset();
    }

    void Jitter::Reset()
    {
        num2_      = biRandGab();
        init_flag_ = true;
        phs_       = 0;
    }

    void Jitter::Init(float sample_rate)
    {
        sample_rate_ = sample_rate;
        amp_         = 0.5;
        cps_min_     = 0.5;
        cps_max_     = 4;
        Reset();
    }

    float Jitter::Process()
    {
        float out;
        if(init_flag_)
        {
            init_flag_ = false;
            out        = num2_ * amp_;
            cps_       = randGab() * (cps_max_ - cps_min_) + cps_min_;
            phs_ &= FT_PHMASK;
            num1_    = num2_;
            num2_    = biRandGab();
            dfd_max_ = 1.0 * (num2_ - num1_) / FT_MAXLEN;
            return out;
        }

        out = (num1_ + (float)phs_ * dfd_max_) * amp_;
        phs_ += (int32_t)(cps_ * (float)(FT_MAXLEN / sample_rate_));

        if(phs_ >= FT_MAXLEN)
        {
            cps_ = randGab() * (cps_max_ - cps_min_) + cps_min_;
            phs_ &= FT_PHMASK;
            num1_    = num2_;
            num2_    = biRandGab();
            dfd_max_ = 1.0 * (num2_ - num1_) / FT_MAXLEN;
        }

        return out;
    }


    /** Multimode audio looper
    * 
    * Modes are:
    *  - Normal
    *  - Onetime Dub
    *  - Replace
    *  - Frippertronics
    *
    * Read more about the looper modes in the mode enum documentation.
    */
    class Looper
    {
    public:
        Looper() {}
        ~Looper() {}

        /** 
         ** Normal Mode: Input is added to the existing loop infinitely while recording 
        ** 
        ** Onetime Dub Mode: Recording starts at the first sample of the buffer and is added
        **     to the existing buffer contents. Recording automatically stops after one full loop.
        **
        ** Replace Mode: Audio in the buffer is replaced while recording is on.
        **
        ** Frippertronics Mode: infinite looping recording with fixed decay on each loop. The module acts like tape-delay set up.
        */
        enum class Mode
        {
            NORMAL,
            ONETIME_DUB,
            REPLACE,
            FRIPPERTRONICS,
        };

        void Init(float *mem, size_t size)
        {
            buffer_size_ = size;
            buff_        = mem;

            InitBuff();
            state_      = State::EMPTY;
            mode_       = Mode::NORMAL;
            half_speed_ = false;
            reverse_    = false;
            rec_queue_  = false;
            win_idx_    = 0;
        }

        /** Handles reading/writing to the Buffer depending on the mode. */
        float Process(const float input)
        {
            float sig = 0.f;
            float inc;
            bool  hitloop = false;
            // Record forward at normal speed during the first loop no matter what.
            inc = state_ == State::EMPTY || state_ == State::REC_FIRST
                    ? 1.f
                    : GetIncrementSize();
            win_ = WindowVal(win_idx_ * kWindowFactor);
            switch(state_)
            {
                case State::EMPTY: sig = 0.0f; break;
                case State::REC_FIRST:
                    sig = 0.f;
                    Write(pos_, input * win_);
                    if(win_idx_ < kWindowSamps - 1)
                        win_idx_ += 1;
                    recsize_ = pos_;
                    pos_ += inc;
                    if(pos_ > buffer_size_ - 1)
                    {
                        state_   = State::PLAYING;
                        recsize_ = pos_ - 1;
                        pos_     = 0;
                    }
                    break;
                case State::PLAYING:
                    sig = Read(pos_);
                    /** This is a way of 'seamless looping'
                     ** The first N samps after recording is done are recorded with the input faded out. 
                    */
                    if(win_idx_ < kWindowSamps - 1)
                    {
                        Write(pos_, sig + input * (1.f - win_));
                        win_idx_ += 1;
                    }

                    pos_ += inc;
                    if(pos_ > recsize_ - 1)
                    {
                        pos_    = 0;
                        hitloop = true;
                    }
                    else if(pos_ < 0)
                    {
                        pos_    = recsize_ - 1;
                        hitloop = true;
                    }
                    if(hitloop)

                    {
                        if(rec_queue_ && mode_ == Mode::ONETIME_DUB)
                        {
                            rec_queue_ = false;
                            state_     = State::REC_DUB;
                            win_idx_   = 0;
                        }
                    }
                    break;
                case State::REC_DUB:
                    sig = Read(pos_);
                    switch(mode_)
                    {
                        case Mode::REPLACE: Write(pos_, input * win_); break;
                        case Mode::FRIPPERTRONICS:
                            Write(pos_, (input * win_) + (sig * kFripDecayVal));
                            break;
                        case Mode::NORMAL:
                        case Mode::ONETIME_DUB:
                        default: Write(pos_, (input * win_) + sig); break;
                    }
                    if(win_idx_ < kWindowSamps - 1)
                        win_idx_ += 1;
                    pos_ += inc;
                    if(pos_ > recsize_ - 1)
                    {
                        pos_    = 0;
                        hitloop = true;
                    }
                    else if(pos_ < 0)
                    {
                        pos_    = recsize_ - 1;
                        hitloop = true;
                    }
                    if(hitloop && mode_ == Mode::ONETIME_DUB)
                    {
                        state_   = State::PLAYING;
                        win_idx_ = 0;
                    }

                    break;
                default: break;
            }
            near_beginning_ = state_ != State::EMPTY && !Recording() && pos_ < 4800
                                ? true
                                : false;

            return sig;
        }

        /** Effectively erases the buffer 
         ** Note: This does not actually change what is in the buffer  */
        inline void Clear() { state_ = State::EMPTY; }

        /** Engages/Disengages the recording, depending on Mode.
         ** In all modes, the first time this is triggered a new loop will be started.
        ** The second trigger will set the loop size, and begin playback of the loop. 
        */
        inline void TrigRecord()
        {
            switch(state_)
            {
                case State::EMPTY:
                    pos_        = 0;
                    recsize_    = 0;
                    state_      = State::REC_FIRST;
                    half_speed_ = false;
                    reverse_    = false;
                    break;
                case State::REC_FIRST:
                case State::REC_DUB: state_ = State::PLAYING; break;
                case State::PLAYING:
                    if(mode_ == Mode::ONETIME_DUB)
                        rec_queue_ = true;
                    else
                        state_ = State::REC_DUB;
                    break;
                default: state_ = State::EMPTY; break;
            }
            if(!rec_queue_)
                win_idx_ = 0;
        }

        /** Returns true if the looper is currently being written to. */
        inline const bool Recording() const
        {
            return state_ == State::REC_DUB || state_ == State::REC_FIRST;
        }

        inline const bool RecordingQueued() const { return rec_queue_; }

        /** Increments the Mode by one step useful for buttons, etc. that need to step through the Looper modes. */
        inline void IncrementMode()
        {
            int m = static_cast<int>(mode_);
            m     = m + 1;
            if(m > kNumModes - 1)
                m = 0;
            mode_ = static_cast<Mode>(m);
        }

        /** Sets the recording mode to the specified Mode. */
        inline void SetMode(Mode mode) { mode_ = mode; }

        /** Returns the specific recording mode that is currently set. */
        inline const Mode GetMode() const { return mode_; }

        inline void ToggleReverse() { reverse_ = !reverse_; }
        inline void SetReverse(bool state) { reverse_ = state; }
        inline bool GetReverse() const { return reverse_; }

        inline void ToggleHalfSpeed() { half_speed_ = !half_speed_; }
        inline void SetHalfSpeed(bool state) { half_speed_ = state; }
        inline bool GetHalfSpeed() const { return half_speed_; }

        inline bool IsNearBeginning() { return near_beginning_; }

    private:
        /** Constants */

        /** Decay value for frippertronics mode is sin(PI / 4) */
        static constexpr float kFripDecayVal      = 0.7071067811865476f;
        static constexpr int   kNumModes          = 4;
        static constexpr int   kNumPlaybackSpeeds = 3;
        static constexpr int   kWindowSamps       = 1200;
        static constexpr float kWindowFactor      = (1.f / kWindowSamps);

        /** Private Member Functions */
        float GetIncrementSize()
        {
            float inc = 1.f;
            if(half_speed_)
                inc = 0.5f;
            return reverse_ ? -inc : inc;
        }

        /** Initialize the buffer */
        void InitBuff() { std::fill(&buff_[0], &buff_[buffer_size_ - 1], 0); }

        /** Get a floating point sample from the buffer */
        inline const float Read(size_t pos) const { return buff_[pos]; }

        /** Reads from a specified point in the delay line using linear interpolation */
        float ReadF(float pos)
        {
            float    a, b, frac;
            uint32_t i_idx = static_cast<uint32_t>(pos);
            frac           = pos - i_idx;
            a              = buff_[i_idx];
            b              = buff_[(i_idx + 1) % buffer_size_];
            return a + (b - a) * frac;
        }

        /** Write to a known location in the buffer */
        inline void Write(size_t pos, float val) { buff_[pos] = val; }

        /** Linear to Constpower approximation for windowing*/
        float WindowVal(float in) { return sin(HALFPI_F * in); }

        // Private Enums

        /** Internal looper state */
        enum class State
        {
            EMPTY,
            REC_FIRST,
            PLAYING,
            REC_DUB,
        };

        /** Private Member Variables */
        Mode   mode_;
        State  state_;
        float *buff_;
        size_t buffer_size_;
        float  pos_, win_;
        size_t win_idx_;
        bool   half_speed_;
        bool   reverse_;
        size_t recsize_;
        bool   rec_queue_;
        bool   near_beginning_;
    };
    

    /** Probabilistic trigger  module

    Original author(s) : Paul Batchelor

    Ported from soundpipe by Ben Sergentanis, May 2020
    */
    class Maytrig
    {
    public:
        Maytrig() {}
        ~Maytrig() {}
        /** probabilistically generates triggers

            \param prob (1 always returns true, 0 always false)

            \return given a probability 0 to 1, returns true or false.
        */
        inline float Process(float prob)
        {
            return ((float)rand() / (float)RAND_MAX) <= prob ? true : false;
        }

    private:
    };

    /** Creates a clock signal at a specific frequency.
    */
    class Metro
    {
    public:
        Metro() {}
        ~Metro() {}
        /** Initializes Metro module.
            Arguments:
            - freq: frequency at which new clock signals will be generated
                Input Range: 
            - sample_rate: sample rate of audio engine
                Input range: 
        */
        void Init(float freq, float sample_rate);

        /** checks current state of Metro object and updates state if necesary.
        */
        uint8_t Process();

        /** resets phase to 0
        */
        inline void Reset() { phs_ = 0.0f; }
        /** Sets frequency at which Metro module will run at.
        */
        void SetFreq(float freq);

        /** Returns current value for frequency.
        */
        inline float GetFreq() { return freq_; }

    private:
        float freq_;
        float phs_, sample_rate_, phs_inc_;
    };

    void Metro::Init(float freq, float sample_rate)
    {
        freq_        = freq;
        phs_         = 0.0f;
        sample_rate_ = sample_rate;
        phs_inc_     = (TWOPI_F * freq_) / sample_rate_;
    }

    uint8_t Metro::Process()
    {
        phs_ += phs_inc_;
        if(phs_ >= TWOPI_F)
        {
            phs_ -= TWOPI_F;
            return 1;
        }
        return 0;
    }

    void Metro::SetFreq(float freq)
    {
        freq_    = freq;
        phs_inc_ = (TWOPI_F * freq_) / sample_rate_;
    }

    /** Applies portamento to an input signal. 

    At each new step value, the input is low-pass filtered to 
    move towards that value at a rate determined by ihtim. ihtim is the half-time of the 
    function (in seconds), during which the curve will traverse half the distance towards the new value, 
    then half as much again, etc., theoretically never reaching its asymptote.

    This code has been ported from Soundpipe to DaisySP by Paul Batchelor. 

    The Soundpipe module was extracted from the Csound opcode "portk".

    Original Author(s): Robbin Whittle, John ffitch

    Year: 1995, 1998

    Location: Opcodes/biquad.c
    */
    class Port
    {
    public:
        Port() {}
        ~Port() {}
        /** Initializes Port module

            \param sample_rate: sample rate of audio engine
            \param htime: half-time of the function, in seconds.
        */

        void Init(float sample_rate, float htime);

        /** Applies portamento to input signal and returns processed signal. 
            \return slewed output signal
        */
        float Process(float in);


        /** Sets htime
        */
        inline void SetHtime(float htime) { htime_ = htime; }
        /** returns current value of htime
        */
        inline float GetHtime() { return htime_; }

    private:
        float htime_;
        float c1_, c2_, yt1_, prvhtim_;
        float sample_rate_, onedsr_;
    };

    

    
    
    void Port::Init(float sample_rate, float htime)
    {
        yt1_     = 0;
        prvhtim_ = -100.0;
        htime_   = htime;

        sample_rate_ = sample_rate;
        onedsr_      = 1.0 / sample_rate_;
    }

    float Port::Process(float in)
    {
        if(prvhtim_ != htime_)
        {
            c2_      = pow(0.5, onedsr_ / htime_);
            c1_      = 1.0 - c2_;
            prvhtim_ = htime_;
        }

        return yt1_ = c1_ * in + c2_ * yt1_;
    }


    /** Dual track and hold / Sample and hold module. \n 
        Ported from soundpipe by Ben Sergentanis, June 2020.
        @author Paul Batchelor
        @date   2015
    */
    class SampleHold
    {
    public:
        SampleHold() {}
        ~SampleHold() {}

        enum Mode
        {
            MODE_SAMPLE_HOLD,
            MODE_TRACK_HOLD,
            MODE_LAST,
        };

        /** Process the next sample. Both sample and track and hold are run in parallel
        \param trigger Trigger the sample/track and hold
        \param input   Signal to be sampled/tracked and held
        \param mode    Whether to output the tracked or sampled values.
        */
        inline float
        Process(bool trigger, float input, Mode mode = MODE_SAMPLE_HOLD)
        {
            Update(trigger, input);
            return mode == MODE_SAMPLE_HOLD ? sample_ : track_;
        }

    private:
        float track_    = 0;
        float sample_   = 0;
        bool  previous_ = false;


        inline void Update(bool trigger, float input)
        {
            if(trigger)
            {
                if(!previous_)
                {
                    sample_ = input;
                }
                track_ = input;
            }
            previous_ = trigger;
        }
    };    


    /**  
         @brief Smooth random generator for internal modulation. \n
        @author Ported by Ben Sergentanis 
        @date Jan 2021 
        Ported from pichenettes/eurorack/plaits/dsp/noise/smooth_random_generator.h \n
        to an independent module. \n
        Original code written by Emilie Gillet in 2016. \n
    */
    class SmoothRandomGenerator
    {
    public:
        SmoothRandomGenerator() {}
        ~SmoothRandomGenerator() {}

        /** Initialize the module
            \param sample_rate Audio engine sample rate.
        */
        void Init(float sample_rate)
        {
            sample_rate_ = sample_rate;

            SetFreq(1.f);
            phase_    = 0.0f;
            from_     = 0.0f;
            interval_ = 0.0f;
        }

        /** Get the next float. Ranges from -1 to 1. */
        float Process()
        {
            phase_ += frequency_;
            if(phase_ >= 1.0f)
            {
                phase_ -= 1.0f;
                from_ += interval_;
                interval_ = rand() * kRandFrac * 2.0f - 1.0f - from_;
            }
            float t = phase_ * phase_ * (3.0f - 2.0f * phase_);
            return from_ + interval_ * t;
        }

        /** How often to slew to a new random value
            \param freq Rate in Hz
        */
        void SetFreq(float freq)
        {
            freq       = freq / sample_rate_;
            frequency_ = fclamp(freq, 0.f, 1.f);
        }

    private:
        float frequency_;
        float phase_;
        float from_;
        float interval_;

        float sample_rate_;

        static constexpr float kRandFrac = 1.f / (float)RAND_MAX;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // ADEnv
    ////////////////////////////////////////////////////////////////////////////////////////////

    enum AdEnvSegment
    {
        /** located at phase location 0, and not currently running */
        ADENV_SEG_IDLE,
        /** First segment of envelope where phase moves from MIN value to MAX value */
        ADENV_SEG_ATTACK,
        /** Second segment of envelope where phase moves from MAX to MIN value */
        ADENV_SEG_DECAY,
        /** The final segment of the envelope (currently decay) */
        ADENV_SEG_LAST,
    };

    /** Trigger-able envelope with adjustable min/max, and independent per-segment time control.

        \author shensley
        \todo - Add Cycling
        \todo - Implement Curve (its only linear for now).
        \todo - Maybe make this an ADsr_ that has AD/AR/Asr_ modes.
    */
    class AdEnv
    {
    public:
        AdEnv() {}
        ~AdEnv() {}
        /** Initializes the ad envelope.

            Defaults:
                - current segment = idle
                - curve = linear
                - phase = 0
                - min = 0
                - max = 1

            \param sample_rate      sample rate of the audio engine being run
        */
        void Init(float sample_rate);

        /** Processes the current sample of the envelope. This should be called once
            per sample period.
            \return the current envelope value.
        */
        float Process();

        /** Starts or retriggers the envelope.*/
        inline void Trigger() { trigger_ = 1; }
        /** Sets the length of time (in seconds) for a specific segment. */
        inline void SetTime(uint8_t seg, float time) { segment_time_[seg] = time; }
        /** Sets the amount of curve applied. A positve value will create a log
            curve. Input range: -100 to 100.  (or more)
        */
        inline void SetCurve(float scalar) { curve_scalar_ = scalar; }
        /** Sets the minimum value of the envelope output.
            Input range: -FLTmax_, to FLTmax_
        */
        inline void SetMin(float min) { min_ = min; }
        /** Sets the maximum value of the envelope output.
            Input range: -FLTmax_, to FLTmax_
        */
        inline void SetMax(float max) { max_ = max; }
        /** Returns the current output value without processing the next sample */
        inline float GetValue() const { return (output_ * (max_ - min_)) + min_; }
        /** Returns the segment of the envelope that the phase is currently located
            in.
        */
        inline uint8_t GetCurrentSegment() { return current_segment_; }
        /** Returns true if the envelope is currently in any stage apart from idle.
        */
        inline bool IsRunning() const { return current_segment_ != ADENV_SEG_IDLE; }

    private:
        uint8_t  current_segment_, prev_segment_;
        float    segment_time_[ADENV_SEG_LAST];
        float    sample_rate_, min_, max_, output_, curve_scalar_;
        float    c_inc_, curve_x_, retrig_val_;
        uint32_t phase_;
        uint8_t  trigger_;
    };


    
    // Private Functions
    inline void AdEnv::Init(float sample_rate)
    {
        sample_rate_     = sample_rate;
        current_segment_ = ADENV_SEG_IDLE;
        curve_scalar_    = 0.0f; // full linear
        phase_           = 0;
        min_             = 0.0f;
        max_             = 1.0f;
        output_          = 0.0001f;
        for(uint8_t i = 0; i < ADENV_SEG_LAST; i++)
        {
            segment_time_[i] = 0.05f;
        }
    }

    inline float AdEnv::Process()
    {
        uint32_t time_samps;
        float    val, out, end, beg, inc;

        // Handle Retriggering
        if(trigger_)
        {
            trigger_         = 0;
            current_segment_ = ADENV_SEG_ATTACK;
            phase_           = 0;
            curve_x_         = 0.0f;
            retrig_val_      = output_;
        }

        time_samps = (uint32_t)(segment_time_[current_segment_] * sample_rate_);

        // Fixed for now, but we could always make this a more flexible multi-segment envelope
        switch(current_segment_)
        {
            case ADENV_SEG_ATTACK:
                beg = retrig_val_;
                end = 1.0f;
                break;
            case ADENV_SEG_DECAY:
                beg = 1.0f;
                end = 0.0f;
                break;
            case ADENV_SEG_IDLE:
            default:
                beg = 0;
                end = 0;
                break;
        }

        if(prev_segment_ != current_segment_)
        {
            //Reset at segment beginning
            curve_x_ = 0;
            phase_   = 0;
        }

        //recalculate increment value
        if(curve_scalar_ == 0.0f)
        {
            c_inc_ = (end - beg) / time_samps;
        }
        else
        {
            c_inc_ = (end - beg) / (1.0f - EXPF(curve_scalar_));
        }


        // update output
        val = output_;
        inc = c_inc_;
        out = val;
        if(curve_scalar_ == 0.0f)
        {
            val += inc;
        }
        else
        {
            curve_x_ += (curve_scalar_ / time_samps);
            val = beg + inc * (1.0f - EXPF(curve_x_));
            if(val != val)
                val = 0.0f; // NaN check
        }

        // Update Segment
        phase_ += 1;
        prev_segment_ = current_segment_;
        if(current_segment_ != ADENV_SEG_IDLE)
        {
            if((out >= 1.f && current_segment_ == ADENV_SEG_ATTACK)
            || (out <= 0.f && current_segment_ == ADENV_SEG_DECAY))
            {
                // Advance segment
                current_segment_++;
                // TODO: Add Cycling feature here.
                if(current_segment_ > ADENV_SEG_DECAY)
                {
                    current_segment_ = ADENV_SEG_IDLE;
                }
            }
        }
        if(current_segment_ == ADENV_SEG_IDLE)
        {
            val = out = 0.0f;
        }
        output_ = val;

        return out * (max_ - min_) + min_;
    }

    /** Distinct stages that the phase of the envelope can be located in.
    - IDLE   = located at phase location 0, and not currently running
    - ATTACK  = First segment of envelope where phase moves from 0 to 1
    - DECAY   = Second segment of envelope where phase moves from 1 to SUSTAIN value
    - RELEASE =     Fourth segment of envelop where phase moves from SUSTAIN to 0
    */
    enum
    {
        ADSR_SEG_IDLE    = 0,
        ADSR_SEG_ATTACK  = 1,
        ADSR_SEG_DECAY   = 2,
        ADSR_SEG_RELEASE = 4
    };


    /** adsr envelope module

    Original author(s) : Paul Batchelor

    Ported from Soundpipe by Ben Sergentanis, May 2020
    
    Remake by Steffan DIedrichsen, May 2021
    */
    class Adsr
    {
    public:
        Adsr() {}
        ~Adsr() {}
        /** Initializes the Adsr module.
            \param sample_rate - The sample rate of the audio engine being run. 
        */
        void Init(float sample_rate, int blockSize = 1);
        /**
         \function Retrigger forces the envelope back to attack phase
        \param hard  resets the history to zero, results in a click.
        */
        void Retrigger(bool hard);
        /** Processes one sample through the filter and returns one sample.
            \param gate - trigger the envelope, hold it to sustain 
        */
        float Process(bool gate);
        /** Sets time
            Set time per segment in seconds
        */
        void SetTime(int seg, float time);
        void SetAttackTime(float timeInS, float shape = 0.0f);
        void SetDecayTime(float timeInS);
        void SetReleaseTime(float timeInS);

    private:
        void SetTimeConstant(float timeInS, float& time, float& coeff);

    public:
        /** Sustain level
            \param sus_level - sets sustain level, 0...1.0
        */
        inline void SetSustainLevel(float sus_level)
        {
            sus_level = (sus_level <= 0.f) ? -0.01f // forces envelope into idle
                                        : (sus_level > 1.f) ? 1.f : sus_level;
            sus_level_ = sus_level;
        }
        /** get the current envelope segment
            \return the segment of the envelope that the phase is currently located in.
        */
        inline uint8_t GetCurrentSegment() { return mode_; }
        /** Tells whether envelope is active
            \return true if the envelope is currently in any stage apart from idle.
        */
        inline bool IsRunning() const { return mode_ != ADSR_SEG_IDLE; }

    private:
        float   sus_level_{0.f};
        float   x_{0.f};
        float   attackShape_{-1.f};
        float   attackTarget_{0.0f};
        float   attackTime_{-1.0f};
        float   decayTime_{-1.0f};
        float   releaseTime_{-1.0f};
        float   attackD0_{0.f};
        float   decayD0_{0.f};
        float   releaseD0_{0.f};
        int     sample_rate_;
        uint8_t mode_{ADSR_SEG_IDLE};
        bool    gate_{false};
    };

    inline void Adsr::Init(float sample_rate, int blockSize)
    {
        sample_rate_  = sample_rate / blockSize;
        attackShape_  = -1.f;
        attackTarget_ = 0.0f;
        attackTime_   = -1.f;
        decayTime_    = -1.f;
        releaseTime_  = -1.f;
        sus_level_    = 0.7f;
        x_            = 0.0f;
        gate_         = false;
        mode_         = ADSR_SEG_IDLE;

        SetTime(ADSR_SEG_ATTACK, 0.1f);
        SetTime(ADSR_SEG_DECAY, 0.1f);
        SetTime(ADSR_SEG_RELEASE, 0.1f);
    }

    inline void Adsr::Retrigger(bool hard)
    {
        mode_ = ADSR_SEG_ATTACK;
        if(hard)
            x_ = 0.f;
    }

    inline void Adsr::SetTime(int seg, float time)
    {
        switch(seg)
        {
            case ADSR_SEG_ATTACK: SetAttackTime(time, 0.0f); break;
            case ADSR_SEG_DECAY:
            {
                SetTimeConstant(time, decayTime_, decayD0_);
            }
            break;
            case ADSR_SEG_RELEASE:
            {
                SetTimeConstant(time, releaseTime_, releaseD0_);
            }
            break;
            default: return;
        }
    }

    inline void Adsr::SetAttackTime(float timeInS, float shape)
    {
        if((timeInS != attackTime_) || (shape != attackShape_))
        {
            attackTime_  = timeInS;
            attackShape_ = shape;
            if(timeInS > 0.f)
            {
                float x         = shape;
                float target    = 9.f * powf(x, 10.f) + 0.3f * x + 1.01f;
                attackTarget_   = target;
                float logTarget = logf(1.f - (1.f / target)); // -1 for decay
                attackD0_       = 1.f - expf(logTarget / (timeInS * sample_rate_));
            }
            else
                attackD0_ = 1.f; // instant change
        }
    }
    inline void Adsr::SetDecayTime(float timeInS)
    {
        SetTimeConstant(timeInS, decayTime_, decayD0_);
    }
    inline void Adsr::SetReleaseTime(float timeInS)
    {
        SetTimeConstant(timeInS, releaseTime_, releaseD0_);
    }
    inline void Adsr::SetTimeConstant(float timeInS, float& time, float& coeff)
    {
        if(timeInS != time)
        {
            time = timeInS;
            if(time > 0.f)
            {
                const float target = logf(1. / M_E);
                coeff              = 1.f - expf(target / (time * sample_rate_));
            }
            else
                coeff = 1.f; // instant change
        }
    }
    inline float Adsr::Process(bool gate)
    {
        float out = 0.0f;

        if(gate && !gate_) // rising edge
            mode_ = ADSR_SEG_ATTACK;
        else if(!gate && gate_) // falling edge
            mode_ = ADSR_SEG_RELEASE;
        gate_ = gate;

        float D0(attackD0_);
        if(mode_ == ADSR_SEG_DECAY)
            D0 = decayD0_;
        else if(mode_ == ADSR_SEG_RELEASE)
            D0 = releaseD0_;

        float target = mode_ == ADSR_SEG_DECAY ? sus_level_ : -0.01f;
        switch(mode_)
        {
            case ADSR_SEG_IDLE: out = 0.0f; break;
            case ADSR_SEG_ATTACK:
                x_ += D0 * (attackTarget_ - x_);
                out = x_;
                if(out > 1.f)
                {
                    x_ = out = 1.f;
                    mode_    = ADSR_SEG_DECAY;
                }
                break;
            case ADSR_SEG_DECAY:
            case ADSR_SEG_RELEASE:
                x_ += D0 * (target - x_);
                out = x_;
                if(out < 0.0f)
                {
                    x_ = out = 0.f;
                    mode_    = ADSR_SEG_IDLE;
                }
            default: break;
        }
        return out;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Ramp
    ////////////////////////////////////////////////////////////////////////////////////////////

    class Line
    {
    public:
        Line() {}
        ~Line() {}
        /** Initializes Line module.
        */
        void Init(float sample_rate);

        /** Processes Line segment. Returns one sample.
            value of finished will be updated to a 1, upon completion of the Line's trajectory.
        */
        float Process(uint8_t *finished);

        /** Begin creation of Line. 
            \param start - beginning value
            \param end - ending value
            \param dur - duration in seconds of Line segment
        */
        void Start(float start, float end, float dur);

    private:
        float   start_, end_, dur_;
        float   inc_, val_, sample_rate_;
        uint8_t finished_;
    };

    inline void Line::Init(float sample_rate)
    {
        sample_rate_ = sample_rate;
        dur_         = 0.5f;
        end_         = 0.0f;
        start_       = 1.0f;
        val_         = 1.0f;
    }

    inline void Line::Start(float start, float end, float dur)
    {
        start_    = start;
        end_      = end;
        dur_      = dur;
        inc_      = (end - start) / ((sample_rate_ * dur_));
        val_      = start_;
        finished_ = 0;
    }

    inline float Line::Process(uint8_t *finished)
    {
        float out;
        out = val_;

        if((end_ > start_ && out >= end_) || (end_ < start_ && out <= end_))
        {
            finished_ = 1;
            val_      = end_;
            out       = end_;
        }
        else
        {
            val_ += inc_;
        }
        *finished = finished_;
        return out;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////
    // Phasor
    ////////////////////////////////////////////////////////////////////////////////////////////

    class Phasor
    {
    public:
        Phasor() {}
        ~Phasor() {}
        /** Initializes the Phasor module
        sample rate, and freq are in Hz
        initial phase is in radians
        Additional Init functions have defaults when arg is not specified:
        - phs = 0.0f
        - freq = 1.0f
        */
        inline void Init(float sample_rate, float freq, float initial_phase)
        {
            sample_rate_ = sample_rate;
            phs_         = initial_phase;
            SetFreq(freq);
        }

        /** Initialize phasor with samplerate and freq
        */
        inline void Init(float sample_rate, float freq)
        {
            Init(sample_rate, freq, 0.0f);
        }

        /** Initialize phasor with samplerate
        */
        inline void Init(float sample_rate) { Init(sample_rate, 1.0f, 0.0f); }
        /** processes Phasor and returns current value
        */
        float Process();


        /** Sets frequency of the Phasor in Hz
        */
        void SetFreq(float freq);


        /** Returns current frequency value in Hz
        */
        inline float GetFreq() { return freq_; }

    private:
        float freq_;
        float sample_rate_, inc_, phs_;
    };

    inline void Phasor::SetFreq(float freq)
    {
        freq_ = freq;
        inc_  = (TWOPI_F * freq_) / sample_rate_;
    }

    inline float Phasor::Process()
    {
        float out;
        out = phs_ / TWOPI_F;
        phs_ += inc_;
        if(phs_ > TWOPI_F)
        {
            phs_ -= TWOPI_F;
        }
        if(phs_ < 0.0f)
        {
            phs_ = 0.0f;
        }
        return out;
    }    


    /** Synthesis of several waveforms, including polyBLEP bandlimited waveforms.
    */
    class Oscillator
    {
    public:
        Oscillator() {}
        ~Oscillator() {}
        /** Choices for output waveforms, POLYBLEP are appropriately labeled. Others are naive forms.
        */
        enum
        {
            WAVE_SIN,
            WAVE_TRI,
            WAVE_SAW,
            WAVE_RAMP,
            WAVE_SQUARE,
            WAVE_POLYBLEP_TRI,
            WAVE_POLYBLEP_SAW,
            WAVE_POLYBLEP_SQUARE,
            WAVE_LAST,
        };


        /** Initializes the Oscillator

            \param sample_rate - sample rate of the audio engine being run, and the frequency that the Process function will be called.

            Defaults:
            - freq_ = 100 Hz
            - amp_ = 0.5
            - waveform_ = sine wave.
        */
        void Init(float sample_rate)
        {
            sr_        = sample_rate;
            sr_recip_  = 1.0f / sample_rate;
            freq_      = 100.0f;
            amp_       = 0.5f;
            phase_     = 0.0f;
            phase_inc_ = CalcPhaseInc(freq_);
            waveform_  = WAVE_SIN;
            eoc_       = true;
            eor_       = true;
        }


        /** Changes the frequency of the Oscillator, and recalculates phase increment.
        */
        inline void SetFreq(const float f)
        {
            freq_      = f;
            phase_inc_ = CalcPhaseInc(f);
        }


        /** Sets the amplitude of the waveform.
        */
        inline void SetAmp(const float a) { amp_ = a; }
        /** Sets the waveform to be synthesized by the Process() function.
        */
        inline void SetWaveform(const uint8_t wf)
        {
            waveform_ = wf < WAVE_LAST ? wf : WAVE_SIN;
        }

        /** Returns true if cycle is at end of rise. Set during call to Process.
        */
        inline bool IsEOR() { return eor_; }

        /** Returns true if cycle is at end of cycle. Set during call to Process.
        */
        inline bool IsEOC() { return eoc_; }

        /** Returns true if cycle rising.
        */
        inline bool IsRising() { return phase_ < PI_F; }

        /** Returns true if cycle falling.
        */
        inline bool IsFalling() { return phase_ >= PI_F; }

        /** Processes the waveform to be generated, returning one sample. This should be called once per sample period.
        */
        float Process();


        /** Adds a value 0.0-1.0 (mapped to 0.0-TWO_PI) to the current phase. Useful for PM and "FM" synthesis.
        */
        void PhaseAdd(float _phase) { phase_ += (_phase * TWOPI_F); }
        /** Resets the phase to the input argument. If no argumeNt is present, it will reset phase to 0.0;
        */
        void Reset(float _phase = 0.0f) { phase_ = _phase; }

    private:
        float   CalcPhaseInc(float f);
        uint8_t waveform_;
        float   amp_, freq_;
        float   sr_, sr_recip_, phase_, phase_inc_;
        float   last_out_, last_freq_;
        bool    eor_, eoc_;
    };


    constexpr float TWO_PI_RECIP = 1.0f / TWOPI_F;

    float Polyblep(float phase_inc, float t)
    {
        float dt = phase_inc * TWO_PI_RECIP;
        if(t < dt)
        {
            t /= dt;
            return t + t - t * t - 1.0f;
        }
        else if(t > 1.0f - dt)
        {
            t = (t - 1.0f) / dt;
            return t * t + t + t + 1.0f;
        }
        else
        {
            return 0.0f;
        }
    }    

    float Oscillator::Process()
    {
        float out, t;
        switch(waveform_)
        {
            case WAVE_SIN: out = sinf(phase_); break;
            case WAVE_TRI:
                t   = -1.0f + (2.0f * phase_ * TWO_PI_RECIP);
                out = 2.0f * (fabsf(t) - 0.5f);
                break;
            case WAVE_SAW:
                out = -1.0f * (((phase_ * TWO_PI_RECIP * 2.0f)) - 1.0f);
                break;
            case WAVE_RAMP: out = ((phase_ * TWO_PI_RECIP * 2.0f)) - 1.0f; break;
            case WAVE_SQUARE: out = phase_ < PI_F ? (1.0f) : -1.0f; break;
            case WAVE_POLYBLEP_TRI:
                t   = phase_ * TWO_PI_RECIP;
                out = phase_ < PI_F ? 1.0f : -1.0f;
                out += Polyblep(phase_inc_, t);
                out -= Polyblep(phase_inc_, fmodf(t + 0.5f, 1.0f));
                // Leaky Integrator:
                // y[n] = A + x[n] + (1 - A) * y[n-1]
                out       = phase_inc_ * out + (1.0f - phase_inc_) * last_out_;
                last_out_ = out;
                break;
            case WAVE_POLYBLEP_SAW:
                t   = phase_ * TWO_PI_RECIP;
                out = (2.0f * t) - 1.0f;
                out -= Polyblep(phase_inc_, t);
                out *= -1.0f;
                break;
            case WAVE_POLYBLEP_SQUARE:
                t   = phase_ * TWO_PI_RECIP;
                out = phase_ < PI_F ? 1.0f : -1.0f;
                out += Polyblep(phase_inc_, t);
                out -= Polyblep(phase_inc_, fmodf(t + 0.5f, 1.0f));
                out *= 0.707f; // ?
                break;
            default: out = 0.0f; break;
        }
        phase_ += phase_inc_;
        if(phase_ > TWOPI_F)
        {
            phase_ -= TWOPI_F;
            eoc_ = true;
        }
        else
        {
            eoc_ = false;
        }
        eor_ = (phase_ - phase_inc_ < PI_F && phase_ >= PI_F);

        return out * amp_;
    }

    float Oscillator::CalcPhaseInc(float f)
    {
        return (TWOPI_F * f) * sr_recip_;
    }

    /**  
         @brief Oscillator Bank module.
        @author Ben Sergentanis
        @date Dec 2020 
        A mixture of 7 sawtooth and square waveforms in the style of divide-down organs \n \n
        Ported from pichenettes/eurorack/plaits/dsp/oscillator/string_synth_oscillator.h \n 
        \n to an independent module. \n
        Original code written by Emilie Gillet in 2016. \n
    */

    class OscillatorBank
    {
    public:
        OscillatorBank() {}
        ~OscillatorBank() {}

        /** Init string synth module
            \param sample_rate Audio engine sample rate
        */
        void Init(float sample_rate);

        /** Get next floating point sample 
        */
        float Process();

        /** Set oscillator frequency (8' oscillator)
            \param freq Frequency in Hz
        */
        void SetFreq(float freq);

        /** Set amplitudes of 7 oscillators. 0-6 are Saw 8', Square 8', Saw 4', Square 4', Saw 2', Square 2', Saw 1'
            \param amplitudes array of 7 floating point amplitudes. Must sum to 1.
        */
        void SetAmplitudes(const float* amplitudes);

        /** Set a single amplitude
            \param amp Amplitude to set.
            \param idx Which wave's amp to set
        */
        void SetSingleAmp(float amp, int idx);

        /** Set overall gain.
            \param gain Gain to set. 0-1.
        */
        void SetGain(float gain);

    private:
        // Oscillator state.
        float phase_;
        float next_sample_;
        int   segment_;
        float gain_;
        float registration_[7];
        float unshifted_registration_[7];

        float frequency_;
        float saw_8_gain_;
        float saw_4_gain_;
        float saw_2_gain_;
        float saw_1_gain_;

        float sample_rate_;
        bool  recalc_, recalc_gain_;

        bool cmp(float a, float b) { return fabsf(a - b) > .0000001; }
    };
    void OscillatorBank::Init(float sample_rate)
    {
        sample_rate_ = sample_rate;

        phase_       = 0.0f;
        next_sample_ = 0.0f;
        segment_     = 0.0f;

        frequency_  = 0.f;
        saw_8_gain_ = 0.0f;
        saw_4_gain_ = 0.0f;
        saw_2_gain_ = 0.0f;
        saw_1_gain_ = 0.0f;

        recalc_ = recalc_gain_ = true;
        SetGain(1.f);

        for(int i = 0; i < 7; i++)
        {
            registration_[i]           = 0.f;
            unshifted_registration_[i] = 0.f;
        }
        SetSingleAmp(1.f, 0);
        SetFreq(440.f);
    }

    float OscillatorBank::Process()
    {
        if(recalc_)
        {
            recalc_ = false;
            frequency_ *= 8.0f;

            // Deal with very high frequencies by shifting everything 1 or 2 octave
            // down: Instead of playing the 1st harmonic of a 8kHz wave, we play the
            // second harmonic of a 4kHz wave.
            size_t shift = 0;
            while(frequency_ > 0.5f)
            {
                shift += 2;
                frequency_ *= 0.5f;
            }

            for(int i = 0; i < 7; i++)
            {
                registration_[i] = 0.f;
            }

            for(size_t i = 0; i < 7 - shift; i++)
            {
                registration_[i + shift] = unshifted_registration_[i];
            }
        }

        if(recalc_gain_ || recalc_)
        {
            saw_8_gain_ = (registration_[0] + 2.0f * registration_[1]) * gain_;
            saw_4_gain_
                = (registration_[2] - registration_[1] + 2.0f * registration_[3])
                * gain_;
            saw_2_gain_
                = (registration_[4] - registration_[3] + 2.0f * registration_[5])
                * gain_;
            saw_1_gain_ = (registration_[6] - registration_[5]) * gain_;
        }

        float this_sample_ = next_sample_;
        next_sample_       = 0.0f;

        phase_ += frequency_;
        int next_segment_ = static_cast<int>(phase_);
        if(next_segment_ != segment_)
        {
            float discontinuity = 0.0f;
            if(next_segment_ == 8)
            {
                phase_ -= 8.0f;
                next_segment_ -= 8;
                discontinuity -= saw_8_gain_;
            }
            if((next_segment_ & 3) == 0)
            {
                discontinuity -= saw_4_gain_;
            }
            if((next_segment_ & 1) == 0)
            {
                discontinuity -= saw_2_gain_;
            }
            discontinuity -= saw_1_gain_;
            if(discontinuity != 0.0f)
            {
                float fraction = phase_ - static_cast<float>(next_segment_);
                float t        = fraction / frequency_;
                this_sample_ += ThisBlepSample(t) * discontinuity;
                next_sample_ += NextBlepSample(t) * discontinuity;
            }
        }
        segment_ = next_segment_;

        next_sample_ += (phase_ - 4.0f) * saw_8_gain_ * 0.125f;
        next_sample_ += (phase_ - float(segment_ & 4) - 2.0f) * saw_4_gain_ * 0.25f;
        next_sample_ += (phase_ - float(segment_ & 6) - 1.0f) * saw_2_gain_ * 0.5f;
        next_sample_ += (phase_ - float(segment_ & 7) - 0.5f) * saw_1_gain_;

        return 2.0f * this_sample_;
    }

    void OscillatorBank::SetFreq(float freq)
    {
        freq       = freq / sample_rate_;
        freq       = freq > 0.5f ? 0.5f : freq;
        recalc_    = cmp(freq, frequency_) || recalc_;
        frequency_ = freq;
    }

    void OscillatorBank::SetAmplitudes(const float* amplitudes)
    {
        for(int i = 0; i < 7; i++)
        {
            recalc_ = cmp(unshifted_registration_[i], amplitudes[i]) || recalc_;
            unshifted_registration_[i] = amplitudes[i];
        }
    }

    void OscillatorBank::SetSingleAmp(float amp, int idx)
    {
        if(idx < 0 || idx > 6)
        {
            return;
        }
        recalc_ = cmp(unshifted_registration_[idx], amp) || recalc_;
        unshifted_registration_[idx] = amp;
    }

    void OscillatorBank::SetGain(float gain)
    {
        gain         = gain > 1.f ? 1.f : gain;
        gain         = gain < 0.f ? 0.f : gain;
        recalc_gain_ = cmp(gain, gain_) || recalc_gain_;
        gain_        = gain;
    }

    /** Band Limited Oscillator

    Based on bltriangle, blsaw, blsquare from soundpipe

    Original Author(s): Paul Batchelor, saw2 Faust by Julius Smith

    Ported by Ben Sergentanis, May 2020
    */
    class BlOsc
    {
    public:
        BlOsc() {}
        ~BlOsc() {}
        /** Bl Waveforms
    */
        enum Waveforms
        {
            WAVE_TRIANGLE,
            WAVE_SAW,
            WAVE_SQUARE,
            WAVE_OFF,
        };


        /** -Initialize oscillator.
            -Defaults to: 440Hz, .5 amplitude, .5 pw, Triangle.
        */
        void Init(float sample_rate);


        /** - Get next floating point oscillator sample.
        */
        float Process();


        /** - Float freq: Set oscillator frequency in Hz.
        */
        inline void SetFreq(float freq) { freq_ = freq; };
        /** - Float amp: Set oscillator amplitude, 0 to 1.
        */
        inline void SetAmp(float amp) { amp_ = amp; };
        /** - Float pw: Set square osc pulsewidth, 0 to 1. (no thru 0 at the moment)
        */
        inline void SetPw(float pw) { pw_ = 1 - pw; };
        /** - uint8_t waveform: select between waveforms from enum above.
            - i.e. SetWaveform(BL_WAVEFORM_SAW); to set waveform to saw
        */
        inline void SetWaveform(uint8_t waveform) { mode_ = waveform; }

        /** - reset the phase of the oscillator.
        */
        void Reset();

    private:
        float rec0_[2], rec1_[2], vec0_[2], vec1_[2], vec2_[4096], freq_, amp_, pw_,
            sampling_freq_, half_sr_, quarter_sr_, sec_per_sample_, two_over_sr_,
            four_over_sr_;
        uint8_t mode_;
        int     iota_;

        float ProcessSquare();
        float ProcessTriangle();
        float ProcessSaw();
    };

    void BlOsc::Init(float sample_rate)
    {
        sampling_freq_  = sample_rate;
        half_sr_        = 0.5 * sampling_freq_;
        quarter_sr_     = sampling_freq_ * 0.25;
        sec_per_sample_ = 1.0 / sampling_freq_;
        two_over_sr_    = (float)(2.0 / sampling_freq_);
        four_over_sr_   = 4.0 / sampling_freq_;

        freq_ = 440;
        amp_  = 0.5;
        pw_   = 0.5;
        iota_ = 0;
        mode_ = WAVE_TRIANGLE;

        for(int i = 0; i < 2; i++)
        {
            rec0_[i] = rec1_[i] = vec0_[i] = vec1_[i] = 0.0;
        }

        for(int i = 0; i < 4096; i++)
        {
            vec2_[i] = 0.0;
        }
    }

    float BlOsc::ProcessSquare()
    {
        float out;
        float fSlow2 = fmin(2047.0, sampling_freq_ * (pw_ / freq_));
        float fSlow5 = (float)((int)fSlow2 + 1) - fSlow2;
        float fSlow6 = (quarter_sr_ / freq_);
        float fSlow7 = (sec_per_sample_ * freq_);
        float fSlow8 = fSlow2 - (int)fSlow2;

        rec0_[0]     = fmodf(rec0_[1] + fSlow7, 1.0);
        float fTemp0 = 2.0 * rec0_[0] - 1.0;
        fTemp0 *= fTemp0; //mult faster than fpow for squaring?
        vec1_[0]            = fTemp0;
        float fTemp1        = (fSlow6 * ((fTemp0 - vec1_[1])));
        vec2_[iota_ & 4095] = fTemp1;

        out = amp_
            * (0.0
                - ((fSlow5 * vec2_[(iota_ - (int)fSlow2) & 4095]
                    + fSlow8 * vec2_[(iota_ - ((int)fSlow2 + 1)) & 4095])
                    - fTemp1));
        rec0_[1] = rec0_[0];
        vec1_[1] = vec1_[0];
        iota_++;

        return out;
    }

    float BlOsc::ProcessTriangle()
    {
        float out;
        float fSlow1 = four_over_sr_ * (amp_ * freq_);
        float fSlow3 = half_sr_ / freq_;
        int   iSlow4 = (int)fSlow3;
        int   iSlow5 = 1 + iSlow4;
        float fSlow6 = iSlow5 - fSlow3;
        float fSlow7 = quarter_sr_ / freq_;
        float fSlow8 = sec_per_sample_ * freq_;
        float fSlow9 = fSlow3 - iSlow4; //decimal portion


        rec1_[0]     = fmodf((fSlow8 + rec1_[1]), 1.0);
        float fTemp0 = 2.0 * rec1_[0] - 1.0;
        fTemp0 *= fTemp0; //mult faster than fpow for squaring?
        vec1_[0]            = fTemp0;
        float fTemp1        = fSlow7 * (fTemp0 - vec1_[1]);
        vec2_[iota_ & 4095] = fTemp1;
        rec0_[0]            = 0.0
                - ((fSlow6 * vec2_[(iota_ - iSlow4) & 4095]
                    + fSlow9 * vec2_[(iota_ - iSlow5) & 4095])
                    - (.999 * rec0_[1] + fTemp1));

        out      = (float)(fSlow1 * rec0_[0]);
        rec1_[1] = rec1_[0];
        rec0_[1] = rec0_[0];
        vec1_[1] = vec1_[0];
        iota_++;

        return out;
    }

    float BlOsc::ProcessSaw()
    {
        float out;
        //fSlow0 = Slider1 = freq
        float fSlow1 = sampling_freq_ * (amp_ / freq_);
        float fSlow2 = (two_over_sr_ * freq_);
        float fSlow3 = (sampling_freq_ / freq_);

        rec0_[0]     = fmod((1.0 + rec0_[1]), fSlow3);
        float fTemp0 = fSlow2 * rec0_[0] - 1.0;
        fTemp0 *= fTemp0; //mult faster than fpow for squaring?
        vec0_[0] = fTemp0;
        vec1_[0] = 0.25;
        out      = (float)(fSlow1 * ((fTemp0 - vec0_[1]) * vec1_[1]));
        rec0_[1] = rec0_[0];
        vec0_[1] = vec0_[0];
        vec1_[1] = vec1_[0];

        return out;
    }

    void BlOsc::Reset()
    {
        iota_ = 0;

        for(int i = 0; i < 2; i++)
        {
            rec0_[i] = rec1_[i] = vec0_[i] = vec1_[i] = 0.0;
        }

        for(int i = 0; i < 4096; i++)
        {
            vec2_[i] = 0.0;
        }
    }

    float BlOsc::Process()
    {
        switch(mode_)
        {
            case WAVE_TRIANGLE: return ProcessTriangle();
            case WAVE_SAW: return ProcessSaw();
            case WAVE_SQUARE: return ProcessSquare();
            default: break;
        }

        return 0.0;
    }    

/** Simple 2 operator FM synth voice.
     
    Date: November, 2020

    Author: Ben Sergentanis
*/

class Fm2
{
  public:
    Fm2() {}
    ~Fm2() {}

    /** Initializes the FM2 module.
        \param samplerate - The sample rate of the audio engine being run. 
    */
    void Init(float samplerate);


    /**  Returns the next sample
    */
    float Process();

    /** Carrier freq. setter
        \param freq Carrier frequency in Hz
    */
    void SetFrequency(float freq);

    /** Set modulator freq. relative to carrier
        \param ratio New modulator freq = carrier freq. * ratio
    */
    void SetRatio(float ratio);

    /** Index setter
      \param FM depth, 5 = 2PI rads
  */
    void SetIndex(float index);

    /** Returns the current FM index. */
    float GetIndex();

    /** Resets both oscillators */
    void Reset();

  private:
    static constexpr float kIdxScalar      = 0.2f;
    static constexpr float kIdxScalarRecip = 1.f / kIdxScalar;

    Oscillator mod_, car_;
    float      idx_;
    float      freq_, lfreq_, ratio_, lratio_;
};


void Fm2::Init(float samplerate)
{
    //init oscillators
    car_.Init(samplerate);
    mod_.Init(samplerate);

    //set some reasonable values
    lfreq_  = 440.f;
    lratio_ = 2.f;
    SetFrequency(lfreq_);
    SetRatio(lratio_);

    car_.SetAmp(1.f);
    mod_.SetAmp(1.f);

    car_.SetWaveform(Oscillator::WAVE_SIN);
    mod_.SetWaveform(Oscillator::WAVE_SIN);

    idx_ = 1.f;
}

float Fm2::Process()
{
    if(lratio_ != ratio_ || lfreq_ != freq_)
    {
        lratio_ = ratio_;
        lfreq_  = freq_;
        car_.SetFreq(lfreq_);
        mod_.SetFreq(lfreq_ * lratio_);
    }

    float modval = mod_.Process();
    car_.PhaseAdd(modval * idx_);
    return car_.Process();
}

void Fm2::SetFrequency(float freq)
{
    freq_ = fabsf(freq);
}

void Fm2::SetRatio(float ratio)
{
    ratio_ = fabsf(ratio);
}

void Fm2::SetIndex(float index)
{
    idx_ = index * kIdxScalar;
}

float Fm2::GetIndex()
{
    return idx_ * kIdxScalarRecip;
}

void Fm2::Reset()
{
    car_.Reset();
    mod_.Reset();
}


/**  
       @brief Harmonic Oscillator Module based on Chebyshev polynomials. 
       @author Ben Sergentanis
       @date Dec 2020 
       Harmonic Oscillator Module based on Chebyshev polynomials \n 
       Works well for a small number of harmonics. For the higher order harmonics. \n
       We need to reinitialize the recurrence by computing two high harmonics. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/harmonic_oscillator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
template <int num_harmonics = 16>
class HarmonicOscillator
{
  public:
    HarmonicOscillator() {}
    ~HarmonicOscillator() {}

    /** Initialize harmonic oscillator
        \param sample_rate Audio engine samplerate
    */
    void Init(float sample_rate)
    {
        sample_rate_ = sample_rate;
        phase_       = 0.0f;

        for(int i = 0; i < num_harmonics; ++i)
        {
            amplitude_[i]    = 0.0f;
            newamplitude_[i] = 0.f;
        }
        amplitude_[0]    = 1.f;
        newamplitude_[0] = 1.f;

        SetFirstHarmIdx(1);
        SetFreq(440.f);

        recalc_ = false;
    }

    /** Get the next floating point sample */
    float Process()
    {
        if(recalc_)
        {
            recalc_ = false;
            for(int i = 0; i < num_harmonics; ++i)
            {
                float f = frequency_
                          * static_cast<float>(first_harmonic_index_ + i);
                if(f >= 0.5f)
                {
                    f = 0.5f;
                }
                amplitude_[i] = newamplitude_[i] * (1.0f - f * 2.0f);
            }
        }

        phase_ += frequency_;
        if(phase_ >= 1.0f)
        {
            phase_ -= 1.0f;
        }
        const float two_x = 2.0f * sinf(phase_ * TWOPI_F);
        float       previous, current;
        if(first_harmonic_index_ == 1)
        {
            previous = 1.0f;
            current  = two_x * 0.5f;
        }
        else
        {
            const float k = first_harmonic_index_;
            previous      = sinf((phase_ * (k - 1.0f) + 0.25f) * TWOPI_F);
            current       = sinf((phase_ * k) * TWOPI_F);
        }

        float sum = 0.0f;
        for(int i = 0; i < num_harmonics; ++i)
        {
            sum += amplitude_[i] * current;
            float temp = current;
            current    = two_x * current - previous;
            previous   = temp;
        }

        return sum;
    }

    /** Set the main frequency 
        \param freq Freq to be set in Hz.
    */
    void SetFreq(float freq)
    {
        //convert from Hz to phase_inc / sample
        freq       = freq / sample_rate_;
        freq       = freq >= .5f ? .5f : freq;
        freq       = freq <= -.5f ? -.5f : freq;
        recalc_    = cmp(freq, frequency_) || recalc_;
        frequency_ = freq;
    }

    /** Offset the set of harmonics. Passing in 3 means "harmonic 0" is the 3rd harm., 1 is the 4th, etc.
        \param idx Default behavior is 1. Values < 0 default to 1.
    */
    void SetFirstHarmIdx(int idx)
    {
        idx                   = idx < 1 ? 1 : idx;
        recalc_               = cmp(idx, first_harmonic_index_) || recalc_;
        first_harmonic_index_ = idx;
    }

    /** Set the amplitudes of each harmonic of the root. 
        \param amplitudes Amplitudes to set. Sum of all amplitudes must be < 1. The array referenced must be at least as large as num_harmonics. 
    */
    void SetAmplitudes(const float* amplitudes)
    {
        for(int i = 0; i < num_harmonics; i++)
        {
            recalc_          = cmp(newamplitude_[i], amplitudes[i]) || recalc_;
            newamplitude_[i] = amplitudes[i];
        }
    }

    /** Sets one amplitude. Does nothing if idx out of range.
        \param amp Amplitude to set
        \param idx Which harmonic to set.
    */
    void SetSingleAmp(const float amp, int idx)
    {
        if(idx < 0 || idx >= num_harmonics)
        {
            return;
        }
        recalc_            = cmp(amplitude_[idx], amp) || recalc_;
        newamplitude_[idx] = amp;
    }


  private:
    bool cmp(float a, float b) { return fabsf(a - b) > .000001f; }

    float sample_rate_;
    float phase_;
    float frequency_;
    float amplitude_[num_harmonics];
    float newamplitude_[num_harmonics];
    bool  recalc_;

    int first_harmonic_index_;
};



/**  
       @brief Formant Oscillator Module. 
       @author Ben Sergentanis
       @date Dec 2020 
       Sinewave with aliasing-free phase reset. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/formant_oscillator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class FormantOscillator
{
  public:
    FormantOscillator() {}
    ~FormantOscillator() {}

    /** Initializes the FormantOscillator module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);

    /** Get the next sample
    */
    float Process();

    /** Set the formant frequency.
        \param freq Frequency in Hz
    */
    void SetFormantFreq(float freq);

    /** Set the carrier frequency. This is the "main" frequency.
        \param freq Frequency in Hz
    */
    void SetCarrierFreq(float freq);

    /** Set the amount of phase shift
        \param ps Typically 0-1. Works with other values though, including negative.
    */
    void SetPhaseShift(float ps);

  private:
    inline float Sine(float phase);
    inline float ThisBlepSample(float t);
    inline float NextBlepSample(float t);

    // Oscillator state.
    float carrier_phase_;
    float formant_phase_;
    float next_sample_;

    // For interpolation of parameters.
    float carrier_frequency_;
    float formant_frequency_;
    float phase_shift_;
    float ps_inc_;

    float sample_rate_;

    //DISALLOW_COPY_AND_ASSIGN(FormantOscillator);
};

void FormantOscillator::Init(float sample_rate)
{
    carrier_phase_ = 0.0f;
    formant_phase_ = 0.0f;
    next_sample_   = 0.0f;

    carrier_frequency_ = 0.0f;
    formant_frequency_ = 100.f;
    phase_shift_       = 0.0f;

    sample_rate_ = sample_rate;
}

float FormantOscillator::Process()
{
    float this_sample = next_sample_;
    float next_sample = 0.0f;
    carrier_phase_ += carrier_frequency_;

    if(carrier_phase_ >= 1.0f)
    {
        carrier_phase_ -= 1.0f;
        float reset_time = carrier_phase_ / carrier_frequency_;

        float formant_phase_at_reset
            = formant_phase_ + (1.0f - reset_time) * formant_frequency_;
        float before        = Sine(formant_phase_at_reset + phase_shift_
                            + (ps_inc_ * (1.0f - reset_time)));
        float after         = Sine(0.0f + phase_shift_ + ps_inc_);
        float discontinuity = after - before;
        this_sample += discontinuity * ThisBlepSample(reset_time);
        next_sample += discontinuity * NextBlepSample(reset_time);
        formant_phase_ = reset_time * formant_frequency_;
    }
    else
    {
        formant_phase_ += formant_frequency_;
        if(formant_phase_ >= 1.0f)
        {
            formant_phase_ -= 1.0f;
        }
    }

    phase_shift_ += ps_inc_;
    ps_inc_ = 0.f;

    next_sample += Sine(formant_phase_ + phase_shift_);

    next_sample_ = next_sample;
    return this_sample;
}

void FormantOscillator::SetFormantFreq(float freq)
{
    //convert from Hz to phase_inc / sample
    formant_frequency_ = freq / sample_rate_;
    formant_frequency_ = formant_frequency_ >= .25f ? .25f : formant_frequency_;
    formant_frequency_
        = formant_frequency_ <= -.25f ? -.25f : formant_frequency_;
}

void FormantOscillator::SetCarrierFreq(float freq)
{
    //convert from Hz to phase_inc / sample
    carrier_frequency_ = freq / sample_rate_;
    carrier_frequency_ = carrier_frequency_ >= .25f ? .25f : carrier_frequency_;
    carrier_frequency_
        = carrier_frequency_ <= -.25f ? -.25f : carrier_frequency_;
}

void FormantOscillator::SetPhaseShift(float ps)
{
    ps_inc_ = ps - phase_shift_;
}

inline float FormantOscillator::Sine(float phase)
{
    return sinf(phase * TWOPI_F);
}

inline float FormantOscillator::ThisBlepSample(float t)
{
    return 0.5f * t * t;
}

inline float FormantOscillator::NextBlepSample(float t)
{
    t = 1.0f - t;
    return -0.5f * t * t;
}


/**  
     @brief Variable Saw Oscillator. 
    @author Ben Sergentanis
    @date Dec 2020 
    Saw with variable slope or notch. \n \n
    Ported from pichenettes/eurorack/plaits/dsp/oscillator/variable_saw_oscillator.h \n 
    \n to an independent module. \n
    Original code written by Emilie Gillet in 2016. \n
*/
class VariableSawOscillator
{
  public:
    VariableSawOscillator() {}
    ~VariableSawOscillator() {}

    void Init(float sample_rate);

    /** Get the next sample */
    float Process();

    /** Set master freq.
        \param frequency Freq in Hz.
    */
    void SetFreq(float frequency);

    /** Adjust the wave depending on the shape
        \param pw Notch or slope. Works best -1 to 1.
    */
    void SetPW(float pw);

    /** Slope or notch
        \param waveshape 0 = notch, 1 = slope
    */
    void SetWaveshape(float waveshape);


  private:
    float ComputeNaiveSample(float phase,
                             float pw,
                             float slope_up,
                             float slope_down,
                             float triangle_amount,
                             float notch_amount);

    float sample_rate_;

    // Oscillator state.
    float phase_;
    float next_sample_;
    float previous_pw_;
    bool  high_;

    const float kVariableSawNotchDepth = 0.2f;

    // For interpolation of parameters.
    float frequency_;
    float pw_;
    float waveshape_;
};

void VariableSawOscillator::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    phase_       = 0.0f;
    next_sample_ = 0.0f;
    previous_pw_ = 0.5f;
    high_        = false;

    SetFreq(220.f);
    SetPW(0.f);
    SetWaveshape(1.f);
}

float VariableSawOscillator::Process()
{
    float next_sample = next_sample_;

    float this_sample = next_sample;
    next_sample       = 0.0f;

    const float triangle_amount = waveshape_;
    const float notch_amount    = 1.0f - waveshape_;
    const float slope_up        = 1.0f / (pw_);
    const float slope_down      = 1.0f / (1.0f - pw_);

    phase_ += frequency_;

    if(!high_ && phase_ >= pw_)
    {
        const float triangle_step
            = (slope_up + slope_down) * frequency_ * triangle_amount;
        const float notch
            = (kVariableSawNotchDepth + 1.0f - pw_) * notch_amount;
        const float t = (phase_ - pw_) / (previous_pw_ - pw_ + frequency_);
        this_sample += notch * ThisBlepSample(t);
        next_sample += notch * NextBlepSample(t);
        this_sample -= triangle_step * ThisIntegratedBlepSample(t);
        next_sample -= triangle_step * NextIntegratedBlepSample(t);
        high_ = true;
    }
    else if(phase_ >= 1.0f)
    {
        phase_ -= 1.0f;
        const float triangle_step
            = (slope_up + slope_down) * frequency_ * triangle_amount;
        const float notch = (kVariableSawNotchDepth + 1.0f) * notch_amount;
        const float t     = phase_ / frequency_;
        this_sample -= notch * ThisBlepSample(t);
        next_sample -= notch * NextBlepSample(t);
        this_sample += triangle_step * ThisIntegratedBlepSample(t);
        next_sample += triangle_step * NextIntegratedBlepSample(t);
        high_ = false;
    }

    next_sample += ComputeNaiveSample(
        phase_, pw_, slope_up, slope_down, triangle_amount, notch_amount);
    previous_pw_ = pw_;

    next_sample_ = next_sample;
    return (2.0f * this_sample - 1.0f) / (1.0f + kVariableSawNotchDepth);
}

void VariableSawOscillator::SetFreq(float frequency)
{
    frequency  = frequency / sample_rate_;
    frequency  = frequency >= .25f ? .25f : frequency;
    pw_        = frequency >= .25f ? .5f : pw_;
    frequency_ = frequency;
}

void VariableSawOscillator::SetPW(float pw)
{
    if(frequency_ >= .25f)
    {
        pw_ = .5f;
    }
    else
    {
        pw_ = fclamp(pw, frequency_ * 2.0f, 1.0f - 2.0f * frequency_);
    }
}

void VariableSawOscillator::SetWaveshape(float waveshape)
{
    waveshape_ = waveshape;
}

float VariableSawOscillator::ComputeNaiveSample(float phase,
                                                float pw,
                                                float slope_up,
                                                float slope_down,
                                                float triangle_amount,
                                                float notch_amount)
{
    float notch_saw = phase < pw ? phase : 1.0f + kVariableSawNotchDepth;
    float triangle
        = phase < pw ? phase * slope_up : 1.0f - (phase - pw) * slope_down;
    return notch_saw * notch_amount + triangle * triangle_amount;
}


/**         
       @brief Variable Waveshape Oscillator. 
       @author Ben Sergentanis
       @date Dec 2020 
       Continuously variable waveform. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/variable_shape_oscillator.h \n 
       \n to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class VariableShapeOscillator
{
  public:
    VariableShapeOscillator() {}
    ~VariableShapeOscillator() {}

    /** Initialize the oscillator
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get next sample
    */
    float Process();

    /** Set master freq.
        \param frequency Freq in Hz.
    */
    void SetFreq(float frequency);

    /** Set pulse width / saw, ramp, tri.
        \param pw PW when shape is square. Saw, ramp, tri otherwise.
    */
    void SetPW(float pw);

    /** Switch from saw/ramp/tri to square.
        \param waveshape 0 is saw/ramp/tri, 1 is square.
    */
    void SetWaveshape(float waveshape);

    /** Whether or not to sync to the sync oscillator
        \param enable_sync True to turn sync on.
    */
    void SetSync(bool enable_sync);

    /** Set sync oscillator freq.
        \param frequency Freq in Hz.
    */
    void SetSyncFreq(float frequency);

  private:
    float ComputeNaiveSample(float phase,
                             float pw,
                             float slope_up,
                             float slope_down,
                             float triangle_amount,
                             float square_amount);

    float sample_rate_;
    bool  enable_sync_;

    // Oscillator state.
    float master_phase_;
    float slave_phase_;
    float next_sample_;
    float previous_pw_;
    bool  high_;

    // For interpolation of parameters.
    float master_frequency_;
    float slave_frequency_;
    float pw_;
    float waveshape_;
};

void VariableShapeOscillator::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    master_phase_ = 0.0f;
    slave_phase_  = 0.0f;
    next_sample_  = 0.0f;
    previous_pw_  = 0.5f;
    high_         = false;

    SetFreq(440.f);
    SetWaveshape(0.f);
    SetPW(0.f);
    SetSync(false);
    SetSyncFreq(220.f);
}

float VariableShapeOscillator::Process()
{
    float next_sample = next_sample_;

    bool  reset                   = false;
    bool  transition_during_reset = false;
    float reset_time              = 0.0f;

    float this_sample = next_sample;
    next_sample       = 0.0f;

    const float square_amount   = fmax(waveshape_ - 0.5f, 0.0f) * 2.0f;
    const float triangle_amount = fmax(1.0f - waveshape_ * 2.0f, 0.0f);
    const float slope_up        = 1.0f / (pw_);
    const float slope_down      = 1.0f / (1.0f - pw_);

    if(enable_sync_)
    {
        master_phase_ += master_frequency_;
        if(master_phase_ >= 1.0f)
        {
            master_phase_ -= 1.0f;
            reset_time = master_phase_ / master_frequency_;

            float slave_phase_at_reset
                = slave_phase_ + (1.0f - reset_time) * slave_frequency_;
            reset = true;
            if(slave_phase_at_reset >= 1.0f)
            {
                slave_phase_at_reset -= 1.0f;
                transition_during_reset = true;
            }
            if(!high_ && slave_phase_at_reset >= pw_)
            {
                transition_during_reset = true;
            }
            float value = ComputeNaiveSample(slave_phase_at_reset,
                                             pw_,
                                             slope_up,
                                             slope_down,
                                             triangle_amount,
                                             square_amount);
            this_sample -= value * ThisBlepSample(reset_time);
            next_sample -= value * NextBlepSample(reset_time);
        }
    }

    slave_phase_ += slave_frequency_;
    while(transition_during_reset || !reset)
    {
        if(!high_)
        {
            if(slave_phase_ < pw_)
            {
                break;
            }
            float t = (slave_phase_ - pw_)
                      / (previous_pw_ - pw_ + slave_frequency_);
            float triangle_step = (slope_up + slope_down) * slave_frequency_;
            triangle_step *= triangle_amount;

            this_sample += square_amount * ThisBlepSample(t);
            next_sample += square_amount * NextBlepSample(t);
            this_sample -= triangle_step * ThisIntegratedBlepSample(t);
            next_sample -= triangle_step * NextIntegratedBlepSample(t);
            high_ = true;
        }

        if(high_)
        {
            if(slave_phase_ < 1.0f)
            {
                break;
            }
            slave_phase_ -= 1.0f;
            float t             = slave_phase_ / slave_frequency_;
            float triangle_step = (slope_up + slope_down) * slave_frequency_;
            triangle_step *= triangle_amount;

            this_sample -= (1.0f - triangle_amount) * ThisBlepSample(t);
            next_sample -= (1.0f - triangle_amount) * NextBlepSample(t);
            this_sample += triangle_step * ThisIntegratedBlepSample(t);
            next_sample += triangle_step * NextIntegratedBlepSample(t);
            high_ = false;
        }
    }

    if(enable_sync_ && reset)
    {
        slave_phase_ = reset_time * slave_frequency_;
        high_        = false;
    }

    next_sample += ComputeNaiveSample(slave_phase_,
                                      pw_,
                                      slope_up,
                                      slope_down,
                                      triangle_amount,
                                      square_amount);
    previous_pw_ = pw_;


    next_sample_ = next_sample;
    return (2.0f * this_sample - 1.0f);
}

void VariableShapeOscillator::SetFreq(float frequency)
{
    frequency         = frequency / sample_rate_;
    frequency         = frequency >= .25f ? .25f : frequency;
    master_frequency_ = frequency;
}

void VariableShapeOscillator::SetPW(float pw)
{
    if(slave_frequency_ >= .25f)
    {
        pw_ = .5f;
    }
    else
    {
        pw_ = fclamp(
            pw, slave_frequency_ * 2.0f, 1.0f - 2.0f * slave_frequency_);
    }
}

void VariableShapeOscillator::SetWaveshape(float waveshape)
{
    waveshape_ = waveshape;
}

void VariableShapeOscillator::SetSync(bool enable_sync)
{
    enable_sync_ = enable_sync;
}

void VariableShapeOscillator::SetSyncFreq(float frequency)
{
    frequency        = frequency / sample_rate_;
    pw_              = frequency >= .25f ? .5f : pw_;
    frequency        = frequency >= .25f ? .25f : frequency;
    slave_frequency_ = frequency;
}

float VariableShapeOscillator::ComputeNaiveSample(float phase,
                                                  float pw,
                                                  float slope_up,
                                                  float slope_down,
                                                  float triangle_amount,
                                                  float square_amount)
{
    float saw    = phase;
    float square = phase < pw ? 0.0f : 1.0f;
    float triangle
        = phase < pw ? phase * slope_up : 1.0f - (phase - pw) * slope_down;
    saw += (square - saw) * square_amount;
    saw += (triangle - saw) * triangle_amount;
    return saw;
}

/**  
         @brief Vosim Oscillator Module \n 
       @author Ben Sergentanis
       @date Dec 2020 
       Two sinewaves multiplied by and sync'ed to a carrier. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/vosim_oscillator.h \n 
       \n to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class VosimOscillator
{
  public:
    VosimOscillator() {}
    ~VosimOscillator() {}

    /** Initializes the FormantOscillator module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);

    /** Get the next sample
    */
    float Process();

    /** Set carrier frequency.
        \param freq Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Set formant 1 frequency.
        \param freq Frequency in Hz.
    */
    void SetForm1Freq(float freq);

    /** Set formant 2 frequency.
        \param freq Frequency in Hz.
    */
    void SetForm2Freq(float freq);

    /** Waveshaping
        \param shape Shape to set. Works -1 to 1
    */
    void SetShape(float shape);

  private:
    float Sine(float phase);

    float sample_rate_;

    // Oscillator state.
    float carrier_phase_;
    float formant_1_phase_;
    float formant_2_phase_;

    // For interpolation of parameters.
    float carrier_frequency_;
    float formant_1_frequency_;
    float formant_2_frequency_;
    float carrier_shape_;
};

void VosimOscillator::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    carrier_phase_   = 0.0f;
    formant_1_phase_ = 0.0f;
    formant_2_phase_ = 0.0f;

    SetFreq(105.f);
    SetForm1Freq(1390.f);
    SetForm2Freq(817.f);
    SetShape(.5f);
}

float VosimOscillator::Process()
{
    carrier_phase_ += carrier_frequency_;
    if(carrier_phase_ >= 1.0f)
    {
        carrier_phase_ -= 1.0f;
        float reset_time = carrier_phase_ / carrier_frequency_;
        formant_1_phase_ = reset_time * formant_1_frequency_;
        formant_2_phase_ = reset_time * formant_2_frequency_;
    }
    else
    {
        formant_1_phase_ += formant_1_frequency_;
        if(formant_1_phase_ >= 1.0f)
        {
            formant_1_phase_ -= 1.0f;
        }
        formant_2_phase_ += formant_2_frequency_;
        if(formant_2_phase_ >= 1.0f)
        {
            formant_2_phase_ -= 1.0f;
        }
    }

    float carrier         = Sine(carrier_phase_ * 0.5f + 0.25f) + 1.0f;
    float reset_phase     = 0.75f - 0.25f * carrier_shape_;
    float reset_amplitude = Sine(reset_phase);
    float formant_0 = Sine(formant_1_phase_ + reset_phase) - reset_amplitude;
    float formant_1 = Sine(formant_2_phase_ + reset_phase) - reset_amplitude;
    return carrier * (formant_0 + formant_1) * 0.25f + reset_amplitude;
}

void VosimOscillator::SetFreq(float freq)
{
    carrier_frequency_ = freq / sample_rate_;
    carrier_frequency_ = carrier_frequency_ > .25f ? .25f : carrier_frequency_;
}

void VosimOscillator::SetForm1Freq(float freq)
{
    formant_1_frequency_ = freq / sample_rate_;
    formant_1_frequency_
        = formant_1_frequency_ > .25f ? .25f : formant_1_frequency_;
}

void VosimOscillator::SetForm2Freq(float freq)
{
    formant_2_frequency_ = freq / sample_rate_;
    formant_2_frequency_
        = formant_2_frequency_ > .25f ? .25f : formant_2_frequency_;
}

void VosimOscillator::SetShape(float shape)
{
    carrier_shape_ = shape;
}

float VosimOscillator::Sine(float phase)
{
    return sinf(TWOPI_F * phase);
}

/**  
    @brief     ZOscillator Module \n 
    @author Ben Sergentanis
    @date Dec 2020 
    Sinewave multiplied by and sync'ed to a carrier. \n \n
    Ported from pichenettes/eurorack/plaits/dsp/oscillator/z_oscillator.h \n 
    \n to an independent module. \n
    Original code written by Emilie Gillet in 2016. \n
*/
class ZOscillator
{
  public:
    ZOscillator() {}
    ~ZOscillator() {}

    /** Init ZOscillator module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get next sample
    */
    float Process();

    /** Set the carrier frequency
        \param freq Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Set the formant osc. freq
        \param freq Frequency in Hz.
    */
    void SetFormantFreq(float freq);

    /** Adjust the contour of the waveform.
        \param shape Waveshape to set. Works best 0-1.
    */
    void SetShape(float shape);

    /** Set the offset amount and phase shift. \n
        < 1/3 is just phase shift, > 2/3 is just offset, and between them is both. \n
        \param mode Mode to set. Works best -1 to 1
    */
    void SetMode(float mode);

  private:
    inline float Sine(float phase);
    inline float Z(float c, float d, float f, float shape, float mode);

    float sample_rate_;

    // Oscillator state.
    float carrier_phase_;
    float discontinuity_phase_;
    float formant_phase_;
    float next_sample_;

    // For interpolation of parameters.
    float carrier_frequency_;
    float formant_frequency_;
    float carrier_shape_, shape_new_;
    float mode_, mode_new_;
};

void ZOscillator::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    carrier_phase_       = 0.0f;
    discontinuity_phase_ = 0.0f;
    formant_phase_       = 0.0f;
    next_sample_         = 0.0f;

    SetFreq(220.f);
    SetFormantFreq(550.f);
    SetMode(0.f);
    SetShape(1.f);
}

float ZOscillator::Process()
{
    float next_sample = next_sample_;

    bool  reset      = false;
    float reset_time = 0.0f;

    float this_sample = next_sample;
    next_sample       = 0.0f;

    discontinuity_phase_ += 2.0f * carrier_frequency_;
    carrier_phase_ += carrier_frequency_;
    reset = discontinuity_phase_ >= 1.0f;

    if(reset)
    {
        discontinuity_phase_ -= 1.0f;
        reset_time = discontinuity_phase_ / (2.0f * carrier_frequency_);

        float carrier_phase_before = carrier_phase_ >= 1.0f ? 1.0f : 0.5f;
        float carrier_phase_after  = carrier_phase_ >= 1.0f ? 0.0f : 0.5f;

        float mode_sub  = mode_ + (1.f - reset_time) * (mode_ - mode_new_);
        float shape_sub = carrier_shape_
                          + (1.0f - reset_time) * (carrier_shape_ - shape_new_);
        float before
            = Z(carrier_phase_before,
                1.0f,
                formant_phase_ + (1.0f - reset_time) * formant_frequency_,
                shape_sub,
                mode_sub);

        float after = Z(carrier_phase_after, 0.0f, 0.0f, shape_new_, mode_new_);

        float discontinuity = after - before;
        this_sample += discontinuity * ThisBlepSample(reset_time);
        next_sample += discontinuity * NextBlepSample(reset_time);
        formant_phase_ = reset_time * formant_frequency_;

        if(carrier_phase_ > 1.0f)
        {
            carrier_phase_ = discontinuity_phase_ * 0.5f;
        }
    }
    else
    {
        formant_phase_ += formant_frequency_;
        if(formant_phase_ >= 1.0f)
        {
            formant_phase_ -= 1.0f;
        }
    }

    if(carrier_phase_ >= 1.0f)
    {
        carrier_phase_ -= 1.0f;
    }

    carrier_shape_ = shape_new_;
    mode_          = mode_new_;
    next_sample += Z(carrier_phase_,
                     discontinuity_phase_,
                     formant_phase_,
                     carrier_shape_,
                     mode_);

    next_sample_ = next_sample;
    return this_sample;
}

inline float ZOscillator::Sine(float phase)
{
    return sinf(phase * TWOPI_F);
}

void ZOscillator::SetFreq(float freq)
{
    //convert from Hz to phase_inc / sample
    carrier_frequency_ = freq / sample_rate_;
    carrier_frequency_ = carrier_frequency_ >= .25f ? .25f : carrier_frequency_;
}

void ZOscillator::SetFormantFreq(float freq)
{
    //convert from Hz to phase_inc / sample
    formant_frequency_ = freq / sample_rate_;
    formant_frequency_ = formant_frequency_ >= .25f ? .25f : formant_frequency_;
}

void ZOscillator::SetShape(float shape)
{
    shape_new_ = shape;
}

void ZOscillator::SetMode(float mode)
{
    mode_new_ = mode;
}

inline float ZOscillator::Z(float c, float d, float f, float shape, float mode)
{
    float ramp_down = 0.5f * (1.0f + Sine(0.5f * d + 0.25f));

    float offset;
    float phase_shift;
    if(mode < 0.333f)
    {
        offset      = 1.0f;
        phase_shift = 0.25f + mode * 1.50f;
    }
    else if(mode < 0.666f)
    {
        phase_shift = 0.7495f - (mode - 0.33f) * 0.75f;
        offset      = -Sine(phase_shift);
    }
    else
    {
        phase_shift = 0.7495f - (mode - 0.33f) * 0.75f;
        offset      = 0.001f;
    }

    float discontinuity = Sine(f + phase_shift);
    float contour;
    if(shape < 0.5f)
    {
        shape *= 2.0f;
        if(c >= 0.5f)
        {
            ramp_down *= shape;
        }
        contour = 1.0f + (Sine(c + 0.25f) - 1.0f) * shape;
    }
    else
    {
        contour = Sine(c + shape * 0.5f);
    }
    return (ramp_down * (offset + discontinuity) - offset) * contour;
}

/**  
       Allpass filter module \n 
       Passes all frequencies at their original levels, with a phase shift. \n 
       Ported from soundpipe by Ben Sergentanis, May 2020 
       @author Barry Vercoe, John ffitch
       @date 1991
*/
class Allpass
{
  public:
    Allpass() {}
    ~Allpass() {}

    /**         
        Initializes the allpass module.
        \param sample_rate  The sample rate of the audio engine being run.
    \param buff Buffer for allpass to use.
    \param size Size of buff.
    */
    void Init(float sample_rate, float* buff, size_t size);

    /** 
     \param in Input sample.
     \return Next floating point sample.
    */
    float Process(float in);

    /**
       Sets the filter frequency (Implemented by delay time).
       \param looptime Filter looptime in seconds.
    */
    void SetFreq(float looptime);

    /**
        \param revtime Reverb time in seconds.
    */
    inline void SetRevTime(float revtime) { rev_time_ = revtime; }


  private:
    float  sample_rate_, rev_time_, loop_time_, prvt_, coef_, max_loop_time_;
    float* buf_;
    int    buf_pos_, mod_;
};

void Allpass::Init(float sample_rate, float* buff, size_t size)
{
    sample_rate_   = sample_rate;
    rev_time_      = 3.5;
    max_loop_time_ = ((float)size / sample_rate_) - .01;
    loop_time_     = max_loop_time_;
    mod_           = (int)(loop_time_ * sample_rate_);
    buf_           = buff;
    prvt_          = 0.0;
    coef_          = 0.0;
    buf_pos_       = 0;
}

float Allpass::Process(float in)
{
    float y, z, out;
    if(prvt_ != rev_time_)
    {
        prvt_ = rev_time_;
        coef_ = expf(-6.9078 * loop_time_ / prvt_);
    }

    y              = buf_[buf_pos_];
    z              = coef_ * y + in;
    buf_[buf_pos_] = z;
    out            = y - coef_ * z;

    buf_pos_++;
    buf_pos_ %= mod_;
    return out;
}

void Allpass::SetFreq(float freq)
{
    loop_time_ = fmaxf(fminf(freq, max_loop_time_), .0001);
    mod_       = fmaxf(loop_time_ * sample_rate_, 0);
}

/** A first-order recursive high-pass filter with variable frequency response.
     Original Author(s): Barry Vercoe, John FFitch, Gabriel Maldonado

     Year: 1991

     Original Location: Csound -- OOps/ugens5.c

     Ported from soundpipe by Ben Sergentanis, May 2020
    */
class ATone
{
  public:
    ATone() {}
    ~ATone() {}
    /** Initializes the ATone module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Processes one sample through the filter and returns one sample.
        \param in - input signal 
    */
    float Process(float &in);

    /** Sets the cutoff frequency or half-way point of the filter.
        \param freq - frequency value in Hz. Range: Any positive value.
    */
    inline void SetFreq(float &freq)
    {
        freq_ = freq;
        CalculateCoefficients();
    }

    /** get current frequency
        \return the current value for the cutoff frequency or half-way point of the filter.
    */
    inline float GetFreq() { return freq_; }

  private:
    void  CalculateCoefficients();
    float out_, prevout_, in_, freq_, c2_, sample_rate_;
};

void ATone::Init(float sample_rate)
{
    prevout_     = 0.0f;
    freq_        = 1000.0f;
    c2_          = 0.5f;
    sample_rate_ = sample_rate;
}

float ATone::Process(float &in)
{
    float out;

    out      = c2_ * (prevout_ + in);
    prevout_ = out - in;

    return out;
}

void ATone::CalculateCoefficients()
{
    float b, c2;

    b   = 2.0f - cosf(TWOPI_F * freq_ / sample_rate_);
    c2  = b - sqrtf(b * b - 1.0f);
    c2_ = c2;
}

/** Two pole recursive filter

    Original author(s) : Hans Mikelson

    Year: 1998

    Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Biquad
{
  public:
    Biquad() {}
    ~Biquad() {}
    /** Initializes the biquad module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Filters the input signal
        \return filtered output
    */
    float Process(float in);


    /** Sets resonance amount
        \param res : Set filter resonance.
    */
    inline void SetRes(float res)
    {
        res_ = res;
        Reset();
    }


    /** Sets filter cutoff in Hz
        \param cutoff : Set filter cutoff.
    */
    inline void SetCutoff(float cutoff)
    {
        cutoff_ = cutoff;
        Reset();
    }

  private:
    float sample_rate_, cutoff_, res_, b0_, b1_, b2_, a0_, a1_, a2_,
        two_pi_d_sr_, xnm1_, xnm2_, ynm1_, ynm2_;
    void Reset();
};
void Biquad::Reset()
{
    float con   = cutoff_ * two_pi_d_sr_;
    float alpha = 1.0f - 2.0f * res_ * cosf(con) * cosf(con)
                  + res_ * res_ * cosf(2 * con);
    float beta  = 1.0f + cosf(con);
    float gamma = 1 + cosf(con);
    float m1    = alpha * gamma + beta * sinf(con);
    float m2    = alpha * gamma - beta * sinf(con);
    float den   = sqrtf(m1 * m1 + m2 * m2);

    b0_ = 1.5f * (alpha * alpha + beta * beta) / den;
    b1_ = b0_;
    b2_ = 0.0f;
    a0_ = 1.0f;
    a1_ = -2.0 * res_ * cosf(con);
    a2_ = res_ * res_;
}

void Biquad::Init(float sample_rate)
{
    sample_rate_ = sample_rate;
    two_pi_d_sr_ = TWOPI_F / sample_rate_;

    cutoff_ = 500;
    res_    = 0.7;

    Reset();

    xnm1_ = xnm2_ = ynm1_ = ynm2_ = 0.0f;
}

float Biquad::Process(float in)
{
    float xn, yn;
    float a0 = a0_, a1 = a1_, a2 = a2_;
    float b0 = b0_, b1 = b1_, b2 = b2_;
    float xnm1 = xnm1_, xnm2 = xnm2_, ynm1 = ynm1_, ynm2 = ynm2_;

    xn   = in;
    yn   = (b0 * xn + b1 * xnm1 + b2 * xnm2 - a1 * ynm1 - a2 * ynm2) / a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;

    xnm1_ = xnm1;
    xnm2_ = xnm2;
    ynm1_ = ynm1;
    ynm2_ = ynm2;

    return yn;
}

/** Comb filter module

    Original author(s) :

    Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Comb
{
  public:
    Comb() {}
    ~Comb() {}

    /** Initializes the Comb module.
        \param sample_rate - The sample rate of the audio engine being run. 
        \param buff - input buffer, kept in either main() or global space
        \param size - size of buff
    */
    void Init(float sample_rate, float* buff, size_t size);


    /** processes the comb filter
    */
    float Process(float in);


    /** Sets the period of the comb filter in seconds
    */
    void SetPeriod(float looptime);

    /** Sets the frequency of the comb filter in Hz
    */
    inline void SetFreq(float freq)
    {
        if(freq > 0)
        {
            SetPeriod(1.f / freq);
        }
    }

    /** Sets the decay time of the comb filter
    */
    inline void SetRevTime(float revtime) { rev_time_ = revtime; }

  private:
    float  sample_rate_, rev_time_, loop_time_, prvt_, coef_, max_loop_time_;
    float* buf_;
    size_t buf_pos_, mod_, max_size_;
};

static float log001 = -6.9078f; // log .001

void Comb::Init(float sample_rate, float* buff, size_t size)
{
    sample_rate_   = sample_rate;
    rev_time_      = 3.5;
    max_size_      = size;
    max_loop_time_ = ((float)size / sample_rate_) - .01;
    loop_time_     = max_loop_time_;
    mod_           = sample_rate_ * loop_time_;
    buf_           = buff;
    prvt_          = 0.0f;
    coef_          = 0.0f;
    buf_pos_       = 0;
}

float Comb::Process(float in)
{
    float tmp     = 0;
    float coef    = coef_;
    float outsamp = 0;

    if(prvt_ != rev_time_)
    {
        prvt_         = rev_time_;
        float exp_arg = (float)(log001 * loop_time_ / prvt_);
        if(exp_arg < -36.8413615)
        {
            coef = coef_ = 0;
        }
        else
        {
            coef = coef_ = expf(exp_arg);
        }
    }

    // internal delay line
    outsamp                = buf_[(buf_pos_ + mod_) % max_size_];
    tmp                    = (outsamp * coef) + in;
    buf_[(size_t)buf_pos_] = tmp;
    buf_pos_               = (buf_pos_ - 1 + max_size_) % max_size_;

    return outsamp;
}

void Comb::SetPeriod(float looptime)
{
    if(looptime > 0)
    {
        loop_time_ = fminf(looptime, max_loop_time_);
        mod_       = loop_time_ * sample_rate_;
        if(mod_ > max_size_)
        {
            mod_ = max_size_ - 1;
        }
    }
}

/** Resonant Modal Filter

    Extracted from soundpipe to work as a Daisy Module,

    originally extracted from csound by Paul Batchelor.

    Original Author(s): Francois Blanc, Steven Yi

    Year: 2001

    Location: Opcodes/biquad.c (csound)
*/
class Mode
{
  public:
    Mode() {}
    ~Mode() {}
    /** Initializes the instance of the module.
        sample_rate: frequency of the audio engine in Hz
    */
    void Init(float sample_rate);

    /** Processes one input sample through the filter, and returns the output.
    */
    float Process(float in);

    /** Clears the filter, returning the output to 0.0
    */
    void Clear();

    /** Sets the resonant frequency of the modal filter.
        Range: Any frequency such that sample_rate / freq < PI (about 15.2kHz at 48kHz)
    */
    inline void SetFreq(float freq) { freq_ = freq; }
    /** Sets the quality factor of the filter.
        Range: Positive Numbers (Good values range from 70 to 1400)
    */
    inline void SetQ(float q) { q_ = q; }

  private:
    float freq_, q_;
    float xnm1_, ynm1_, ynm2_, a0_, a1_, a2_;
    float d_, lfq_, lq_, sr_;
};





void Mode::Init(float sample_rate)
{
    freq_ = 500.0f;
    q_    = 50;

    xnm1_ = ynm1_ = ynm2_ = 0.0f;
    a0_ = a1_ = a2_ = d_ = 0.0f;
    lfq_ = lq_ = -1.0f;
    sr_        = sample_rate;
}

void Mode::Clear()
{
    xnm1_ = ynm1_ = ynm2_ = 0.0f;
    a0_ = a1_ = a2_ = 0.0f;
    d_              = 0.0f;
    lfq_            = -1.0f;
    lq_             = -1.0f;
}

float Mode::Process(float in)
{
    float out;
    float lfq = lfq_, lq = lq_;
    float xn, yn, a0 = a0_, a1 = a1_, a2 = a2_, d = d_;
    float xnm1 = xnm1_, ynm1 = ynm1_, ynm2 = ynm2_;
    float kfq = freq_;
    float kq  = q_;

    if(lfq != kfq || lq != kq)
    {
        float kfreq  = kfq * (2.0f * (float)M_PI);
        float kalpha = (sr_ / kfreq);
        float kbeta  = kalpha * kalpha;
        d            = 0.5f * kalpha;
        lq           = kq;
        lfq          = kfq;
        a0           = 1.0f / (kbeta + d / kfreq);
        a1           = a0 * (1.0f - 2.0f * kbeta);
        a2           = a0 * (kbeta - d / kq);
    }
    xn = in;

    yn = a0 * xnm1 - a1 * ynm1 - a2 * ynm2;

    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;

    yn = yn * d;

    out   = yn;
    xnm1_ = xnm1;
    ynm1_ = ynm1;
    ynm2_ = ynm2;
    lfq_  = lfq;
    lq_   = lq;
    d_    = d;
    a0_   = a0;
    a1_   = a1;
    a2_   = a2;
    return out;
}

/** Moog ladder filter module

Ported from soundpipe

Original author(s) : Victor Lazzarini, John ffitch (fast tanh), Bob Moog

*/
class MoogLadder
{
  public:
    MoogLadder() {}
    ~MoogLadder() {}
    /** Initializes the MoogLadder module.
        sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Processes the lowpass filter
    */
    float Process(float in);

    /** 
        Sets the cutoff frequency or half-way point of the filter.
        Arguments
        - freq - frequency value in Hz. Range: Any positive value.
    */
    inline void SetFreq(float freq) { freq_ = freq; }
    /** 
        Sets the resonance of the filter.
    */
    inline void SetRes(float res) { res_ = res; }

  private:
    float istor_, res_, freq_, delay_[6], tanhstg_[3], old_freq_, old_res_,
        sample_rate_, old_acr_, old_tune_;
    float my_tanh(float x);
};

float MoogLadder::my_tanh(float x)
{
    int sign = 1;
    if(x < 0)
    {
        sign = -1;
        x    = -x;
        return x * sign;
    }
    else if(x >= 4.0f)
    {
        return sign;
    }
    else if(x < 0.5f)
    {
        return x * sign;
    }
    return sign * tanhf(x);
}

void MoogLadder::Init(float sample_rate)
{
    sample_rate_ = sample_rate;
    istor_       = 0.0f;
    res_         = 0.4f;
    freq_        = 1000.0f;

    for(int i = 0; i < 6; i++)
    {
        delay_[i]       = 0.0;
        tanhstg_[i % 3] = 0.0;
    }

    old_freq_ = 0.0f;
    old_res_  = -1.0f;
}

float MoogLadder::Process(float in)
{
    float  freq = freq_;
    float  res  = res_;
    float  res4;
    float* delay   = delay_;
    float* tanhstg = tanhstg_;
    float  stg[4];
    float  acr, tune;

    float THERMAL = 0.000025;

    if(res < 0)
    {
        res = 0;
    }

    if(old_freq_ != freq || old_res_ != res)
    {
        float f, fc, fc2, fc3, fcr;
        old_freq_ = freq;
        fc        = (freq / sample_rate_);
        f         = 0.5f * fc;
        fc2       = fc * fc;
        fc3       = fc2 * fc2;

        fcr  = 1.8730f * fc3 + 0.4955f * fc2 - 0.6490f * fc + 0.9988f;
        acr  = -3.9364f * fc2 + 1.8409f * fc + 0.9968f;
        tune = (1.0f - expf(-((2 * PI_F) * f * fcr))) / THERMAL;

        old_res_  = res;
        old_acr_  = acr;
        old_tune_ = tune;
    }
    else
    {
        res  = old_res_;
        acr  = old_acr_;
        tune = old_tune_;
    }

    res4 = 4.0f * res * acr;

    for(int j = 0; j < 2; j++)
    {
        in -= res4 * delay[5];
        delay[0] = stg[0]
            = delay[0] + tune * (my_tanh(in * THERMAL) - tanhstg[0]);
        for(int k = 1; k < 4; k++)
        {
            in     = stg[k - 1];
            stg[k] = delay[k]
                     + tune
                           * ((tanhstg[k - 1] = my_tanh(in * THERMAL))
                              - (k != 3 ? tanhstg[k]
                                        : my_tanh(delay[k] * THERMAL)));
            delay[k] = stg[k];
        }
        delay[5] = (stg[3] + delay[4]) * 0.5f;
        delay[4] = stg[3];
    }
    return delay[5];
}


/** Non-linear filter
    
    port by: Stephen Hensley, December 2019
    
    The four 5-coefficients: a, b, d, C, and L are used to configure different filter types.

    Structure for Dobson/Fitch nonlinear filter 

    Revised Formula from Risto Holopainen 12 Mar 2004

    `Y{n} =tanh(a Y{n-1} + b Y{n-2} + d Y^2{n-L} + X{n} - C)`

    Though traditional filter types can be made, 
    the effect will always respond differently to different input.

    This Source is a heavily modified version of the original
    source from Csound.

    \todo make this work on a single sample instead of just on blocks at a time.
    */
class NlFilt
{
  public:
    /** Initializes the NlFilt object.
        */
    void Init();

    /** Process the array pointed to by \*in and updates the output to \*out;
        This works on a block of audio at once, the size of which is set with the size. 
        */
    void ProcessBlock(float *in, float *out, size_t size);


    /** inputs these are the five coefficients for the filter.
        */
    inline void SetCoefficients(float a, float b, float d, float C, float L)
    {
        a_ = a;
        b_ = b;
        d_ = d;
        C_ = C;
        L_ = L;
    }


    /** Set Coefficient a
        */
    inline void SetA(float a) { a_ = a; }
    /** Set Coefficient b
        */
    inline void SetB(float b) { b_ = b; }
    /** Set Coefficient d
        */
    inline void SetD(float d) { d_ = d; }
    /** Set Coefficient C
        */
    inline void SetC(float C) { C_ = C; }
    /** Set Coefficient L
        */
    inline void SetL(float L) { L_ = L; }

  private:
    int32_t Set();

    float   in_, a_, b_, d_, C_, L_;
    float   delay_[DSY_NLFILT_MAX_DELAY];
    int32_t point_;
};

/**      Double Sampled, Stable State Variable Filter

Credit to Andrew Simper from musicdsp.org

This is his "State Variable Filter (Double Sampled, Stable)"

Additional thanks to Laurent de Soras for stability limit, and 
Stefan Diedrichsen for the correct notch output

Ported by: Stephen Hensley
*/
class Svf
{
  public:
    Svf() {}
    ~Svf() {}
    /** Initializes the filter
        float sample_rate - sample rate of the audio engine being run, and the frequency that the Process function will be called.
    */
    void Init(float sample_rate);


    /** 
        Process the input signal, updating all of the outputs.
    */
    void Process(float in);


    /** sets the frequency of the cutoff frequency. 
        f must be between 0.0 and sample_rate / 3
    */
    void SetFreq(float f);

    /** sets the resonance of the filter.
        Must be between 0.0 and 1.0 to ensure stability.
    */
    void SetRes(float r);

    /** sets the drive of the filter 
        affects the response of the resonance of the filter
    */
    void SetDrive(float d);
    /** lowpass output
        \return low pass output of the filter
    */
    inline float Low() { return out_low_; }
    /** highpass output
        \return high pass output of the filter
    */
    inline float High() { return out_high_; }
    /** bandpass output
        \return band pass output of the filter
    */
    inline float Band() { return out_band_; }
    /** notchpass output
        \return notch pass output of the filter
    */
    inline float Notch() { return out_notch_; }
    /** peak output
        \return peak output of the filter
    */
    inline float Peak() { return out_peak_; }

  private:
    float sr_, fc_, res_, drive_, freq_, damp_;
    float notch_, low_, high_, band_, peak_;
    float input_;
    float out_low_, out_high_, out_band_, out_peak_, out_notch_;
    float pre_drive_, fc_max_;
};



void Svf::Init(float sample_rate)
{
    sr_        = sample_rate;
    fc_        = 200.0f;
    res_       = 0.5f;
    drive_     = 0.5f;
    pre_drive_ = 0.5f;
    freq_      = 0.25f;
    damp_      = 0.0f;
    notch_     = 0.0f;
    low_       = 0.0f;
    high_      = 0.0f;
    band_      = 0.0f;
    peak_      = 0.0f;
    input_     = 0.0f;
    out_notch_ = 0.0f;
    out_low_   = 0.0f;
    out_high_  = 0.0f;
    out_peak_  = 0.0f;
    out_band_  = 0.0f;
    fc_max_    = sr_ / 3.f;
}

void Svf::Process(float in)
{
    input_ = in;
    // first pass
    notch_ = input_ - damp_ * band_;
    low_   = low_ + freq_ * band_;
    high_  = notch_ - low_;
    band_  = freq_ * high_ + band_ - drive_ * band_ * band_ * band_;
    // take first sample of output
    out_low_   = 0.5f * low_;
    out_high_  = 0.5f * high_;
    out_band_  = 0.5f * band_;
    out_peak_  = 0.5f * (low_ - high_);
    out_notch_ = 0.5f * notch_;
    // second pass
    notch_ = input_ - damp_ * band_;
    low_   = low_ + freq_ * band_;
    high_  = notch_ - low_;
    band_  = freq_ * high_ + band_ - drive_ * band_ * band_ * band_;
    // average second pass outputs
    out_low_ += 0.5f * low_;
    out_high_ += 0.5f * high_;
    out_band_ += 0.5f * band_;
    out_peak_ += 0.5f * (low_ - high_);
    out_notch_ += 0.5f * notch_;
}

void Svf::SetFreq(float f)
{
    fc_ = fclamp(f, 1.0e-6, fc_max_);
    // Set Internal Frequency for fc_
    freq_ = 2.0f
            * sinf(PI_F
                   * MIN(0.25f,
                         fc_ / (sr_ * 2.0f))); // fs*2 because double sampled
    // recalculate damp
    damp_ = MIN(2.0f * (1.0f - powf(res_, 0.25f)),
                MIN(2.0f, 2.0f / freq_ - freq_ * 0.5f));
}

void Svf::SetRes(float r)
{
    float res = fclamp(r, 0.f, 1.f);
    res_      = res;
    // recalculate damp
    damp_  = MIN(2.0f * (1.0f - powf(res_, 0.25f)),
                MIN(2.0f, 2.0f / freq_ - freq_ * 0.5f));
    drive_ = pre_drive_ * res_;
}

void Svf::SetDrive(float d)
{
    float drv  = fclamp(d * 0.1f, 0.f, 1.f);
    pre_drive_ = drv;
    drive_     = pre_drive_ * res_;
}

/** A first-order recursive low-pass filter with variable frequency response.
*/
class Tone
{
  public:
    Tone() {}
    ~Tone() {}
    /** Initializes the Tone module.
        sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Processes one sample through the filter and returns one sample.
        in - input signal 
    */
    float Process(float &in);

    /** Sets the cutoff frequency or half-way point of the filter.

        \param freq - frequency value in Hz. Range: Any positive value.
    */
    inline void SetFreq(float &freq)
    {
        freq_ = freq;
        CalculateCoefficients();
    }

    /** 
        \return the current value for the cutoff frequency or half-way point of the filter.
    */
    inline float GetFreq() { return freq_; }

  private:
    void  CalculateCoefficients();
    float out_, prevout_, in_, freq_, c1_, c2_, sample_rate_;
};


/**  
       @brief 808 bass drum model, revisited.
       @author Ben Sergentanis
       @date Jan 2021 
       Ported from pichenettes/eurorack/plaits/dsp/drums/analog_bass_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class AnalogBassDrum
{
  public:
    AnalogBassDrum() {}
    ~AnalogBassDrum() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next sample.
        \param trigger True strikes the drum. Defaults to false.
    */
    float Process(bool trigger = false);

    /** Strikes the drum. */
    void Trig();

    /** Set the bassdrum to play infinitely
        \param sustain True = infinite length
    */
    void SetSustain(bool sustain);

    /** Set a small accent.
        \param accent Works 0-1
    */
    void SetAccent(float accent);

    /** Set the drum's root frequency
        \param f0 Frequency in Hz
    */
    void SetFreq(float f0);

    /** Set the amount of click.
        \param tone Works 0-1.
    */
    void SetTone(float tone);

    /** Set the decay length of the drum.
        \param decay Works best 0-1.
    */
    void SetDecay(float decay);

    /** Set the amount of fm attack. Works together with self fm.
        \param attack_fm_amount Works best 0-1.
    */
    void SetAttackFmAmount(float attack_fm_amount);

    /**Set the amount of self fm. Also affects fm attack, and volume decay.
        \param self_fm_amount Works best 0-1.
    */
    void SetSelfFmAmount(float self_fm_amount);

  private:
    inline float Diode(float x);

    float sample_rate_;

    float accent_, f0_, tone_, decay_;
    float attack_fm_amount_, self_fm_amount_;

    bool trig_, sustain_;

    int   pulse_remaining_samples_;
    int   fm_pulse_remaining_samples_;
    float pulse_;
    float pulse_height_;
    float pulse_lp_;
    float fm_pulse_lp_;
    float retrig_pulse_;
    float lp_out_;
    float tone_lp_;
    float sustain_gain_;

    Svf resonator_;

    //for use in sin + cos osc. in sustain mode
    float phase_;
};
void AnalogBassDrum::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    trig_ = false;

    pulse_remaining_samples_    = 0;
    fm_pulse_remaining_samples_ = 0;
    pulse_                      = 0.0f;
    pulse_height_               = 0.0f;
    pulse_lp_                   = 0.0f;
    fm_pulse_lp_                = 0.0f;
    retrig_pulse_               = 0.0f;
    lp_out_                     = 0.0f;
    tone_lp_                    = 0.0f;
    sustain_gain_               = 0.0f;
    phase_                      = 0.f;


    SetSustain(false);
    SetAccent(.1f);
    SetFreq(50.f);
    SetTone(.1f);
    SetDecay(.3f);
    SetSelfFmAmount(1.f);
    SetAttackFmAmount(.5f);

    resonator_.Init(sample_rate_);
}

inline float AnalogBassDrum::Diode(float x)
{
    if(x >= 0.0f)
    {
        return x;
    }
    else
    {
        x *= 2.0f;
        return 0.7f * x / (1.0f + fabsf(x));
    }
}

float AnalogBassDrum::Process(bool trigger)
{
    const int kTriggerPulseDuration  = static_cast<int>(1.0e-3f * sample_rate_);
    const int kFMPulseDuration       = static_cast<int>(6.0e-3f * sample_rate_);
    const float kPulseDecayTime      = 0.2e-3f * sample_rate_;
    const float kPulseFilterTime     = 0.1e-3f * sample_rate_;
    const float kRetrigPulseDuration = 0.05f * sample_rate_;

    const float scale = 0.001f / f0_;
    const float q     = 1500.0f * powf(2.f, kOneTwelfth * decay_ * 80.0f);
    const float tone_f
        = fmin(4.0f * f0_ * powf(2.f, kOneTwelfth * tone_ * 108.0f), 1.0f);
    const float exciter_leak = 0.08f * (tone_ + 0.25f);


    if(trigger || trig_)
    {
        trig_ = false;

        pulse_remaining_samples_    = kTriggerPulseDuration;
        fm_pulse_remaining_samples_ = kFMPulseDuration;
        pulse_height_               = 3.0f + 7.0f * accent_;
        lp_out_                     = 0.0f;
    }

    // Q39 / Q40
    float pulse = 0.0f;
    if(pulse_remaining_samples_)
    {
        --pulse_remaining_samples_;
        pulse = pulse_remaining_samples_ ? pulse_height_ : pulse_height_ - 1.0f;
        pulse_ = pulse;
    }
    else
    {
        pulse_ *= 1.0f - 1.0f / kPulseDecayTime;
        pulse = pulse_;
    }
    if(sustain_)
    {
        pulse = 0.0f;
    }

    // C40 / R163 / R162 / D83
    fonepole(pulse_lp_, pulse, 1.0f / kPulseFilterTime);
    pulse = Diode((pulse - pulse_lp_) + pulse * 0.044f);

    // Q41 / Q42
    float fm_pulse = 0.0f;
    if(fm_pulse_remaining_samples_)
    {
        --fm_pulse_remaining_samples_;
        fm_pulse = 1.0f;
        // C39 / C52
        retrig_pulse_ = fm_pulse_remaining_samples_ ? 0.0f : -0.8f;
    }
    else
    {
        // C39 / R161
        retrig_pulse_ *= 1.0f - 1.0f / kRetrigPulseDuration;
    }
    if(sustain_)
    {
        fm_pulse = 0.0f;
    }
    fonepole(fm_pulse_lp_, fm_pulse, 1.0f / kPulseFilterTime);

    // Q43 and R170 leakage
    float punch = 0.7f + Diode(10.0f * lp_out_ - 1.0f);

    // Q43 / R165
    float attack_fm = fm_pulse_lp_ * 1.7f * attack_fm_amount_;
    float self_fm   = punch * 0.08f * self_fm_amount_;
    float f         = f0_ * (1.0f + attack_fm + self_fm);
    f               = fclamp(f, 0.0f, 0.4f);

    float resonator_out;
    if(sustain_)
    {
        sustain_gain_ = accent_ * decay_;
        phase_ += f;
        phase_ = phase_ >= 1.f ? phase_ - 1.f : phase_;

        resonator_out = sin(TWOPI_F * phase_) * sustain_gain_;
        lp_out_       = cos(TWOPI_F * phase_) * sustain_gain_;
    }
    else
    {
        resonator_.SetFreq(f * sample_rate_);
        //resonator_.SetRes(1.0f + q * f);
        resonator_.SetRes(.4f * q * f);

        resonator_.Process((pulse - retrig_pulse_ * 0.2f) * scale);
        resonator_out = resonator_.Band();
        lp_out_       = resonator_.Low();
    }

    fonepole(tone_lp_, pulse * exciter_leak + resonator_out, tone_f);

    return tone_lp_;
}

void AnalogBassDrum::Trig()
{
    trig_ = true;
}

void AnalogBassDrum::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void AnalogBassDrum::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void AnalogBassDrum::SetFreq(float f0)
{
    f0 /= sample_rate_;
    f0_ = fclamp(f0, 0.f, .5f);
}

void AnalogBassDrum::SetTone(float tone)
{
    tone_ = fclamp(tone, 0.f, 1.f);
}

void AnalogBassDrum::SetDecay(float decay)
{
    decay_ = decay * .1f;
    decay_ -= .1f;
}

void AnalogBassDrum::SetAttackFmAmount(float attack_fm_amount)
{
    attack_fm_amount_ = attack_fm_amount * 50.f;
}

void AnalogBassDrum::SetSelfFmAmount(float self_fm_amount)
{
    self_fm_amount_ = self_fm_amount * 50.f;
}


/**  
       @brief 808 snare drum model, revisited.
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/analog_snare_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class AnalogSnareDrum
{
  public:
    AnalogSnareDrum() {}
    ~AnalogSnareDrum() {}

    static const int kNumModes = 5;

    /** Init the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next sample
        \param trigger Hit the drum with true. Defaults to false.
    */
    float Process(bool trigger = false);

    /** Trigger the drum */
    void Trig();

    /** Init the module
        \param sample_rate Audio engine sample rate
    */
    void SetSustain(bool sustain);

    /** Set how much accent to use
        \param accent Works 0-1.
    */
    void SetAccent(float accent);

    /** Set the drum's root frequency
        \param f0 Freq in Hz
    */
    void SetFreq(float f0);

    /** Set the brightness of the drum tone.
        \param tone Works 0-1. 1 = bright, 0 = dark.
    */
    void SetTone(float tone);

    /** Set the length of the drum decay
        \param decay Works with positive numbers
    */
    void SetDecay(float decay);

    /** Sets the mix between snare and drum.
        \param snappy 1 = just snare. 0 = just drum.
    */
    void SetSnappy(float snappy);

  private:
    float sample_rate_;

    float f0_, tone_, accent_, snappy_, decay_;
    bool  sustain_;
    bool  trig_;

    inline float SoftLimit(float x);
    inline float SoftClip(float x);

    int   pulse_remaining_samples_;
    float pulse_;
    float pulse_height_;
    float pulse_lp_;
    float noise_envelope_;
    float sustain_gain_;

    Svf resonator_[kNumModes];
    Svf noise_filter_;

    // Replace the resonators in "free running" (sustain) mode.
    float phase_[kNumModes];
};

static const int kNumModes = 5;

void AnalogSnareDrum::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    trig_ = false;

    pulse_remaining_samples_ = 0;
    pulse_                   = 0.0f;
    pulse_height_            = 0.0f;
    pulse_lp_                = 0.0f;
    noise_envelope_          = 0.0f;
    sustain_gain_            = 0.0f;

    SetSustain(false);
    SetAccent(.6f);
    SetFreq(200.f);
    SetDecay(.3f);
    SetSnappy(.7f);
    SetTone(.5f);

    for(int i = 0; i < kNumModes; ++i)
    {
        resonator_[i].Init(sample_rate_);
        phase_[i] = 0.f;
    }
    noise_filter_.Init(sample_rate_);
}

/** Trigger the drum */
void AnalogSnareDrum::Trig()
{
    trig_ = true;
}

void AnalogSnareDrum::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void AnalogSnareDrum::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void AnalogSnareDrum::SetFreq(float f0)
{
    f0  = f0 / sample_rate_;
    f0_ = fclamp(f0, 0.f, .4f);
}

void AnalogSnareDrum::SetTone(float tone)
{
    tone_ = fclamp(tone, 0.f, 1.f);
    tone_ *= 2.f;
}

void AnalogSnareDrum::SetDecay(float decay)
{
    decay_ = decay;
    return;
    decay_ = fmax(decay, 0.f);
}

void AnalogSnareDrum::SetSnappy(float snappy)
{
    snappy_ = fclamp(snappy, 0.f, 1.f);
}

float AnalogSnareDrum::Process(bool trigger)
{
    const float decay_xt = decay_ * (1.0f + decay_ * (decay_ - 1.0f));
    const int   kTriggerPulseDuration = 1.0e-3 * sample_rate_;
    const float kPulseDecayTime       = 0.1e-3 * sample_rate_;
    const float q = 2000.0f * powf(2.f, kOneTwelfth * decay_xt * 84.0f);
    const float noise_envelope_decay
        = 1.0f
          - 0.0017f
                * powf(2.f,
                       kOneTwelfth * (-decay_ * (50.0f + snappy_ * 10.0f)));
    const float exciter_leak = snappy_ * (2.0f - snappy_) * 0.1f;

    float snappy = snappy_ * 1.1f - 0.05f;
    snappy       = fclamp(snappy, 0.0f, 1.0f);

    float tone = tone_;

    if(trigger || trig_)
    {
        trig_                    = false;
        pulse_remaining_samples_ = kTriggerPulseDuration;
        pulse_height_            = 3.0f + 7.0f * accent_;
        noise_envelope_          = 2.0f;
    }

    static const float kModeFrequencies[kNumModes]
        = {1.00f, 2.00f, 3.18f, 4.16f, 5.62f};

    float f[kNumModes];
    float gain[kNumModes];

    for(int i = 0; i < kNumModes; ++i)
    {
        f[i] = fmin(f0_ * kModeFrequencies[i], 0.499f);
        resonator_[i].SetFreq(f[i] * sample_rate_);
        //        resonator_[i].SetRes(1.0f + f[i] * (i == 0 ? q : q * 0.25f));
        resonator_[i].SetRes((f[i] * (i == 0 ? q : q * 0.25f)) * .2);
    }

    if(tone < 0.666667f)
    {
        // 808-style (2 modes)
        tone *= 1.5f;
        gain[0] = 1.5f + (1.0f - tone) * (1.0f - tone) * 4.5f;
        gain[1] = 2.0f * tone + 0.15f;
        for(int i = 2; i < kNumModes; i++)
        {
            gain[i] = 0.f;
        }
    }
    else
    {
        // What the 808 could have been if there were extra modes!
        tone    = (tone - 0.666667f) * 3.0f;
        gain[0] = 1.5f - tone * 0.5f;
        gain[1] = 2.15f - tone * 0.7f;
        for(int i = 2; i < kNumModes; ++i)
        {
            gain[i] = tone;
            tone *= tone;
        }
    }

    float f_noise = f0_ * 16.0f;
    fclamp(f_noise, 0.0f, 0.499f);
    noise_filter_.SetFreq(f_noise * sample_rate_);
    //noise_filter_.SetRes(1.0f + f_noise * 1.5f);
    noise_filter_.SetRes(f_noise * 1.5f);

    // Q45 / Q46
    float pulse = 0.0f;
    if(pulse_remaining_samples_)
    {
        --pulse_remaining_samples_;
        pulse = pulse_remaining_samples_ ? pulse_height_ : pulse_height_ - 1.0f;
        pulse_ = pulse;
    }
    else
    {
        pulse_ *= 1.0f - 1.0f / kPulseDecayTime;
        pulse = pulse_;
    }

    float sustain_gain_value = sustain_gain_ = accent_ * decay_;

    // R189 / C57 / R190 + C58 / C59 / R197 / R196 / IC14
    pulse_lp_ = fclamp(pulse_lp_, pulse, 0.75f);

    float shell = 0.0f;
    for(int i = 0; i < kNumModes; ++i)
    {
        float excitation
            = i == 0 ? (pulse - pulse_lp_) + 0.006f * pulse : 0.026f * pulse;

        phase_[i] += f[i];
        phase_[i] = phase_[i] >= 1.f ? phase_[i] - 1.f : phase_[i];

        resonator_[i].Process(excitation);

        shell += gain[i]
                 * (sustain_
                        ? sin(phase_[i] * TWOPI_F) * sustain_gain_value * 0.25f
                        : resonator_[i].Band() + excitation * exciter_leak);
    }
    shell = SoftClip(shell);

    // C56 / R194 / Q48 / C54 / R188 / D54
    float noise = 2.0f * rand() * kRandFrac - 1.0f;
    if(noise < 0.0f)
        noise = 0.0f;
    noise_envelope_ *= noise_envelope_decay;
    noise *= (sustain_ ? sustain_gain_value : noise_envelope_) * snappy * 2.0f;

    // C66 / R201 / C67 / R202 / R203 / Q49
    noise_filter_.Process(noise);
    noise = noise_filter_.Band();

    // IC13
    return noise + shell * (1.0f - snappy);
}

inline float AnalogSnareDrum::SoftLimit(float x)
{
    return x * (27.0f + x * x) / (27.0f + 9.0f * x * x);
}

inline float AnalogSnareDrum::SoftClip(float x)
{
    if(x < -3.0f)
    {
        return -1.0f;
    }
    else if(x > 3.0f)
    {
        return 1.0f;
    }
    else
    {
        return SoftLimit(x);
    }
}


/**  
       @brief 808 style "metallic noise" with 6 square oscillators.
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/hihat.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SquareNoise
{
  public:
    SquareNoise() {}
    ~SquareNoise() {}

    void Init(float sample_rate);

    float Process(float f0);

  private:
    uint32_t phase_[6];
};

/**  
       @brief Ring mod style metallic noise generator.
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/hihat.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class RingModNoise
{
  public:
    RingModNoise() {}
    ~RingModNoise() {}

    void Init(float sample_rate);

    float Process(float f0);

  private:
    float      ProcessPair(Oscillator* osc, float f1, float f2);
    Oscillator oscillator_[6];

    float sample_rate_;
};

/**  
       @brief Swing type VCA
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/hihat.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SwingVCA
{
  public:
    float operator()(float s, float gain)
    {
        s *= s > 0.0f ? 10.0f : 0.1f;
        s = s / (1.0f + fabsf(s));
        return (s + 1.0f) * gain;
    }
};

/**  
       @brief Linear type VCA
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/hihat.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class LinearVCA
{
  public:
    float operator()(float s, float gain) { return s * gain; }
};

/**  
       @brief 808 HH, with a few extra parameters to push things to the CY territory...
       @author Ben Sergentanis
       @date Jan 2021
       The template parameter MetallicNoiseSource allows another kind of "metallic \n
       noise" to be used, for results which are more similar to KR-55 or FM hi-hats. \n \n 
       Ported from pichenettes/eurorack/plaits/dsp/drums/hihat.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
template <typename MetallicNoiseSource = SquareNoise,
          typename VCA                 = LinearVCA,
          bool resonance               = true>
class HiHat
{
  public:
    HiHat() {}
    ~HiHat() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate)
    {
        sample_rate_ = sample_rate;

        trig_ = false;

        envelope_     = 0.0f;
        noise_clock_  = 0.0f;
        noise_sample_ = 0.0f;
        sustain_gain_ = 0.0f;

        SetFreq(3000.f);
        SetTone(.5f);
        SetDecay(.2f);
        SetNoisiness(.8f);
        SetAccent(.8f);
        SetSustain(false);

        metallic_noise_.Init(sample_rate_);
        noise_coloration_svf_.Init(sample_rate_);
        hpf_.Init(sample_rate_);
    }

    /** Get the next sample
        \param trigger Hit the hihat with true. Defaults to false.
    */
    float Process(bool trigger = false)
    {
        const float envelope_decay
            = 1.0f - 0.003f * SemitonesToRatio(-decay_ * 84.0f);
        const float cut_decay
            = 1.0f - 0.0025f * SemitonesToRatio(-decay_ * 36.0f);

        if(trigger || trig_)
        {
            trig_ = false;

            envelope_
                = (1.5f + 0.5f * (1.0f - decay_)) * (0.3f + 0.7f * accent_);
        }

        // Process the metallic noise.
        float out = metallic_noise_.Process(2.0f * f0_);

        // Apply BPF on the metallic noise.
        float cutoff = 150.0f / sample_rate_ * SemitonesToRatio(tone_ * 72.0f);

        cutoff = fclamp(cutoff, 0.0f, 16000.0f / sample_rate_);


        noise_coloration_svf_.SetFreq(cutoff * sample_rate_);
        noise_coloration_svf_.SetRes(resonance ? 3.0f + 6.0f * tone_ : 1.0f);

        noise_coloration_svf_.Process(out);
        out = noise_coloration_svf_.Band();

        // This is not at all part of the 808 circuit! But to add more variety, we
        // add a variable amount of clocked noise to the output of the 6 schmitt
        // trigger oscillators.
        float noise_f = f0_ * (16.0f + 16.0f * (1.0f - noisiness_));
        noise_f       = fclamp(noise_f, 0.0f, 0.5f);

        noise_clock_ += noise_f;
        if(noise_clock_ >= 1.0f)
        {
            noise_clock_ -= 1.0f;
            noise_sample_ = rand() * kRandFrac - 0.5f;
        }
        out += noisiness_ * (noise_sample_ - out);

        // Apply VCA.
        sustain_gain_ = accent_ * decay_;
        VCA vca;
        envelope_ *= envelope_ > 0.5f ? envelope_decay : cut_decay;
        out = vca(out, sustain_ ? sustain_gain_ : envelope_);

        hpf_.SetFreq(cutoff * sample_rate_);
        hpf_.SetRes(.5f);
        hpf_.Process(out);
        out = hpf_.High();

        return out;
    }

    /** Trigger the hihat */
    void Trig() { trig_ = true; }

    /** Make the hihat ring out infinitely.
        \param sustain True = infinite sustain.
    */
    void SetSustain(bool sustain) { sustain_ = sustain; }

    /** Set how much accent to use
        \param accent Works 0-1.
    */
    void SetAccent(float accent) { accent_ = fclamp(accent, 0.f, 1.f); }

    /** Set the hihat tone's root frequency
        \param f0 Freq in Hz
    */
    void SetFreq(float f0)
    {
        f0 /= sample_rate_;
        f0_ = fclamp(f0, 0.f, 1.f);
    }

    /** Set the overall brightness / darkness of the hihat.
        \param tone Works from 0-1.
    */
    void SetTone(float tone) { tone_ = fclamp(tone, 0.f, 1.f); }

    /** Set the length of the hihat decay
        \param decay Works > 0. Tuned for 0-1.
    */
    void SetDecay(float decay)
    {
        decay_ = fmax(decay, 0.f);
        decay_ *= 1.7;
        decay_ -= 1.2;
    }

    /** Sets the mix between tone and noise
        \param snappy 1 = just noise. 0 = just tone.
    */
    void SetNoisiness(float noisiness)
    {
        noisiness_ = fclamp(noisiness, 0.f, 1.f);
        noisiness_ *= noisiness_;
    }


  private:
    float sample_rate_;

    float accent_, f0_, tone_, decay_, noisiness_;
    bool  sustain_;
    bool  trig_;

    float SemitonesToRatio(float in) { return powf(2.f, in * kOneTwelfth); }

    float envelope_;
    float noise_clock_;
    float noise_sample_;
    float sustain_gain_;

    MetallicNoiseSource metallic_noise_;
    Svf                 noise_coloration_svf_;
    Svf                 hpf_;
};


void SquareNoise::Init(float sample_rate)
{
    for(int i = 0; i < 6; i++)
    {
        phase_[i] = 0;
    }
}

float SquareNoise::Process(float f0)
{
    const float ratios[6] = {// Nominal f0: 414 Hz
                             1.0f,
                             1.304f,
                             1.466f,
                             1.787f,
                             1.932f,
                             2.536f};

    uint32_t increment[6];
    uint32_t phase[6];
    for(int i = 0; i < 6; ++i)
    {
        float f = f0 * ratios[i];
        if(f >= 0.499f)
            f = 0.499f;
        increment[i] = static_cast<uint32_t>(f * 4294967296.0f);
        phase[i]     = phase_[i];
    }

    phase[0] += increment[0];
    phase[1] += increment[1];
    phase[2] += increment[2];
    phase[3] += increment[3];
    phase[4] += increment[4];
    phase[5] += increment[5];
    uint32_t noise = 0;
    noise += (phase[0] >> 31);
    noise += (phase[1] >> 31);
    noise += (phase[2] >> 31);
    noise += (phase[3] >> 31);
    noise += (phase[4] >> 31);
    noise += (phase[5] >> 31);

    for(int i = 0; i < 6; ++i)
    {
        phase_[i] = phase[i];
    }

    return 0.33f * static_cast<float>(noise) - 1.0f;
}

void RingModNoise::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    for(int i = 0; i < 6; ++i)
    {
        oscillator_[i].Init(sample_rate_);
    }
}

float RingModNoise::Process(float f0)
{
    const float ratio = f0 / (0.01f + f0);
    const float f1a   = 200.0f / sample_rate_ * ratio;
    const float f1b   = 7530.0f / sample_rate_ * ratio;
    const float f2a   = 510.0f / sample_rate_ * ratio;
    const float f2b   = 8075.0f / sample_rate_ * ratio;
    const float f3a   = 730.0f / sample_rate_ * ratio;
    const float f3b   = 10500.0f / sample_rate_ * ratio;

    float out = ProcessPair(&oscillator_[0], f1a, f1b);
    out += ProcessPair(&oscillator_[2], f2a, f2b);
    out += ProcessPair(&oscillator_[4], f3a, f3b);

    return out;
}

float RingModNoise::ProcessPair(Oscillator* osc, float f1, float f2)
{
    osc[0].SetWaveform(Oscillator::WAVE_SQUARE);
    osc[0].SetFreq(f1 * sample_rate_);
    float temp_1 = osc[0].Process();

    osc[1].SetWaveform(Oscillator::WAVE_SAW);
    osc[1].SetFreq(f2 * sample_rate_);
    float temp_2 = osc[1].Process();

    return temp_1 * temp_2;
}


/**  
       @brief Click noise for SyntheticBassDrum
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/synthetic_bass_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SyntheticBassDrumClick
{
  public:
    SyntheticBassDrumClick() {}
    ~SyntheticBassDrumClick() {}

    /** Init the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get the next sample.
        \param in Trigger the click.
    */
    float Process(float in);

  private:
    float lp_;
    float hp_;
    Svf   filter_;
};

/**  
       @brief Attack Noise generator for SyntheticBassDrum. 
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/drums/synthetic_bass_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SyntheticBassDrumAttackNoise
{
  public:
    SyntheticBassDrumAttackNoise() {}
    ~SyntheticBassDrumAttackNoise() {}

    /** Init the module */
    void Init();

    /** Get the next sample. */
    float Process();

  private:
    float lp_;
    float hp_;
};

/**  
       @brief Naive bass drum model (modulated oscillator with FM + envelope).
       @author Ben Sergentanis
       @date Jan 2021
       Inadvertently 909-ish. \n \n 
       Ported from pichenettes/eurorack/plaits/dsp/drums/synthetic_bass_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SyntheticBassDrum
{
  public:
    SyntheticBassDrum() {}
    ~SyntheticBassDrum() {}

    /** Init the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Generates a distorted sine wave */
    inline float DistortedSine(float phase, float phase_noise, float dirtiness);

    /** Transistor VCA simulation.
        \param s Input sample.
        \param gain VCA gain.
    */
    inline float TransistorVCA(float s, float gain);

    /** Get the next sample.
        \param trigger True triggers the BD. This is optional.
    */
    float Process(bool trigger = false);

    /** Trigger the drum */
    void Trig();

    /** Allows the drum to play continuously
        \param sustain True sets the drum on infinite sustain.
    */
    void SetSustain(bool sustain);

    /** Sets the amount of accent.
        \param accent Works 0-1.
    */
    void SetAccent(float accent);

    /** Set the bass drum's root frequency.
        \param Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Sets the overall bright / darkness of the drum.
        \param tone Works 0-1.
    */
    void SetTone(float tone);

    /** Sets how long the drum's volume takes to decay.
        \param Works 0-1.
    */
    void SetDecay(float decay);

    /** Makes things grimy
        \param dirtiness Works 0-1.
    */
    void SetDirtiness(float dirtiness);

    /** Sets how much of a pitch sweep the drum experiences when triggered.
        \param fm_envelope_amount Works 0-1.
    */
    void SetFmEnvelopeAmount(float fm_envelope_amount);

    /** Sets how long the initial pitch sweep takes.
        \param fm_envelope_decay Works 0-1.
    */
    void SetFmEnvelopeDecay(float fm_envelope_decay);

  private:
    float sample_rate_;

    bool  trig_;
    bool  sustain_;
    float accent_, new_f0_, tone_, decay_;
    float dirtiness_, fm_envelope_amount_, fm_envelope_decay_;

    float f0_;
    float phase_;
    float phase_noise_;

    float fm_;
    float fm_lp_;
    float body_env_;
    float body_env_lp_;
    float transient_env_;
    float transient_env_lp_;

    float sustain_gain_;

    float tone_lp_;

    SyntheticBassDrumClick       click_;
    SyntheticBassDrumAttackNoise noise_;

    int body_env_pulse_width_;
    int fm_pulse_width_;
};


void SyntheticBassDrumClick::Init(float sample_rate)
{
    lp_ = 0.0f;
    hp_ = 0.0f;
    filter_.Init(sample_rate);
    filter_.SetFreq(5000.0f);
    filter_.SetRes(1.f); //2.f
}

float SyntheticBassDrumClick::Process(float in)
{
    //SLOPE(lp_, in, 0.5f, 0.1f);
    float error = in - lp_;
    lp_ += (error > 0 ? .5f : .1f) * error;

    fonepole(hp_, lp_, 0.04f);
    filter_.Process(lp_ - hp_);
    return filter_.Low();
}

void SyntheticBassDrumAttackNoise::Init()
{
    lp_ = 0.0f;
    hp_ = 0.0f;
}

float SyntheticBassDrumAttackNoise::Process()
{
    float sample = rand() * kRandFrac;
    fonepole(lp_, sample, 0.05f);
    fonepole(hp_, lp_, 0.005f);
    return lp_ - hp_;
}

void SyntheticBassDrum::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    trig_ = false;

    phase_                = 0.0f;
    phase_noise_          = 0.0f;
    f0_                   = 0.0f;
    fm_                   = 0.0f;
    fm_lp_                = 0.0f;
    body_env_lp_          = 0.0f;
    body_env_             = 0.0f;
    body_env_pulse_width_ = 0;
    fm_pulse_width_       = 0;
    tone_lp_              = 0.0f;
    sustain_gain_         = 0.0f;

    SetFreq(100.f);
    SetSustain(false);
    SetAccent(.2f);
    SetTone(.6f);
    SetDecay(.7f);
    SetDirtiness(.3f);
    SetFmEnvelopeAmount(.6);
    SetFmEnvelopeDecay(.3);

    click_.Init(sample_rate);
    noise_.Init();
}

inline float SyntheticBassDrum::DistortedSine(float phase,
                                              float phase_noise,
                                              float dirtiness)
{
    phase += phase_noise * dirtiness;

    //MAKE_INTEGRAL_FRACTIONAL(phase);
    int32_t phase_integral   = static_cast<int32_t>(phase);
    float   phase_fractional = phase - static_cast<float>(phase_integral);

    phase            = phase_fractional;
    float triangle   = (phase < 0.5f ? phase : 1.0f - phase) * 4.0f - 1.0f;
    float sine       = 2.0f * triangle / (1.0f + fabsf(triangle));
    float clean_sine = sinf(TWOPI_F * (phase + 0.75f));
    return sine + (1.0f - dirtiness) * (clean_sine - sine);
}

inline float SyntheticBassDrum::TransistorVCA(float s, float gain)
{
    s = (s - 0.6f) * gain;
    return 3.0f * s / (2.0f + fabsf(s)) + gain * 0.3f;
}

float SyntheticBassDrum::Process(bool trigger)
{
    float dirtiness = dirtiness_;
    dirtiness *= fmax(1.0f - 8.0f * new_f0_, 0.0f);

    const float fm_decay
        = 1.0f
          - 1.0f / (0.008f * (1.0f + fm_envelope_decay_ * 4.0f) * sample_rate_);

    const float body_env_decay
        = 1.0f
          - 1.0f / (0.02f * sample_rate_)
                * powf(2.f, (-decay_ * 60.0f) * kOneTwelfth);
    const float transient_env_decay = 1.0f - 1.0f / (0.005f * sample_rate_);
    const float tone_f              = fmin(
        4.0f * new_f0_ * powf(2.f, (tone_ * 108.0f) * kOneTwelfth), 1.0f);
    const float transient_level = tone_;

    if(trigger || trig_)
    {
        trig_     = false;
        fm_       = 1.0f;
        body_env_ = transient_env_ = 0.3f + 0.7f * accent_;
        body_env_pulse_width_      = sample_rate_ * 0.001f;
        fm_pulse_width_            = sample_rate_ * 0.0013f;
    }

    sustain_gain_ = accent_ * decay_;

    fonepole(phase_noise_, rand() * kRandFrac - 0.5f, 0.002f);

    float mix = 0.0f;

    if(sustain_)
    {
        f0_ = new_f0_;
        phase_ += f0_;
        if(phase_ >= 1.0f)
        {
            phase_ -= 1.0f;
        }
        float body = DistortedSine(phase_, phase_noise_, dirtiness);
        mix -= TransistorVCA(body, sustain_gain_);
    }
    else
    {
        if(fm_pulse_width_)
        {
            --fm_pulse_width_;
            phase_ = 0.25f;
        }
        else
        {
            fm_ *= fm_decay;
            float fm = 1.0f + fm_envelope_amount_ * 3.5f * fm_lp_;
            f0_      = new_f0_;
            phase_ += fmin(f0_ * fm, 0.5f);
            if(phase_ >= 1.0f)
            {
                phase_ -= 1.0f;
            }
        }

        if(body_env_pulse_width_)
        {
            --body_env_pulse_width_;
        }
        else
        {
            body_env_ *= body_env_decay;
            transient_env_ *= transient_env_decay;
        }

        const float envelope_lp_f = 0.1f;
        fonepole(body_env_lp_, body_env_, envelope_lp_f);
        fonepole(transient_env_lp_, transient_env_, envelope_lp_f);
        fonepole(fm_lp_, fm_, envelope_lp_f);

        float body      = DistortedSine(phase_, phase_noise_, dirtiness);
        float transient = click_.Process(body_env_pulse_width_ ? 0.0f : 1.0f)
                          + noise_.Process();

        mix -= TransistorVCA(body, body_env_lp_);
        mix -= transient * transient_env_lp_ * transient_level;
    }

    fonepole(tone_lp_, mix, tone_f);
    return tone_lp_;
}

void SyntheticBassDrum::Trig()
{
    trig_ = true;
}

void SyntheticBassDrum::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void SyntheticBassDrum::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void SyntheticBassDrum::SetFreq(float freq)
{
    freq /= sample_rate_;
    new_f0_ = fclamp(freq, 0.f, 1.f);
}

void SyntheticBassDrum::SetTone(float tone)
{
    tone_ = fclamp(tone, 0.f, 1.f);
}

void SyntheticBassDrum::SetDecay(float decay)
{
    decay  = fclamp(decay, 0.f, 1.f);
    decay_ = decay * decay;
}

void SyntheticBassDrum::SetDirtiness(float dirtiness)
{
    dirtiness_ = fclamp(dirtiness, 0.f, 1.f);
}

void SyntheticBassDrum::SetFmEnvelopeAmount(float fm_envelope_amount)
{
    fm_envelope_amount_ = fclamp(fm_envelope_amount, 0.f, 1.f);
}

void SyntheticBassDrum::SetFmEnvelopeDecay(float fm_envelope_decay)
{
    fm_envelope_decay  = fclamp(fm_envelope_decay, 0.f, 1.f);
    fm_envelope_decay_ = fm_envelope_decay * fm_envelope_decay;
}


/**  
       @brief Naive snare drum model (two modulated oscillators + filtered noise).
       @author Ben Sergentanis
       @date Jan 2021
       Uses a few magic numbers taken from the 909 schematics: \n 
       - Ratio between the two modes of the drum set to 1.47. \n
       - Funky coupling between the two modes. \n
       - Noise coloration filters and envelope shapes for the snare. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/drums/synthetic_snare_drum.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class SyntheticSnareDrum
{
  public:
    SyntheticSnareDrum() {}
    ~SyntheticSnareDrum() {}

    /** Init the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next sample.
        \param trigger True = hit the drum. This argument is optional.
    */
    float Process(bool trigger = false);

    /** Trigger the drum */
    void Trig();

    /** Make the drum ring out infinitely.
        \param sustain True = infinite sustain.
    */
    void SetSustain(bool sustain);

    /** Set how much accent to use
        \param accent Works 0-1.
    */
    void SetAccent(float accent);

    /** Set the drum's root frequency
        \param f0 Freq in Hz
    */
    void SetFreq(float f0);

    /** Set the amount of fm sweep.
        \param fm_amount Works from 0 - 1.
    */
    void SetFmAmount(float fm_amount);

    /** Set the length of the drum decay
        \param decay Works with positive numbers
    */
    void SetDecay(float decay);

    /** Sets the mix between snare and drum.
        \param snappy 1 = just snare. 0 = just drum.
    */
    void SetSnappy(float snappy);

  private:
    inline float DistortedSine(float phase);

    float sample_rate_;

    bool  trig_;
    bool  sustain_;
    float accent_, f0_, fm_amount_, decay_, snappy_;

    float phase_[2];
    float drum_amplitude_;
    float snare_amplitude_;
    float fm_;
    float sustain_gain_;
    int   hold_counter_;

    Svf drum_lp_;
    Svf snare_hp_;
    Svf snare_lp_;
};

void SyntheticSnareDrum::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    phase_[0]        = 0.0f;
    phase_[1]        = 0.0f;
    drum_amplitude_  = 0.0f;
    snare_amplitude_ = 0.0f;
    fm_              = 0.0f;
    hold_counter_    = 0;
    sustain_gain_    = 0.0f;

    SetSustain(false);
    SetAccent(.6f);
    SetFreq(200.f);
    SetFmAmount(.1f);
    SetDecay(.3f);
    SetSnappy(.7f);

    trig_ = false;

    drum_lp_.Init(sample_rate_);
    snare_hp_.Init(sample_rate_);
    snare_lp_.Init(sample_rate_);
}

inline float SyntheticSnareDrum::DistortedSine(float phase)
{
    float triangle = (phase < 0.5f ? phase : 1.0f - phase) * 4.0f - 1.3f;
    return 2.0f * triangle / (1.0f + fabsf(triangle));
}

bool even = true;

float SyntheticSnareDrum::Process(bool trigger)
{
    const float decay_xt = decay_ * (1.0f + decay_ * (decay_ - 1.0f));
    const float drum_decay
        = 1.0f
          - 1.0f / (0.015f * sample_rate_)
                * powf(2.f,
                       kOneTwelfth
                           * (-decay_xt * 72.0f - fm_amount_ * 12.0f
                              + snappy_ * 7.0f));

    const float snare_decay
        = 1.0f
          - 1.0f / (0.01f * sample_rate_)
                * powf(2.f, kOneTwelfth * (-decay_ * 60.0f - snappy_ * 7.0f));
    const float fm_decay = 1.0f - 1.0f / (0.007f * sample_rate_);

    float snappy = snappy_ * 1.1f - 0.05f;
    snappy       = fclamp(snappy, 0.0f, 1.0f);

    const float drum_level  = sqrtf(1.0f - snappy);
    const float snare_level = sqrtf(snappy);

    const float snare_f_min = fmin(10.0f * f0_, 0.5f);
    const float snare_f_max = fmin(35.0f * f0_, 0.5f);

    snare_hp_.SetFreq(snare_f_min * sample_rate_);
    snare_lp_.SetFreq(snare_f_max * sample_rate_);
    snare_lp_.SetRes(0.5f + 2.0f * snappy);

    drum_lp_.SetFreq(3.0f * f0_ * sample_rate_);

    if(trigger || trig_)
    {
        trig_            = false;
        snare_amplitude_ = drum_amplitude_ = 0.3f + 0.7f * accent_;
        fm_                                = 1.0f;
        phase_[0] = phase_[1] = 0.0f;
        hold_counter_
            = static_cast<int>((0.04f + decay_ * 0.03f) * sample_rate_);
    }

    even = !even;
    if(sustain_)
    {
        sustain_gain_ = snare_amplitude_ = accent_ * decay_;
        drum_amplitude_                  = snare_amplitude_;
        fm_                              = 0.0f;
    }
    else
    {
        // Compute all D envelopes.
        // The envelope for the drum has a very long tail.
        // The envelope for the snare has a "hold" stage which lasts between
        // 40 and 70 ms
        drum_amplitude_
            *= (drum_amplitude_ > 0.03f || even) ? drum_decay : 1.0f;
        if(hold_counter_)
        {
            --hold_counter_;
        }
        else
        {
            snare_amplitude_ *= snare_decay;
        }
        fm_ *= fm_decay;
    }

    // The 909 circuit has a funny kind of oscillator coupling - the signal
    // leaving Q40's collector and resetting all oscillators allow some
    // intermodulation.
    float reset_noise        = 0.0f;
    float reset_noise_amount = (0.125f - f0_) * 8.0f;
    reset_noise_amount       = fclamp(reset_noise_amount, 0.0f, 1.0f);
    reset_noise_amount *= reset_noise_amount;
    reset_noise_amount *= fm_amount_;
    reset_noise += phase_[0] > 0.5f ? -1.0f : 1.0f;
    reset_noise += phase_[1] > 0.5f ? -1.0f : 1.0f;
    reset_noise *= reset_noise_amount * 0.025f;

    float f = f0_ * (1.0f + fm_amount_ * (4.0f * fm_));
    phase_[0] += f;
    phase_[1] += f * 1.47f;
    if(reset_noise_amount > 0.1f)
    {
        if(phase_[0] >= 1.0f + reset_noise)
        {
            phase_[0] = 1.0f - phase_[0];
        }
        if(phase_[1] >= 1.0f + reset_noise)
        {
            phase_[1] = 1.0f - phase_[1];
        }
    }
    else
    {
        if(phase_[0] >= 1.0f)
        {
            phase_[0] -= 1.0f;
        }
        if(phase_[1] >= 1.0f)
        {
            phase_[1] -= 1.0f;
        }
    }

    float drum = -0.1f;
    drum += DistortedSine(phase_[0]) * 0.60f;
    drum += DistortedSine(phase_[1]) * 0.25f;
    drum *= drum_amplitude_ * drum_level;

    drum_lp_.Process(drum);
    drum = drum_lp_.Low();

    float noise = rand() * kRandFrac;
    snare_lp_.Process(noise);
    float snare = snare_lp_.Low();
    snare_hp_.Process(snare);
    snare = snare_hp_.High();
    snare = (snare + 0.1f) * (snare_amplitude_ + fm_) * snare_level;

    return snare + drum; // It's a snare, it's a drum, it's a snare drum.
}

void SyntheticSnareDrum::Trig()
{
    trig_ = true;
}

void SyntheticSnareDrum::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void SyntheticSnareDrum::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void SyntheticSnareDrum::SetFreq(float f0)
{
    f0 /= sample_rate_;
    f0_ = fclamp(f0, 0.f, 1.f);
}

void SyntheticSnareDrum::SetFmAmount(float fm_amount)
{
    fm_amount  = fclamp(fm_amount, 0.f, 1.f);
    fm_amount_ = fm_amount * fm_amount;
}

void SyntheticSnareDrum::SetDecay(float decay)
{
    decay_ = fmax(decay, 0.f);
}

void SyntheticSnareDrum::SetSnappy(float snappy)
{
    snappy_ = fclamp(snappy, 0.f, 1.f);
}



/** Balances two sound sources. Sig is boosted to the level of comp.

    *Original author(s) : Barry Vercoe, john ffitch, Gabriel Maldonado

    *Year: 1991

    *Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Balance
{
  public:
    Balance() {}
    ~Balance() {}
    /** Initializes the balance module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** adjust sig level to level of comp
    */
    float Process(float sig, float comp);


    /** adjusts the rate at which level compensation happens
        \param cutoff : Sets half power point of special internal cutoff filter.

        defaults to 10
    */
    inline void SetCutoff(float cutoff) { ihp_ = cutoff; }

  private:
    float sample_rate_, ihp_, c2_, c1_, prvq_, prvr_, prva_;
};


void Balance::Init(float sample_rate)
{
    float b;
    sample_rate_ = sample_rate;
    ihp_         = 10.0f;
    b            = 2.0f - cosf(ihp_ * (TWOPI_F / sample_rate_));
    c2_          = b - sqrtf(b * b - 1.0f);
    c1_          = 1.0f - c2_;
    prvq_ = prvr_ = prva_ = 0.0f;
}

float Balance::Process(float sig, float comp)
{
    float q, r, a, diff, out;
    float c1 = c1_;
    float c2 = c2_;

    q        = prvq_;
    r        = prvr_;
    float as = sig;
    float cs = comp;

    q = c1 * as * as + c2 * q;
    r = c1 * cs * cs + c2 * r;

    prvq_ = q;
    prvr_ = r;

    if(q != 0.0f)
    {
        a = sqrtf(r / q);
    }
    else
    {
        a = sqrtf(r);
    }

    if((diff = a - prva_) != 0.0f)
    {
        out = sig * prva_;
    }
    else
    {
        out = sig * a;
    }

    prva_ = a;

    return out;
}

/** dynamics compressor

influenced by compressor in soundpipe (from faust).

Modifications made to do:
- Less calculations during each process loop (coefficients recalculated on parameter change).
- C++-ified
- added sidechain support
- pulled gain apart for monitoring and multichannel support
- improved readability
- improved makeup-gain calculations
- changing controls now costs a lot less
- a lot less expensive

by: shensley, improved upon by AvAars
\todo Add soft/hard knee settings
*/
class Compressor
{
  public:
    Compressor() {}
    ~Compressor() {}
    /** Initializes compressor
        \param sample_rate rate at which samples will be produced by the audio engine.
    */
    void Init(float sample_rate);

    /** Compress the audio input signal, saves the calculated gain
        \param in audio input signal
    */
    float Process(float in);

    /** Compresses the audio input signal, keyed by a secondary input.
        \param in audio input signal (to be compressed)
        \param key audio input that will be used to side-chain the compressor
    */
    float Process(float in, float key)
    {
        Process(key);
        return Apply(in);
    }

    /** Apply compression to the audio signal, based on the previously calculated gain
        \param in audio input signal
    */
    float Apply(float in) { return gain_ * in; }

    /** Compresses a block of audio
        \param in audio input signal
        \param out audio output signal
        \param size the size of the block
    */
    void ProcessBlock(float *in, float *out, size_t size)
    {
        ProcessBlock(in, out, in, size);
    }

    /** Compresses a block of audio, keyed by a secondary input
        \param in audio input signal (to be compressed)
        \param out audio output signal
        \param key audio input that will be used to side-chain the compressor
        \param size the size of the block
    */
    void ProcessBlock(float *in, float *out, float *key, size_t size);

    /** Compresses a block of multiple channels of audio, keyed by a secondary input
        \param in audio input signals (to be compressed)
        \param out audio output signals
        \param key audio input that will be used to side-chain the compressor
        \param channels the number of audio channels
        \param size the size of the block
    */
    void ProcessBlock(float **in,
                      float **out,
                      float * key,
                      size_t  channels,
                      size_t  size);

    /** Gets the amount of gain reduction */
    float GetRatio() { return ratio_; }

    /** Sets the amount of gain reduction applied to compressed signals
     \param ratio Expects 1.0 -> 40. (untested with values < 1.0)
    */
    void SetRatio(float ratio)
    {
        ratio_ = ratio;
        RecalculateRatio();
    }

    /** Gets the threshold in dB */
    float GetThreshold() { return thresh_; }

    /** Sets the threshold in dB at which compression will be applied
     \param threshold Expects 0.0 -> -80.
    */
    void SetThreshold(float threshold)
    {
        thresh_ = threshold;
        RecalculateMakeup();
    }

    /** Gets the envelope time for onset of compression */
    float GetAttack() { return atk_; }

    /** Sets the envelope time for onset of compression for signals above the threshold.
        \param attack Expects 0.001 -> 10
    */
    void SetAttack(float attack)
    {
        atk_ = attack;
        RecalculateAttack();
    }

    /** Gets the envelope time for release of compression */
    float GetRelease() { return rel_; }

    /** Sets the envelope time for release of compression as input signal falls below threshold.
        \param release Expects 0.001 -> 10
    */
    void SetRelease(float release)
    {
        rel_ = release;
        RecalculateRelease();
    }

    /** Gets the additional gain to make up for the compression */
    float GetMakeup() { return makeup_gain_; }

    /** Manually sets the additional gain to make up for the compression
        \param gain Expects 0.0 -> 80
    */
    void SetMakeup(float gain) { makeup_gain_ = gain; }

    /** Enables or disables the automatic makeup gain. Disabling sets the makeup gain to 0.0
        \param enable true to enable, false to disable
    */
    void AutoMakeup(bool enable)
    {
        makeup_auto_ = enable;
        makeup_gain_ = 0.0f;
        RecalculateMakeup();
    }

    /** Gets the gain reduction in dB
    */
    float GetGain() { return fastlog10f(gain_) * 20.0f; }

  private:
    float ratio_, thresh_, atk_, rel_;
    float makeup_gain_;
    float gain_;

    // Recorded slope and gain, used in next sample
    float slope_rec_, gain_rec_;

    // Internals from faust
    float atk_slo2_, ratio_mul_, atk_slo_, rel_slo_;

    int   sample_rate_;
    float sample_rate_inv2_, sample_rate_inv_;

    // Auto makeup gain enable
    bool makeup_auto_;

    // Methods for recalculating internals
    void RecalculateRatio()
    {
        ratio_mul_ = ((1.0f - atk_slo2_) * ((1.0f / ratio_) - 1.0f));
    }

    void RecalculateAttack()
    {
        atk_slo_  = expf(-(sample_rate_inv_ / atk_));
        atk_slo2_ = expf(-(sample_rate_inv2_ / atk_));

        RecalculateRatio();
    }

    void RecalculateRelease() { rel_slo_ = expf((-(sample_rate_inv_ / rel_))); }

    void RecalculateMakeup()
    {
        if(makeup_auto_)
            makeup_gain_ = fabsf(thresh_ - thresh_ / ratio_) * 0.5f;
    }
};


void Compressor::Init(float sample_rate)
{
    sample_rate_      = min(192000, max(1, sample_rate));
    sample_rate_inv_  = 1.0f / (float)sample_rate_;
    sample_rate_inv2_ = 2.0f / (float)sample_rate_;

    // Initializing the params in this order to avoid dividing by zero

    SetRatio(2.0f);
    SetAttack(0.1f);
    SetRelease(0.1f);
    SetThreshold(-12.0f);
    AutoMakeup(true);

    gain_rec_  = 0.1f;
    slope_rec_ = 0.1f;
}

float Compressor::Process(float in)
{
    float inAbs   = fabsf(in);
    float cur_slo = ((slope_rec_ > inAbs) ? rel_slo_ : atk_slo_);
    slope_rec_    = ((slope_rec_ * cur_slo) + ((1.0f - cur_slo) * inAbs));
    gain_rec_     = ((atk_slo2_ * gain_rec_)
                 + (ratio_mul_
                    * fmax(((20.f * fastlog10f(slope_rec_)) - thresh_), 0.f)));
    gain_         = pow10f(0.05f * (gain_rec_ + makeup_gain_));

    return gain_ * in;
}

void Compressor::ProcessBlock(float *in, float *out, float *key, size_t size)
{
    for(size_t i = 0; i < size; i++)
    {
        Process(key[i]);
        out[i] = Apply(in[i]);
    }
}

// Multi-channel block processing
void Compressor::ProcessBlock(float **in,
                              float **out,
                              float * key,
                              size_t  channels,
                              size_t  size)
{
    for(size_t i = 0; i < size; i++)
    {
        Process(key[i]);
        for(size_t c = 0; c < channels; c++)
        {
            out[c][i] = Apply(in[c][i]);
        }
    }
}


/** Curve applied to the CrossFade
- LIN = linear
- CPOW = constant power
- LOG = logarithmic
- EXP  exponential
- LAST = end of enum (used for array indexing)
*/
enum
{
    CROSSFADE_LIN,
    CROSSFADE_CPOW,
    CROSSFADE_LOG,
    CROSSFADE_EXP,
    CROSSFADE_LAST,
};

/** Performs a CrossFade between two signals

    Original author: Paul Batchelor

    Ported from Soundpipe by Andrew Ikenberry

    added curve option for constant power, etc.
*/
class CrossFade
{
  public:
    CrossFade() {}
    ~CrossFade() {}
    /** Initializes CrossFade module
        Defaults
        - current position = .5
        - curve = linear
    */
    inline void Init(int curve)
    {
        pos_   = 0.5f;
        curve_ = curve < CROSSFADE_LAST ? curve : CROSSFADE_LIN;
    }

    /** Initialize with default linear curve 
    */
    inline void Init() { Init(CROSSFADE_LIN); }
    /** processes CrossFade and returns single sample
    */
    float Process(float &in1, float &in2);


    /** Sets position of CrossFade between two input signals
        Input range: 0 to 1
    */
    inline void SetPos(float pos) { pos_ = pos; }
    /** Sets current curve applied to CrossFade 
    Expected input: See [Curve Options](##curve-options)
    */
    inline void SetCurve(uint8_t curve) { curve_ = curve; }
    /** Returns current position
    */
    inline float GetPos(float pos) { return pos_; }
    /** Returns current curve
    */
    inline uint8_t GetCurve(uint8_t curve) { return curve_; }

  private:
    float   pos_;
    uint8_t curve_;
};



/** Simple Peak Limiter

This was extracted from pichenettes/stmlib.

Credit to pichenettes/Mutable Instruments
*/
class Limiter
{
  public:
    Limiter() {}
    ~Limiter() {}
    /** Initializes the Limiter instance. 
    */
    void Init();

    /** Processes a block of audio through the limiter.
        \param in - pointer to a block of audio samples to be processed. The buffer is operated on directly.
        \param size - size of the buffer "in"
        \param pre_gain - amount of pre_gain applied to the signal.
    */
    void ProcessBlock(float *in, size_t size, float pre_gain);

  private:
    float peak_;
};


void Limiter::Init()
{
    peak_ = 0.5f;
}

void Limiter::ProcessBlock(float *in, size_t size, float pre_gain)
{
    while(size--)
    {
        float pre  = *in * pre_gain;
        float peak = fabsf(pre);
        SLOPE(peak_, peak, 0.05f, 0.00002f);
        float gain = (peak_ <= 1.0f ? 1.0f : 1.0f / peak_);
        *in++      = SoftLimit(pre * gain * 0.7f);
    }
}


/** Autowah module

    Original author(s) :

    Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Autowah
{
  public:
    Autowah() {}
    ~Autowah() {}
    /** Initializes the Autowah module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Initializes the Autowah module.
        \param in - input signal to be wah'd
    */
    float Process(float in);


    /** sets wah
        \param wah : set wah amount, , 0...1.0
    */
    inline void SetWah(float wah) { wah_ = wah; }
    /** sets mix amount
        \param drywet : set effect dry/wet, 0...100.0
    */
    inline void SetDryWet(float drywet) { wet_dry_ = drywet; }
    /** sets wah level
        \param level : set wah level, 0...1.0
    */
    inline void SetLevel(float level) { level_ = level; }

  private:
    float sampling_freq_, const1_, const2_, const4_, wah_, level_, wet_dry_,
        rec0_[3], rec1_[2], rec2_[2], rec3_[2], rec4_[2], rec5_[2];
};


void Autowah::Init(float sample_rate)
{
    sampling_freq_ = sample_rate;
    const1_        = 1413.72f / sampling_freq_;
    const2_        = expf(0.0f - (100.0f / sampling_freq_));
    const4_        = expf(0.0f - (10.0f / sampling_freq_));

    wet_dry_ = 100.0f;
    level_   = 0.1f;
    wah_     = 0.0;


    for(int i = 0; i < 2; i++)
    {
        rec1_[i] = rec2_[i] = rec3_[i] = rec4_[i] = rec5_[i] = 0.0f;
    }

    for(int i = 0; i < 3; i++)
    {
        rec0_[i] = 0.0f;
    }
}

float Autowah::Process(float in)
{
    float out;
    float fSlow2 = (0.01f * (wet_dry_ * level_));
    float fSlow3 = (1.0f - 0.01f * wet_dry_) + (1.f - wah_);

    float fTemp1 = fabs(in);
    rec3_[0]
        = fmaxf(fTemp1, (const4_ * rec3_[1]) + ((1.0f - const4_) * fTemp1));
    rec2_[0]     = (const2_ * rec2_[1]) + ((1.0f - const2_) * rec3_[0]);
    float fTemp2 = fminf(1.0f, rec2_[0]);
    float fTemp3 = powf(2.0f, (2.3f * fTemp2));
    float fTemp4
        = 1.0f
          - (const1_ * fTemp3 / powf(2.0f, (1.0f + 2.0f * (1.0f - fTemp2))));
    rec1_[0]
        = ((0.999f * rec1_[1])
           + (0.001f
              * (0.0f - (2.0f * (fTemp4 * cosf((const1_ * 2 * fTemp3)))))));
    rec4_[0] = ((0.999f * rec4_[1]) + (0.001f * fTemp4 * fTemp4));
    rec5_[0] = ((0.999f * rec5_[1]) + (0.0001f * powf(4.0f, fTemp2)));
    rec0_[0] = (0.0f
                - (((rec1_[0] * rec0_[1]) + (rec4_[0] * rec0_[2]))
                   - (fSlow2 * (rec5_[0] * in))));

    out      = ((wah_ * (rec0_[0] - rec0_[1])) + (fSlow3 * in));
    rec3_[1] = rec3_[0];
    rec2_[1] = rec2_[0];
    rec1_[1] = rec1_[0];
    rec4_[1] = rec4_[0];
    rec5_[1] = rec5_[0];
    rec0_[2] = rec0_[1];
    rec0_[1] = rec0_[0];

    return out;
}



/**  
    @brief Single Chorus engine. Used in Chorus.
    @author Ben Sergentanis
*/
class ChorusEngine
{
  public:
    ChorusEngine() {}
    ~ChorusEngine() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get the next sample
        \param in Sample to process
    */
    float Process(float in);

    /** How much to modulate the delay by.
        \param depth Works 0-1.
    */
    void SetLfoDepth(float depth);

    /** Set lfo frequency.
        \param freq Frequency in Hz
    */
    void SetLfoFreq(float freq);

    /** Set the internal delay rate. 
        \param delay Tuned for 0-1. Maps to .1 to 50 ms.
    */
    void SetDelay(float delay);

    /** Set the delay time in ms.
        \param ms Delay time in ms, 0 to 50 ms.
    */
    void SetDelayMs(float ms);

    /** Set the feedback amount.
        \param feedback Amount from 0-1.
    */
    void SetFeedback(float feedback);

  private:
    float                    sample_rate_;
    static constexpr int32_t kDelayLength
        = 2400; // 50 ms at 48kHz = .05 * 48000

    //triangle lfos
    float lfo_phase_;
    float lfo_freq_;
    float lfo_amp_;

    float feedback_;

    float delay_;

    DelayLine<float, kDelayLength> del_;

    float ProcessLfo();
};

//wraps up all of the chorus engines
/**  
    @brief Chorus Effect.
    @author Ben Sergentanis
    @date Jan 2021
    Based on https://www.izotope.com/en/learn/understanding-chorus-flangers-and-phasers-in-audio-production.html \n
    and https://www.researchgate.net/publication/236629475_Implementing_Professional_Audio_Effects_with_DSPs \n
*/
class Chorus
{
  public:
    Chorus() {}
    ~Chorus() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the net floating point sample. Defaults to left channel.
        \param in Sample to process
    */
    float Process(float in);

    /** Get the left channel's last sample */
    float GetLeft();

    /** Get the right channel's last sample */
    float GetRight();

    /** Pan both channels individually.
        \param panl Pan the left channel. 0 is left, 1 is right.
        \param panr Pan the right channel.
    */
    void SetPan(float panl, float panr);

    /** Pan both channels.
        \param pan Where to pan both channels to. 0 is left, 1 is right.
    */
    void SetPan(float pan);

    /** Set both lfo depths individually.
        \param depthl Left channel lfo depth. Works 0-1.
        \param depthr Right channel lfo depth.
    */
    void SetLfoDepth(float depthl, float depthr);

    /** Set both lfo depths.
        \param depth Both channels lfo depth. Works 0-1.
    */
    void SetLfoDepth(float depth);

    /** Set both lfo frequencies individually.
        \param depthl Left channel lfo freq in Hz.
        \param depthr Right channel lfo freq in Hz.
    */
    void SetLfoFreq(float freql, float freqr);

    /** Set both lfo frequencies.
        \param depth Both channel lfo freqs in Hz.
    */
    void SetLfoFreq(float freq);

    /** Set both channel delay amounts individually.
        \param delayl Left channel delay amount. Works 0-1.
        \param delayr Right channel delay amount.
    */
    void SetDelay(float delayl, float delayr);

    /** Set both channel delay amounts.
        \param delay Both channel delay amount. Works 0-1.
    */
    void SetDelay(float delay);

    /** Set both channel delay individually.
        \param msl Left channel delay in ms.
        \param msr Right channel delay in ms.
    */
    void SetDelayMs(float msl, float msr);

    /** Set both channel delay in ms.
        \param ms Both channel delay amounts in ms.
    */
    void SetDelayMs(float ms);

    /** Set both channels feedback individually.
        \param feedbackl Left channel feedback. Works 0-1.
        \param feedbackr Right channel feedback.
    */
    void SetFeedback(float feedbackl, float feedbackr);

    /** Set both channels feedback.
        \param feedback Both channel feedback. Works 0-1.
    */
    void SetFeedback(float feedback);

  private:
    ChorusEngine engines_[2];
    float        gain_frac_;
    float        pan_[2];

    float sigl_, sigr_;
};

//ChorusEngine stuff
void ChorusEngine::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    del_.Init();
    lfo_amp_  = 0.f;
    feedback_ = .2f;
    SetDelay(.75);

    lfo_phase_ = 0.f;
    SetLfoFreq(.3f);
    SetLfoDepth(.9f);
}

float ChorusEngine::Process(float in)
{
    float lfo_sig = ProcessLfo();
    del_.SetDelay(lfo_sig + delay_);

    float out = del_.Read();
    del_.Write(in + out * feedback_);

    return (in + out) * .5f; //equal mix
}

void ChorusEngine::SetLfoDepth(float depth)
{
    depth    = fclamp(depth, 0.f, .93f);
    lfo_amp_ = depth * delay_;
}

void ChorusEngine::SetLfoFreq(float freq)
{
    freq = 4.f * freq / sample_rate_;
    freq *= lfo_freq_ < 0.f ? -1.f : 1.f;  //if we're headed down, keep going
    lfo_freq_ = fclamp(freq, -.25f, .25f); //clip at +/- .125 * sr
}

void ChorusEngine::SetDelay(float delay)
{
    delay = (.1f + delay * 7.9f); //.1 to 8 ms
    SetDelayMs(delay);
}

void ChorusEngine::SetDelayMs(float ms)
{
    ms     = fmax(.1f, ms);
    delay_ = ms * .001f * sample_rate_; //ms to samples

    lfo_amp_ = fmin(lfo_amp_, delay_); //clip this if needed
}

void ChorusEngine::SetFeedback(float feedback)
{
    feedback_ = fclamp(feedback, 0.f, 1.f);
}

float ChorusEngine::ProcessLfo()
{
    lfo_phase_ += lfo_freq_;

    //wrap around and flip direction
    if(lfo_phase_ > 1.f)
    {
        lfo_phase_ = 1.f - (lfo_phase_ - 1.f);
        lfo_freq_ *= -1.f;
    }
    else if(lfo_phase_ < -1.f)
    {
        lfo_phase_ = -1.f - (lfo_phase_ + 1.f);
        lfo_freq_ *= -1.f;
    }

    return lfo_phase_ * lfo_amp_;
}

//Chorus Stuff
void Chorus::Init(float sample_rate)
{
    engines_[0].Init(sample_rate);
    engines_[1].Init(sample_rate);
    SetPan(.25f, .75f);

    gain_frac_ = .5f;
    sigl_ = sigr_ = 0.f;
}

float Chorus::Process(float in)
{
    sigl_ = 0.f;
    sigr_ = 0.f;

    for(int i = 0; i < 2; i++)
    {
        float sig = engines_[i].Process(in);
        sigl_ += (1.f - pan_[i]) * sig;
        sigr_ += pan_[i] * sig;
    }

    sigl_ *= gain_frac_;
    sigr_ *= gain_frac_;

    return sigl_;
}

float Chorus::GetLeft()
{
    return sigl_;
}

float Chorus::GetRight()
{
    return sigr_;
}

void Chorus::SetPan(float panl, float panr)
{
    pan_[0] = fclamp(panl, 0.f, 1.f);
    pan_[1] = fclamp(panr, 0.f, 1.f);
}

void Chorus::SetPan(float pan)
{
    SetPan(pan, pan);
}

void Chorus::SetLfoDepth(float depthl, float depthr)
{
    engines_[0].SetLfoDepth(depthl);
    engines_[1].SetLfoDepth(depthr);
}

void Chorus::SetLfoDepth(float depth)
{
    SetLfoDepth(depth, depth);
}

void Chorus::SetLfoFreq(float freql, float freqr)
{
    engines_[0].SetLfoFreq(freql);
    engines_[1].SetLfoFreq(freqr);
}

void Chorus::SetLfoFreq(float freq)
{
    SetLfoFreq(freq, freq);
}

void Chorus::SetDelay(float delayl, float delayr)
{
    engines_[0].SetDelay(delayl);
    engines_[1].SetDelay(delayr);
}

void Chorus::SetDelay(float delay)
{
    SetDelay(delay, delay);
}

void Chorus::SetDelayMs(float msl, float msr)
{
    engines_[0].SetDelayMs(msl);
    engines_[1].SetDelayMs(msr);
}

void Chorus::SetDelayMs(float ms)
{
    SetDelayMs(ms, ms);
}

void Chorus::SetFeedback(float feedbackl, float feedbackr)
{
    engines_[0].SetFeedback(feedbackl);
    engines_[1].SetFeedback(feedbackr);
}

void Chorus::SetFeedback(float feedback)
{
    SetFeedback(feedback, feedback);
}



/** Performs downsampling and bitcrush effects
*/
class Decimator
{
  public:
    Decimator() {}
    ~Decimator() {}
    /** Initializes downsample module
    */
    void Init();

    /** Applies downsample and bitcrush effects to input signal.
        \return one sample. This should be called once per sample period. 
    */
    float Process(float input);


    /** Sets amount of downsample 
        Input range: 
    */
    inline void SetDownsampleFactor(float downsample_factor)
    {
        downsample_factor_ = downsample_factor;
    }

    /** Sets amount of bitcrushing 
        Input range: 0...1.0
    */
    inline void SetBitcrushFactor(float bitcrush_factor)
    {
        //            bitcrush_factor_ = bitcrush_factor;
        bits_to_crush_ = (uint32_t)(bitcrush_factor * kMaxBitsToCrush);
    }

    /** Sets the exact number of bits to crush
        0-16 bits
    */
    inline void SetBitsToCrush(const uint8_t &bits)
    {
        bits_to_crush_ = bits <= kMaxBitsToCrush ? bits : kMaxBitsToCrush;
    }


    /** Returns current setting of downsample
    */
    inline float GetDownsampleFactor() { return downsample_factor_; }
    /** Returns current setting of bitcrush
    */
    inline float GetBitcrushFactor() { return bitcrush_factor_; }

  private:
    const uint8_t kMaxBitsToCrush = 16;
    float         downsample_factor_, bitcrush_factor_;
    uint32_t      bits_to_crush_;
    float         downsampled_, bitcrushed_;
    uint32_t      inc_, threshold_;
};


void Decimator::Init()
{
    downsample_factor_ = 1.0f;
    bitcrush_factor_   = 0.0f;
    downsampled_       = 0.0f;
    bitcrushed_        = 0.0f;
    inc_               = 0;
    threshold_         = 0;
}

float Decimator::Process(float input)
{
    int32_t temp;
    //downsample
    threshold_ = (uint32_t)((downsample_factor_ * downsample_factor_) * 96.0f);
    inc_ += 1;
    if(inc_ > threshold_)
    {
        inc_         = 0;
        downsampled_ = input;
    }
    //bitcrush
    temp = (int32_t)(downsampled_ * 65536.0f);
    temp >>= bits_to_crush_; // shift off
    temp <<= bits_to_crush_; // move back with zeros
    bitcrushed_ = (float)temp / 65536.0f;
    return bitcrushed_;
}


/** @brief Flanging Audio Effect 
 *
 * Generates a modulating phase shifted copy of a signal, and recombines
 * with the original to create a 'flanging' sound effect.
 */
class Flanger
{
  public:
    /** Initialize the modules
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get the next sample
        \param in Sample to process
    */
    float Process(float in);

    /** How much of the signal to feedback into the delay line.
        \param feedback Works 0-1.
    */
    void SetFeedback(float feedback);

    /** How much to modulate the delay by.
        \param depth Works 0-1.
    */
    void SetLfoDepth(float depth);

    /** Set lfo frequency.
        \param freq Frequency in Hz
    */
    void SetLfoFreq(float freq);

    /** Set the internal delay rate. 
        \param delay Tuned for 0-1. Maps to .1 to 7 ms.
    */
    void SetDelay(float delay);

    /** Set the delay time in ms.
        \param ms Delay time in ms, .1 to 7 ms.
    */
    void SetDelayMs(float ms);

  private:
    float                    sample_rate_;
    static constexpr int32_t kDelayLength = 960; // 20 ms at 48kHz = .02 * 48000

    float feedback_;

    //triangle lfos
    float lfo_phase_;
    float lfo_freq_;
    float lfo_amp_;

    float delay_;

    DelayLine<float, kDelayLength> del_;

    float ProcessLfo();
};


void Flanger::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    SetFeedback(.2f);

    del_.Init();
    lfo_amp_ = 0.f;
    SetDelay(.75);

    lfo_phase_ = 0.f;
    SetLfoFreq(.3);
    SetLfoDepth(.9);
}

float Flanger::Process(float in)
{
    float lfo_sig = ProcessLfo();
    del_.SetDelay(1.f + lfo_sig + delay_);

    float out = del_.Read();
    del_.Write(in + out * feedback_);

    return (in + out) * .5f; //equal mix
}

void Flanger::SetFeedback(float feedback)
{
    feedback_ = fclamp(feedback, 0.f, 1.f);
    feedback_ *= .97f;
}

void Flanger::SetLfoDepth(float depth)
{
    depth    = fclamp(depth, 0.f, .93f);
    lfo_amp_ = depth * delay_;
}

void Flanger::SetLfoFreq(float freq)
{
    freq = 4.f * freq / sample_rate_;
    freq *= lfo_freq_ < 0.f ? -1.f : 1.f;  //if we're headed down, keep going
    lfo_freq_ = fclamp(freq, -.25f, .25f); //clip at +/- .125 * sr
}

void Flanger::SetDelay(float delay)
{
    delay = (.1f + delay * 6.9); //.1 to 7 ms
    SetDelayMs(delay);
}

void Flanger::SetDelayMs(float ms)
{
    ms     = fmax(.1, ms);
    delay_ = ms * .001f * sample_rate_; //ms to samples

    lfo_amp_ = fmin(lfo_amp_, delay_); //clip this if needed
}

float Flanger::ProcessLfo()
{
    lfo_phase_ += lfo_freq_;

    //wrap around and flip direction
    if(lfo_phase_ > 1.f)
    {
        lfo_phase_ = 1.f - (lfo_phase_ - 1.f);
        lfo_freq_ *= -1.f;
    }
    else if(lfo_phase_ < -1.f)
    {
        lfo_phase_ = -1.f - (lfo_phase_ + 1.f);
        lfo_freq_ *= -1.f;
    }

    return lfo_phase_ * lfo_amp_;
}


/** fold module

Original author(s) : John FFitch, Gabriel Maldonado

Year : 1998

Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Fold
{
  public:
    Fold() {}
    ~Fold() {}
    /** Initializes the fold module.
    */
    void Init();


    /** applies foldover distortion to input 
    */
    float Process(float in);


    /** 
        \param incr : set fold increment
    */
    inline void SetIncrement(float incr) { incr_ = incr; }

  private:
    float incr_, index_, value_;
    int   sample_index_;
};


void Fold::Init()
{
    incr_         = 1000.f;
    sample_index_ = 0;
    index_        = 0.0f;
    value_        = 0.0f;
}

float Fold::Process(float in)
{
    float out;

    if(index_ < sample_index_)
    {
        index_ += incr_;
        out = value_ = in;
    }
    else
    {
        out = value_;
    }

    sample_index_++;
    return out;
}


/** bitcrush module

Original author(s) : Paul Batchelor, 

Ported from soundpipe by Ben Sergentanis, May 2020
*/
class Bitcrush
{
  public:
    Bitcrush() {}
    ~Bitcrush() {}
    /** Initializes the bitcrush module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** bit crushes and downsamples the input
    */
    float Process(float in);


    /** adjusts bitdepth
        \param bitdepth : Sets bit depth, 0...16
    */
    inline void SetBitDepth(int bitdepth) { bit_depth_ = bitdepth; }
    /** adjusts the downsampling frequency
        \param crushrate : Sets rate to downsample to, 0...SampleRate
    */
    inline void SetCrushRate(float crushrate) { crush_rate_ = crushrate; }

  private:
    float sample_rate_, crush_rate_;
    int   bit_depth_;
};


static Fold fold;

void Bitcrush::Init(float sample_rate)
{
    bit_depth_   = 8;
    crush_rate_  = 10000;
    sample_rate_ = sample_rate;
    fold.Init();
}

float Bitcrush::Process(float in)
{
    float bits    = pow(2, bit_depth_);
    float foldamt = sample_rate_ / crush_rate_;
    float out;

    out = in * 65536.0f;
    out += 32768;
    out *= (bits / 65536.0f);
    out = floor(out);
    out *= (65536.0f / bits) - 32768;

    fold.SetIncrement(foldamt);
    out = fold.Process(out);
    out /= 65536.0;

    return out;
}


/**  
       @brief Distortion / Overdrive Module
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Ported from pichenettes/eurorack/plaits/dsp/fx/overdrive.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2014. \n
*/
class Overdrive
{
  public:
    Overdrive() {}
    ~Overdrive() {}

    /** Initializes the module with 0 gain */
    void Init();

    /** Get the next sample
      \param in Input to be overdriven
    */
    float Process(float in);

    /** Set the amount of drive
          \param drive Works from 0-1
      */
    void SetDrive(float drive);

  private:
    float drive_;
    float pre_gain_;
    float post_gain_;
};


void Overdrive::Init()
{
    SetDrive(.5f);
}

float Overdrive::Process(float in)
{
    float pre = pre_gain_ * in;
    return SoftClip(pre) * post_gain_;
}

void Overdrive::SetDrive(float drive)
{
    drive  = fclamp(drive, 0.f, 1.f);
    drive_ = 2.f * drive;

    const float drive_2    = drive_ * drive_;
    const float pre_gain_a = drive_ * 0.5f;
    const float pre_gain_b = drive_2 * drive_2 * drive_ * 24.0f;
    pre_gain_              = pre_gain_a + (pre_gain_b - pre_gain_a) * drive_2;

    const float drive_squashed = drive_ * (2.0f - drive_);
    post_gain_ = 1.0f / SoftClip(0.33f + drive_squashed * (pre_gain_ - 0.33f));
}



/**  
    @brief Single Phaser engine. Used in Phaser.
    @author Ben Sergentanis
*/
class PhaserEngine
{
  public:
    PhaserEngine() {}
    ~PhaserEngine() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get the next sample
        \param in Sample to process
    */
    float Process(float in);

    /** How much to modulate the allpass filter by.
        \param depth Works 0-1.
    */
    void SetLfoDepth(float depth);

    /** Set lfo frequency.
        \param lfo_freq Frequency in Hz
    */
    void SetLfoFreq(float lfo_freq);

    /** Set the allpass frequency
        \param ap_freq Frequency in Hz.
    */
    void SetFreq(float ap_freq);

    /** Set the feedback amount.
        \param feedback Amount from 0-1.
    */
    void SetFeedback(float feedback);

  private:
    float                    sample_rate_;
    static constexpr int32_t kDelayLength
        = 2400; // 50 ms at 48kHz = .05 * 48000

    //triangle lfo
    float lfo_phase_;
    float lfo_freq_;
    float lfo_amp_;

    float os_;

    float feedback_;
    float ap_freq_;

    float deltime_;
    float last_sample_;

    DelayLine<float, kDelayLength> del_;

    float ProcessLfo();
};

//wraps up all of the phaser engines
/**  
    @brief Phaser Effect.
    @author Ben Sergentanis
    @date Jan 2021
*/
class Phaser
{
  public:
    Phaser() {}
    ~Phaser() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next floating point sample.
        \param in Sample to process
    */
    float Process(float in);

    /** Number of allpass stages.
        \param poles Works 1 to 8.
    */
    void SetPoles(int poles);

    /** Set all lfo depths
        \param depth Works 0-1.
    */
    void SetLfoDepth(float depth);

    /** Set all lfo frequencies.
        \param lfo_freq Lfo freq in Hz.
    */
    void SetLfoFreq(float lfo_freq);

    /** Set all channel allpass freq in Hz.
        \param ap_freq Frequency in Hz.
    */
    void SetFreq(float ap_freq);

    /** Set all channels feedback.
        \param feedback Works 0-1.
    */
    void SetFeedback(float feedback);

  private:
    static constexpr int kMaxPoles = 8;
    PhaserEngine         engines_[kMaxPoles];
    float                gain_frac_;
    int                  poles_;
};


//PhaserEngine stuff
void PhaserEngine::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    del_.Init();
    lfo_amp_  = 0.f;
    feedback_ = .2f;
    SetFreq(200.f);

    del_.SetDelay(0.f);

    os_ = 30.f; //30 hertz frequency offset, lower than this introduces crunch
    deltime_ = 0.f;

    last_sample_ = 0.f;
    lfo_phase_   = 0.f;
    SetLfoFreq(.3);
    SetLfoDepth(.9);
}

float PhaserEngine::Process(float in)
{
    float lfo_sig = ProcessLfo();
    fonepole(deltime_, sample_rate_ / (lfo_sig + ap_freq_ + os_), .0001f);

    last_sample_ = del_.Allpass(in + feedback_ * last_sample_, deltime_, .3f);

    return (in + last_sample_) * .5f; //equal mix
}

void PhaserEngine::SetLfoDepth(float depth)
{
    lfo_amp_ = fclamp(depth, 0.f, 1.f);
}

void PhaserEngine::SetLfoFreq(float lfo_freq)
{
    lfo_freq = 4.f * lfo_freq / sample_rate_;
    lfo_freq *= lfo_freq_ < 0.f ? -1.f : 1.f; //if we're headed down, keep going
    lfo_freq_ = fclamp(lfo_freq, -.25f, .25f); //clip at +/- .125 * sr
}

void PhaserEngine::SetFreq(float ap_freq)
{
    ap_freq_ = fclamp(ap_freq, 0.f, 20000.f); //0 - 20kHz
}

void PhaserEngine::SetFeedback(float feedback)
{
    feedback_ = fclamp(feedback, 0.f, .75f);
}

float PhaserEngine::ProcessLfo()
{
    lfo_phase_ += lfo_freq_;

    //wrap around and flip direction
    if(lfo_phase_ > 1.f)
    {
        lfo_phase_ = 1.f - (lfo_phase_ - 1.f);
        lfo_freq_ *= -1.f;
    }
    else if(lfo_phase_ < -1.f)
    {
        lfo_phase_ = -1.f - (lfo_phase_ + 1.f);
        lfo_freq_ *= -1.f;
    }

    return lfo_phase_ * lfo_amp_ * ap_freq_;
}

//Phaser Stuff
void Phaser::Init(float sample_rate)
{
    for(size_t i = 0; i < kMaxPoles; i++)
    {
        engines_[i].Init(sample_rate);
    }

    poles_     = 4;
    gain_frac_ = .5f;
}

float Phaser::Process(float in)
{
    float sig = 0.f;

    for(int i = 0; i < poles_; i++)
    {
        sig += engines_[i].Process(in);
    }

    return sig;
}

void Phaser::SetPoles(int poles)
{
    poles_ = DSY_CLAMP(poles, 1, 8);
}

void Phaser::SetLfoDepth(float depth)
{
    for(int i = 0; i < kMaxPoles; i++)
    {
        engines_[i].SetLfoDepth(depth);
    }
}

void Phaser::SetLfoFreq(float lfo_freq)
{
    for(int i = 0; i < kMaxPoles; i++)
    {
        engines_[i].SetLfoFreq(lfo_freq);
    }
}

void Phaser::SetFreq(float ap_freq)
{
    for(int i = 0; i < kMaxPoles; i++)
    {
        engines_[i].SetFreq(ap_freq);
    }
}

void Phaser::SetFeedback(float feedback)
{
    for(int i = 0; i < kMaxPoles; i++)
    {
        engines_[i].SetFeedback(feedback);
    }
}

static inline uint32_t hash_xs32(uint32_t x)
{
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

inline uint32_t myrand()
{
    static uint32_t seed = 1;
    seed                 = hash_xs32(seed);
    return seed;
}

/**  time-domain pitchshifter

Author: shensley

Based on "Pitch Shifting" from ucsd.edu 

t = 1 - ((s *f) / R)

where:
    s is the size of the delay
    f is the frequency of the lfo
    r is the sample_rate

solving for t = 12.0
f = (12 - 1) * 48000 / SHIFT_BUFFER_SIZE;

\todo - move hash_xs32 and myrand to dsp.h and give appropriate names
*/
class PitchShifter
{
  public:
    PitchShifter() {}
    ~PitchShifter() {}
    /** Initialize pitch shifter
    */
    void Init(float sr)
    {
        force_recalc_ = false;
        sr_           = sr;
        mod_freq_     = 5.0f;
        SetSemitones();
        for(uint8_t i = 0; i < 2; i++)
        {
            gain_[i] = 0.0f;
            d_[i].Init();
            phs_[i].Init(sr, 50, i == 0 ? 0 : PI_F);
        }
        shift_up_ = true;
        del_size_ = SHIFT_BUFFER_SIZE;
        SetDelSize(del_size_);
        fun_ = 0.0f;
    }

    /** process pitch shifter
    */
    float Process(float &in)
    {
        float val, fade1, fade2;
        // First Process delay mod/crossfade
        fade1 = phs_[0].Process();
        fade2 = phs_[1].Process();
        if(prev_phs_a_ > fade1)
        {
            mod_a_amt_ = fun_ * ((float)(myrand() % 255) / 255.0f)
                         * (del_size_ * 0.5f);
            mod_coeff_[0]
                = 0.0002f + (((float)(myrand() % 255) / 255.0f) * 0.001f);
        }
        if(prev_phs_b_ > fade2)
        {
            mod_b_amt_ = fun_ * ((float)(myrand() % 255) / 255.0f)
                         * (del_size_ * 0.5f);
            mod_coeff_[1]
                = 0.0002f + (((float)(myrand() % 255) / 255.0f) * 0.001f);
        }
        slewed_mod_[0] += mod_coeff_[0] * (mod_a_amt_ - slewed_mod_[0]);
        slewed_mod_[1] += mod_coeff_[1] * (mod_b_amt_ - slewed_mod_[1]);
        prev_phs_a_ = fade1;
        prev_phs_b_ = fade2;
        if(shift_up_)
        {
            fade1 = 1.0f - fade1;
            fade2 = 1.0f - fade2;
        }
        mod_[0] = fade1 * (del_size_ - 1);
        mod_[1] = fade2 * (del_size_ - 1);
#ifdef USE_ARM_DSP
        gain_[0] = arm_sin_f32(fade1 * (float)M_PI);
        gain_[1] = arm_sin_f32(fade2 * (float)M_PI);
#else
        gain_[0] = sinf(fade1 * PI_F);
        gain_[1] = sinf(fade2 * PI_F);
#endif

        // Handle Delay Writing
        d_[0].Write(in);
        d_[1].Write(in);
        // Modulate Delay Lines
        //mod_a_amt = mod_b_amt = 0.0f;
        d_[0].SetDelay(mod_[0] + mod_a_amt_);
        d_[1].SetDelay(mod_[1] + mod_b_amt_);
        d_[0].SetDelay(mod_[0] + slewed_mod_[0]);
        d_[1].SetDelay(mod_[1] + slewed_mod_[1]);
        val = 0.0f;
        val += (d_[0].Read() * gain_[0]);
        val += (d_[1].Read() * gain_[1]);
        return val;
    }

    /** sets transposition in semitones
    */
    void SetTransposition(const float &transpose)
    {
        float   ratio;
        uint8_t idx;
        if(transpose_ != transpose || force_recalc_)
        {
            transpose_ = transpose;
            idx        = (uint8_t)fabsf(transpose);
            ratio      = semitone_ratios_[idx % 12];
            ratio *= (uint8_t)(fabsf(transpose) / 12) + 1;
            if(transpose > 0.0f)
            {
                shift_up_ = true;
            }
            else
            {
                shift_up_ = false;
            }
            mod_freq_ = ((ratio - 1.0f) * sr_) / del_size_;
            if(mod_freq_ < 0.0f)
            {
                mod_freq_ = 0.0f;
            }
            phs_[0].SetFreq(mod_freq_);
            phs_[1].SetFreq(mod_freq_);
            if(force_recalc_)
            {
                force_recalc_ = false;
            }
        }
    }

    /** sets delay size changing the timbre of the pitchshifting 
    */
    void SetDelSize(uint32_t size)
    {
        del_size_     = size < SHIFT_BUFFER_SIZE ? size : SHIFT_BUFFER_SIZE;
        force_recalc_ = true;
        SetTransposition(transpose_);
    }

    /** sets an amount of internal random modulation, kind of sounds like tape-flutter
    */
    inline void SetFun(float f) { fun_ = f; }

  private:
    inline void SetSemitones()
    {
        for(size_t i = 0; i < 12; i++)
        {
            semitone_ratios_[i] = powf(2.0f, (float)i / 12);
        }
    }
    typedef DelayLine<float, SHIFT_BUFFER_SIZE> ShiftDelay;
    ShiftDelay                                  d_[2];
    float                                       pitch_shift_, mod_freq_;
    uint32_t                                    del_size_;
    /** lfo stuff
*/
    bool   force_recalc_;
    float  sr_;
    bool   shift_up_;
    Phasor phs_[2];
    float  gain_[2], mod_[2], transpose_;
    float  fun_, mod_a_amt_, mod_b_amt_, prev_phs_a_, prev_phs_b_;
    float  slewed_mod_[2], mod_coeff_[2];
    /** pitch stuff
*/
    float semitone_ratios_[12];
};



/**Delay line for internal reverb use
*/
typedef struct
{
    int    write_pos;         /**< write position */
    int    buffer_size;       /**< buffer size */
    int    read_pos;          /**< read position */
    int    read_pos_frac;     /**< fractional component of read pos */
    int    read_pos_frac_inc; /**< increment for fractional */
    int    dummy;             /**<  dummy var */
    int    seed_val;          /**< randseed */
    int    rand_line_cnt;     /**< number of random lines */
    float  filter_state;      /**< state of filter */
    float *buf;               /**< buffer ptr */
} ReverbScDl;

/** Stereo Reverb

Reverb SC:               Ported from csound/soundpipe

Original author(s):        Sean Costello, Istvan Varga

Year:                    1999, 2005

Ported to soundpipe by:  Paul Batchelor

Ported by:                Stephen Hensley
*/
class ReverbSc
{
  public:
    ReverbSc() {}
    ~ReverbSc() {}
    /** Initializes the reverb module, and sets the sample_rate at which the Process function will be called.
        Returns 0 if all good, or 1 if it runs out of delay times exceed maximum allowed.
    */
    int Init(float sample_rate);

    /** Process the input through the reverb, and updates values of out1, and out2 with the new processed signal.
    */
    int Process(const float &in1, const float &in2, float *out1, float *out2);

    /** controls the reverb time. reverb tail becomes infinite when set to 1.0
        \param fb - sets reverb time. range: 0.0 to 1.0
    */
    inline void SetFeedback(const float &fb) { feedback_ = fb; }
    /** controls the internal dampening filter's cutoff frequency.
        \param freq - low pass frequency. range: 0.0 to sample_rate / 2
    */
    inline void SetLpFreq(const float &freq) { lpfreq_ = freq; }

  private:
    void       NextRandomLineseg(ReverbScDl *lp, int n);
    int        InitDelayLine(ReverbScDl *lp, int n);
    float      feedback_, lpfreq_;
    float      i_sample_rate_, i_pitch_mod_, i_skip_init_;
    float      sample_rate_;
    float      damp_fact_;
    float      prv_lpfreq_;
    int        init_done_;
    ReverbScDl delay_lines_[8];
    float      aux_[DSY_REVERBSC_MAX_SIZE];
};




/* kReverbParams[n][0] = delay time (in seconds)                     */
/* kReverbParams[n][1] = random variation in delay time (in seconds) */
/* kReverbParams[n][2] = random variation frequency (in 1/sec)       */
/* kReverbParams[n][3] = random seed (0 - 32767)                     */

static const float kReverbParams[8][4]
    = {{(2473.0 / DEFAULT_SRATE), 0.0010, 3.100, 1966.0},
       {(2767.0 / DEFAULT_SRATE), 0.0011, 3.500, 29491.0},
       {(3217.0 / DEFAULT_SRATE), 0.0017, 1.110, 22937.0},
       {(3557.0 / DEFAULT_SRATE), 0.0006, 3.973, 9830.0},
       {(3907.0 / DEFAULT_SRATE), 0.0010, 2.341, 20643.0},
       {(4127.0 / DEFAULT_SRATE), 0.0011, 1.897, 22937.0},
       {(2143.0 / DEFAULT_SRATE), 0.0017, 0.891, 29491.0},
       {(1933.0 / DEFAULT_SRATE), 0.0006, 3.221, 14417.0}};

static int DelayLineMaxSamples(float sr, float i_pitch_mod, int n);
//static int InitDelayLine(dsy_reverbsc_dl *lp, int n);
static int         DelayLineBytesAlloc(float sr, float i_pitch_mod, int n);
static const float kOutputGain = 0.35;
static const float kJpScale    = 0.25;

int ReverbSc::Init(float sr)
{
    i_sample_rate_ = sr;
    sample_rate_   = sr;
    feedback_      = 0.97;
    lpfreq_        = 10000;
    i_pitch_mod_   = 1;
    i_skip_init_   = 0;
    damp_fact_     = 1.0;
    prv_lpfreq_    = 0.0;
    init_done_     = 1;
    int i, n_bytes = 0;
    n_bytes = 0;
    for(i = 0; i < 8; i++)
    {
        if(n_bytes > DSY_REVERBSC_MAX_SIZE)
            return 1;
        delay_lines_[i].buf = (aux_) + n_bytes;
        InitDelayLine(&delay_lines_[i], i);
        n_bytes += DelayLineBytesAlloc(sr, 1, i);
    }
    return 0;
}

static int DelayLineMaxSamples(float sr, float i_pitch_mod, int n)
{
    float max_del;

    max_del = kReverbParams[n][0];
    max_del += (kReverbParams[n][1] * (float)i_pitch_mod * 1.125);
    return (int)(max_del * sr + 16.5);
}

static int DelayLineBytesAlloc(float sr, float i_pitch_mod, int n)
{
    int n_bytes = 0;

    n_bytes += (DelayLineMaxSamples(sr, i_pitch_mod, n) * (int)sizeof(float));
    return n_bytes;
}

void ReverbSc::NextRandomLineseg(ReverbScDl *lp, int n)
{
    float prv_del, nxt_del, phs_inc_val;

    /* update random seed */
    if(lp->seed_val < 0)
        lp->seed_val += 0x10000;
    lp->seed_val = (lp->seed_val * 15625 + 1) & 0xFFFF;
    if(lp->seed_val >= 0x8000)
        lp->seed_val -= 0x10000;
    /* length of next segment in samples */
    lp->rand_line_cnt = (int)((sample_rate_ / kReverbParams[n][2]) + 0.5);
    prv_del           = (float)lp->write_pos;
    prv_del -= ((float)lp->read_pos
                + ((float)lp->read_pos_frac / (float)DELAYPOS_SCALE));
    while(prv_del < 0.0)
        prv_del += lp->buffer_size;
    prv_del = prv_del / sample_rate_; /* previous delay time in seconds */
    nxt_del = (float)lp->seed_val * kReverbParams[n][1] / 32768.0;
    /* next delay time in seconds */
    nxt_del = kReverbParams[n][0] + (nxt_del * (float)i_pitch_mod_);
    /* calculate phase increment per sample */
    phs_inc_val           = (prv_del - nxt_del) / (float)lp->rand_line_cnt;
    phs_inc_val           = phs_inc_val * sample_rate_ + 1.0;
    lp->read_pos_frac_inc = (int)(phs_inc_val * DELAYPOS_SCALE + 0.5);
}

int ReverbSc::InitDelayLine(ReverbScDl *lp, int n)
{
    float read_pos;
    /* int     i; */

    /* calculate length of delay line */
    lp->buffer_size = DelayLineMaxSamples(sample_rate_, 1, n);
    lp->dummy       = 0;
    lp->write_pos   = 0;
    /* set random seed */
    lp->seed_val = (int)(kReverbParams[n][3] + 0.5);
    /* set initial delay time */
    read_pos     = (float)lp->seed_val * kReverbParams[n][1] / 32768;
    read_pos     = kReverbParams[n][0] + (read_pos * (float)i_pitch_mod_);
    read_pos     = (float)lp->buffer_size - (read_pos * sample_rate_);
    lp->read_pos = (int)read_pos;
    read_pos     = (read_pos - (float)lp->read_pos) * (float)DELAYPOS_SCALE;
    lp->read_pos_frac = (int)(read_pos + 0.5);
    /* initialise first random line segment */
    NextRandomLineseg(lp, n);
    /* clear delay line to zero */
    lp->filter_state = 0.0;
    for(int i = 0; i < lp->buffer_size; i++)
    {
        lp->buf[i] = 0;
    }
    return REVSC_OK;
}

int ReverbSc::Process(const float &in1,
                      const float &in2,
                      float *      out1,
                      float *      out2)
{
    float       a_in_l, a_in_r, a_out_l, a_out_r;
    float       vm1, v0, v1, v2, am1, a0, a1, a2, frac;
    ReverbScDl *lp;
    int         read_pos;
    uint32_t    n;
    int         buffer_size; /* Local copy */
    float       damp_fact = damp_fact_;

    //if (init_done_ <= 0) return REVSC_NOT_OK;
    if(init_done_ <= 0)
        return REVSC_NOT_OK;

    /* calculate tone filter coefficient if frequency changed */
    if(lpfreq_ != prv_lpfreq_)
    {
        prv_lpfreq_ = lpfreq_;
        damp_fact
            = 2.0f - cosf(prv_lpfreq_ * (2.0f * (float)M_PI) / sample_rate_);
        damp_fact = damp_fact_
            = damp_fact - sqrtf(damp_fact * damp_fact - 1.0f);
    }

    /* calculate "resultant junction pressure" and mix to input signals */

    a_in_l = a_out_l = a_out_r = 0.0;
    for(n = 0; n < 8; n++)
    {
        a_in_l += delay_lines_[n].filter_state;
    }
    a_in_l *= kJpScale;
    a_in_r = a_in_l + in2;
    a_in_l = a_in_l + in1;

    /* loop through all delay lines */

    for(n = 0; n < 8; n++)
    {
        lp          = &delay_lines_[n];
        buffer_size = lp->buffer_size;

        /* send input signal and feedback to delay line */

        lp->buf[lp->write_pos]
            = (float)((n & 1 ? a_in_r : a_in_l) - lp->filter_state);
        if(++lp->write_pos >= buffer_size)
        {
            lp->write_pos -= buffer_size;
        }

        /* read from delay line with cubic interpolation */

        if(lp->read_pos_frac >= DELAYPOS_SCALE)
        {
            lp->read_pos += (lp->read_pos_frac >> DELAYPOS_SHIFT);
            lp->read_pos_frac &= DELAYPOS_MASK;
        }
        if(lp->read_pos >= buffer_size)
            lp->read_pos -= buffer_size;
        read_pos = lp->read_pos;
        frac     = (float)lp->read_pos_frac * (1.0 / (float)DELAYPOS_SCALE);

        /* calculate interpolation coefficients */

        a2 = frac * frac;
        a2 -= 1.0;
        a2 *= (1.0 / 6.0);
        a1 = frac;
        a1 += 1.0;
        a1 *= 0.5;
        am1 = a1 - 1.0;
        a0  = 3.0 * a2;
        a1 -= a0;
        am1 -= a2;
        a0 -= frac;

        /* read four samples for interpolation */

        if(read_pos > 0 && read_pos < (buffer_size - 2))
        {
            vm1 = (float)(lp->buf[read_pos - 1]);
            v0  = (float)(lp->buf[read_pos]);
            v1  = (float)(lp->buf[read_pos + 1]);
            v2  = (float)(lp->buf[read_pos + 2]);
        }
        else
        {
            /* at buffer wrap-around, need to check index */

            if(--read_pos < 0)
                read_pos += buffer_size;
            vm1 = (float)lp->buf[read_pos];
            if(++read_pos >= buffer_size)
                read_pos -= buffer_size;
            v0 = (float)lp->buf[read_pos];
            if(++read_pos >= buffer_size)
                read_pos -= buffer_size;
            v1 = (float)lp->buf[read_pos];
            if(++read_pos >= buffer_size)
                read_pos -= buffer_size;
            v2 = (float)lp->buf[read_pos];
        }
        v0 = (am1 * vm1 + a0 * v0 + a1 * v1 + a2 * v2) * frac + v0;

        /* update buffer read position */

        lp->read_pos_frac += lp->read_pos_frac_inc;

        /* apply feedback gain and lowpass filter */

        v0 *= (float)feedback_;
        v0               = (lp->filter_state - v0) * damp_fact + v0;
        lp->filter_state = v0;

        /* mix to output */

        if(n & 1)
        {
            a_out_r += v0;
        }
        else
        {
            a_out_l += v0;
        }

        /* start next random line segment if current one has reached endpoint */

        if(--(lp->rand_line_cnt) <= 0)
        {
            NextRandomLineseg(lp, n);
        }
    }
    /* someday, use a_out_r for multimono out */

    *out1 = a_out_l * kOutputGain;
    *out2 = a_out_r * kOutputGain;
    return REVSC_OK;
}


/** @file sampleratereducer.h */

/** 
       @brief Sample rate reducer.
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/fx/sample_rate_reducer.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2014. \n
*/
class SampleRateReducer
{
  public:
    SampleRateReducer() {}
    ~SampleRateReducer() {}

    /** Initialize the module */
    void Init();

    /** Get the next floating point sample 
        \param in Sample to be processed.
    */
    float Process(float in);

    /** Set the new sample rate.
        \param Works over 0-1. 1 is full quality, .5 is half sample rate, etc.
    */
    void SetFreq(float frequency);

  private:
    float frequency_;
    float phase_;
    float sample_;
    float previous_sample_;
    float next_sample_;
};

void SampleRateReducer::Init()
{
    frequency_       = .2f;
    phase_           = 0.0f;
    sample_          = 0.0f;
    next_sample_     = 0.0f;
    previous_sample_ = 0.0f;
}

float SampleRateReducer::Process(float in)
{
    float this_sample = next_sample_;
    next_sample_      = 0.f;
    phase_ += frequency_;
    if(phase_ >= 1.0f)
    {
        phase_ -= 1.0f;
        float t = phase_ / frequency_;
        // t = 0: the transition occurred right at this sample.
        // t = 1: the transition occurred at the previous sample.
        // Use linear interpolation to recover the fractional sample.
        float new_sample
            = previous_sample_ + (in - previous_sample_) * (1.0f - t);
        float discontinuity = new_sample - sample_;
        this_sample += discontinuity * ThisBlepSample(t);
        next_sample_ = discontinuity * NextBlepSample(t);
        sample_      = new_sample;
    }
    next_sample_ += sample_;
    previous_sample_ = in;

    return this_sample;
}

void SampleRateReducer::SetFreq(float frequency)
{
    frequency_ = fclamp(frequency, 0.f, 1.f);
}



/**  
    @brief Tremolo effect.
    @author Ben Sergentanis
    @date Jan 2021
    Based on https://christianfloisand.wordpress.com/2012/04/18/coding-some-tremolo/ \n
*/
class Tremolo
{
  public:
    Tremolo() {}
    ~Tremolo() {}

    /** Initializes the module
        \param sample_rate  The sample rate of the audio engine being run.
    */
    void Init(float sample_rate);

    /** 
     \param in Input sample.
     \return Next floating point sample.
    */
    float Process(float in);

    /** Sets the tremolo rate.
       \param freq Tremolo freq in Hz.
    */
    void SetFreq(float freq);

    /** Shape of the modulating lfo 
        \param waveform Oscillator waveform. Use Oscillator::WAVE_SIN for example.
    */
    void SetWaveform(int waveform);

    /** How much to modulate your volume.
        \param depth Works 0-1.
    */
    void SetDepth(float depth);


  private:
    float      sample_rate_, dc_os_;
    Oscillator osc_;
};


void Tremolo::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    osc_.Init(sample_rate_);
    SetDepth(1.f);
    SetFreq(1.f);
}

float Tremolo::Process(float in)
{
    float modsig = dc_os_ + osc_.Process();
    return in * modsig;
}

void Tremolo::SetFreq(float freq)
{
    osc_.SetFreq(freq);
}

void Tremolo::SetWaveform(int waveform)
{
    osc_.SetWaveform(waveform);
}
void Tremolo::SetDepth(float depth)
{
    depth = fclamp(depth, 0.f, 1.f);
    depth *= .5f;
    osc_.SetAmp(depth);
    dc_os_ = 1.f - depth;
}


/**  
           @brief Clocked Noise Module
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Noise processed by a sample and hold running at a target frequency. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/noise/clocked_noise.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class ClockedNoise
{
  public:
    ClockedNoise() {}
    ~ClockedNoise() {}

    /** Initialize module    
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next floating point sample */
    float Process();

    /** Set the frequency at which the next sample is generated.
        \param freq Frequency in Hz
    */
    void SetFreq(float freq);

    /** Calling this forces another random float to be generated */
    void Sync();

  private:
    // Oscillator state.
    float phase_;
    float sample_;
    float next_sample_;

    // For interpolation of parameters.
    float frequency_;

    float sample_rate_;

    static constexpr float kRandFrac = 1.f / (float)RAND_MAX;
};

void ClockedNoise::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    phase_       = 0.0f;
    sample_      = 0.0f;
    next_sample_ = 0.0f;
    frequency_   = 0.001f;
}

float ClockedNoise::Process()
{
    float next_sample = next_sample_;
    float sample      = sample_;

    float this_sample = next_sample;
    next_sample       = 0.0f;

    const float raw_sample = rand() * kRandFrac * 2.0f - 1.0f;
    float       raw_amount = 4.0f * (frequency_ - 0.25f);
    raw_amount             = fclamp(raw_amount, 0.0f, 1.0f);

    phase_ += frequency_;

    if(phase_ >= 1.0f)
    {
        phase_ -= 1.0f;
        float t             = phase_ / frequency_;
        float new_sample    = raw_sample;
        float discontinuity = new_sample - sample;
        this_sample += discontinuity * ThisBlepSample(t);
        next_sample += discontinuity * NextBlepSample(t);
        sample = new_sample;
    }

    next_sample += sample;
    next_sample_ = next_sample;
    sample_      = sample;

    return this_sample + raw_amount * (raw_sample - this_sample);
}

void ClockedNoise::SetFreq(float freq)
{
    freq       = freq / sample_rate_;
    freq       = fclamp(freq, 0.0f, 1.0f);
    frequency_ = freq;
}

void ClockedNoise::Sync()
{
    phase_ = 1.0f;
}


/**  
       @brief Dust Module
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Randomly Clocked Samples \n \n
       Ported from pichenettes/eurorack/plaits/dsp/noise/dust.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n

*/
class Dust
{
  public:
    Dust() {}
    ~Dust() {}

    void Init() { SetDensity(.5f); }

    float Process()
    {
        float inv_density = 1.0f / density_;
        float u           = rand() * kRandFrac;
        if(u < density_)
        {
            return u * inv_density;
        }
        return 0.0f;
    }

    void SetDensity(float density)
    {
        density_ = fclamp(density, 0.f, 1.f);
        density_ = density_ * .3f;
    }

  private:
    float                  density_;
    static constexpr float kRandFrac = 1.f / (float)RAND_MAX;
};


/**  
       @brief Fractal Noise, stacks octaves of a noise source.
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       T is the noise source to use. T must have SetFreq() and Init(sample_rate) functions. \n
       Order is the number of noise sources to stack. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/noise/fractal_random_generator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
       
*/
template <typename T, int order>
class FractalRandomGenerator
{
  public:
    FractalRandomGenerator() {}
    ~FractalRandomGenerator() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate)
    {
        sample_rate_ = sample_rate;

        SetColor(.5f);
        SetFreq(440.f);
        for(int i = 0; i < order; ++i)
        {
            generator_[i].Init(sample_rate_);
        }
    }

    /** Get the next sample. */
    float Process()
    {
        float gain      = 0.5f;
        float sum       = 0.0f;
        float frequency = frequency_;

        for(int i = 0; i < order; ++i)
        {
            generator_[i].SetFreq(frequency);
            sum += generator_[i].Process() * gain;
            gain *= decay_;
            frequency *= 2.0f;
        }

        return sum;
    }

    /** Set the lowest noise frequency.
        \param freq Frequency of the lowest noise source in Hz.
    */
    void SetFreq(float freq) { frequency_ = fclamp(freq, 0.f, sample_rate_); }

    /** Sets the amount of high frequency noise.
        \** Works 0-1. 1 is the brightest, and 0 is the darkest.
    */
    void SetColor(float color) { decay_ = fclamp(color, 0.f, 1.f); }

  private:
    float sample_rate_;
    float frequency_;
    float decay_;

    T generator_[order];
};


/**  
       @brief Granular Oscillator Module. 
       @author Ben Sergentanis
       @date Dec 2020 
       A phase-distorted single cycle sine * another continuously running sine, \n
       the whole thing synced to a main oscillator. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/grainlet_oscillator.h \n 
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/

class GrainletOscillator
{
  public:
    GrainletOscillator() {}
    ~GrainletOscillator() {}

    /** Initialize the oscillator
        \param sample_rate Sample rate of audio engine
    */
    void Init(float sample_rate);

    /** Get the next sample */
    float Process();

    /** Sets the carrier frequency
        \param freq Frequency in Hz
    */
    void SetFreq(float freq);

    /** Sets the formant frequency
        \param freq Frequency in Hz
    */
    void SetFormantFreq(float freq);

    /** Waveshaping
        \param shape Shapes differently from 0-1, 1-2, and > 2.
    */
    void SetShape(float shape);

    /** Sets the amount of formant to bleed through
        \param bleed Works best 0-1
    */
    void SetBleed(float bleed);

  private:
    float Sine(float phase);

    float Carrier(float phase, float shape);

    float Grainlet(float carrier_phase,
                   float formant_phase,
                   float shape,
                   float bleed);

    // Oscillator state.
    float carrier_phase_;
    float formant_phase_;
    float next_sample_;

    // For interpolation of parameters.
    float carrier_frequency_;
    float formant_frequency_;
    float carrier_shape_;
    float carrier_bleed_;

    float new_carrier_shape_;
    float new_carrier_bleed_;

    float sample_rate_;
};


void GrainletOscillator::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    carrier_phase_ = 0.0f;
    formant_phase_ = 0.0f;
    next_sample_   = 0.0f;

    carrier_shape_ = 0.f;
    carrier_bleed_ = 0.f;

    SetFreq(440.f);
    SetFormantFreq(220.f);
    SetShape(.5f);
    SetBleed(.5f);
}

float GrainletOscillator::Process()
{
    float this_sample = next_sample_;
    float next_sample = 0.0f;

    carrier_phase_ += carrier_frequency_;

    if(carrier_phase_ >= 1.0f)
    {
        carrier_phase_ -= 1.0f;
        float reset_time = carrier_phase_ / carrier_frequency_;

        float shape_inc = new_carrier_shape_ - carrier_shape_;
        float bleed_inc = new_carrier_bleed_ - carrier_bleed_;

        float before = Grainlet(
            1.0f,
            formant_phase_ + (1.0f - reset_time) * formant_frequency_,
            new_carrier_shape_ + shape_inc * (1.0f - reset_time),
            new_carrier_bleed_ + bleed_inc * (1.0f - reset_time));

        float after
            = Grainlet(0.0f, 0.0f, new_carrier_shape_, new_carrier_bleed_);

        float discontinuity = after - before;
        this_sample += discontinuity * ThisBlepSample(reset_time);
        next_sample += discontinuity * NextBlepSample(reset_time);
        formant_phase_ = reset_time * formant_frequency_;
    }
    else
    {
        formant_phase_ += formant_frequency_;
        if(formant_phase_ >= 1.0f)
        {
            formant_phase_ -= 1.0f;
        }
    }

    carrier_bleed_ = new_carrier_bleed_;
    carrier_shape_ = new_carrier_shape_;
    next_sample += Grainlet(
        carrier_phase_, formant_phase_, carrier_shape_, carrier_bleed_);
    next_sample_ = next_sample;
    return this_sample;
}

void GrainletOscillator::SetFreq(float freq)
{
    carrier_frequency_ = freq / sample_rate_;
    carrier_frequency_ = carrier_frequency_ > 0.5f ? 0.5f : carrier_frequency_;
}

void GrainletOscillator::SetFormantFreq(float freq)
{
    formant_frequency_ = freq / sample_rate_;
    formant_frequency_ = formant_frequency_ > 0.5f ? 0.5f : formant_frequency_;
}

void GrainletOscillator::SetShape(float shape)
{
    new_carrier_shape_ = shape;
}

void GrainletOscillator::SetBleed(float bleed)
{
    new_carrier_bleed_ = bleed;
}


float GrainletOscillator::Sine(float phase)
{
    return sinf(phase * TWOPI_F);
}

float GrainletOscillator::Carrier(float phase, float shape)
{
    shape *= 3.0f;
    int   shape_integral   = static_cast<int>(shape);
    float shape_fractional = shape - static_cast<float>(shape_integral);

    float t = 1.0f - shape_fractional;

    if(shape_integral == 0)
    {
        phase = phase * (1.0f + t * t * t * 15.0f);
        if(phase >= 1.0f)
        {
            phase = 1.0f;
        }
        phase += 0.75f;
    }
    else if(shape_integral == 1)
    {
        float breakpoint = 0.001f + 0.499f * t * t * t;
        if(phase < breakpoint)
        {
            phase *= (0.5f / breakpoint);
        }
        else
        {
            phase = 0.5f + (phase - breakpoint) * 0.5f / (1.0f - breakpoint);
        }
        phase += 0.75f;
    }
    else
    {
        t     = 1.0f - t;
        phase = 0.25f + phase * (0.5f + t * t * t * 14.5f);
        if(phase >= 0.75f)
            phase = 0.75f;
    }
    return (Sine(phase) + 1.0f) * 0.25f;
}

float GrainletOscillator::Grainlet(float carrier_phase,
                                   float formant_phase,
                                   float shape,
                                   float bleed)
{
    float carrier = Carrier(carrier_phase, shape);
    float formant = Sine(formant_phase);
    return carrier * (formant + bleed) / (1.0f + bleed);
}



/**  
       @brief Random impulse train processed by a resonant filter.
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Noise processed by a sample and hold running at a target frequency. \n \n
       Ported from pichenettes/eurorack/plaits/dsp/noise/particle.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class Particle
{
  public:
    Particle() {}
    ~Particle() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate.
    */
    void Init(float sample_rate);

    /** Get the next sample */
    float Process();

    /** Get the raw noise output. Must call Process() first. */
    float GetNoise();

    /** Set the resonant filter frequency
        \param freq Frequency in Hz
    */
    void SetFreq(float frequency);

    /** Set the filter resonance
        \param resonance Works 0-1
    */
    void SetResonance(float resonance);

    /** How often to randomize filter frequency
        \param freq Frequency in Hz.
    */
    void SetRandomFreq(float freq);

    /** Noise density
        \param Works 0-1.
    */
    void SetDensity(float density);

    /** Overall module gain
        \param Works 0-1.
    */
    void SetGain(float gain);

    /** How much to randomize the set filter frequency.
        \param spread Works over positive numbers.
    */
    void SetSpread(float spread);

    /** Force randomize the frequency.
        \param sync True to randomize freq.
    */
    void SetSync(bool sync);

  private:
    static constexpr float kRandFrac  = 1.f / (float)RAND_MAX;
    static constexpr float kRatioFrac = 1.f / 12.f;
    float                  sample_rate_;
    float aux_, frequency_, density_, gain_, spread_, resonance_;
    bool  sync_;


    float rand_phase_;
    float rand_freq_;


    float pre_gain_;
    Svf   filter_;
};


void Particle::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    sync_ = false;
    aux_  = 0.f;
    SetFreq(440.f);
    resonance_ = .9f;
    density_   = .5f;
    gain_      = 1.f;
    spread_    = 1.f;

    SetRandomFreq(sample_rate_ / 48.f); //48 is the default block size
    rand_phase_ = 0.f;

    pre_gain_ = 0.0f;
    filter_.Init(sample_rate_);
    filter_.SetDrive(.7f);
}

float Particle::Process()
{
    float u = rand() * kRandFrac;
    float s = 0.0f;

    if(u <= density_ || sync_)
    {
        s = u <= density_ ? u * gain_ : s;
        rand_phase_ += rand_freq_;

        if(rand_phase_ >= 1.f || sync_)
        {
            rand_phase_ = rand_phase_ >= 1.f ? rand_phase_ - 1.f : rand_phase_;

            const float u = 2.0f * rand() * kRandFrac - 1.0f;
            const float f
                = fmin(powf(2.f, kRatioFrac * spread_ * u) * frequency_, .25f);
            pre_gain_ = 0.5f / sqrtf(resonance_ * f * sqrtf(density_));
            filter_.SetFreq(f * sample_rate_);
            filter_.SetRes(resonance_);
        }
    }
    aux_ = s;

    filter_.Process(pre_gain_ * s);
    return filter_.Band();
}

float Particle::GetNoise()
{
    return aux_;
}

void Particle::SetFreq(float freq)
{
    freq /= sample_rate_;
    frequency_ = fclamp(freq, 0.f, 1.f);
}

void Particle::SetResonance(float resonance)
{
    resonance_ = fclamp(resonance, 0.f, 1.f);
}

void Particle::SetRandomFreq(float freq)
{
    freq /= sample_rate_;
    rand_freq_ = fclamp(freq, 0.f, 1.f);
}

void Particle::SetDensity(float density)
{
    density_ = fclamp(density * .3f, 0.f, 1.f);
}

void Particle::SetGain(float gain)
{
    gain_ = fclamp(gain, 0.f, 1.f);
}

void Particle::SetSpread(float spread)
{
    spread_ = spread < 0.f ? 0.f : spread;
}

void Particle::SetSync(bool sync)
{
    sync_ = sync;
}


/** fast white noise generator

    I think this came from musicdsp.org at some point
*/
class WhiteNoise
{
  public:
    WhiteNoise() {}
    ~WhiteNoise() {}
    /** Initializes the WhiteNoise object
    */
    void Init()
    {
        amp_      = 1.0f;
        randseed_ = 1;
    }

    /** sets the amplitude of the noise output
    */
    inline void SetAmp(float a) { amp_ = a; }
    /** returns a new sample of noise in the range of -amp_ to amp_
    */
    inline float Process()
    {
        randseed_ *= 16807;
        return (randseed_ * coeff_) * amp_;
    }

  private:
    static constexpr float coeff_ = 4.6566129e-010f;
    float                  amp_;
    int32_t                randseed_;
};


/** 
    Imitates the sound of dripping water via Physical Modeling Synthesis. \n 
    Ported from soundpipe by Ben Sergentanis, May 2020 
    @author Perry Cook  
    @date 2000  
*/
class Drip
{
  public:
    Drip() {}
    ~Drip() {}

    /** 
        Initializes the Drip module. 
        \param sample_rate The sample rate of the audio engine being run.
    \param dettack The period of time over which all sound is stopped.
    */
    void Init(float sample_rate, float dettack);

    /** 
        Process the next floating point sample.
        \param trig If true, begins a new drip.
        \return Next sample.
    */
    float Process(bool trig);

  private:
    float gains0_, gains1_, gains2_, kloop_, dettack_, num_tubes_, damp_,
        shake_max_, freq_, freq1_, freq2_, amp_, snd_level_, outputs00_,
        outputs01_, outputs10_, outputs11_, outputs20_, outputs21_,
        total_energy_, center_freqs0_, center_freqs1_, center_freqs2_,
        num_objects_save_, sound_decay_, system_decay_, finalZ0_, finalZ1_,
        finalZ2_, coeffs01_, coeffs00_, coeffs11_, coeffs10_, coeffs21_,
        coeffs20_, shake_energy_, shake_damp_, shake_max_save_, num_objects_,
        sample_rate_, res_freq0_, res_freq1_, res_freq2_, inputs1_, inputs2_;

    int   my_random(int max);
    float noise_tick();
};



int Drip::my_random(int max)
{
    return (rand() % (max + 1));
}

float Drip::noise_tick()
{
    float temp;
    temp = 1.0f * rand() - 1073741823.5f;
    return temp * (1.0f / 1073741823.0f);
}

void Drip::Init(float sample_rate, float dettack)
{
    sample_rate_ = sample_rate;
    float temp;
    dettack_   = dettack;
    num_tubes_ = 10;
    damp_      = 0.2f;
    shake_max_ = 0.0f;
    freq_      = 450.0f;
    freq1_     = 600.0f;
    freq2_     = 720.0f;
    amp_       = 0.3f;

    snd_level_   = 0.0;
    float tpidsr = 2.0 * PI_F / sample_rate_;

    kloop_     = (sample_rate_ * dettack_);
    outputs00_ = 0.0f;
    outputs01_ = 0.0f;
    outputs10_ = 0.0f;
    outputs11_ = 0.0f;
    outputs20_ = 0.0f;
    outputs21_ = 0.0f;

    total_energy_ = 0.0f;

    center_freqs0_ = res_freq0_ = WUTR_CENTER_FREQ0;
    center_freqs1_ = res_freq1_ = WUTR_CENTER_FREQ1;
    center_freqs2_ = res_freq2_ = WUTR_CENTER_FREQ2;
    num_objects_save_ = num_objects_ = WUTR_NUM_SOURCES;
    sound_decay_                     = WUTR_SOUND_DECAY;
    system_decay_                    = WUTR_SYSTEM_DECAY;
    temp    = logf(WUTR_NUM_SOURCES) * WUTR_GAIN / WUTR_NUM_SOURCES;
    gains0_ = gains1_ = gains2_ = temp;
    coeffs01_                   = WUTR_RESON * WUTR_RESON;
    coeffs00_ = -WUTR_RESON * 2.0f * cosf(WUTR_CENTER_FREQ0 * tpidsr);
    coeffs11_ = WUTR_RESON * WUTR_RESON;
    coeffs10_ = -WUTR_RESON * 2.0f * cosf(WUTR_CENTER_FREQ1 * tpidsr);
    coeffs21_ = WUTR_RESON * WUTR_RESON;
    coeffs20_ = -WUTR_RESON * 2.0f * cosf(WUTR_CENTER_FREQ2 * tpidsr);

    shake_energy_ = amp_ * 1.0f * MAX_SHAKE * 0.1f;
    shake_damp_   = 0.0f;
    if(shake_energy_ > MAX_SHAKE)
        shake_energy_ = MAX_SHAKE;
    shake_max_save_ = 0.0f;
    num_objects_    = 10.0f;
    finalZ0_ = finalZ1_ = finalZ2_ = 0.0f;
}

float Drip::Process(bool trig)
{
    float data;
    float lastOutput;

    float tpidsr = 2.0f * PI_F / sample_rate_;

    if(trig)
    {
        Init(sample_rate_, dettack_);
    }
    if(num_tubes_ != 0.0f && num_tubes_ != num_objects_)
    {
        num_objects_ = num_tubes_;
        if(num_objects_ < 1.0f)
        {
            num_objects_ = 1.0f;
        }
    }
    if(freq_ != 0.0f && freq_ != res_freq0_)
    {
        res_freq0_ = freq_;
        coeffs00_  = -WUTR_RESON * 2.0f * cosf(res_freq0_ * tpidsr);
    }
    if(damp_ != 0.0f && damp_ != shake_damp_)
    {
        shake_damp_   = damp_;
        system_decay_ = WUTR_SYSTEM_DECAY + (shake_damp_ * 0.002f);
    }
    if(shake_max_ != 0.0f && shake_max_ != shake_max_save_)
    {
        shake_max_save_ = shake_max_;
        shake_energy_ += shake_max_save_ * MAX_SHAKE * 0.1f;
        if(shake_energy_ > MAX_SHAKE)
            shake_energy_ = MAX_SHAKE;
    }
    if(freq1_ != 0.0f && freq1_ != res_freq1_)
    {
        res_freq1_ = freq1_;
        coeffs10_  = -WUTR_RESON * 2.0f * cosf(res_freq1_ * tpidsr);
    }
    if(freq2_ != 0.0f && freq2_ != res_freq2_)
    {
        res_freq2_ = freq2_;
        coeffs20_  = -WUTR_RESON * 2.0f * cosf(res_freq2_ * tpidsr);
    }
    if((--kloop_) == 0.0f)
    {
        shake_energy_ = 0.0f;
    }

    float shakeEnergy = shake_energy_;
    float systemDecay = system_decay_;
    float sndLevel    = snd_level_;
    float num_objects = num_objects_;
    float soundDecay  = sound_decay_;
    float inputs0, inputs1, inputs2;

    shakeEnergy *= systemDecay; /* Exponential system decay */

    sndLevel = shakeEnergy;
    if(my_random(32767) < num_objects)
    {
        int j;
        j = my_random(3);
        if(j == 0)
        {
            center_freqs0_ = res_freq1_ * (0.75f + (0.25f * noise_tick()));
            gains0_        = fabsf(noise_tick());
        }
        else if(j == 1)
        {
            center_freqs1_ = res_freq1_ * (1.0f + (0.25f * noise_tick()));
            gains1_        = fabsf(noise_tick());
        }
        else
        {
            center_freqs2_ = res_freq1_ * (1.25f + (0.25f * noise_tick()));
            gains2_        = fabsf(noise_tick());
        }
    }

    gains0_ *= WUTR_RESON;
    if(gains0_ > 0.001f)
    {
        center_freqs0_ *= WUTR_FREQ_SWEEP;
        coeffs00_ = -WUTR_RESON * 2.0f * cosf(center_freqs0_ * tpidsr);
    }
    gains1_ *= WUTR_RESON;
    if(gains1_ > 0.00f)
    {
        center_freqs1_ *= WUTR_FREQ_SWEEP;
        coeffs10_ = -WUTR_RESON * 2.0f * cosf(center_freqs1_ * tpidsr);
    }
    gains2_ *= WUTR_RESON;
    if(gains2_ > 0.001f)
    {
        center_freqs2_ *= WUTR_FREQ_SWEEP;
        coeffs20_ = -WUTR_RESON * 2.0f * cosf(center_freqs2_ * tpidsr);
    }

    sndLevel *= soundDecay;
    inputs0 = sndLevel;
    inputs0 *= noise_tick();
    inputs1 = inputs0 * gains1_;
    inputs2 = inputs0 * gains2_;
    inputs0 *= gains0_;
    inputs0 -= outputs00_ * coeffs00_;
    inputs0 -= outputs01_ * coeffs01_;
    outputs01_ = outputs00_;
    outputs00_ = inputs0;
    data       = gains0_ * outputs00_;
    inputs1 -= outputs10_ * coeffs10_;
    inputs1 -= outputs11_ * coeffs11_;
    outputs11_ = outputs10_;
    outputs10_ = inputs1_;
    data += gains1_ * outputs10_;
    inputs2 -= outputs20_ * coeffs20_;
    inputs2 -= outputs21_ * coeffs21_;
    outputs21_ = outputs20_;
    outputs20_ = inputs2_;
    data += gains2_ * outputs20_;

    finalZ2_ = finalZ1_;
    finalZ1_ = finalZ0_;
    finalZ0_ = data * 4.0f;

    lastOutput = finalZ2_ - finalZ0_;
    lastOutput *= 0.005f;
    shake_energy_ = shakeEnergy;
    snd_level_    = sndLevel;
    return lastOutput;
}


/**  
       @brief Comb filter / KS string.
       @author Ben Sergentanis
       @date Jan 2021 
       "Lite" version of the implementation used in Rings \n \n 
       Ported from pichenettes/eurorack/plaits/dsp/oscillator/formant_oscillator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class String
{
  public:
    String() {}
    ~String() {}

    /** Initialize the module.
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Clear the delay line */
    void Reset();

    /** Get the next floating point sample
        \param in Signal to excite the string.
    */
    float Process(const float in);

    /** Set the string frequency.
        \param freq Frequency in Hz
    */
    void SetFreq(float freq);

    /** Set the string's behavior.
        \param -1 to 0 is curved bridge, 0 to 1 is dispersion.
    */
    void SetNonLinearity(float non_linearity_amount);

    /** Set the string's overall brightness
        \param Works 0-1.
    */
    void SetBrightness(float brightness);

    /** Set the string's decay time.
        \param damping Works 0-1.
    */
    void SetDamping(float damping);


  private:
    static constexpr size_t kDelayLineSize = 1024;

    enum StringNonLinearity
    {
        NON_LINEARITY_CURVED_BRIDGE,
        NON_LINEARITY_DISPERSION
    };

    template <String::StringNonLinearity non_linearity>
    float ProcessInternal(const float in);

    DelayLine<float, kDelayLineSize>     string_;
    DelayLine<float, kDelayLineSize / 4> stretch_;

    float frequency_, non_linearity_amount_, brightness_, damping_;

    float sample_rate_;

    Tone iir_damping_filter_;

    DcBlock dc_blocker_;

    CrossFade crossfade_;

    float dispersion_noise_;
    float curved_bridge_;

    // Very crappy linear interpolation upsampler used for low pitches that
    // do not fit the delay line. Rarely used.
    float src_phase_;
    float out_sample_[2];
};


// We render 4 modes simultaneously since there are enough registers to hold
// all state variables.
/**  
       @brief SVF for use in the Resonator Class \n 
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Ported from pichenettes/eurorack/plaits/dsp/physical_modelling/resonator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
template <int batch_size>
class ResonatorSvf
{
  public:
    enum FilterMode
    {
        LOW_PASS,
        BAND_PASS,
        BAND_PASS_NORMALIZED,
        HIGH_PASS
    };

    ResonatorSvf() {}
    ~ResonatorSvf() {}

    void Init()
    {
        for(int i = 0; i < batch_size; ++i)
        {
            state_1_[i] = state_2_[i] = 0.0f;
        }
    }

    template <FilterMode mode, bool add>
    void Process(const float* f,
                 const float* q,
                 const float* gain,
                 const float  in,
                 float*       out)
    {
        float g[batch_size];
        float r[batch_size];
        float r_plus_g[batch_size];
        float h[batch_size];
        float state_1[batch_size];
        float state_2[batch_size];
        float gains[batch_size];
        for(int i = 0; i < batch_size; ++i)
        {
            g[i]        = fasttan(f[i]);
            r[i]        = 1.0f / q[i];
            h[i]        = 1.0f / (1.0f + r[i] * g[i] + g[i] * g[i]);
            r_plus_g[i] = r[i] + g[i];
            state_1[i]  = state_1_[i];
            state_2[i]  = state_2_[i];
            gains[i]    = gain[i];
        }

        float s_in  = in;
        float s_out = 0.0f;
        for(int i = 0; i < batch_size; ++i)
        {
            const float hp
                = (s_in - r_plus_g[i] * state_1[i] - state_2[i]) * h[i];
            const float bp = g[i] * hp + state_1[i];
            state_1[i]     = g[i] * hp + bp;
            const float lp = g[i] * bp + state_2[i];
            state_2[i]     = g[i] * bp + lp;
            s_out += gains[i] * ((mode == LOW_PASS) ? lp : bp);
        }
        if(add)
        {
            *out++ += s_out;
        }
        else
        {
            *out++ = s_out;
        }

        for(int i = 0; i < batch_size; ++i)
        {
            state_1_[i] = state_1[i];
            state_2_[i] = state_2[i];
        }
    }

  private:
    static constexpr float kPiPow3 = PI_F * PI_F * PI_F;
    static constexpr float kPiPow5 = kPiPow3 * PI_F * PI_F;
    static inline float    fasttan(float f)
    {
        const float a  = 3.260e-01 * kPiPow3;
        const float b  = 1.823e-01 * kPiPow5;
        float       f2 = f * f;
        return f * (PI_F + f2 * (a + b * f2));
    }

    float state_1_[batch_size];
    float state_2_[batch_size];
};


/**         
       @brief Resonant Body Simulation
       @author Ported by Ben Sergentanis 
       @date Jan 2021 
       Ported from pichenettes/eurorack/plaits/dsp/physical_modelling/resonator.h \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n 
*/
class Resonator
{
  public:
    Resonator() {}
    ~Resonator() {}

    /** Initialize the module
        \param position    Offset the phase of the amplitudes. 0-1
        \param resolution Quality vs speed scalar
        \param sample_rate Samplerate of the audio engine being run.
    */
    void Init(float position, int resolution, float sample_rate);

    /** Get the next sample_rate
        \param in The signal to excited the resonant body
    */
    float Process(const float in);

    /** Resonator frequency.
        \param freq Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Changes the general charater of the resonator (stiffness, brightness)
        \param structure Works best from 0-1
    */
    void SetStructure(float structure);

    /** Set the brighness of the resonator
        \param brightness Works best 0-1
    */
    void SetBrightness(float brightness);

    /** How long the resonant body takes to decay.
        \param damping Works best 0-1
    */
    void SetDamping(float damping);

  private:
    int   resolution_;
    float frequency_, brightness_, structure_, damping_;

    static constexpr int   kMaxNumModes   = 24;
    static constexpr int   kModeBatchSize = 4;
    static constexpr float ratiofrac_     = 1.f / 12.f;
    static constexpr float stiff_frac_    = 1.f / 64.f;
    static constexpr float stiff_frac_2   = 1.f / .6f;

    float sample_rate_;

    float CalcStiff(float sig);

    float                        mode_amplitude_[kMaxNumModes];
    ResonatorSvf<kModeBatchSize> mode_filters_[kMaxNumModes / kModeBatchSize];
};


void Resonator::Init(float position, int resolution, float sample_rate)
{
    sample_rate_ = sample_rate;

    SetFreq(440.f);
    SetStructure(.5f);
    SetBrightness(.5f);
    SetDamping(.5f);

    resolution_ = fmin(resolution, kMaxNumModes);

    for(int i = 0; i < resolution; ++i)
    {
        mode_amplitude_[i] = cos(position * TWOPI_F) * 0.25f;
    }

    for(int i = 0; i < kMaxNumModes / kModeBatchSize; ++i)
    {
        mode_filters_[i].Init();
    }
}

inline float NthHarmonicCompensation(int n, float stiffness)
{
    float stretch_factor = 1.0f;
    for(int i = 0; i < n - 1; ++i)
    {
        stretch_factor += stiffness;
        if(stiffness < 0.0f)
        {
            stiffness *= 0.93f;
        }
        else
        {
            stiffness *= 0.98f;
        }
    }
    return 1.0f / stretch_factor;
}

float Resonator::Process(const float in)
{
    //convert Hz to cycles / sample
    float out = 0.f;

    float stiffness  = CalcStiff(structure_);
    float f0         = frequency_ * NthHarmonicCompensation(3, stiffness);
    float brightness = brightness_;

    float harmonic       = f0;
    float stretch_factor = 1.0f;

    float input  = damping_ * 79.7f;
    float q_sqrt = powf(2.f, input * ratiofrac_);

    float q = 500.0f * q_sqrt * q_sqrt;
    brightness *= 1.0f - structure_ * 0.3f;
    brightness *= 1.0f - damping_ * 0.3f;
    float q_loss = brightness * (2.0f - brightness) * 0.85f + 0.15f;

    float mode_q[kModeBatchSize];
    float mode_f[kModeBatchSize];
    float mode_a[kModeBatchSize];
    int   batch_counter = 0;

    ResonatorSvf<kModeBatchSize>* batch_processor = &mode_filters_[0];

    for(int i = 0; i < resolution_; ++i)
    {
        float mode_frequency = harmonic * stretch_factor;
        if(mode_frequency >= 0.499f)
        {
            mode_frequency = 0.499f;
        }
        const float mode_attenuation = 1.0f - mode_frequency * 2.0f;

        mode_f[batch_counter] = mode_frequency;
        mode_q[batch_counter] = 1.0f + mode_frequency * q;
        mode_a[batch_counter] = mode_amplitude_[i] * mode_attenuation;
        ++batch_counter;

        if(batch_counter == kModeBatchSize)
        {
            batch_counter = 0;
            batch_processor
                ->Process<ResonatorSvf<kModeBatchSize>::BAND_PASS, true>(
                    mode_f, mode_q, mode_a, in, &out);
            ++batch_processor;
        }

        stretch_factor += stiffness;
        if(stiffness < 0.0f)
        {
            // Make sure that the partials do not fold back into negative frequencies.
            stiffness *= 0.93f;
        }
        else
        {
            // This helps adding a few extra partials in the highest frequencies.
            stiffness *= 0.98f;
        }
        harmonic += f0;
        q *= q_loss;
    }

    return out;
}

void Resonator::SetFreq(float freq)
{
    frequency_ = freq / sample_rate_;
}

void Resonator::SetStructure(float structure)
{
    structure_ = fmax(fmin(structure, 1.f), 0.f);
}

void Resonator::SetBrightness(float brightness)
{
    brightness_ = fmax(fmin(brightness, 1.f), 0.f);
}

void Resonator::SetDamping(float damping)
{
    damping_ = fmax(fmin(damping, 1.f), 0.f);
}

float Resonator::CalcStiff(float sig)
{
    if(sig < .25f)
    {
        sig = .25 - sig;
        sig = -sig * .25;
    }
    else if(sig < .3f)
    {
        sig = 0.f;
    }
    else if(sig < .9f)
    {
        sig -= .3f;
        sig *= stiff_frac_2;
    }
    else
    {
        sig -= .9f;
        sig *= 10; // div by .1
        sig *= sig;
        sig = 1.5 - cos(sig * PI_F) * .5f;
    }
    return sig;
}



/**  
       @brief Simple modal synthesis voice with a mallet exciter: click -> LPF -> resonator.
       @author Ben Sergentanis
       @date Jan 2021
       The click can be replaced by continuous white noise. \n \n 
       Ported from pichenettes/eurorack/plaits/dsp/physical_modelling/modal_voice.h \n
       and pichenettes/eurorack/plaits/dsp/physical_modelling/modal_voice.cc \n 
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class ModalVoice
{
  public:
    ModalVoice() {}
    ~ModalVoice() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Get the next sample
        \param trigger Strike the resonator. Defaults to false.
    */
    float Process(bool trigger = false);

    /** Continually excite the resonator with noise.
        \param sustain True turns on the noise.
    */
    void SetSustain(bool sustain);

    /** Strike the resonator. */
    void Trig();

    /** Set the resonator root frequency.
        \param freq Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Hit the resonator a bit harder.
        \param accent Works 0-1.
    */
    void SetAccent(float accent);

    /** Changes the general charater of the resonator (stiffness, brightness)
        \param structure Works best from 0-1
    */
    void SetStructure(float structure);

    /** Set the brighness of the resonator, and the noise density.
        \param brightness Works best 0-1
    */
    void SetBrightness(float brightness);

    /** How long the resonant body takes to decay.
        \param damping Works best 0-1
    */
    void SetDamping(float damping);

    /** Get the raw excitation signal. Must call Process() first. */
    float GetAux();

  private:
    float sample_rate_;

    bool  sustain_, trig_;
    float f0_, structure_, brightness_, damping_;
    float density_, accent_;
    float aux_;

    ResonatorSvf<1> excitation_filter_;
    Resonator       resonator_;
    Dust            dust_;
};



void ModalVoice::Init(float sample_rate)
{
    sample_rate_ = sample_rate;
    aux_         = 0.f;

    excitation_filter_.Init();
    resonator_.Init(0.015f, 24, sample_rate_);
    excitation_filter_.Init();
    dust_.Init();

    SetSustain(false);
    SetFreq(440.f);
    SetAccent(.3f);
    SetStructure(.6f);
    SetBrightness(.8f);
    SetDamping(.6f);
}

void ModalVoice::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void ModalVoice::Trig()
{
    trig_ = true;
}

void ModalVoice::SetFreq(float freq)
{
    resonator_.SetFreq(freq);
    f0_ = freq / sample_rate_;
    f0_ = fclamp(f0_, 0.f, .25f);
}

void ModalVoice::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void ModalVoice::SetStructure(float structure)
{
    resonator_.SetStructure(structure);
}

void ModalVoice::SetBrightness(float brightness)
{
    brightness_ = fclamp(brightness, 0.f, 1.f);
    density_    = brightness_ * brightness_;
}

void ModalVoice::SetDamping(float damping)
{
    damping_ = fclamp(damping, 0.f, 1.f);
}

float ModalVoice::GetAux()
{
    return aux_;
}

float ModalVoice::Process(bool trigger)
{
    float brightness = brightness_ + 0.25f * accent_ * (1.0f - brightness_);
    float damping    = damping_ + 0.25f * accent_ * (1.0f - damping_);

    const float range  = sustain_ ? 36.0f : 60.0f;
    const float f      = sustain_ ? 4.0f * f0_ : 2.0f * f0_;
    const float cutoff = fmin(
        f
            * powf(2.f,
                   kOneTwelfth
                       * ((brightness * (2.0f - brightness) - 0.5f) * range)),
        0.499f);
    const float q = sustain_ ? 0.7f : 1.5f;

    float temp = 0.f;
    // Synthesize excitation signal.
    if(sustain_)
    {
        const float dust_f = 0.00005f + 0.99995f * density_ * density_;
        dust_.SetDensity(dust_f);
        temp = dust_.Process() * (4.0f - dust_f * 3.0f) * accent_;
    }
    else if(trigger || trig_)
    {
        const float attenuation = 1.0f - damping * 0.5f;
        const float amplitude   = (0.12f + 0.08f * accent_) * attenuation;
        temp = amplitude * powf(2.f, kOneTwelfth * (cutoff * cutoff * 24.0f))
               / cutoff;
        trig_ = false;
    }

    const float one = 1.0f;
    excitation_filter_.Process<ResonatorSvf<1>::LOW_PASS, false>(
        &cutoff, &q, &one, temp, &temp);

    aux_ = temp;

    resonator_.SetBrightness(brightness);
    resonator_.SetDamping(damping);

    return resonator_.Process(temp);
}


/** The method of natural decay that the algorithm will use.
    - RECURSIVE: 1st order recursive filter, with coefs .5.
    - WEIGHTED_AVERAGE: weighted averaging. 
    */
enum
{
    PLUCK_MODE_RECURSIVE,
    PLUCK_MODE_WEIGHTED_AVERAGE,
    PLUCK_LAST,
};

/** Produces a naturally decaying plucked string or drum sound based on the Karplus-Strong algorithms.

    Ported from soundpipe to DaisySP
    
    This code was originally extracted from the Csound opcode "pluck"

    Original Author(s): Barry Vercoe, John ffitch
    Year: 1991

    Location: OOps/ugens4.c
    */
class Pluck
{
  public:
    Pluck() {}
    ~Pluck() {}
    /** Initializes the Pluck module.

            \param sample_rate: Sample rate of the audio engine being run.
            \param buf: buffer used as an impulse when triggering the Pluck algorithm
            \param npt: number of elementes in buf.
            \param mode: Sets the mode of the algorithm.
    */
    void Init(float sample_rate, float *buf, int32_t npt, int32_t mode);


    /** Processes the waveform to be generated, returning one sample. This should be called once per sample period.
    */
    float Process(float &trig);

    /** 
        Sets the amplitude of the output signal.
        Input range: 0-1?
    */
    inline void SetAmp(float amp) { amp_ = amp; }
    /** Sets the frequency of the output signal in Hz.
        Input range: Any positive value
    */
    inline void SetFreq(float freq) { freq_ = freq; }
    /** Sets the time it takes for a triggered note to end in seconds.
        Input range: 0-1
    */
    inline void SetDecay(float decay) { decay_ = decay; }
    /** Sets the dampening factor applied by the filter (based on PLUCK_MODE)
        Input range: 0-1
    */
    inline void SetDamp(float damp) { damp_ = damp; }
    /** Sets the mode of the algorithm.
    */
    inline void SetMode(int32_t mode) { mode_ = mode; }
    /** Returns the current value for amp.
    */
    inline float GetAmp() { return amp_; }
    /** Returns the current value for freq.
    */
    inline float GetFreq() { return freq_; }
    /** Returns the current value for decay.
    */
    inline float GetDecay() { return decay_; }
    /** Returns the current value for damp.
    */
    inline float GetDamp() { return damp_; }
    /** Returns the current value for mode.
    */
    inline int32_t GetMode() { return mode_; }

  private:
    void    Reinit();
    float   amp_, freq_, decay_, damp_, ifreq_;
    float   sicps_;
    int32_t phs256_, npts_, maxpts_;
    float * buf_;
    float   sample_rate_;
    char    init_;
    int32_t mode_;
};



void Pluck::Reinit()
{
    int    n;
    float  val = 0;
    float *ap  = buf_;
    //npts_ = (int32_t)roundf(decay_ * (float)(maxpts_ - PLUKMIN) + PLUKMIN);
    npts_ = (int32_t)(decay_ * (float)(maxpts_ - PLUKMIN) + PLUKMIN);
    //sicps_ = ((float)npts_ * INTERPFACTOR + INTERPFACTOR/2.0f) * (1.0f / _sr);
    sicps_ = ((float)npts_ * 256.0f + 128.0f) * (1.0f / sample_rate_);
    for(n = npts_; n--;)
    {
        val   = (float)((float)rand() / (float)RAND_MAX);
        *ap++ = (val * 2.0f) - 1.0f;
    }
    phs256_ = 0;
}

void Pluck::Init(float sample_rate, float *buf, int32_t npts, int32_t mode)
{
    amp_         = 0.5f;
    freq_        = 300;
    decay_       = 1.0f;
    sample_rate_ = sample_rate;
    mode_        = mode;

    maxpts_ = npts;
    npts_   = npts;
    buf_    = buf;

    Reinit();
    /* tuned pitch convt */
    sicps_ = (npts * 256.0f + 128.0f) * (1.0f / sample_rate_);
    init_  = 1;
}

float Pluck::Process(float &trig)
{
    float * fp, out;
    int32_t phs256, phsinc, ltwopi, offset;
    float   coeff;

    // unused variable
    // float inv_coeff;

    float frac, diff;
    float dampmin = 0.42f;

    if(trig != 0)
    {
        init_ = 0;
        Reinit();
    }

    if(init_)
    {
        return 0;
    }
    // Set Coeff for mode.
    switch(mode_)
    {
        case PLUCK_MODE_RECURSIVE:
            coeff = ((0.5f - dampmin) * damp_) + dampmin;
            break;
        case PLUCK_MODE_WEIGHTED_AVERAGE:
            coeff = 0.05f + (damp_ * 0.90f);
            break;
        default: coeff = 0.5f; break;
    }

    // variable set but not used
    //inv_coeff = 1.0f - coeff;

    phsinc = (int32_t)(freq_ * sicps_);
    phs256 = phs256_;
    ltwopi = npts_ << 8;
    offset = phs256 >> 8;
    fp     = (float *)buf_ + offset; /* lookup position   */
    diff   = fp[1] - fp[0];
    frac   = (float)(phs256 & 255) / 256.0f; /*  w. interpolation */
    out    = (fp[0] + diff * frac) * amp_;   /*  gives output val */
    if((phs256 += phsinc) >= ltwopi)
    {
        int   nn;
        float preval;
        phs256 -= ltwopi;
        fp     = buf_;
        preval = fp[0];
        fp[0]  = fp[npts_];
        fp++;
        nn = npts_;
        do
        {
            /* 1st order recursive filter*/
            //preval = (*fp + preval) * coeff;
            /* weighted average - stretches decay times */
            switch(mode_)
            {
                case PLUCK_MODE_RECURSIVE:
                    preval = (*fp + preval) * coeff;
                    break;
                case PLUCK_MODE_WEIGHTED_AVERAGE:
                    preval = (*fp * coeff) + (preval * (1.0f - coeff));
                    break;
                default: break;
            }
            *fp++ = preval;
        } while(--nn);
    }
    phs256_ = phs256;
    return out;
}


/** Simplified Pseudo-Polyphonic Pluck Voice

Template Based Pluck Voice, with configurable number of voices and simple pseudo-polyphony.

DC Blocking included to prevent biases from causing unwanted saturation distortion.

**Author**: shensley

**Date Added**: March 2020
*/
template <size_t num_voices>
class PolyPluck
{
  public:
    /** Initializes the PolyPluck instance.
        \param sample_rate: rate in Hz that the Process() function will be called.
    */
    void Init(float sample_rate)
    {
        active_voice_ = 0;
        p_damp_       = 0.95f;
        p_decay_      = 0.75f;
        for(size_t i = 0; i < num_voices; i++)
        {
            plk_[i].Init(sample_rate, plkbuff_[i], 256, PLUCK_MODE_RECURSIVE);
            plk_[i].SetDamp(0.85f);
            plk_[i].SetAmp(0.18f);
            plk_[i].SetDecay(0.85f);
        }
        blk_.Init(sample_rate);
    }

    /** Process function, synthesizes and sums the output of all voices,
        triggering a new voice with frequency of MIDI note number when trig > 0.

        \param trig: value by reference of trig. When trig > 0 a the next voice will be triggered, and trig will be set to 0.
        \param note: MIDI note number for the active_voice.
    */
    float Process(float &trig, float note)
    {
        float sig, tval;
        sig = 0.0f;
        if(trig > 0.0f)
        {
            // increment active voice
            active_voice_ = (active_voice_ + 1) % num_voices;
            // set new voice to new note
            plk_[active_voice_].SetDamp(p_damp_);
            plk_[active_voice_].SetDecay(p_decay_);
            plk_[active_voice_].SetAmp(0.25f);
        }
        plk_[active_voice_].SetFreq(mtof(note));

        for(size_t i = 0; i < num_voices; i++)
        {
            tval = (trig > 0.0f && i == active_voice_) ? 1.0f : 0.0f;
            sig += plk_[i].Process(tval);
        }
        if(trig > 0.0f)
            trig = 0.0f;
        return blk_.Process(sig);
    }

    /** Sets the decay coefficients of the pluck voices. 
        \param p expects 0.0-1.0 input.
    */
    void SetDecay(float p) { p_damp_ = p; }

  private:
    DcBlock blk_;
    Pluck   plk_[num_voices];
    float   plkbuff_[num_voices][256];
    float   p_damp_, p_decay_;
    size_t  active_voice_;
};



/**  
       @brief Extended Karplus-Strong, with all the niceties from Rings 
       @author Ben Sergentanis
       @date Jan 2021
       Ported from pichenettes/eurorack/plaits/dsp/physical_modelling/string_voice.h \n
       and pichenettes/eurorack/plaits/dsp/physical_modelling/string_voice.cc \n
       to an independent module. \n
       Original code written by Emilie Gillet in 2016. \n
*/
class StringVoice
{
  public:
    StringVoice() {}
    ~StringVoice() {}

    /** Initialize the module
        \param sample_rate Audio engine sample rate
    */
    void Init(float sample_rate);

    /** Reset the string oscillator */
    void Reset();

    /** Get the next sample
        \param trigger Strike the string. Defaults to false.
    */
    float Process(bool trigger = false);

    /** Continually excite the string with noise.
        \param sustain True turns on the noise.
    */
    void SetSustain(bool sustain);

    /** Strike the string. */
    void Trig();

    /** Set the string root frequency.
        \param freq Frequency in Hz.
    */
    void SetFreq(float freq);

    /** Hit the string a bit harder. Influences brightness and decay.
        \param accent Works 0-1.
    */
    void SetAccent(float accent);

    /** Changes the string's nonlinearity (string type).
        \param structure Works 0-1. 0-.26 is curved bridge, .26-1 is dispersion.
    */
    void SetStructure(float structure);

    /** Set the brighness of the string, and the noise density.
        \param brightness Works best 0-1
    */
    void SetBrightness(float brightness);

    /** How long the resonant body takes to decay relative to the accent level.
        \param damping Works best 0-1. Full damp is only achieved with full accent.
    */
    void SetDamping(float damping);

    /** Get the raw excitation signal. Must call Process() first. */
    float GetAux();

  private:
    float sample_rate_;

    bool  sustain_, trig_;
    float f0_, brightness_, damping_;
    float density_, accent_;
    float aux_;

    Dust   dust_;
    Svf    excitation_filter_;
    String string_;
    size_t remaining_noise_samples_;
};

void StringVoice::Init(float sample_rate)
{
    sample_rate_ = sample_rate;

    excitation_filter_.Init(sample_rate);
    string_.Init(sample_rate_);
    dust_.Init();
    remaining_noise_samples_ = 0;

    SetSustain(false);
    SetFreq(440.f);
    SetAccent(.8f);
    SetStructure(.7f);
    SetBrightness(.2f);
    SetDamping(.7f);
}

void StringVoice::Reset()
{
    string_.Reset();
}

void StringVoice::SetSustain(bool sustain)
{
    sustain_ = sustain;
}

void StringVoice::Trig()
{
    trig_ = true;
}

void StringVoice::SetFreq(float freq)
{
    string_.SetFreq(freq);
    f0_ = freq / sample_rate_;
    f0_ = fclamp(f0_, 0.f, .25f);
}

void StringVoice::SetAccent(float accent)
{
    accent_ = fclamp(accent, 0.f, 1.f);
}

void StringVoice::SetStructure(float structure)
{
    structure = fclamp(structure, 0.f, 1.f);
    const float non_linearity
        = structure < 0.24f
              ? (structure - 0.24f) * 4.166f
              : (structure > 0.26f ? (structure - 0.26f) * 1.35135f : 0.0f);
    string_.SetNonLinearity(non_linearity);
}

void StringVoice::SetBrightness(float brightness)
{
    brightness_ = fclamp(brightness, 0.f, 1.f);
    density_    = brightness_ * brightness_;
}

void StringVoice::SetDamping(float damping)
{
    damping_ = fclamp(damping, 0.f, 1.f);
}

float StringVoice::GetAux()
{
    return aux_;
}

float StringVoice::Process(bool trigger)
{
    const float brightness = brightness_ + .25 * accent_ * (1.f - brightness_);
    const float damping    = damping_ + .25 * accent_ * (1.f - damping_);

    // Synthesize excitation signal.
    if(trigger || trig_ || sustain_)
    {
        trig_              = false;
        const float range  = 72.0f;
        const float f      = 4.0f * f0_;
        const float cutoff = fmin(
            f
                * powf(2.f,
                       kOneTwelfth * (brightness * (2.0f - brightness) - 0.5f)
                           * range),
            0.499f);
        const float q            = sustain_ ? 1.0f : 0.5f;
        remaining_noise_samples_ = static_cast<size_t>(1.0f / f0_);
        excitation_filter_.SetFreq(cutoff * sample_rate_);
        excitation_filter_.SetRes(q);
    }

    float temp = 0.f;

    if(sustain_)
    {
        const float dust_f = 0.00005f + 0.99995f * density_ * density_;
        dust_.SetDensity(dust_f);
        temp = dust_.Process() * (8.0f - dust_f * 6.0f) * accent_;
    }
    else if(remaining_noise_samples_)
    {
        temp = 2.0f * rand() * kRandFrac - 1.0f;
        remaining_noise_samples_--;
        remaining_noise_samples_ = DSY_MAX(remaining_noise_samples_, 0.f);
    }

    excitation_filter_.Process(temp);
    temp = excitation_filter_.Low();

    aux_ = temp;

    string_.SetBrightness(brightness);
    string_.SetDamping(damping);

    return string_.Process(temp);
}

}