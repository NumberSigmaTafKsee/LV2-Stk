#pragma once

#include "synthbase.h"
#include "synthfunctions.h"

namespace SynthLab
{
	/**
	\class AddOscCore
	\ingroup ModuleCores
	\brief
	This is the "blank" core template for compiling your own Cores as dynamic modules.

	Base Class: ModuleCore
	- Overrides the five (5) common functions plus a special getParameters() method to
	return a shared pointer to the parameters structure.
	- NOTE: These functions have identical names as the SynthModules that own them,
	however the arguments are different. ModuleCores use the CoreProcData structure
	for passing arguments into the cores because they are thunk-barrier compliant.
	- This means that the owning SynthModule must prepare this structure and populate it prior to
	function calls. The large majority of this preparation is done in the SynthModule constructor
	and is one-time in nature.

	GUI Parameters: Depends on the type of Core you are implementing
	- GUI parameters are delivered into the core via the thunk-barrier compliant CoreProcData
	argument that is passed into each function identically
	- processInfo.moduleParameters contains a void* version of the GUI parameter structure pointer
	- the Core function casts the GUI parameter pointer prior to usage

	Access to Modulators is done via the thunk-barrier compliant CoreProcData argument
	- processInfo.modulationInputs
	- processInfo.modulationOutputs

	Access to audio buffers (I/O/FM) is done via the thunk-barrier compliant CoreProcData argument
	- processInfo.inputBuffers
	- processInfo.outputBuffers
	- processInfo.fmBuffers

	Construction: Cores follow the same construction pattern
	- set the Module type and name parameters
	- expose the 16 module strings
	- expose the 4 mod knob label strings
	- intialize any internal variables

	Standalone Mode:
	- These objects are designed to be internal members of the outer SynthModule that owns them.
	They may be used in standalone mode without modification, and you will use the CoreProcData
	structure to pass information into the functions.

	Render:
	- Oscillators and Filters:
	1. processes all audio samples in block
	2. processes from input buffer to output buffer using pointers in the CoreProcData argument

	- EGs
	1. renders into the modulation output array that is passed into the function via the
	CoreProcData structure and populates the arrays with index values of:
	2. kEGNormalOutput normal EG output value
	3. kEGBiasedOutput biased EG (sustain level subtracted out)

	- LFOs
	1.renders into the modulation output array that is passed into the function via the
	CoreProcData structure and populates the arrays with index values of:
	2. kLFONormalOutput normal LFO output value
	3. kLFOInvertedOutput 180 degrees out of phase with normal output
	4. kUnipolarFromMax unipolar version of the waveform that sits at +1 when the output amplitude is at 0,
	and descends from +1 downward as the output amplitude is increases; used for tremolo
	5. kUnipolarFromMin unipolar version of the waveform that sits at 0 when the output amplitude is at 0,
	and ascends from 0 upward as the output amplitude is increases; the ordinary unipolar version

	- Filters
	1. renders into the output buffer using pointers in the CoreProcData argument to the render function
	2. processes one block of audio input into one block of audio output per render cycle
	3. processes in mono that is copied to the right channel as dual-mono stereo

	\author Will Pirkle http://www.willpirkle.com
	\remark This object is included and described in further detail in
	Designing Software Synthesizer Plugins in C++ 2nd Ed. by Will Pirkle
	\version Revision : 1.0
	\date Date : 2021 / 04 / 26
	*/
	class AddOscCore : public ModuleCore
	{
	public:
		// --- constructor/destructor
		AddOscCore();					/* C-TOR */
		virtual ~AddOscCore(){}		/* D-TOR */

		virtual bool reset(CoreProcData& processInfo) override;
		virtual bool update(CoreProcData& processInfo) override;
		virtual bool render(CoreProcData& processInfo) override;
		virtual bool doNoteOn(CoreProcData& processInfo) override;
		virtual bool doNoteOff(CoreProcData& processInfo) override;

	protected:
		double sampleRate = 1.0;
		double midiPitch = 0.0;			///< the midi pitch

		// --- timebase
		SynthClock oscClock;
		
		enum { ChebyT2, ChebyT3, ChebyT4, ChebyT5 };
		// --- NOT the most efficient implementation!!!!!
		inline double doChebyWaveshaper(double xn, uint32_t order)
		{
			double xSquared = xn*xn;
			double xCubed = xn*xSquared;
			double x4th = xn*xCubed;
			double x5th = xn*x4th;

			if (order == ChebyT2)
				return (2.0 * xSquared) - 1.0;
			else if (order == ChebyT3)
				return (4.0 * xCubed) - (3.0 * xn);
			else if (order == ChebyT4)
				return (8.0 * x4th) - (8.0 * xSquared) + 1.0;
			else if (order == ChebyT5)
				return (16.0 * x5th) - (20.0 * xCubed) + (5.0 * xn);

			return xn; // not found
		}


	};



} // namespace

