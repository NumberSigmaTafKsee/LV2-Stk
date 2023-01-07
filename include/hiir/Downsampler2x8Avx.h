/*****************************************************************************

        Downsampler2x8Avx.h
        Ported  Downsampler2x4Sse.h from SSE to AVX by Dario Mambro
        Downsampler2x4Sse.h by Laurent de Soras

Downsamples vectors of 8 float by a factor 2 the input signal, using AVX
instruction set.

This object must be aligned on a 32-byte boundary!

Template parameters:
	- NC: number of coefficients, > 0

--- Legal stuff ---

This program is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it
and/or modify it under the terms of the Do What The Fuck You Want
To Public License, Version 2, as published by Sam Hocevar. See
http://sam.zoy.org/wtfpl/COPYING for more details.

*Tab=3***********************************************************************/



#pragma once
#if ! defined (hiir_Downsampler2x8Avx_HEADER_INCLUDED)
#define hiir_Downsampler2x8Avx_HEADER_INCLUDED

#if defined (_MSC_VER)
	#pragma warning (4 : 4250)
#endif



/*\\\ INCLUDE FILES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

#include "hiir/def.h"
#include "hiir/StageDataAvx.h"

#include <xmmintrin.h>

#include <array>



namespace hiir
{



template <int NC>
class Downsampler2x8Avx
{

	static_assert ((NC > 0), "Number of coefficient must be positive.");

/*\\\ PUBLIC \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

public:

	typedef float DataType;
	static constexpr int _nbr_chn  = 8;
	static constexpr int NBR_COEFS = NC;

	               Downsampler2x8Avx ();
	               Downsampler2x8Avx (const Downsampler2x8Avx <NC> &other) = default;
	               Downsampler2x8Avx (Downsampler2x8Avx <NC> &&other) = default;
	               ~Downsampler2x8Avx ()                            = default;

	Downsampler2x8Avx <NC> &
	               operator = (const Downsampler2x8Avx <NC> &other) = default;
	Downsampler2x8Avx <NC> &
	               operator = (Downsampler2x8Avx <NC> &&other)      = default;

	void           set_coefs (const double coef_arr []);

	hiir_FORCEINLINE __m256
	               process_sample (const float in_ptr [_nbr_chn * 2]);
	hiir_FORCEINLINE __m256
	               process_sample (__m256 in_0, __m256 in_1);
	void           process_block (float out_ptr [], const float in_ptr [], long nbr_spl);

	hiir_FORCEINLINE void
	               process_sample_split (__m256 &low, __m256 &high, const float in_ptr [_nbr_chn * 2]);
	hiir_FORCEINLINE void
	               process_sample_split (__m256 &low, __m256 &high, __m256 in_0, __m256 in_1);
	void           process_block_split (float out_l_ptr [], float out_h_ptr [], const float in_ptr [], long nbr_spl);

	void           clear_buffers ();



/*\\\ PROTECTED \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

protected:



/*\\\ PRIVATE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

	// Stages 0 and 1 contain only input memories
	typedef std::array <StageDataAvx, NBR_COEFS + 2> Filter;

	Filter         _filter; // Should be the first member (thus easier to align)



/*\\\ FORBIDDEN MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

	bool           operator == (const Downsampler2x8Avx <NC> &other) const = delete;
	bool           operator != (const Downsampler2x8Avx <NC> &other) const = delete;

}; // class Downsampler2x8Avx



}  // namespace hiir



#include "hiir/Downsampler2x8Avx.hpp"



#endif   // hiir_Downsampler2x8Avx_HEADER_INCLUDED



/*\\\ EOF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
