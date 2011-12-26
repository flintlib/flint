/* mul_fft -- radix 2 fft routines for MPIR.

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

*/

/******************************************************************************

    Copyright (C) 2009, 2011 William Hart
 
******************************************************************************/

#ifndef FFT_H
#define FFT_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include "mpir.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

void fft_combine_limbs(mp_limb_t * res, mp_limb_t ** poly, long length, 
            mp_size_t coeff_limbs, mp_size_t output_limbs, mp_size_t total_limbs);


void fft_combine_bits(mp_limb_t * res, mp_limb_t ** poly, long length, 
                 mp_bitcnt_t bits, mp_size_t output_limbs, mp_size_t total_limbs);

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_limb_t * limbs, 
            mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs);

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_limb_t * limbs, 
                 mp_size_t total_limbs, mp_bitcnt_t bits, mp_size_t output_limbs);

#ifdef __cplusplus
}
#endif

#endif

