/* 

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

#include "gmp.h"
#include "flint.h"
#include "fft.h"

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_limb_t * limbs, 
                mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs)
{
   mp_size_t i, skip, length = (total_limbs - 1)/coeff_limbs + 1;
   
   for (skip = 0, i = 0; skip + coeff_limbs <= total_limbs; skip += coeff_limbs, i++)
   {
      flint_mpn_zero(poly[i], output_limbs + 1);
      flint_mpn_copyi(poly[i], limbs + skip, coeff_limbs);
   }
   
   if (i < length) 
      flint_mpn_zero(poly[i], output_limbs + 1);
   
   if (total_limbs > skip) 
      flint_mpn_copyi(poly[i], limbs + skip, total_limbs - skip);
   
   return length;
}

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_limb_t * limbs, 
               mp_size_t total_limbs, mp_bitcnt_t bits, mp_size_t output_limbs)
{
   mp_size_t i, coeff_limbs, limbs_left, length = (FLINT_BITS*total_limbs - 1)/bits + 1;
   mp_bitcnt_t shift_bits, top_bits = ((FLINT_BITS - 1) & bits);
   mp_limb_t * limb_ptr;
   mp_limb_t mask;
   
   if (top_bits == 0)
      return fft_split_limbs(poly, limbs, total_limbs, bits/FLINT_BITS, output_limbs);

   coeff_limbs = (bits/FLINT_BITS) + 1;
   mask = (1L<<top_bits) - 1L;
   shift_bits = 0L;
   limb_ptr = limbs;                      
    
   for (i = 0; i < length - 1; i++)
   {
      flint_mpn_zero(poly[i], output_limbs + 1);
      
      if (!shift_bits)
      {
         flint_mpn_copyi(poly[i], limb_ptr, coeff_limbs);
         poly[i][coeff_limbs - 1] &= mask;
         limb_ptr += (coeff_limbs - 1);
         shift_bits += top_bits;
      } else
      {
         mpn_rshift(poly[i], limb_ptr, coeff_limbs, shift_bits);
         limb_ptr += (coeff_limbs - 1);
         shift_bits += top_bits;

         if (shift_bits >= FLINT_BITS)
         {
            limb_ptr++;
            poly[i][coeff_limbs - 1] += (limb_ptr[0] << (FLINT_BITS - (shift_bits - top_bits)));
            shift_bits -= FLINT_BITS; 
         }
         
         poly[i][coeff_limbs - 1] &= mask;
         
      } 
   }
   
   flint_mpn_zero(poly[i], output_limbs + 1);
   
   limbs_left = total_limbs - (limb_ptr - limbs);
   
   if (!shift_bits)
      flint_mpn_copyi(poly[i], limb_ptr, limbs_left);
   else
      mpn_rshift(poly[i], limb_ptr, limbs_left, shift_bits);                   
     
   return length;
}

