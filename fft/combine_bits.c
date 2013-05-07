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

#include "stdlib.h"
#include "gmp.h"
#include "flint.h"
#include "fft.h"

void fft_combine_limbs(mp_limb_t * res, mp_limb_t ** poly, long length, 
            mp_size_t coeff_limbs, mp_size_t output_limbs, mp_size_t total_limbs)
{
   mp_size_t skip, i;
   
   for (skip = 0, i = 0; i < length && skip + output_limbs + 1 <= total_limbs; i++, skip += coeff_limbs) 
      mpn_add(res + skip, res + skip, output_limbs + 1, poly[i], output_limbs); 

   while ((skip < total_limbs) && (i < length))
   {
      mpn_add(res + skip, res + skip, total_limbs - skip, poly[i], FLINT_MIN(total_limbs - skip, output_limbs));
      
      i++;
      
      skip += coeff_limbs;
   }  
}

void fft_combine_bits(mp_limb_t * res, mp_limb_t ** poly, long length, 
                  mp_bitcnt_t bits, mp_size_t output_limbs, mp_size_t total_limbs)
{
   mp_bitcnt_t shift_bits, top_bits = ((FLINT_BITS - 1) & bits);
   mp_size_t coeff_limbs, i;
   mp_limb_t * temp, * limb_ptr, * end;
   
   if (top_bits == 0)
   {
      fft_combine_limbs(res, poly, length, bits/FLINT_BITS, output_limbs, total_limbs);
      return;
   }
   
   coeff_limbs = (bits/FLINT_BITS) + 1;
   temp = flint_malloc((output_limbs + 1)*sizeof(mp_limb_t));
   shift_bits = 0;
   limb_ptr = res;
   end = res + total_limbs;
   
   for (i = 0; i < length && limb_ptr + output_limbs + 1 < end; i++)
   { 
      if (shift_bits)
      {
         mpn_lshift(temp, poly[i], output_limbs + 1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, output_limbs + 1);
      } else
         mpn_add(limb_ptr, limb_ptr, output_limbs + 1, poly[i], output_limbs);

      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);

      if (shift_bits >= FLINT_BITS)
      {
         limb_ptr++;
         shift_bits -= FLINT_BITS;
      }      
   } 

   while (limb_ptr < end && i < length)
   {
      if (shift_bits)
      {
         mpn_lshift(temp, poly[i], output_limbs + 1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, end - limb_ptr);
      } else
         mpn_add_n(limb_ptr, limb_ptr, poly[i], end - limb_ptr);

      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);

      if (shift_bits >= FLINT_BITS)
      {
         limb_ptr++;
         shift_bits -= FLINT_BITS;
      }  

      i++;    
   }
   
   flint_free(temp);     
}
