/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_srcptr limbs, 
                mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs)
{
   mp_size_t i, skip, length = (total_limbs - 1)/coeff_limbs + 1;
   mp_size_t num = total_limbs/coeff_limbs;

#pragma omp parallel for private(i, skip)
   for (i = 0; i < num; i++)
   {
      skip = i*coeff_limbs;

      flint_mpn_zero(poly[i], output_limbs + 1);
      flint_mpn_copyi(poly[i], limbs + skip, coeff_limbs);
   }

   i = num;
   skip = i*coeff_limbs;
   
   if (i < length) 
      flint_mpn_zero(poly[i], output_limbs + 1);
   
   if (total_limbs > skip) 
      flint_mpn_copyi(poly[i], limbs + skip, total_limbs - skip);
   
   return length;
}

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs, 
               mp_size_t total_limbs, mp_bitcnt_t bits, mp_size_t output_limbs)
{
   mp_size_t i, coeff_limbs, limbs_left, length = (FLINT_BITS*total_limbs - 1)/bits + 1;
   mp_bitcnt_t shift_bits, top_bits = ((FLINT_BITS - 1) & bits);
   mp_srcptr limb_ptr;
   mp_limb_t mask;
   
   if (top_bits == 0)
      return fft_split_limbs(poly, limbs, total_limbs, bits/FLINT_BITS, output_limbs);

   coeff_limbs = (bits/FLINT_BITS) + 1;
   mask = (WORD(1)<<top_bits) - WORD(1);
   shift_bits = WORD(0);
   limb_ptr = limbs;                      
    
#pragma omp parallel for private(i, limb_ptr, shift_bits)
   for (i = 0; i < length - 1; i++)
   {
      flint_mpn_zero(poly[i], output_limbs + 1);
      
      limb_ptr = limbs + i*(coeff_limbs - 1) + (i*top_bits)/FLINT_BITS;
      shift_bits = (i*top_bits) % FLINT_BITS;

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
   
   i = length - 1;
   limb_ptr = limbs + i*(coeff_limbs - 1) + (i*top_bits)/FLINT_BITS;
   shift_bits = (i*top_bits) % FLINT_BITS;

   flint_mpn_zero(poly[i], output_limbs + 1);
   
   limbs_left = total_limbs - (limb_ptr - limbs);
   
   if (!shift_bits)
      flint_mpn_copyi(poly[i], limb_ptr, limbs_left);
   else
      mpn_rshift(poly[i], limb_ptr, limbs_left, shift_bits);                   
     
   return length;
}

