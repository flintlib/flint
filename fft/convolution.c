/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2008-2011 William Hart
    
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft.h"

void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, long depth, 
                              long limbs, long trunc, mp_limb_t ** t1, 
                          mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t * tt)
{
   long n = (1L<<depth), j;
   long w = (limbs*FLINT_BITS)/n;
   long sqrt = (1L<<(depth/2));
   
   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);
      
      fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);
   
      if (ii != jj)
         fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);
         
         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, tt);
      }

      ifft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   } else
   {
      trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt));
      
      fft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
      
      if (ii != jj)
         fft_mfa_truncate_sqrt2_outer(jj, n, w, t1, t2, s1, sqrt, trunc);
      
      fft_mfa_truncate_sqrt2_inner(ii, jj, n, w, t1, t2, s1, sqrt, trunc, tt);
      
      ifft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
   }
}
