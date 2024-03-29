/*
    Copyright (C) 2008-2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fft.h"

void fft_precache(mp_limb_t ** jj, slong depth, slong limbs, slong trunc,
                             mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1)
{
   slong n = (WORD(1)<<depth);
   slong w = (limbs*FLINT_BITS)/n;
   slong sqrt = (WORD(1)<<(depth/2));
   slong j, s, t, u, trunc2;

   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);

      fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
         mpn_normmod_2expp1(jj[j], limbs);
   } else
   {
      trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt));

      fft_mfa_truncate_sqrt2(jj, n, w, t1, t2, s1, sqrt, trunc);

      for (j = 0; j < 2*n; j++)
         mpn_normmod_2expp1(jj[j], limbs);

      trunc2 = (trunc - 2*n)/sqrt;

      for (s = 0; s < trunc2; s++)
      {
         t = n_revbin(s, depth - depth/2 + 1);

         for (u = 0; u < sqrt; u++)
         {
            j = 2*n + t*sqrt + u;
            mpn_normmod_2expp1(jj[j], limbs);
         }
      }
   }
}

void fft_convolution_precache(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
                              slong limbs, slong trunc, mp_limb_t ** t1,
                          mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)
{
   slong n = (WORD(1)<<depth), j, s, t, u, trunc2;
   slong w = (limbs*FLINT_BITS)/n;
   slong sqrt = (WORD(1)<<(depth/2));

   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);

      fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
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

      fft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

      for (j = 0; j < 2*n; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
      }

      trunc2 = (trunc - 2*n)/sqrt;

      for (s = 0; s < trunc2; s++)
      {
         t = n_revbin(s, depth - depth/2 + 1);

         for (u = 0; u < sqrt; u++)
         {
            j = 2*n + t*sqrt + u;
            mpn_normmod_2expp1(ii[j], limbs);

            fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
         }
      }

      ifft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   }
}
