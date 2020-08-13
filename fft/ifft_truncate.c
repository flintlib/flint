/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
      
void ifft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, 
                    mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w*n)/FLINT_BITS;
    
    if (trunc == 2*n)
        ifft_radix2(ii, n, w, t1, t2);
    else if (trunc <= n)
    {
        for (i = trunc; i < n; i++)
        {
            mpn_add_n(ii[i], ii[i], ii[i+n], limbs + 1);
            mpn_div_2expmod_2expp1(ii[i], ii[i], limbs, 1);
        }
      
        ifft_truncate1(ii, n/2, 2*w, t1, t2, trunc);

        for (i = 0; i < trunc; i++)
        {
#if HAVE_ADDSUB_N
            mpn_addsub_n(ii[i], ii[i], ii[i], ii[n+i], limbs + 1);
#else
            mpn_add_n(ii[i], ii[i], ii[i], limbs + 1);
            mpn_sub_n(ii[i], ii[i], ii[n+i], limbs + 1);
#endif
        }
    } else
    {
        ifft_radix2(ii, n/2, 2*w, t1, t2);

        for (i = trunc - n; i < n; i++)
        {
            mpn_sub_n(ii[i+n], ii[i], ii[i+n], limbs + 1);
            fft_adjust(*t1, ii[i+n], i, limbs, w);
            mpn_add_n(ii[i], ii[i], ii[i+n], limbs + 1);
            SWAP_PTRS(ii[i+n], *t1);
        }
   
        ifft_truncate1(ii+n, n/2, 2*w, t1, t2, trunc - n);

        for (i = 0; i < trunc - n; i++) 
        {   
            ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
   
            SWAP_PTRS(ii[i],   *t1);
            SWAP_PTRS(ii[n+i], *t2);
        }
    }
}

void ifft_truncate(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, 
                   mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w*n)/GMP_LIMB_BITS;
   
    if (trunc == 2*n)
        ifft_radix2(ii, n, w, t1, t2);
    else if (trunc <= n)
    {
        ifft_truncate(ii, n/2, 2*w, t1, t2, trunc);

        for (i = 0; i < trunc; i++)
            mpn_add_n(ii[i], ii[i], ii[i], limbs + 1);
    } else
    {
        ifft_radix2(ii, n/2, 2*w, t1, t2);

        for (i = trunc - n; i < n; i++)
            fft_adjust(ii[i+n], ii[i], i, limbs, w);
        
        ifft_truncate1(ii+n, n/2, 2*w, t1, t2, trunc - n);

        for (i = 0; i < trunc - n; i++) 
        {   
            ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
   
            SWAP_PTRS(ii[i],   *t1);
            SWAP_PTRS(ii[n+i], *t2);
        }

        for (i = trunc - n; i < n; i++)
            mpn_add_n(ii[i], ii[i], ii[i], limbs + 1);
    }
}
