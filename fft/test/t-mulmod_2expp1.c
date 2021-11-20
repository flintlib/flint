/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fft.h"

int
main(void)
{
    flint_bitcnt_t depth, w;
    int iters;

    FLINT_TEST_INIT(state);

    flint_printf("mulmod_2expp1....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    for (iters = 0; iters < 100; iters++)
    {
        for (depth = 6; depth <= 18; depth++)
        {
            for (w = 1; w <= 2; w++)
            {
                mp_size_t n = (UWORD(1)<<depth);
                flint_bitcnt_t bits = n*w;
                mp_size_t int_limbs = bits/FLINT_BITS;
                mp_size_t j;
                mp_limb_t c, * i1, * i2, * r1, * r2, * tt;
        
                i1 = flint_malloc(6*(int_limbs+1)*sizeof(mp_limb_t));
                i2 = i1 + int_limbs + 1;
                r1 = i2 + int_limbs + 1;
                r2 = r1 + int_limbs + 1;
                tt = r2 + int_limbs + 1;

                random_fermat(i1, state, int_limbs);
                random_fermat(i2, state, int_limbs);
                mpn_normmod_2expp1(i1, int_limbs);
                mpn_normmod_2expp1(i2, int_limbs);

                fft_mulmod_2expp1(r2, i1, i2, n, w, tt);
                c = 2*i1[int_limbs] + i2[int_limbs];
                c = flint_mpn_mulmod_2expp1_basecase(r1, i1, i2, c, int_limbs*FLINT_BITS, tt);
            
                for (j = 0; j < int_limbs; j++)
                {
                    if (r1[j] != r2[j]) 
                    {
                        flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                if (c != r2[int_limbs])
                {
                    flint_printf("error in limb %wd, %wx != %wx\n", j, c, r2[j]);
                    fflush(stdout);
                    flint_abort();
                }

                flint_free(i1);
            }
        }
    }
    
    /* test squaring */
    for (iters = 0; iters < 100; iters++)
    {
        for (depth = 6; depth <= 18; depth++)
        {
            for (w = 1; w <= 2; w++)
            {
                mp_size_t n = (UWORD(1)<<depth);
                flint_bitcnt_t bits = n*w;
                mp_size_t int_limbs = bits/FLINT_BITS;
                mp_size_t j;
                mp_limb_t c, * i1, * r1, * r2, * tt;
        
                i1 = flint_malloc(5*(int_limbs+1)*sizeof(mp_limb_t));
                r1 = i1 + int_limbs + 1;
                r2 = r1 + int_limbs + 1;
                tt = r2 + int_limbs + 1;

                random_fermat(i1, state, int_limbs);
                mpn_normmod_2expp1(i1, int_limbs);
                
                fft_mulmod_2expp1(r2, i1, i1, n, w, tt);
                c = i1[int_limbs] + 2*i1[int_limbs];
                c = flint_mpn_mulmod_2expp1_basecase(r1, i1, i1, c, int_limbs*FLINT_BITS, tt);
            
                for (j = 0; j < int_limbs; j++)
                {
                    if (r1[j] != r2[j]) 
                    {
                        flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                if (c != r2[int_limbs])
                {
                    flint_printf("error in limb %wd, %wx != %wx\n", j, c, r2[j]);
                    fflush(stdout);
                    flint_abort();
                }

                flint_free(i1);
            }
        }
    }
    
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
