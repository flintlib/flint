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
    
    FLINT_TEST_INIT(state);

    flint_printf("mul_mfa_truncate_sqrt2....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    for (depth = 6; depth <= 13; depth++)
    {
        for (w = 1; w <= 3 - (depth >= 12); w++)
        {
            mp_size_t n = (UWORD(1)<<depth);
            flint_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
            mp_size_t trunc = 2*n + 2*n_randint(state, n) + 2; /* trunc is even */
            flint_bitcnt_t bits = (trunc/2)*bits1;
            mp_size_t int_limbs = (bits - 1)/FLINT_BITS + 1;
            mp_size_t j;
            mp_limb_t * i1, *i2, *r1, *r2;
        
            i1 = flint_malloc(6*int_limbs*sizeof(mp_limb_t));
            i2 = i1 + int_limbs;
            r1 = i2 + int_limbs;
            r2 = r1 + 2*int_limbs;
   
            random_fermat(i1, state, int_limbs);
            random_fermat(i2, state, int_limbs);
            
            mpn_mul(r2, i1, int_limbs, i2, int_limbs);
            mul_mfa_truncate_sqrt2(r1, i1, int_limbs, i2, int_limbs, depth, w);
            
            for (j = 0; j < 2*int_limbs; j++)
            {
                if (r1[j] != r2[j]) 
                {
                    flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                    abort();
                }
            }

            flint_free(i1);
        }
    }

    /* test squaring */
    for (depth = 6; depth <= 13; depth++)
    {
        for (w = 1; w <= 3 - (depth >= 12); w++)
        {
            mp_size_t n = (UWORD(1)<<depth);
            flint_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
            mp_size_t trunc = 2*n + 2*n_randint(state, n) + 2; /* trunc is even */
            flint_bitcnt_t bits = (trunc/2)*bits1;
            mp_size_t int_limbs = (bits - 1)/FLINT_BITS + 1;
            mp_size_t j;
            mp_limb_t * i1, *r1, *r2;
        
            i1 = flint_malloc(5*int_limbs*sizeof(mp_limb_t));
            r1 = i1 + int_limbs;
            r2 = r1 + 2*int_limbs;
   
            random_fermat(i1, state, int_limbs);
            
            mpn_mul(r2, i1, int_limbs, i1, int_limbs);
            mul_mfa_truncate_sqrt2(r1, i1, int_limbs, i1, int_limbs, depth, w);
            
            for (j = 0; j < 2*int_limbs; j++)
            {
                if (r1[j] != r2[j]) 
                {
                    flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                    abort();
                }
            }

            flint_free(i1);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
