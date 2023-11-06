/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fft.h"

TEST_FUNCTION_START(fft_ifft_truncate, state)
{
    flint_bitcnt_t depth, w, maxdepth;

    _flint_rand_init_gmp(state);

    maxdepth = (flint_test_multiplier() > 10) ? 12 :
               (flint_test_multiplier() > 1)  ? 11 : 10;

    for (depth = 6; depth <= maxdepth; depth++)
    {
        for (w = 1; w <= 5; w++)
        {
            mp_size_t n = (UWORD(1)<<depth);
            mp_size_t trunc = n_randint(state, 2*n) + 1;
            mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
            mp_size_t size = limbs + 1;
            mp_size_t i;
            mp_limb_t * ptr;
            mp_limb_t ** ii, ** jj, *t1, *t2;

            trunc = 2*((trunc + 1)/2);

            ii = flint_malloc((2*(n + n*size) + 2*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size)
            {
                ii[i] = ptr;
                random_fermat(ii[i], state, limbs);
            }
            t1 = ptr;
            t2 = t1 + size;

            for (i = 0; i < 2*n; i++)
               mpn_normmod_2expp1(ii[i], limbs);

            jj = flint_malloc(2*(n + n*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size)
            {
                jj[i] = ptr;
                flint_mpn_copyi(jj[i], ii[i], size);
            }

            fft_truncate(ii, n, w, &t1, &t2, trunc);
            ifft_truncate(ii, n, w, &t1, &t2, trunc);
            for (i = 0; i < trunc; i++)
            {
                mpn_div_2expmod_2expp1(ii[i], ii[i], limbs, depth + 1);
                mpn_normmod_2expp1(ii[i], limbs);
            }

            for (i = 0; i < trunc; i++)
            {
                if (mpn_cmp(ii[i], jj[i], size) != 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Error in entry %wd\n", i);
                    fflush(stdout);
                    flint_abort();
                }
            }

            flint_free(ii);
            flint_free(jj);
        }
    }

    TEST_FUNCTION_END(state);
}
