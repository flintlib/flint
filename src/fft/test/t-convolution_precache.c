/*
    Copyright (C) 2009, 2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fft.h"

TEST_FUNCTION_START(fft_convolution_precache, state)
{
    flint_bitcnt_t depth, w, maxdepth;

    _flint_rand_init_gmp(state);

    maxdepth = (flint_test_multiplier() > 10) ? 13 :
               (flint_test_multiplier() > 1)  ? 12 : 11;

    for (depth = 6; depth <= maxdepth; depth++)
    {
        for (w = 1; w <= 5; w++)
        {
            mp_size_t n = (UWORD(1)<<depth);
            mp_size_t trunc = 2*n + n_randint(state, 2*n) + 1;
            mp_size_t n1 = (UWORD(1)<<(depth/2));
            mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
            mp_size_t size = limbs + 1;
            mp_size_t i, len1, len2;
            mp_limb_t * ptr;
            mp_limb_t ** ii, ** jj, ** ii2, ** jj2, * t1, * t2, * s1, * tt;

            trunc = 2*n1*((trunc + 2*n1 - 1)/(2*n1));
            len1 = n_randint(state, trunc);
            len2 = trunc - len1 + 1;

            ii = flint_malloc((4*(n + n*size) + 5*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size)
            {
                ii[i] = ptr;
                if (i < len1)
                    random_fermat(ii[i], state, limbs);
                else
                    flint_mpn_zero(ii[i], size);
            }
            t1 = ptr;
            t2 = t1 + size;
            s1 = t2 + size;
            tt = s1 + size;

            jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size)
            {
                jj[i] = ptr;

                if (i < len2)
                    random_fermat(jj[i], state, limbs);
                else
                    flint_mpn_zero(jj[i], size);
            }

            for (i = 0; i < 4*n; i++)
            {
               mpn_normmod_2expp1(ii[i], limbs);
               mpn_normmod_2expp1(jj[i], limbs);
            }

            ii2 = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) ii2 + 4*n; i < 4*n; i++, ptr += size)
            {
                ii2[i] = ptr;
                flint_mpn_copyi(ii2[i], ii[i], size);
            }

            jj2 = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
            for (i = 0, ptr = (mp_limb_t *) jj2 + 4*n; i < 4*n; i++, ptr += size)
            {
                jj2[i] = ptr;
                flint_mpn_copyi(jj2[i], jj[i], size);
            }

            fft_precache(jj, depth, limbs, trunc, &t1, &t2, &s1);
            fft_convolution_precache(ii, jj, depth, limbs, trunc, &t1, &t2, &s1, &tt);
            fft_convolution_basic(ii2, jj2, depth, limbs, trunc, &t1, &t2, &s1, &tt);

            for (i = 0; i < trunc; i++)
            {
                if (mpn_cmp(ii[i], ii2[i], size) != 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("n = %wd, trunc = %wd\n", n, trunc);
                    flint_printf("Error in entry %wd\n", i);
                    fflush(stdout);
                    flint_abort();
                }
            }

            flint_free(ii);
            flint_free(jj);
            flint_free(ii2);
            flint_free(jj2);
        }
    }

    TEST_FUNCTION_END(state);
}
