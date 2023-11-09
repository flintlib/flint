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

TEST_FUNCTION_START(flint_mpn_mul_fft_main, state)
{
    flint_bitcnt_t depth, w, maxdepth;

    _flint_rand_init_gmp(state);

    maxdepth = (flint_test_multiplier() > 10) ? 12 :
               (flint_test_multiplier() > 1)  ? 11 : 10;

    for (depth = 6; depth <= maxdepth; depth++)
    {
        for (w = 1; w <= 3 - (depth >= 12); w++)
        {
            int iter = 1 + 200*(depth <= 8) + 80*(depth <= 9) + 10*(depth <= 10), i;

            for (i = 0; i < iter; i++)
            {
               mp_size_t n = (UWORD(1)<<depth);
               flint_bitcnt_t bits1 = (n*w - (depth + 1))/2;
               mp_size_t len1 = 2*n + n_randint(state, 2*n) + 1;
               mp_size_t len2 = 2*n + 2 - len1 + n_randint(state, 2*n);

               flint_bitcnt_t b1 = len1*bits1, b2;
               mp_size_t n1, n2;
               mp_size_t j;
               mp_limb_t * i1, *i2, *r1, *r2;

               if (len2 <= 0)
                  len2 = 2*n + n_randint(state, 2*n) + 1;

               flint_set_num_threads(1 + n_randint(state, 4));

               b2 = len2*bits1;

               n1 = (b1 - 1)/FLINT_BITS + 1;
               n2 = (b2 - 1)/FLINT_BITS + 1;

               if (n1 < n2) /* ensure b1 >= b2 */
               {
                  mp_size_t t = n1;
                  flint_bitcnt_t tb = b1;
                  n1 = n2;
                  b1 = b2;
                  n2 = t;
                  b2 = tb;
               }

               i1 = flint_malloc(3*(n1 + n2)*sizeof(mp_limb_t));
               i2 = i1 + n1;
               r1 = i2 + n2;
               r2 = r1 + n1 + n2;

               flint_mpn_urandomb(i1, state->gmp_state, b1);
               flint_mpn_urandomb(i2, state->gmp_state, b2);

               mpn_mul(r2, i1, n1, i2, n2);
               flint_mpn_mul_fft_main(r1, i1, n1, i2, n2);

               for (j = 0; j < n1 + n2; j++)
               {
                   if (r1[j] != r2[j])
                   {
                       flint_printf("error in limb %wd, %wx != %wx\n", j, r1[j], r2[j]);
                       fflush(stdout);
                       flint_abort();
                   }
               }

               flint_free(i1);
            }
        }
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
