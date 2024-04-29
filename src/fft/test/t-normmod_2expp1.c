/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fft.h"

/* Defined in t-adjust.c, t-adjust_sqrt2.c, t-butterfly.c, t-butterfly_lshB.c,
 * t-butterfly_rshB.c, t-butterfly_sqrt2.c, t-butterfly_twiddle.c,
 * t-div_2expmod_2expp1.c, t-mul_2expmod_2expp1.c, t-negmod_2expp1.c,
 * t-normmod_2expp1.c */
#ifndef set_p
#define set_p set_p
/* set p = 2^wn + 1 */
void set_p(mpz_t p, mp_size_t n, flint_bitcnt_t w)
{
   flint_mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, n*w);
   flint_mpz_add_ui(p, p, 1);
}
#endif

TEST_FUNCTION_START(mpn_normmod_2expp1, state)
{
    flint_bitcnt_t bits;
    mp_size_t j, k, n, w, limbs;
    mp_limb_t * nn;
    mpz_t p, m1, m2;

    mpz_init(m1);
    mpz_init(m2);
    mpz_init(p);

    /* normalisation mod p = 2^wn + 1 where B divides nw and n is a power of 2 */
    for (bits = FLINT_BITS; bits < 32*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 32; j++)
        {
            for (k = 1; k <= GMP_NUMB_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;
                limbs = (n*w)/GMP_LIMB_BITS;

                nn = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                random_fermat(nn, state, limbs);
                fermat_to_mpz(m1, nn, limbs);
                set_p(p, n, w);

                mpn_normmod_2expp1(nn, limbs);
                fermat_to_mpz(m2, nn, limbs);
                mpz_mod(m1, m1, p);

                if (mpz_cmp(m1, m2) != 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("mpn_normmod_2expp1 error\n");
                    gmp_printf("want %Zx\n\n", m1);
                    gmp_printf("got  %Zx\n", m2);
                    fflush(stdout);
                    flint_abort();
                }

                flint_free(nn);
            }
        }
    }

    mpz_clear(m2);
    mpz_clear(m1);
    mpz_clear(p);

    TEST_FUNCTION_END(state);
}
