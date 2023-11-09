/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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

TEST_FUNCTION_START(mpn_negmod_2expp1, state)
{
    flint_bitcnt_t bits;
    mp_size_t j, k, n, w, limbs;
    mp_limb_t * a, * z;
    mpz_t p, z1, z2;

    _flint_rand_init_gmp(state);

    mpz_init(z1);
    mpz_init(z2);
    mpz_init(p);

    /* normalisation mod p = 2^wn + 1 where B divides nw and n is a power of 2 */
    for (bits = FLINT_BITS; bits < 16*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 32; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;
                limbs = (n*w)/GMP_LIMB_BITS;

                set_p(p, n, w);

                a = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                z = flint_malloc((limbs + 1)*sizeof(mp_limb_t));

                random_fermat(a, state, limbs);
                mpn_normmod_2expp1(a, limbs);

                mpn_negmod_2expp1(z, a, limbs);

                fermat_to_mpz(z1, a, limbs);

                mpz_neg(z1, z1);
                mpz_mod(z1, z1, p);

                fermat_to_mpz(z2, z, limbs);

                if (mpz_cmp(z1, z2) != 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("mpn_negmod_2expp1 error\n");
                    gmp_printf("want %Zx\n\n", z1);
                    gmp_printf("got  %Zx\n", z2);
                    fflush(stdout);
                    flint_abort();
                }

                flint_free(a);
                flint_free(z);
            }
        }
    }

    mpz_clear(p);
    mpz_clear(z2);
    mpz_clear(z1);

    TEST_FUNCTION_END(state);
}
