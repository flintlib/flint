/*
    Copyright (C) 2009, 2011 William Hart

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

void ref_butterfly_rshB(mpz_t t, mpz_t u, mpz_t i1, mpz_t i2,
                                 mpz_t p, mp_size_t x, mp_size_t y)
{
   mpz_t mult1, mult2;

   mpz_init(mult1);
   mpz_init(mult2);

   flint_mpz_set_ui(mult1, 1);
   mpz_mul_2exp(mult1, mult1, x*FLINT_BITS);
   mpz_invert(mult1, mult1, p);
   flint_mpz_set_ui(mult2, 1);
   mpz_mul_2exp(mult2, mult2, y*FLINT_BITS);
   mpz_invert(mult2, mult2, p);
   mpz_mul(mult1, mult1, i1);
   mpz_mul(mult2, mult2, i2);
   mpz_add(t, mult1, mult2);
   mpz_sub(u, mult1, mult2);
   mpz_mod(t, t, p);
   mpz_mod(u, u, p);

   mpz_clear(mult1);
   mpz_clear(mult2);
}

TEST_FUNCTION_START(butterfly_rshB, state)
{
    mp_size_t c, bits, j, k, x, y, n, w, limbs;
    mpz_t p, ma, mb, m2a, m2b, mn1, mn2;
    mp_limb_t * nn1, * nn2, * r1, * r2;

    _flint_rand_init_gmp(state);

    mpz_init(p);
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(m2a);
    mpz_init(m2b);
    mpz_init(mn1);
    mpz_init(mn2);

    for (bits = FLINT_BITS; bits < 20*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 10; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;

                limbs = (n*w)/FLINT_BITS;

                for (c = 0; c < limbs; c++)
                {
                    if (n_randint(state, 100) > 2.0 + flint_test_multiplier() * 10)
                        continue;

                    x = n_randint(state, limbs);
                    y = n_randint(state, limbs);

                    nn1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    nn2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    flint_mpn_rrandom(nn1, state->gmp_state, limbs);
                    random_fermat(nn1, state, limbs);
                    random_fermat(nn2, state, limbs);

                    fermat_to_mpz(mn1, nn1, limbs);
                    fermat_to_mpz(mn2, nn2, limbs);
                    set_p(p, n, w);

                    butterfly_rshB(r1, r2, nn1, nn2, limbs, x, y);
                    fermat_to_mpz(m2a, r1, limbs);
                    fermat_to_mpz(m2b, r2, limbs);

                    mpz_mod(m2a, m2a, p);
                    mpz_mod(m2b, m2b, p);
                    ref_butterfly_rshB(ma, mb, mn1, mn2, p, x, y);

                    if (mpz_cmp(ma, m2a) != 0)
                    {
                        flint_printf("FAIL:\n");
                        flint_printf("butterfly_rshB error a\n");
                        flint_printf("x = %wd, y = %wd, limbs = %wd\n", x, y, limbs);
                        flint_printf("n = %wd, w = %wd, k = %wd, c = %wd\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", ma);
                        gmp_printf("got  %Zx\n", m2a);
                        fflush(stdout);
                        flint_abort();
                    }
                    if (mpz_cmp(mb, m2b) != 0)
                    {
                        flint_printf("FAIL:\n");
                        flint_printf("butterfly_rshB error b\n");
                        flint_printf("x = %wd, y = %wd, limbs = %wd\n", x, y, limbs);
                        flint_printf("n = %wd, w = %wd, k = %wd, c = %wd\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", mb);
                        gmp_printf("got  %Zx\n", m2b);
                        fflush(stdout);
                        flint_abort();
                    }

                    flint_free(nn1);
                    flint_free(nn2);
                    flint_free(r1);
                    flint_free(r2);
                }
            }
        }
    }

    mpz_clear(p);
    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(m2a);
    mpz_clear(m2b);
    mpz_clear(mn1);
    mpz_clear(mn2);

    TEST_FUNCTION_END(state);
}
