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

void ref_adjust_sqrt2(mpz_t r, mpz_t i1, mpz_t p, mp_size_t i, mp_size_t limbs, mp_size_t w)
{
   mpz_mul_2exp(r, i1, (w/2)*i + i/2);
   if (i & 1)
   {
       mpz_mul_2exp(i1, r, 3*limbs*FLINT_BITS/4);
       mpz_mul_2exp(r, r, limbs*FLINT_BITS/4);
       mpz_sub(r, i1, r);
   }
   mpz_mod(r, r, p);
}

TEST_FUNCTION_START(fft_adjust_sqrt2, state)
{
    mp_size_t c, bits, j, k, n, w, limbs;
    mpz_t p, m2a, m2b, mn1;
    mp_limb_t * nn1, * r1, * temp;

    _flint_rand_init_gmp(state);

    mpz_init(p);
    mpz_init(m2a);
    mpz_init(m2b);
    mpz_init(mn1);

    for (bits = FLINT_BITS; bits < 20*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 10; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;

                limbs = (n*w)/FLINT_BITS;

                for (c = 1; c < 2*n; c+=2)
                {
                    if (n_randint(state, 100) > 2.0 + flint_test_multiplier() * 10)
                        continue;

                    set_p(p, n, w);

                    nn1 = flint_malloc((limbs+1)*sizeof(mp_limb_t));
                    r1 = flint_malloc((limbs+1)*sizeof(mp_limb_t));
                    temp = flint_malloc((limbs+1)*sizeof(mp_limb_t));

                    random_fermat(nn1, state, limbs);
                    fermat_to_mpz(mn1, nn1, limbs);
                    ref_adjust_sqrt2(m2a, mn1, p, c, limbs, w);

                    fft_adjust_sqrt2(r1, nn1, c, limbs, w, temp);
                    fermat_to_mpz(m2b, r1, limbs);

                    mpz_mod(m2a, m2a, p);
                    mpz_mod(m2b, m2b, p);

                    if (mpz_cmp(m2a, m2b) != 0)
                    {
                        flint_printf("FAIL:\n");
                        flint_printf("adjust_sqrt2 error a\n");
                        flint_printf("limbs = %wd\n", limbs);
                        flint_printf("n = %wd, w = %wd, c = %wd\n", n, w, c);
                        gmp_printf("want %Zx\n\n", m2a);
                        gmp_printf("got  %Zx\n", m2b);
                        fflush(stdout);
                        flint_abort();
                    }

                    flint_free(temp);
                    flint_free(nn1);
                    flint_free(r1);
                }
            }
        }
    }

    mpz_clear(p);
    mpz_clear(m2a);
    mpz_clear(m2b);
    mpz_clear(mn1);

    TEST_FUNCTION_END(state);
}
