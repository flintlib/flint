/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("derivative....");
    fflush(stdout);

    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        ordering_t ord;
        mp_limb_t modulus;
        slong nvars, len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(fp, ctx);
        nmod_mpoly_init(gp, ctx);
        nmod_mpoly_init(hp, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bits(hp, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(fp, state, len, exp_bits1, ctx);
        nmod_mpoly_randtest_bits(gp, state, len, exp_bits2, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            nmod_mpoly_derivative(hp, h, idx, ctx);
            nmod_mpoly_test(hp, ctx);

            nmod_mpoly_derivative(fp, f, idx, ctx);
            nmod_mpoly_test(fp, ctx);
            nmod_mpoly_derivative(gp, g, idx, ctx);
            nmod_mpoly_test(gp, ctx);

            nmod_mpoly_mul_johnson(t1, f, gp, ctx);
            nmod_mpoly_test(t1, ctx);
            nmod_mpoly_mul_johnson(t2, g, fp, ctx);
            nmod_mpoly_test(t2, ctx);
            nmod_mpoly_add(t1, t1, t2, ctx);
            nmod_mpoly_test(t1, ctx);

            result = nmod_mpoly_equal(hp, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(fp, ctx);
        nmod_mpoly_clear(gp, ctx);
        nmod_mpoly_clear(hp, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, fp, gp, t1, t2;
        ordering_t ord;
        mp_limb_t modulus;
        slong nvars, len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(fp, ctx);
        nmod_mpoly_init(gp, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(fp, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(gp, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);

            nmod_mpoly_derivative(h, h, idx, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_set(fp, f, ctx);
            nmod_mpoly_derivative(fp, fp, idx, ctx);
            nmod_mpoly_test(fp, ctx);
            nmod_mpoly_set(gp, g, ctx);
            nmod_mpoly_derivative(gp, gp, idx, ctx);
            nmod_mpoly_test(gp, ctx);

            nmod_mpoly_mul_johnson(t1, f, gp, ctx);
            nmod_mpoly_test(t1, ctx);
            nmod_mpoly_mul_johnson(t2, g, fp, ctx);
            nmod_mpoly_test(t2, ctx);
            nmod_mpoly_add(t1, t1, t2, ctx);
            nmod_mpoly_test(t1, ctx);

            result = nmod_mpoly_equal(h, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(fp, ctx);
        nmod_mpoly_clear(gp, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

