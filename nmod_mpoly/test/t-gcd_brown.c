/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 50) + 1;
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul_johnson(a, a, t, ctx);
            nmod_mpoly_mul_johnson(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd_brown(g, a, b, ctx);
            if (!res) {
                continue;
            }

            if (nmod_mpoly_is_zero(g, ctx))
            {
                if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (g->coeffs[0] != UWORD(1))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && nmod_mpoly_divides_monagan_pearce(ca, a, g, ctx);
            res = res && nmod_mpoly_divides_monagan_pearce(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = nmod_mpoly_gcd_brown(cg, ca, cb, ctx);

            if (!res)
                continue;

            if (!nmod_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ca, ctx);
        nmod_mpoly_clear(cb, ctx);
        nmod_mpoly_clear(cg, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

