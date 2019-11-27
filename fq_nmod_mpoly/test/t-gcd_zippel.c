/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_zippel....");
    fflush(stdout);

    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;
        int res;
        flint_bitcnt_t pbits;
        slong deg;

        pbits = 1 + n_randint(state, FLINT_BITS);
        pbits = 1 + n_randint(state, pbits);
        deg = 1 + n_randint(state, 4);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, pbits, deg);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(ca, ctx);
        fq_nmod_mpoly_init(cb, ctx);
        fq_nmod_mpoly_init(cg, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        degbound = 50/(2*ctx->minfo->nvars - 1);
        degbounds = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds1 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds2 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds[j] = n_randint(state, degbound + 1) + 1;
            degbounds1[j] = n_randint(state, degbound + 1) + 1;
            degbounds2[j] = n_randint(state, degbound + 1) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            do {
                fq_nmod_mpoly_randtest_bounds(t, state, len, degbounds, ctx);
            } while (t->length == 0);
            fq_nmod_mpoly_randtest_bounds(a, state, len1, degbounds1, ctx);
            fq_nmod_mpoly_randtest_bounds(b, state, len2, degbounds2, ctx);

            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = fq_nmod_mpoly_gcd_zippel(g, a, b, ctx);

            if (!res)
            {
                continue;
                printf("FAIL\n");
                flint_printf("Check that gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
            fq_nmod_mpoly_assert_canonical(g, ctx);

            if (fq_nmod_mpoly_is_zero(g, ctx))
            {
                if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (!fq_nmod_is_one(g->coeffs + 0, ctx->fqctx))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fq_nmod_mpoly_divides(ca, a, g, ctx);
            res = res && fq_nmod_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = fq_nmod_mpoly_gcd_zippel(cg, ca, cb, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check that cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            if (!fq_nmod_mpoly_is_one(cg, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        flint_free(degbounds);
        flint_free(degbounds1);
        flint_free(degbounds2);

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(ca, ctx);
        fq_nmod_mpoly_clear(cb, ctx);
        fq_nmod_mpoly_clear(cg, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

