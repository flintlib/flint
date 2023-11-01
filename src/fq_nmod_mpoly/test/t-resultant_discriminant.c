/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 4, ORD_DEGLEX, 3, 4);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(d, ctx);
        fq_nmod_mpoly_init(d1, ctx);

        fq_nmod_mpoly_set_str_pretty(f, "a^10*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(d1, "b^200 - 4*a^10*c^100000000000000000000", vars, ctx);
        if (!fq_nmod_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fq_nmod_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(d, ctx);
        fq_nmod_mpoly_clear(d1, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

/* fq_nmod_poly_resultant is missing */
#if 0
    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, r;
        fq_nmod_poly_t au, bu;
        fq_nmod_t ru;
        slong len1, len2, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 1, FLINT_BITS, 5);

        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(r, ctx);
        fq_nmod_poly_init(au, ctx->fqctx);
        fq_nmod_poly_init(bu, ctx->fqctx);
        fq_nmod_init(ru, ctx->fqctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        fq_nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fq_nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        fq_nmod_mpoly_get_fq_nmod_poly(au, a, 0, ctx);
        fq_nmod_mpoly_get_fq_nmod_poly(bu, b, 0, ctx);

        fq_nmod_poly_resultant(ru, au, bu);

        if (!fq_nmod_mpoly_resultant(r, a, b, 0, ctx) ||
            !fq_nmod_mpoly_equal_fq_nmod(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(ru, ctx->fqctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(r, ctx);
        fq_nmod_poly_clear(au, ctx->fqctx);
        fq_nmod_poly_clear(bu, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }
#endif

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, c, ab, ra, rb, rab, p;
        slong len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 3, FLINT_BITS, 4);

        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(c, ctx);
        fq_nmod_mpoly_init(ab, ctx);
        fq_nmod_mpoly_init(ra, ctx);
        fq_nmod_mpoly_init(rb, ctx);
        fq_nmod_mpoly_init(rab, ctx);
        fq_nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 12);
        len2 = n_randint(state, 12);
        len3 = n_randint(state, 12);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        fq_nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fq_nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);
        fq_nmod_mpoly_randtest_bound(c, state, len3, exp_bound3, ctx);

        for (j = 0; j < fq_nmod_mpoly_ctx_nvars(ctx); j++)
        {
            fq_nmod_mpoly_mul(ab, a, b, ctx);

            if (!fq_nmod_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            fq_nmod_mpoly_assert_canonical(ra, ctx);

            if (!fq_nmod_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            fq_nmod_mpoly_assert_canonical(rb, ctx);

            if (!fq_nmod_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            fq_nmod_mpoly_assert_canonical(rab, ctx);

            fq_nmod_mpoly_mul(p, ra, rb, ctx);

            if (!fq_nmod_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(c, ctx);
        fq_nmod_mpoly_clear(ab, ctx);
        fq_nmod_mpoly_clear(ra, ctx);
        fq_nmod_mpoly_clear(rb, ctx);
        fq_nmod_mpoly_clear(rab, ctx);
        fq_nmod_mpoly_clear(p, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, ab, r, da, db, dab, p;
        slong len1, len2, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 3, FLINT_BITS, 4);

        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(ab, ctx);
        fq_nmod_mpoly_init(da, ctx);
        fq_nmod_mpoly_init(db, ctx);
        fq_nmod_mpoly_init(dab, ctx);
        fq_nmod_mpoly_init(r, ctx);
        fq_nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 12);
        len2 = n_randint(state, 12);
        exp_bound1 = n_randint(state, 4) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        fq_nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fq_nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        for (j = 0; j < fq_nmod_mpoly_ctx_nvars(ctx); j++)
        {
            if (fq_nmod_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (fq_nmod_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            fq_nmod_mpoly_mul(ab, a, b, ctx);

            if (!fq_nmod_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!fq_nmod_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!fq_nmod_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!fq_nmod_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            fq_nmod_mpoly_mul(p, da, db, ctx);
            fq_nmod_mpoly_mul(p, p, r, ctx);
            fq_nmod_mpoly_mul(p, p, r, ctx);

            if (!fq_nmod_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");

                flint_printf("a: "); fq_nmod_mpoly_print_pretty(a, NULL, ctx); flint_printf("\n");
                flint_printf("b: "); fq_nmod_mpoly_print_pretty(b, NULL, ctx); flint_printf("\n");

                flint_printf("disc(a*b): "); fq_nmod_mpoly_print_pretty(dab, NULL, ctx); flint_printf("\n");
                flint_printf("disc(a): "); fq_nmod_mpoly_print_pretty(da, NULL, ctx); flint_printf("\n");
                flint_printf("disc(b): "); fq_nmod_mpoly_print_pretty(db, NULL, ctx); flint_printf("\n");
                flint_printf("res(a, b): "); fq_nmod_mpoly_print_pretty(r, NULL, ctx); flint_printf("\n");

                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(ab, ctx);
        fq_nmod_mpoly_clear(da, ctx);
        fq_nmod_mpoly_clear(db, ctx);
        fq_nmod_mpoly_clear(dab, ctx);
        fq_nmod_mpoly_clear(r, ctx);
        fq_nmod_mpoly_clear(p, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
