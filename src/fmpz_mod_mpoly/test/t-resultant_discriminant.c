/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        fmpz_t p;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        fmpz_init_set_ui(p, 100003);
        fmpz_mod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, p);
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(d, ctx);
        fmpz_mod_mpoly_init(d1, ctx);

        fmpz_mod_mpoly_set_str_pretty(f, "a^10*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(d1, "b^200 - 4*a^10*c^100000000000000000000", vars, ctx);
        if (!fmpz_mod_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mod_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(d, ctx);
        fmpz_mod_mpoly_clear(d1, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, r;
        fmpz_mod_poly_t au, bu;
        fmpz_t ru;
        slong len1, len2, exp_bound1, exp_bound2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 1, 100);
        if (fmpz_mod_mpoly_ctx_nvars(ctx) < 1)
        {
            fmpz_mod_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(r, ctx);
        fmpz_mod_poly_init(au, ctx->ffinfo);
        fmpz_mod_poly_init(bu, ctx->ffinfo);
        fmpz_init(ru);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        fmpz_mod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fmpz_mod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        fmpz_mod_mpoly_get_fmpz_mod_poly(au, a, 0, ctx);
        fmpz_mod_mpoly_get_fmpz_mod_poly(bu, b, 0, ctx);

        fmpz_mod_poly_resultant(ru, au, bu, ctx->ffinfo);

        if (!fmpz_mod_mpoly_resultant(r, a, b, 0, ctx) ||
            !fmpz_mod_mpoly_equal_fmpz(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(r, ctx);
        fmpz_mod_poly_clear(au, ctx->ffinfo);
        fmpz_mod_poly_clear(bu, ctx->ffinfo);
        fmpz_mod_mpoly_ctx_clear(ctx);
        fmpz_clear(ru);
    }

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, c, ab, ra, rb, rab, p;
        slong len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 3, 200);

        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(c, ctx);
        fmpz_mod_mpoly_init(ab, ctx);
        fmpz_mod_mpoly_init(ra, ctx);
        fmpz_mod_mpoly_init(rb, ctx);
        fmpz_mod_mpoly_init(rab, ctx);
        fmpz_mod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        len3 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        fmpz_mod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fmpz_mod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);
        fmpz_mod_mpoly_randtest_bound(c, state, len3, exp_bound3, ctx);

        for (j = 0; j < fmpz_mod_mpoly_ctx_nvars(ctx); j++)
        {
            fmpz_mod_mpoly_mul(ab, a, b, ctx);

            if (!fmpz_mod_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(ra, ctx);

            if (!fmpz_mod_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(rb, ctx);

            if (!fmpz_mod_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(rab, ctx);

            fmpz_mod_mpoly_mul(p, ra, rb, ctx);

            if (!fmpz_mod_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(c, ctx);
        fmpz_mod_mpoly_clear(ab, ctx);
        fmpz_mod_mpoly_clear(ra, ctx);
        fmpz_mod_mpoly_clear(rb, ctx);
        fmpz_mod_mpoly_clear(rab, ctx);
        fmpz_mod_mpoly_clear(p, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, ab, r, da, db, dab, p;
        slong len1, len2, exp_bound1, exp_bound2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 3, 200);

        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(ab, ctx);
        fmpz_mod_mpoly_init(da, ctx);
        fmpz_mod_mpoly_init(db, ctx);
        fmpz_mod_mpoly_init(dab, ctx);
        fmpz_mod_mpoly_init(r, ctx);
        fmpz_mod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 4) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        fmpz_mod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        fmpz_mod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        for (j = 0; j < fmpz_mod_mpoly_ctx_nvars(ctx); j++)
        {
            if (fmpz_mod_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (fmpz_mod_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            fmpz_mod_mpoly_mul(ab, a, b, ctx);

            if (!fmpz_mod_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!fmpz_mod_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!fmpz_mod_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!fmpz_mod_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            fmpz_mod_mpoly_mul(p, da, db, ctx);
            fmpz_mod_mpoly_mul(p, p, r, ctx);
            fmpz_mod_mpoly_mul(p, p, r, ctx);

            if (!fmpz_mod_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");

                flint_printf("a: "); fmpz_mod_mpoly_print_pretty(a, NULL, ctx); flint_printf("\n");
                flint_printf("b: "); fmpz_mod_mpoly_print_pretty(b, NULL, ctx); flint_printf("\n");

                flint_printf("disc(a*b): "); fmpz_mod_mpoly_print_pretty(dab, NULL, ctx); flint_printf("\n");
                flint_printf("disc(a): "); fmpz_mod_mpoly_print_pretty(da, NULL, ctx); flint_printf("\n");
                flint_printf("disc(b): "); fmpz_mod_mpoly_print_pretty(db, NULL, ctx); flint_printf("\n");
                flint_printf("res(a, b): "); fmpz_mod_mpoly_print_pretty(r, NULL, ctx); flint_printf("\n");

                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(ab, ctx);
        fmpz_mod_mpoly_clear(da, ctx);
        fmpz_mod_mpoly_clear(db, ctx);
        fmpz_mod_mpoly_clear(dab, ctx);
        fmpz_mod_mpoly_clear(r, ctx);
        fmpz_mod_mpoly_clear(p, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
