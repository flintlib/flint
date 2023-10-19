/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        fmpq_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(d, ctx);
        fmpq_mpoly_init(d1, ctx);

        fmpq_mpoly_set_str_pretty(f, "a^10/3*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        fmpq_mpoly_set_str_pretty(d1, "b^200 - 4/3*a^10*c^100000000000000000000", vars, ctx);
        if (!fmpq_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpq_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(d, ctx);
        fmpq_mpoly_clear(d1, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, r;
        fmpq_poly_t au, bu;
        fmpq_t ru;
        slong len1, len2, exp_bound1, exp_bound2;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init(ctx, 1, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(r, ctx);
        fmpq_poly_init(au);
        fmpq_poly_init(bu);
        fmpq_init(ru);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);

        fmpq_mpoly_get_fmpq_poly(au, a, 0, ctx);
        fmpq_mpoly_get_fmpq_poly(bu, b, 0, ctx);

        fmpq_poly_resultant(ru, au, bu);

        if (!fmpq_mpoly_resultant(r, a, b, 0, ctx) ||
            !fmpq_mpoly_equal_fmpq(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(ru);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(r, ctx);
        fmpq_poly_clear(au);
        fmpq_poly_clear(bu);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, c, ab, ra, rb, rab, p;
        slong len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 3);

        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(c, ctx);
        fmpq_mpoly_init(ab, ctx);
        fmpq_mpoly_init(ra, ctx);
        fmpq_mpoly_init(rb, ctx);
        fmpq_mpoly_init(rab, ctx);
        fmpq_mpoly_init(p, ctx);

        len1 = n_randint(state, 12);
        len2 = n_randint(state, 12);
        len3 = n_randint(state, 12);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);
        fmpq_mpoly_randtest_bound(c, state, len3, coeff_bits, exp_bound3, ctx);

        for (j = 0; j < fmpq_mpoly_ctx_nvars(ctx); j++)
        {
            fmpq_mpoly_mul(ab, a, b, ctx);

            if (!fmpq_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            fmpq_mpoly_assert_canonical(ra, ctx);

            if (!fmpq_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            fmpq_mpoly_assert_canonical(rb, ctx);

            if (!fmpq_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            fmpq_mpoly_assert_canonical(rab, ctx);

            fmpq_mpoly_mul(p, ra, rb, ctx);

            if (!fmpq_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(c, ctx);
        fmpq_mpoly_clear(ab, ctx);
        fmpq_mpoly_clear(ra, ctx);
        fmpq_mpoly_clear(rb, ctx);
        fmpq_mpoly_clear(rab, ctx);
        fmpq_mpoly_clear(p, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, ab, r, da, db, dab, p;
        slong len1, len2, exp_bound1, exp_bound2;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 3);

        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(ab, ctx);
        fmpq_mpoly_init(da, ctx);
        fmpq_mpoly_init(db, ctx);
        fmpq_mpoly_init(dab, ctx);
        fmpq_mpoly_init(r, ctx);
        fmpq_mpoly_init(p, ctx);

        len1 = n_randint(state, 12);
        len2 = n_randint(state, 12);
        exp_bound1 = n_randint(state, 4) + 1;
        exp_bound2 = n_randint(state, 4) + 1;
        coeff_bits = n_randint(state, 80) + 1;
        fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);

        for (j = 0; j < fmpq_mpoly_ctx_nvars(ctx); j++)
        {
            if (fmpq_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (fmpq_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            fmpq_mpoly_mul(ab, a, b, ctx);

            if (!fmpq_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!fmpq_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!fmpq_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!fmpq_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            fmpq_mpoly_mul(p, da, db, ctx);
            fmpq_mpoly_mul(p, p, r, ctx);
            fmpq_mpoly_mul(p, p, r, ctx);

            if (!fmpq_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");

                flint_printf("a: "); fmpq_mpoly_print_pretty(a, NULL, ctx); flint_printf("\n");
                flint_printf("b: "); fmpq_mpoly_print_pretty(b, NULL, ctx); flint_printf("\n");

                flint_printf("disc(a*b): "); fmpq_mpoly_print_pretty(dab, NULL, ctx); flint_printf("\n");
                flint_printf("disc(a): "); fmpq_mpoly_print_pretty(da, NULL, ctx); flint_printf("\n");
                flint_printf("disc(b): "); fmpq_mpoly_print_pretty(db, NULL, ctx); flint_printf("\n");
                flint_printf("res(a, b): "); fmpq_mpoly_print_pretty(r, NULL, ctx); flint_printf("\n");

                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(ab, ctx);
        fmpq_mpoly_clear(da, ctx);
        fmpq_mpoly_clear(db, ctx);
        fmpq_mpoly_clear(dab, ctx);
        fmpq_mpoly_clear(r, ctx);
        fmpq_mpoly_clear(p, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
