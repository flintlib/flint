/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, 100003);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(d, ctx);
        nmod_mpoly_init(d1, ctx);

        nmod_mpoly_set_str_pretty(f, "a^10*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        nmod_mpoly_set_str_pretty(d1, "b^200 - 4*a^10*c^100000000000000000000", vars, ctx);
        if (!nmod_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!nmod_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(d, ctx);
        nmod_mpoly_clear(d1, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, r;
        nmod_poly_t au, bu;
        mp_limb_t ru;
        slong len1, len2, exp_bound1, exp_bound2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, 1, mpoly_ordering_randtest(state), modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_poly_init(au, nmod_mpoly_ctx_modulus(ctx));
        nmod_poly_init(bu, nmod_mpoly_ctx_modulus(ctx));

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        nmod_mpoly_get_nmod_poly(au, a, 0, ctx);
        nmod_mpoly_get_nmod_poly(bu, b, 0, ctx);

        ru = nmod_poly_resultant(au, bu);

        if (!nmod_mpoly_resultant(r, a, b, 0, ctx) ||
            !nmod_mpoly_equal_ui(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_poly_clear(au);
        nmod_poly_clear(bu);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, c, ab, ra, rb, rab, p;
        slong len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 3, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(c, ctx);
        nmod_mpoly_init(ab, ctx);
        nmod_mpoly_init(ra, ctx);
        nmod_mpoly_init(rb, ctx);
        nmod_mpoly_init(rab, ctx);
        nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        len3 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);
        nmod_mpoly_randtest_bound(c, state, len3, exp_bound3, ctx);

        for (j = 0; j < nmod_mpoly_ctx_nvars(ctx); j++)
        {
            nmod_mpoly_mul(ab, a, b, ctx);

            if (!nmod_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(ra, ctx);

            if (!nmod_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(rb, ctx);

            if (!nmod_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(rab, ctx);

            nmod_mpoly_mul(p, ra, rb, ctx);

            if (!nmod_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(c, ctx);
        nmod_mpoly_clear(ab, ctx);
        nmod_mpoly_clear(ra, ctx);
        nmod_mpoly_clear(rb, ctx);
        nmod_mpoly_clear(rab, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, ab, r, da, db, dab, p;
        slong len1, len2, exp_bound1, exp_bound2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 3, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ab, ctx);
        nmod_mpoly_init(da, ctx);
        nmod_mpoly_init(db, ctx);
        nmod_mpoly_init(dab, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 4) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        for (j = 0; j < nmod_mpoly_ctx_nvars(ctx); j++)
        {
            if (nmod_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (nmod_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            nmod_mpoly_mul(ab, a, b, ctx);

            if (!nmod_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            nmod_mpoly_mul(p, da, db, ctx);
            nmod_mpoly_mul(p, p, r, ctx);
            nmod_mpoly_mul(p, p, r, ctx);

            if (!nmod_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");

                flint_printf("a: "); nmod_mpoly_print_pretty(a, NULL, ctx); flint_printf("\n");
                flint_printf("b: "); nmod_mpoly_print_pretty(b, NULL, ctx); flint_printf("\n");

                flint_printf("disc(a*b): "); nmod_mpoly_print_pretty(dab, NULL, ctx); flint_printf("\n");
                flint_printf("disc(a): "); nmod_mpoly_print_pretty(da, NULL, ctx); flint_printf("\n");
                flint_printf("disc(b): "); nmod_mpoly_print_pretty(db, NULL, ctx); flint_printf("\n");
                flint_printf("res(a, b): "); nmod_mpoly_print_pretty(r, NULL, ctx); flint_printf("\n");

                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ab, ctx);
        nmod_mpoly_clear(da, ctx);
        nmod_mpoly_clear(db, ctx);
        nmod_mpoly_clear(dab, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
