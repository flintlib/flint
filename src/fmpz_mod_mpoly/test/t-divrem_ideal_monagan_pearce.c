/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_divrem_ideal_monagan_pearce, state)
{
    slong i, j, w;

    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_struct q[4], r[1], g[1], v[4], s[1], t[1];
        fmpz_mod_mpoly_struct* Q[4], * V[4];
        fmpz_t m;

        fmpz_init_set_ui(m, 23);
        fmpz_mod_mpoly_ctx_init(ctx, 8, ORD_LEX, m);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_set_str_pretty(g, "12*x1^83*x2^72*x3^58*x4^35*x5^92*x6^56*x7^87*x8^69 + 22*x1^51*x2^51*x3^42*x4^12*x5^33*x6^61*x7^97*x8^14 + x1^28*x2^60*x3^77*x4^29*x5^4*x6^63*x7^51*x8^92", NULL, ctx);

        fmpz_mod_mpoly_init(v+0, ctx); V[0] = v+0;
        fmpz_mod_mpoly_set_str_pretty(V[0], "12*x1^78*x2^71*x3^53*x4*x5^63*x6^49*x7^49*x8^14 + 17*x1^21*x2^68*x3^59*x4^41*x5^40*x6^43*x7^89*x8^87 + 7*x1^4*x2^22*x3^20*x4^87*x5^98*x6^84*x7^36*x8^4", NULL, ctx);
        fmpz_mod_mpoly_init(v+1, ctx); V[1] = v+1;
        fmpz_mod_mpoly_set_str_pretty(V[1], "4*x1^78*x2^96*x3^85*x4^29*x5^97*x6^79*x7^65*x8^58 + 13*x1^77*x2^86*x3^99*x4^39*x5^41*x6^83*x7^21*x8^31", NULL, ctx);
        fmpz_mod_mpoly_init(v+2, ctx); V[2] = v+2;
        fmpz_mod_mpoly_set_str_pretty(V[2], "6*x1^72*x2^43*x3^56*x4^76*x5^81*x6^94*x7^67*x8^44 + 4*x1^43*x2^11*x3^31*x4^89*x5^55*x6^55*x7^10*x8^78 + 11*x1^38*x2^23*x3^17*x4^87*x5^89*x6^10*x7^87*x8^96", NULL, ctx);
        fmpz_mod_mpoly_init(v+3, ctx); V[3] = v+3;
        fmpz_mod_mpoly_set_str_pretty(V[3], "7*x1^100*x2^94*x3^48*x4^46*x5^81*x6^15*x7^58*x8^83 + 22*x1^84*x2^64*x3^32*x4^43*x5^72*x6^47*x7^62*x8^83 + 5*x1^42*x2^10*x3^22*x4^69*x5^62*x6^36*x7^11*x8^72 + 13*x1^22*x2^62*x3^62*x4^29*x5^97*x6^3*x7^96*x8^40", NULL, ctx);

        fmpz_mod_mpoly_init(s, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        fmpz_mod_mpoly_init(q+0, ctx); Q[0] = q+0;
        fmpz_mod_mpoly_init(q+1, ctx); Q[1] = q+1;
        fmpz_mod_mpoly_init(q+2, ctx); Q[2] = q+2;
        fmpz_mod_mpoly_init(q+3, ctx); Q[3] = q+3;

        fmpz_mod_mpoly_divrem_ideal_monagan_pearce(Q, r, g, V, 4, ctx);

        fmpz_mod_mpoly_mul(t, Q[0], V[0], ctx);
        fmpz_mod_mpoly_add(s, r, t, ctx);
        fmpz_mod_mpoly_mul(t, Q[1], V[1], ctx);
        fmpz_mod_mpoly_add(s, s, t, ctx);
        fmpz_mod_mpoly_mul(t, Q[2], V[2], ctx);
        fmpz_mod_mpoly_add(s, s, t, ctx);
        fmpz_mod_mpoly_mul(t, Q[3], V[3], ctx);
        fmpz_mod_mpoly_add(s, s, t, ctx);

        if (!fmpz_mod_mpoly_equal(s, g, ctx))
        {
            flint_printf("FAIL: Check exponent overflow\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(g, ctx);

        fmpz_mod_mpoly_clear(v+0, ctx);
        fmpz_mod_mpoly_clear(v+1, ctx);
        fmpz_mod_mpoly_clear(v+2, ctx);
        fmpz_mod_mpoly_clear(v+3, ctx);

        fmpz_mod_mpoly_clear(s, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(r, ctx);

        fmpz_mod_mpoly_clear(q+0, ctx);
        fmpz_mod_mpoly_clear(q+1, ctx);
        fmpz_mod_mpoly_clear(q+2, ctx);
        fmpz_mod_mpoly_clear(q+3, ctx);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;
        fmpz_mod_mpoly_struct * qarr[1], * darr[1];

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 100) + 1;
        exp_bits1 = n_randint(state, 100) + 1;
        exp_bits2 = n_randint(state, 100) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2 + 1, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(r, state, len, exp_bits, ctx);

            fmpz_mod_mpoly_mul(h, f, g, ctx);

            qarr[0] = k;
            darr[0] = g;

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);
            fmpz_mod_mpoly_assert_canonical(qarr[0], ctx);
            fmpz_mod_mpoly_assert_canonical(r, ctx);

            if (!fmpz_mod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_clear(r, ctx);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, r, k1, k2;
        fmpz_mod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fmpz_mod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;

        g = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);
        q = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 200);

        for (w = 0; w < num; w++)
        {
            fmpz_mod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fmpz_mod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(k1, ctx);
        fmpz_mod_mpoly_init(k2, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fmpz_mod_mpoly_is_zero(darr[w], ctx))
                    fmpz_mod_mpoly_one(darr[w], ctx);
                fmpz_mod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fmpz_mod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);
            fmpz_mod_mpoly_assert_canonical(r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mod_mpoly_remainder_strongtest(r, darr[w], ctx);
            }

            fmpz_mod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mod_mpoly_add(k2, k2, k1, ctx);
            }
            fmpz_mod_mpoly_add(k2, k2, r, ctx);

            if (!fmpz_mod_mpoly_equal(f, k2, ctx))
            {
                flint_printf("FAIL: Check f = g1*q1 + ... + gn*qn + r\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(darr[w], ctx);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(k1, ctx);
        fmpz_mod_mpoly_clear(k2, ctx);
        fmpz_mod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, r, k1, k2;
        fmpz_mod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fmpz_mod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;

        g = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);
        q = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 200);

        for (w = 0; w < num; w++)
        {
            fmpz_mod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fmpz_mod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(k1, ctx);
        fmpz_mod_mpoly_init(k2, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fmpz_mod_mpoly_is_zero(darr[w], ctx))
                    fmpz_mod_mpoly_one(darr[w], ctx);
                fmpz_mod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fmpz_mod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_set(r, f, ctx);

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mod_mpoly_remainder_strongtest(f, darr[w], ctx);
            }

            fmpz_mod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mod_mpoly_add(k2, k2, k1, ctx);
            }
            fmpz_mod_mpoly_add(k2, k2, f, ctx);

            if (!fmpz_mod_mpoly_equal(r, k2, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(darr[w], ctx);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(k1, ctx);
        fmpz_mod_mpoly_clear(k2, ctx);
        fmpz_mod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
