/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_gcd_zippel, state)
{
    slong i, j;

    /* examples from Zippel's 1979 paper */
    if (1) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t r, d, f, g;
        int success;
        const char* vars[] =
            {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};

        const char * example[][3] =
        {{
            "x1^2 + x1 + 3",
            "2*x1^2 + 2*x1 + 1",
            "x1^2 + 2*x1 + 2"
        }, {
            "2*x1^2*x2^2 + x1*x2 + 2*x1",
            "x2^2 + 2*x1^2*x2 + x1^2 + 1",
            "x1^2*x2^2 + x1^2*x2 + x1*x2 + x1^2 + x1"
        }, {
            "x2^2*x3^2 + x2^2*x3 + 2*x1^2*x2*x3 + x1*x3",
            "x3^2 + x2^2*x3 + x1^2*x2*x3 + x1*x3 + x1^2*x2^2",
            "x2*x3 + 2*x1*x3 + x3 + x1"
        }, {
            "x1^2*x4^2 + x2^2*x3*x4 + x1^2*x2*x4 + x2*x4 + x1^2*x2*x3",
            "x1*x2*x3^2*x4^2 + x1*x3^2*x4^2 + x1*x4^2 + x4^2 + x1*x3*x4",
            "x1*x3^2*x4^2 + x3^2*x4^2 + x4^2 + x1*x2^2*x3*x4 + x1*x2^2"
        }, {
            "x1^3*x2^2*x3^2*x4*x5^2 + x1*x2^2*x5^2 + x1^3*x3*x4^2*x5"
                                    " + x1^3*x2*x3^2*x4*x5 + x1^2*x2*x3^2*x4^2"
            ,
            "x1*x2^2*x5^2 + x1*x2*x3^2*x4*x5 + x1*x2*x3^2*x4^2"
                                                          " + x1*x2^2*x4^2 + 1"
            ,
            "x1*x3^2*x4*x5^2 + x2*x5^2 + x1*x2*x4*x5 + x2*x5 + x1*x2*x3*x4^2"
        }, {
            "x1*x2*x4^2*x5^2*x6^2 + x1*x2^2*x3^2*x4*x5^2*x6^2  + x1^2*x3*x6^2"
                                   " + x1^2*x2*x3^2*x4*x5^2*x6 + x1^2*x3*x5*x6"
            ,
            "x1^2*x2*x4*x5^2*x6^2 + x1*x3*x5^2*x6^2 + x1*x2^2*x6^2"
                                      " + x1^2*x2^2*x3^2*x5*x6 + x1*x3^2*x4*x5"
            ,
            "x2^2*x3^2*x4*x5^2*x6 + x1*x4^2*x5*x6 + x2^2*x3^2*x4*x5*x6"
                                         " + x1*x2^2*x3*x4^2*x6 + x1^2*x3*x5^2"
        }, {
            "x1*x2^2*x4^2*x6^2*x7^2 + x1^2*x3*x4*x6^2*x7^2 + x3^2*x4^2*x7^2"
                                              " + x1^2*x2*x4^2*x6 + x3*x4*x5^2"
            ,
            "x1^2*x2*x4^2*x5*x6^2*x7^2 + x1*x2*x3*x6*x7 + x3*x4^2*x5^2*x7"
                                   " + x1*x4^2*x5^2*x7 + x1^2*x2*x3*x4^2+x5*x6"
            ,
            "x1*x3*x5*x6^2*x7^2 + x2^2*x3^2*x4^2*x5*x6*x7^2 + x4*x6*x7^2"
                                   " + x1^2*x2*x3*x5*x6*x7 + x1^2*x3^2*x4*x5^2"
        }, {
            "x2^2*x4*x5*x6*x7*x8^2 + x1^2*x2*x3^2*x4^2*x6^2*x7^2*x8"
                   " + x1^2*x3*x4^2*x6^2*x7^2 + x1^2*x2^2*x3^2*x4*x5^2*x6*x7^2"
                                                                " + x2^2*x4*x6"
            ,
            "x1^2*x2^2*x3*x4^2*x5*x6^2*x8^2 + x2*x5*x6^2*x8^2"
              " + x1^2*x2^2*x3^2*x4^2*x6^2*x7^2*x8 + x1^2*x3^2*x4*x5^2*x7^2*x8"
                                                      " + x1*x2^2*x3^2*x5^2*x7"
            ,
            "x1*x4^2*x5*x6*x7*x8^2 + x1*x2^2*x4^2*x5^2*x6^2*x8"
                       " + x1^2*x2*x3*x4^2*x6^2*x8 + x1^2*x2^2*x3^2*x4*x5^2*x8"
                                                           " + x1*x2*x4^2*x5^2"
        }, {
            "x1^2*x3^3*x4*x6*x8*x9^2 + x1*x2*x3*x4^2*x5^2*x8*x9"
                      " + x2*x3*x4*x5^2*x8*x9 + x1*x3^3*x4^2*x5^2*x6^2*x7*x8^2"
                                                  " + x2*x3*x4*x5^2*x6*x7*x8^2"
            ,
            "x1^2*x2^2*x3*x7^2*x8*x9 + x2^2*x9 + x1^2*x3*x4^2*x5^2*x6*x7^2"
                                            " + x4^2*x5^2*x7^2 + x3*x4^2*x6*x7"
            ,
            "x1^2*x2*x4*x5*x6*x7^2*x8^2*x9^2 + x1^2*x2*x3*x5*x6^2*x7^2*x8*x9^2"
                              " + x1^2*x3*x4*x6*x7^2*x8*x9 + x1^2*x2^2*x6*x8^2"
                                                        " + x2^2*x4*x5*x6^2*x7"
        }, {
            "x1*x2^2*x4^2*x8*x9^2*x10^2 + x2^2*x4*x5^2*x6*x7*x9*x10^2"
                        " + x1^2*x2*x3*x5^2*x7^2*x9^2 + x1*x3^2*x4^2*x7^2*x9^2"
                                                      " + x1^2*x3*x4*x7^2*x8^2"
            ,
            "x1*x2*x3^2*x4*x6*x7*x8*x9^2*x10^2 + x2^2*x3^2*x4^2*x6^2*x9*x10^2"
                                    " + x1*x2^2*x3^2*x4*x5*x6*x7*x8^2*x9^2*x10"
               " + x1^2*x2*x4^2*x5^2*x8^2*x9^2*x10 + x3*x4^2*x5*x6*x7^2*x9*x10"
            ,
            "x1*x2^2*x3^2*x5^2*x6^2*x7*x8*x9^2*x10^2 + x3*x8*x9^2*x10^2"
                  " + x1*x2^2*x3*x4*x5^2*x6^2*x8^2*x9*x10 + x1*x3*x6*x7*x8*x10"
                                                    " + x4^2*x5^2*x6^2*x7*x9^2"
        }};

        for (i = 1; i <= 10; i++)
        {
            fmpz_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX);
            fmpz_mpoly_init(r, ctx);
            fmpz_mpoly_init(d, ctx);
            fmpz_mpoly_init(f, ctx);
            fmpz_mpoly_init(g, ctx);
            fmpz_mpoly_set_str_pretty(d, example[i - 1][0], vars, ctx);
            fmpz_mpoly_set_str_pretty(f, example[i - 1][1], vars, ctx);
            fmpz_mpoly_set_str_pretty(g, example[i - 1][2], vars, ctx);
            fmpz_mpoly_mul_johnson(f, f, d, ctx);
            fmpz_mpoly_mul_johnson(g, g, d, ctx);
            success = fmpz_mpoly_gcd_zippel(r, f, g, ctx);
            if (!success || !fmpz_mpoly_equal(r, d, ctx))
            {
                flint_printf("FAIL\ncheck example %wd\n",i);
                fflush(stdout);
                flint_abort();
            }
            fmpz_mpoly_clear(r, ctx);
            fmpz_mpoly_clear(d, ctx);
            fmpz_mpoly_clear(f, ctx);
            fmpz_mpoly_clear(g, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }

    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, ca, cb, cg, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds;
        int res;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ca, ctx);
        fmpz_mpoly_init(cb, ctx);
        fmpz_mpoly_init(cg, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);

        degbound = 100/(2*ctx->minfo->nvars - 1);
        degbounds = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        for (j = 0; j < ctx->minfo->nvars; j++)
            degbounds[j] = n_randint(state, degbound + UWORD(1)) + UWORD(1);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bounds(t, state, len, coeff_bits + 1, degbounds, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bounds(a, state, len1, coeff_bits, degbounds, ctx);
            fmpz_mpoly_randtest_bounds(b, state, len2, coeff_bits, degbounds, ctx);
            fmpz_mpoly_mul_johnson(a, a, t, ctx);
            fmpz_mpoly_mul_johnson(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpz_mpoly_gcd_zippel(g, a, b, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check that gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }
                continue;
            }

            if (fmpz_sgn(g->coeffs + 0) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            res = 1;
            res = res && fmpz_mpoly_divides_monagan_pearce(ca, a, g, ctx);
            res = res && fmpz_mpoly_divides_monagan_pearce(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            res = fmpz_mpoly_gcd_zippel(cg, ca, cb, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check that cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            if (!fmpz_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(degbounds);

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ca, ctx);
        fmpz_mpoly_clear(cb, ctx);
        fmpz_mpoly_clear(cg, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
