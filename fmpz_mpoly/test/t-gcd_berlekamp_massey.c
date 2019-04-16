/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
void usleep(ulong);

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_berlekamp_massey....");
    fflush(stdout);

    if (1)
    {
        int res;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t G, A, B;
        const char* vars[] = {"y", "t", "x", "z"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(G, ctx);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(B, ctx);

        fmpz_mpoly_set_str_pretty(A, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(B, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_set_str_pretty(G, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11"
      " - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20"
         " + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);

        fmpz_mpoly_mul(A, A, G, ctx);
        fmpz_mpoly_mul(B, B, G, ctx);
        res = fmpz_mpoly_gcd_berlekamp_massey(G, A, B, ctx);
        if (!res)
        {
            printf("FAIL\n");
            flint_printf("Check example could be computed\n");
            flint_abort();
        }

        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(G, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, ca, cb, cg, t;
        mp_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;
        int res;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ca, ctx);
        fmpz_mpoly_init(cb, ctx);
        fmpz_mpoly_init(cg, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 100/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpz_mpoly_gcd_berlekamp_massey(g, a, b, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            fmpz_mpoly_assert_canonical(g, ctx);

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (fmpz_sgn(g->coeffs + 0) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fmpz_mpoly_divides(ca, a, g, ctx);
            res = res && fmpz_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = fmpz_mpoly_gcd_berlekamp_massey(cg, ca, cb, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            if (!fmpz_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ca, ctx);
        fmpz_mpoly_clear(cb, ctx);
        fmpz_mpoly_clear(cg, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
