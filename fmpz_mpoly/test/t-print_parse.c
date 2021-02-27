/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void check_match(
    const char * sa,
    const char * sb,
    const char ** vars,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t a, b;
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_init(b, ctx);

    if (0 != fmpz_mpoly_set_str_pretty(a, sa, vars, ctx))
    {
        flint_printf("FAIL: set_str failed\n");
        flint_printf("sa = %s\n", sa);
        flint_abort();
    }

    if (0 != fmpz_mpoly_set_str_pretty(b, sb, vars, ctx))
    {
        flint_printf("FAIL: set_str failed\n");
        flint_printf("sb = %s\n", sb);
        flint_abort();
    }

    if (!fmpz_mpoly_equal(a, b, ctx))
    {
        flint_printf("FAIL: results do not match\n");
        flint_printf("sa = %s, sb = %s\n", sa, sb);
        flint_printf("a: ");
        fmpz_mpoly_print_pretty(a, vars, ctx);
        flint_printf("\n");
        flint_printf("b: ");
        fmpz_mpoly_print_pretty(b, vars, ctx);
        flint_printf("\n");
        flint_abort();
    }

    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_clear(b, ctx);
}

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);
    flint_printf("print_parse....");
    fflush(stdout);

    {
        slong len1, exp_bits, coeff_bits;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, f1;
        char * str;
        const char * vars[] = {"x","xy","y","yx","z","zz"};

        {
            fmpz_mpoly_ctx_init(ctx, 6, ORD_LEX);
            check_match("-(1+x)^2", "-1-2*x-x^2", vars, ctx);
            check_match("x+--y", "x+y", vars, ctx);
            check_match("-x+-x*-y", "-x+x*y", vars, ctx);
            check_match("-x+-x^2*-y^3/-x*-y-y", "-x+x*y^4-y", vars, ctx);
            check_match("1+-x+x^2+-x^3+x^4+-x^5+x^6+-x^7", "1-x+x^2-x^3+x^4-x^5+x^6-x^7", vars, ctx);

            fmpz_mpoly_ctx_clear(ctx);
        }

        /* check that parsing inverts printing */
        for (i = 0; i < flint_test_multiplier(); i++)
        {
            fmpz_mpoly_ctx_init_rand(ctx, state, 6);
            fmpz_mpoly_init(f, ctx);
            fmpz_mpoly_init(f1, ctx);

            for (len1 = 3; len1 < 4000; len1 += len1/2)
            {
                coeff_bits = 100;
                exp_bits = 100;
                fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);
                str = fmpz_mpoly_get_str_pretty(f, vars, ctx);
                if (0 != fmpz_mpoly_set_str_pretty(f1, str, vars, ctx))
                {
                    flint_printf("FAIL: set_str failed\n");
                    flint_printf("i = %wd, len1 = %wd\n", i ,len1);
                    flint_abort();
                }

                fmpz_mpoly_assert_canonical(f1, ctx);
                flint_free(str);

                if (!fmpz_mpoly_equal(f, f1, ctx))
                {
                    flint_printf("FAIL: check that parsing inverts printing\n");
                    flint_printf("i = %wd, len1 = %wd\n", i ,len1);
                    flint_abort();
                }
            }

            fmpz_mpoly_clear(f, ctx);
            fmpz_mpoly_clear(f1, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }
    }

    flint_printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

