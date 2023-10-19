/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_divides, state)
{
    int i, result;

    /* Random polynomials */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, q, prod;
        fmpz_mod_ctx_t ctx;
        int divides;
        fmpz_t n;

        fmpz_init(n);

        do fmpz_randtest_unsigned(n, state, 2 * FLINT_BITS);
        while (!fmpz_is_probabprime(n));
        fmpz_mod_ctx_init(ctx, n);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(prod, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 300), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 300), ctx);

        divides = fmpz_mod_poly_divides(q, a, b, ctx);
        fmpz_mod_poly_mul(prod, q, b, ctx);

        result = ((divides && fmpz_mod_poly_equal(a, prod, ctx)) ||
                (!divides && fmpz_mod_poly_is_zero(q, ctx) &&
                 !fmpz_mod_poly_equal(prod, a, ctx)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides = %d\n", divides);
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(prod, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_ctx_clear(ctx);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(prod, ctx);
    }

    /* Random divisible polynomials */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, q, prod;
        fmpz_mod_ctx_t ctx;
        int divides;
        fmpz_t n;

        fmpz_init(n);

        do fmpz_randtest_unsigned(n, state, 2 * FLINT_BITS);
        while (!fmpz_is_probabprime(n));
        fmpz_mod_ctx_init(ctx, n);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(prod, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 300), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 300), ctx);
        fmpz_mod_poly_mul(a, b, a, ctx);

        divides = fmpz_mod_poly_divides(q, a, b, ctx);
        fmpz_mod_poly_mul(prod, q, b, ctx);

        result = (divides && fmpz_mod_poly_equal(a, prod, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides = %d\n", divides);
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(prod, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_ctx_clear(ctx);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(prod, ctx);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, q;
        fmpz_mod_ctx_t ctx;
        int divides1, divides2;
        fmpz_t n;

        fmpz_init(n);

        do fmpz_randtest_unsigned(n, state, 2 * FLINT_BITS);
        while (!fmpz_is_probabprime(n));
        fmpz_mod_ctx_init(ctx, n);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        divides1 = fmpz_mod_poly_divides(q, a, b, ctx);
        divides2 = fmpz_mod_poly_divides(a, a, b, ctx);

        result = (divides1 == divides2 && fmpz_mod_poly_equal(q, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_ctx_clear(ctx);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, q;
        fmpz_mod_ctx_t ctx;
        int divides1, divides2;
        fmpz_t n;

        fmpz_init(n);

        do fmpz_randtest_unsigned(n, state, 2 * FLINT_BITS);
        while (!fmpz_is_probabprime(n));
        fmpz_mod_ctx_init(ctx, n);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        divides1 = fmpz_mod_poly_divides(q, a, b, ctx);
        divides2 = fmpz_mod_poly_divides(b, a, b, ctx);

        result = (divides1 == divides2 && fmpz_mod_poly_equal(q, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_ctx_clear(ctx);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
    }

    TEST_FUNCTION_END(state);
}
