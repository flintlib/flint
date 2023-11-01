/*
    Copyright (C) 2012 Fredrik Johansson
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

TEST_FUNCTION_START(fmpz_mod_poly_sqrt, state)
{
    int i;

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b;
        int square1, square2;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, FLINT_BITS);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);

        fmpz_mod_poly_randtest(a, state, 1 + n_randint(state, 50), ctx);

        if (n_randint(state, 2))
            fmpz_mod_poly_mul(a, a, a, ctx);

        square1 = fmpz_mod_poly_sqrt(b, a, ctx);
        square2 = fmpz_mod_poly_sqrt(a, a, ctx);

        if (square1 != square2 || (square1 && !fmpz_mod_poly_equal(a, b, ctx)))
        {
            flint_printf("FAIL: aliasing:\n");
            flint_printf("square1 = %d, square2 = %d\n\n", square1, square2);
            flint_printf("a: "); fmpz_mod_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); fmpz_mod_poly_print(b, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        int square;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, FLINT_BITS);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        fmpz_mod_poly_randtest(a, state, 1 + n_randint(state, 50), ctx);
        fmpz_mod_poly_mul(b, a, a, ctx);
        square = fmpz_mod_poly_sqrt(c, b, ctx);

        if (!square)
        {
            flint_printf("FAIL: square reported nonsquare:\n");
            flint_printf("a: "); fmpz_mod_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); fmpz_mod_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("c: "); fmpz_mod_poly_print(c, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_mul(c, c, c, ctx);
        if (!fmpz_mod_poly_equal(c, b, ctx))
        {
            flint_printf("FAIL: sqrt(b)^2 != b:\n");
            flint_printf("a: "); fmpz_mod_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); fmpz_mod_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("c: "); fmpz_mod_poly_print(c, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        slong j;
        int square;
        fmpz_mod_ctx_t ctx;
        fmpz_t t;

        fmpz_init(t);

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, FLINT_BITS);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        fmpz_mod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 50), ctx);
        fmpz_mod_poly_mul(b, a, a, ctx);

        j = n_randint(state, fmpz_mod_poly_length(b, ctx));
        fmpz_randm(t, state, ctx->n);

        fmpz_set(b->coeffs + j, t);
        _fmpz_mod_poly_normalise(b);

        square = fmpz_mod_poly_sqrt(c, b, ctx);

        if (square)
        {
            fmpz_mod_poly_mul(c, c, c, ctx);
            if (!fmpz_mod_poly_equal(c, b, ctx))
            {
                flint_printf("FAIL: sqrt(b)^2 != b:\n");
                flint_printf("a: "); fmpz_mod_poly_print(a, ctx); flint_printf("\n\n");
                flint_printf("b: "); fmpz_mod_poly_print(b, ctx); flint_printf("\n\n");
                flint_printf("c: "); fmpz_mod_poly_print(c, ctx); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);

        fmpz_clear(t);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
