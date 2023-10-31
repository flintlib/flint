/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_compose_pow, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    /* Aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;
        slong k;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a, ctx);
        k = n_randint(state, 20) + 1;

        padic_poly_compose_pow(c, b, k, ctx);
        padic_poly_compose_pow(b, b, k, ctx);

        result = (padic_poly_equal(b, c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            padic_poly_print(a, ctx), flint_printf("\n\n");
            padic_poly_print(b, ctx), flint_printf("\n\n");
            padic_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Compare with usual composition */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t f, g, h1, h2;
        slong k;
        padic_t one;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(f,  0, N);
        padic_poly_init2(g,  0, WORD_MAX);  /* TODO:  Check this is OK */
        padic_poly_init2(h1, 0, N);
        padic_poly_init2(h2, 0, N);

        padic_poly_randtest(f, state, n_randint(state, 40), ctx);
        k = n_randint(state, 20) + 1;

        padic_poly_compose_pow(h1, f, k, ctx);

        padic_init(one);
        padic_one(one);
        padic_poly_set_coeff_padic(g, k, one, ctx);
        padic_clear(one);
        padic_poly_compose(h2, f, g, ctx);

        result = (padic_poly_equal(h1, h2) && padic_poly_is_reduced(h1, ctx));
        if (!result)
        {
            flint_printf("FAIL (cmp with composition):\n");
            flint_printf("f  = "), padic_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("g  = "), padic_poly_print(g, ctx), flint_printf("\n\n");
            flint_printf("h1 = "), padic_poly_print(h1, ctx), flint_printf("\n\n");
            flint_printf("h2 = "), padic_poly_print(h2, ctx), flint_printf("\n\n");
            flint_printf("p  = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("N  = %wd\n\n", N);
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(f);
        padic_poly_clear(g);
        padic_poly_clear(h1);
        padic_poly_clear(h2);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
