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

TEST_FUNCTION_START(padic_poly_inv_series, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;
        slong n;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(a, state, n_randint(state, 100) + 1, ctx);
        if (fmpz_is_zero(a->coeffs))
        {
            fmpz_randtest_not_zero(a->coeffs, state, 20);
            fmpz_remove(a->coeffs, a->coeffs, p);
            padic_poly_reduce(a, ctx);
        } else
            fmpz_remove(a->coeffs, a->coeffs, p);

        padic_poly_set(b, a, ctx);
        n = n_randint(state, 100) + 1;

        padic_poly_inv_series(c, b, n, ctx);
        padic_poly_inv_series(b, b, n, ctx);

        result = (padic_poly_equal(b, c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), padic_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), padic_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), padic_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /*
        Check correctness:

        If ord_p(a) = v then we can compute b = a^{-1} mod p^N
        and we will have a b = 1 mod p^{N-|v|}.  Thus, require
        that N - |v| > 0.
     */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;
        slong n, N2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - 1) + 1;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);

        {
            slong i, len = n_randint(state, 10) + 1;
            int alloc;
            fmpz_t pow;

            padic_poly_fit_length(a, len);
            _padic_poly_set_length(a, len);
            a->val = n_randint(state, N);
            if (n_randint(state, 2))
                a->val = - a->val;

            alloc = _padic_ctx_pow_ui(pow, N - a->val, ctx);

            for (i = 0; i < len; i++)
                fmpz_randm(a->coeffs + i, state, pow);
            while (fmpz_is_zero(a->coeffs))
                fmpz_randm(a->coeffs, state, pow);
            fmpz_remove(a->coeffs, a->coeffs, p);
            _padic_poly_normalise(a);

            if (alloc)
                fmpz_clear(pow);
        }

        n = n_randint(state, 100) + 1;

        N2 = N - FLINT_ABS(a->val);
        padic_poly_init2(c, 0, N2);

        padic_poly_inv_series(b, a, n, ctx);
        padic_poly_mul(c, a, b, ctx);
        padic_poly_truncate(c, n, p);

        result = (padic_poly_is_one(c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), padic_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), padic_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), padic_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("N = %wd\n", N);
            flint_printf("N2 = %wd\n", N2);
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
