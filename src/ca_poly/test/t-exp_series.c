/*
    Copyright (C) 2020, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_poly.h"

/* Defined in t-div_series.c, t-exp_series.c, t-inv_series.c, t-mul.c */
#ifndef ca_poly_randtest_same_nf
#define ca_poly_randtest_same_nf ca_poly_randtest_same_nf
void
ca_poly_randtest_same_nf(ca_poly_t res, flint_rand_t state, const ca_t x, slong len, slong bits, slong den_bits, ca_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    ca_poly_fit_length(res, len, ctx);
    for (i = 0; i < len; i++)
        ca_randtest_same_nf(res->coeffs + i, state, x, bits, 1, ctx);
    _ca_poly_set_length(res, len, ctx);
    _ca_poly_normalise(res, ctx);

    fmpz_init(t);
    fmpz_randtest_not_zero(t, state, den_bits);
    ca_poly_div_fmpz(res, res, t, ctx);
    fmpz_clear(t);
}
#endif

/* Defined in t-div_series.c, t-exp_series.c, t-inv_series.c, t-log_series.c */
#ifndef ca_poly_truncate
#define ca_poly_truncate ca_poly_truncate
/* todo */
void
ca_poly_truncate(ca_poly_t poly, slong newlen, ca_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            ca_zero(poly->coeffs + i, ctx);
        poly->length = newlen;
        _ca_poly_normalise(poly, ctx);
    }
}
#endif

TEST_FUNCTION_START(ca_poly_exp_series, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, AB, expA, expB, expAB, expAexpB;
        slong n1, n2, n3;

        /* Test exp(A)*exp(B) = exp(A+B) */
        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(AB, ctx);
        ca_poly_init(expA, ctx);
        ca_poly_init(expB, ctx);
        ca_poly_init(expAB, ctx);
        ca_poly_init(expAexpB, ctx);

        if (n_randint(state, 2))
        {
            ca_poly_randtest(A, state, 4, 2, 10, ctx);
            ca_poly_randtest(B, state, 4, 2, 10, ctx);

            n1 = n_randint(state, 6);
            n2 = n_randint(state, 6);
            n3 = n_randint(state, 6);
        }
        else
        {
            qqbar_t t;
            ca_t x;

            qqbar_init(t);
            ca_init(x, ctx);

            qqbar_randtest(t, state, 8, 10);
            ca_set_qqbar(x, t, ctx);

            ca_poly_randtest_same_nf(A, state, x, 1 + n_randint(state, 20), 2 + n_randint(state, 100), 2 + n_randint(state, 100), ctx);
            ca_poly_randtest_same_nf(B, state, x, 1 + n_randint(state, 20), 2 + n_randint(state, 100), 2 + n_randint(state, 100), ctx);

            n1 = n_randint(state, 50);
            n2 = n_randint(state, 50);
            n3 = n_randint(state, 50);

            qqbar_clear(t);
            ca_clear(x, ctx);
        }

        ca_poly_randtest(expA, state, 4, 2, 10, ctx);
        ca_poly_randtest(expB, state, 4, 2, 10, ctx);
        ca_poly_randtest(expAB, state, 4, 2, 10, ctx);

        ca_poly_add(AB, A, B, ctx);

        ca_poly_exp_series(expA, A, n1, ctx);
        ca_poly_exp_series(expB, B, n2, ctx);
        ca_poly_exp_series(expAB, AB, n3, ctx);

        ca_poly_mullow(expAexpB, expA, expB, FLINT_MIN(FLINT_MIN(n1, n2), n3), ctx);
        ca_poly_truncate(expAB, FLINT_MIN(FLINT_MIN(n1, n2), n3), ctx);

        if (ca_poly_check_equal(expAB, expAexpB, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
            flint_printf("expAexpB = "); ca_poly_print(expAexpB, ctx); flint_printf("\n");
            flint_printf("expAB = "); ca_poly_print(expAB, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(AB, ctx);
        ca_poly_clear(expA, ctx);
        ca_poly_clear(expB, ctx);
        ca_poly_clear(expAB, ctx);
        ca_poly_clear(expAexpB, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
