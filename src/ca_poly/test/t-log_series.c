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

TEST_FUNCTION_START(ca_poly_log_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, C, D;
        slong n1;

        /* Test exp(log(A)) == A */
        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(C, ctx);
        ca_poly_init(D, ctx);

        ca_poly_randtest_rational(A, state, 10, 10, ctx);

        n1 = n_randint(state, 10);

        if (n_randint(state, 2))
        {
            ca_poly_log_series(B, A, n1, ctx);
        }
        else
        {
            ca_poly_set(B, A, ctx);
            ca_poly_log_series(B, B, n1, ctx);
        }

        ca_poly_exp_series(C, B, n1, ctx);
        ca_poly_set(D, A, ctx);
        ca_poly_truncate(D, n1, ctx);

        if (!(A->length == 0 || ca_check_is_zero(A->coeffs, ctx) != T_FALSE))
        {
            if (ca_poly_check_equal(C, D, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_poly_print(C, ctx); flint_printf("\n");
                flint_printf("D = "); ca_poly_print(D, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(C, ctx);
        ca_poly_clear(D, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
