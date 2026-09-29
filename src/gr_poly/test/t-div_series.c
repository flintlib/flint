/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

int
test_div_series(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t A, B, C, D, E;
    int status = GR_SUCCESS;

    gr_ctx_init_random_commutative_ring(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);
    gr_poly_init(D, ctx);
    gr_poly_init(E, ctx);

    if (ctx->methods == _ca_methods)
        n = n_randint(state, 5);
    else
        n = n_randint(state, 20);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(C, state, 1 + n_randint(state, 20), ctx));

    /* todo: randomly make exact multiple? */

    switch (which)
    {
        case 0:
            status |= gr_poly_div_series_basecase(C, A, B, n, ctx);
            break;
        case 1:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_basecase(C, C, B, n, ctx);
            break;
        case 2:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_basecase(C, A, C, n, ctx);
            break;

        case 3:
            status |= gr_poly_div_series_newton(C, A, B, n, n_randint(state, 20), ctx);
            break;
        case 4:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_newton(C, C, B, n, n_randint(state, 20), ctx);
            break;
        case 5:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_newton(C, A, C, n, n_randint(state, 20), ctx);
            break;

        case 6:
            status |= gr_poly_div_series(C, A, B, n, ctx);
            break;
        case 7:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series(C, C, B, n, ctx);
            break;
        case 8:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series(C, A, C, n, ctx);
            break;

        case 9:
            status |= gr_poly_div_series_invmul(C, A, B, n, ctx);
            break;
        case 10:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_invmul(C, C, B, n, ctx);
            break;
        case 11:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_invmul(C, A, C, n, ctx);
            break;

        case 12:
            status |= gr_poly_div_series_divconquer(C, A, B, n, n_randint(state, 20), ctx);
            break;
        case 13:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_divconquer(C, C, B, n, n_randint(state, 20), ctx);
            break;
        case 14:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_divconquer(C, A, C, n, n_randint(state, 20), ctx);
            break;

        default:
            flint_abort();
    }

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_mullow(D, C, B, n, ctx);
        status |= gr_poly_truncate(E, A, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(D, E, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("which = %d, n = %wd\n\n", which, n);
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n\n");
            flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n\n");
            flint_printf("E = "); gr_poly_print(E, ctx); flint_printf("\n\n");
            flint_abort();
        }
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);
    gr_poly_clear(D, ctx);
    gr_poly_clear(E, ctx);

    gr_ctx_clear(ctx);

    return status;
}

/* Exact division by a series whose constant term is not a unit
   (e.g. 2 in Z): the quotient must be recovered exactly. */
static void
test_div_series_exact_nonunit(flint_rand_t state)
{
    gr_ctx_t ctx;
    gr_ptr A, B, Q, Q2;
    slong len, Alen, Blen, sz;
    int status;

    gr_ctx_init_fmpz(ctx);
    sz = ctx->sizeof_elem;
    len = 1 + n_randint(state, 4);
    Blen = 1 + n_randint(state, 4);

    GR_TMP_INIT_VEC(A, len, ctx);
    GR_TMP_INIT_VEC(B, Blen, ctx);
    GR_TMP_INIT_VEC(Q, len, ctx);
    GR_TMP_INIT_VEC(Q2, len, ctx);

    GR_MUST_SUCCEED(_gr_vec_randtest(B, state, Blen, ctx));
    GR_MUST_SUCCEED(_gr_vec_randtest(Q, state, len, ctx));
    GR_MUST_SUCCEED(gr_set_si(B, (2 + (slong) n_randint(state, 5)) * (n_randint(state, 2) ? 1 : -1), ctx));

    /* sometimes make the product shorter than len, with a quotient
       whose coefficients are not divisible by B[0] */
    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(B, 1, sz), Blen - 1, ctx));
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(Q, 1, sz), len - 1, ctx));
        GR_MUST_SUCCEED(gr_set_si(Q, 1 + n_randint(state, 100), ctx));
    }

    GR_MUST_SUCCEED(_gr_poly_mullow(A, B, Blen, Q, len, len, ctx));
    GR_MUST_SUCCEED(_gr_vec_normalise(&Alen, A, len, ctx));

    status = _gr_poly_div_series(Q2, A, Alen, B, Blen, len, ctx);

    if (status != GR_SUCCESS || _gr_vec_equal(Q, Q2, len, ctx) != T_TRUE)
    {
        flint_printf("FAIL (exact division, non-unit constant term)\n");
        flint_printf("len = %wd, Blen = %wd, status = %d\n", len, Blen, status);
        flint_printf("B = "); _gr_vec_print(B, Blen, ctx); flint_printf("\n");
        flint_printf("Q = "); _gr_vec_print(Q, len, ctx); flint_printf("\n");
        flint_printf("Q2 = "); _gr_vec_print(Q2, len, ctx); flint_printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR_VEC(A, len, ctx);
    GR_TMP_CLEAR_VEC(B, Blen, ctx);
    GR_TMP_CLEAR_VEC(Q, len, ctx);
    GR_TMP_CLEAR_VEC(Q2, len, ctx);

    gr_ctx_clear(ctx);
}

TEST_FUNCTION_START(gr_poly_div_series, state)
{
    slong iter;

    for (iter = 0; iter < 10000; iter++)
    {
        test_div_series(state, n_randint(state, 15));
    }

    for (iter = 0; iter < 1000; iter++)
        test_div_series_exact_nonunit(state);

    TEST_FUNCTION_END(state);
}
