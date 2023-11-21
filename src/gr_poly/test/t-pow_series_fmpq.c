/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

/* Defined in t-pow_series_fmpq.c, t-pow_series_ui.c and t-pow_ui.c */
#define test test_pow_series_fmpq
int
test(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    if (n_randint(state, 2))
    {
        /* check (A^(p/q))^q = A^p */
        fmpq_t p, q, pq;
        gr_poly_t A, Apq, Apqq, Ap;

        gr_poly_init(A, ctx);
        gr_poly_init(Apq, ctx);
        gr_poly_init(Apqq, ctx);
        gr_poly_init(Ap, ctx);
        fmpq_init(pq);
        fmpq_init(p);
        fmpq_init(q);

        fmpq_randtest(pq, state, 4);
        fmpq_set_fmpz(p, fmpq_numref(pq));
        fmpq_set_fmpz(q, fmpq_denref(pq));

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 5);
        else
            n = n_randint(state, 10);

        GR_MUST_SUCCEED(gr_poly_randtest(A, state, 10, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(Apq, state, 10, ctx));

        status |= gr_poly_pow_series_fmpq_recurrence(Apq, A, pq, n, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_pow_series_fmpq_recurrence(Apqq, Apq, q, n, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_pow_series_fmpq_recurrence(Ap, A, p, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(Apqq, Ap, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("pq = "); fmpq_print(pq); printf("\n");
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("Apq = "); gr_poly_print(Apq, ctx); flint_printf("\n");
            flint_printf("Apqq = "); gr_poly_print(Apqq, ctx); flint_printf("\n");
            flint_printf("Ap = "); gr_poly_print(Ap, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(Apq, ctx);
        gr_poly_clear(Apqq, ctx);
        gr_poly_clear(Ap, ctx);
        fmpq_clear(pq);
        fmpq_clear(p);
        fmpq_clear(q);
    }
    else
    {
        /* check A^p A^q = A^(p+q) */
        fmpq_t p, q, pq;
        gr_poly_t A, Ap, Aq, Apq, ApAq;

        gr_poly_init(A, ctx);
        gr_poly_init(Ap, ctx);
        gr_poly_init(Aq, ctx);
        gr_poly_init(Apq, ctx);
        gr_poly_init(ApAq, ctx);

        fmpq_init(pq);
        fmpq_init(p);
        fmpq_init(q);

        fmpq_randtest(p, state, 4);
        fmpq_randtest(q, state, 4);
        fmpq_add(pq, p, q);

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 5);
        else
            n = n_randint(state, 10);

        GR_MUST_SUCCEED(gr_poly_randtest(A, state, 10, ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(Aq, state, 10, ctx));

        status |= gr_poly_set(Ap, A, ctx); /* aliasing */
        status |= gr_poly_pow_series_fmpq_recurrence(Ap, Ap, p, n + n_randint(state, 3), ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_pow_series_fmpq_recurrence(Aq, A, q, n + n_randint(state, 3), ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_pow_series_fmpq_recurrence(Apq, A, pq, n, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_mullow(ApAq, Ap, Aq, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(Apq, ApAq, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("p = "); fmpq_print(p); printf("\n");
            flint_printf("q = "); fmpq_print(q); printf("\n");
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("Ap = "); gr_poly_print(Ap, ctx); flint_printf("\n");
            flint_printf("Aq = "); gr_poly_print(Aq, ctx); flint_printf("\n");
            flint_printf("Apq = "); gr_poly_print(Apq, ctx); flint_printf("\n");
            flint_printf("ApAq = "); gr_poly_print(ApAq, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(Ap, ctx);
        gr_poly_clear(Aq, ctx);
        gr_poly_clear(Apq, ctx);
        gr_poly_clear(ApAq, ctx);
        fmpq_clear(pq);
        fmpq_clear(p);
        fmpq_clear(q);

    }

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_pow_series_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test(state, 0);
    }

    TEST_FUNCTION_END(state);
}
#undef test
