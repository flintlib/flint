/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "nfloat.h"

/*
    Tests gr_mat_hessenberg_reduce / gr_mat_hessenberg_expand_q.
    Checks: H is upper Hessenberg, Q is unitary, Q^H A Q = H, and
    Q e_{n-1} = e_{n-1} exactly (the reflectors act on indices < n-1).
*/

static void
_hred_random_float_ctx(gr_ctx_t ctx, slong * prec, flint_rand_t state)
{
    *prec = 64 + n_randint(state, 5) * 64;

    if (n_randint(state, 2))
        gr_ctx_init_complex_float_acf(ctx, *prec);
    else
        GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, *prec, 0));
}

static int
_hred_random_matrix(gr_mat_t A, slong n, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    slong sz = ctx->sizeof_elem;
    gr_ptr t, im;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, im, ctx);
    status |= gr_i(im, ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            gr_ptr e = GR_MAT_ENTRY(A, i, j, sz);
            status |= gr_set_si(e, (slong) n_randint(state, 21) - 10, ctx);
            status |= gr_set_si(t, (slong) n_randint(state, 21) - 10, ctx);
            status |= gr_mul(t, t, im, ctx);
            status |= gr_add(e, e, t, ctx);
            status |= gr_div_ui(e, e, 17, ctx);
        }
    }

    GR_TMP_CLEAR2(t, im, ctx);
    return status;
}

TEST_FUNCTION_START(gr_mat_hessenberg_reduce, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, H, Q, QH, T, U;
        gr_ptr tau, err, tol, amag;
        slong n, prec, i, j, fixed;
        slong sz;
        int status = GR_SUCCESS;

        _hred_random_float_ctx(ctx, &prec, state);
        sz = ctx->sizeof_elem;

        n = n_randint(state, 10);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(H, n, n, ctx);
        gr_mat_init(Q, n, n, ctx);
        gr_mat_init(QH, n, n, ctx);
        gr_mat_init(T, n, n, ctx);
        gr_mat_init(U, n, n, ctx);
        tau = gr_heap_init_vec(n, ctx);
        GR_TMP_INIT3(err, tol, amag, ctx);

        status |= _hred_random_matrix(A, n, state, ctx);
        status |= gr_mat_set(H, A, ctx);

        status |= gr_mat_hessenberg_reduce(H, tau, ctx);
        status |= gr_mat_set(Q, H, ctx);
        status |= gr_mat_hessenberg_expand_q(Q, tau, ctx);
        fixed = n - 1;

        /* discard reflector storage */
        for (i = 0; i < n; i++)
            for (j = i + 2; j < n; j++)
                status |= gr_zero(GR_MAT_ENTRY(H, j, i, sz), ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_abort();
        }

        status |= gr_one(tol, ctx);
        status |= gr_mul_2exp_si(tol, tol, -prec + 12, ctx);
        status |= gr_mul_ui(tol, tol, n + 1, ctx);
        status |= gr_mat_norm_max(amag, A, ctx);

        /* Q e_{n-1} = e_{n-1} exactly */
        if (n > 0)
        {
            for (i = 0; i < n; i++)
            {
                truth_t ok;
                if (i == fixed)
                    ok = gr_is_one(GR_MAT_ENTRY(Q, i, fixed, sz), ctx);
                else
                    ok = gr_is_zero(GR_MAT_ENTRY(Q, i, fixed, sz), ctx);

                if (ok == T_FALSE)
                {
                    flint_printf("FAIL: Q does not fix e_%wd\n", fixed);
                    gr_ctx_println(ctx);
                    flint_printf("Q = "); gr_mat_print(Q, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        /* Q^H Q = I */
        status |= gr_mat_conjugate_transpose(QH, Q, ctx);
        status |= gr_mat_mul(T, QH, Q, ctx);
        status |= gr_mat_one(U, ctx);
        status |= gr_mat_sub(T, T, U, ctx);
        status |= gr_mat_norm_max(err, T, ctx);

        if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
        {
            flint_printf("FAIL: Q not unitary\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("Q = "); gr_mat_print(Q, ctx); flint_printf("\n");
            flint_printf("err = "); gr_println(err, ctx);
            flint_printf("tol = "); gr_println(tol, ctx);
            flint_abort();
        }

        /* Q^H A Q = H */
        status |= gr_mat_mul(T, A, Q, ctx);
        status |= gr_mat_mul(U, QH, T, ctx);
        status |= gr_mat_sub(U, U, H, ctx);
        status |= gr_mat_norm_max(err, U, ctx);
        status |= gr_mul(tol, tol, amag, ctx);

        if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
        {
            flint_printf("FAIL: Q^H A Q != H\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("H = "); gr_mat_print(H, ctx); flint_printf("\n");
            flint_printf("Q = "); gr_mat_print(Q, ctx); flint_printf("\n");
            flint_printf("err = "); gr_println(err, ctx);
            flint_printf("tol = "); gr_println(tol, ctx);
            flint_abort();
        }

        if (gr_mat_is_hessenberg(H, ctx) == T_FALSE)
        {
            flint_printf("FAIL: H not Hessenberg\n");
            gr_ctx_println(ctx);
            flint_printf("H = "); gr_mat_print(H, ctx); flint_printf("\n");
            flint_abort();
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d in verification\n", status);
            gr_ctx_println(ctx);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(H, ctx);
        gr_mat_clear(Q, ctx);
        gr_mat_clear(QH, ctx);
        gr_mat_clear(T, ctx);
        gr_mat_clear(U, ctx);
        gr_heap_clear_vec(tau, n, ctx);
        GR_TMP_CLEAR3(err, tol, amag, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
