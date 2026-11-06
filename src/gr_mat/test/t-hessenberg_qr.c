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
#include "dfloat.h"

/*
    Tests the QR iteration on a Hessenberg matrix: the result must be
    upper triangular, Q unitary, and Q^H H Q equal to the result.
*/

static void
_hqr_random_float_ctx(gr_ctx_t ctx, slong * prec, flint_rand_t state)
{
    if (n_randint(state, 3) == 0 && dfloat_is_supported())
    {
        /* d1-d4 complex, plain or with canonical results */
        int n = 1 + n_randint(state, 4);
        *prec = 53 * n;
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, DFLOAT_COMPLEX
            | (n_randint(state, 4) == 0 ? DFLOAT_STRONG : 0)));
        return;
    }

    *prec = 64 + n_randint(state, 5) * 64;

    if (n_randint(state, 2))
        gr_ctx_init_complex_float_acf(ctx, *prec);
    else
        GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, *prec, 0));
}

static int
_hqr_random_hessenberg(gr_mat_t H, slong n, flint_rand_t state, gr_ctx_t ctx)
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
            gr_ptr e = GR_MAT_ENTRY(H, i, j, sz);

            if (i > j + 1)
            {
                status |= gr_zero(e, ctx);
                continue;
            }

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

static void
_hqr_check(const gr_mat_t H, const gr_mat_t T, const gr_mat_t Q, slong prec, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(H, ctx);
    gr_mat_t QH, W, X;
    gr_ptr err, tol, hmag;
    int status = GR_SUCCESS;

    gr_mat_init(QH, n, n, ctx);
    gr_mat_init(W, n, n, ctx);
    gr_mat_init(X, n, n, ctx);
    GR_TMP_INIT3(err, tol, hmag, ctx);

    status |= gr_one(tol, ctx);
    status |= gr_mul_2exp_si(tol, tol, -prec + 16, ctx);
    status |= gr_mul_ui(tol, tol, n + 1, ctx);
    status |= gr_mat_norm_max(hmag, H, ctx);

    if (gr_mat_is_upper_triangular(T, ctx) == T_FALSE)
    {
        flint_printf("FAIL: result not upper triangular\n");
        gr_ctx_println(ctx);
        flint_printf("H = "); gr_mat_print(H, ctx); flint_printf("\n");
        flint_printf("T = "); gr_mat_print(T, ctx); flint_printf("\n");
        flint_abort();
    }

    status |= gr_mat_conjugate_transpose(QH, Q, ctx);
    status |= gr_mat_mul(W, QH, Q, ctx);
    status |= gr_mat_one(X, ctx);
    status |= gr_mat_sub(W, W, X, ctx);
    status |= gr_mat_norm_max(err, W, ctx);

    if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
    {
        flint_printf("FAIL: Q not unitary\n");
        gr_ctx_println(ctx);
        flint_printf("n = %wd\n", n);
        flint_printf("err = "); gr_println(err, ctx);
        flint_printf("tol = "); gr_println(tol, ctx);
        flint_abort();
    }

    status |= gr_mat_mul(X, H, Q, ctx);
    status |= gr_mat_mul(W, QH, X, ctx);
    status |= gr_mat_sub(W, W, T, ctx);
    status |= gr_mat_norm_max(err, W, ctx);
    status |= gr_mul(tol, tol, hmag, ctx);

    if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
    {
        flint_printf("FAIL: Q^H H Q != T\n");
        gr_ctx_println(ctx);
        flint_printf("n = %wd\n", n);
        flint_printf("err = "); gr_println(err, ctx);
        flint_printf("tol = "); gr_println(tol, ctx);
        flint_abort();
    }

    if (status != GR_SUCCESS)
    {
        flint_printf("FAIL: status = %d in verification\n", status);
        flint_abort();
    }

    gr_mat_clear(QH, ctx);
    gr_mat_clear(W, ctx);
    gr_mat_clear(X, ctx);
    GR_TMP_CLEAR3(err, tol, hmag, ctx);
}

TEST_FUNCTION_START(gr_mat_hessenberg_qr, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t H, T, Q;
        slong n, prec;
        int status = GR_SUCCESS;

        _hqr_random_float_ctx(ctx, &prec, state);
        n = n_randint(state, 14);

        gr_mat_init(H, n, n, ctx);
        gr_mat_init(T, n, n, ctx);
        gr_mat_init(Q, n, n, ctx);

        status |= _hqr_random_hessenberg(H, n, state, ctx);
        status |= gr_mat_set(T, H, ctx);
        status |= gr_mat_one(Q, ctx);
        status |= gr_mat_hessenberg_qr(T, Q, NULL, 0, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("H = "); gr_mat_print(H, ctx); flint_printf("\n");
            flint_abort();
        }

        _hqr_check(H, T, Q, prec, ctx);

        gr_mat_clear(H, ctx);
        gr_mat_clear(T, ctx);
        gr_mat_clear(Q, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
