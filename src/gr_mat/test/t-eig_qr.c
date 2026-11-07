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

/* A random complex floating-point context: acf, nfloat_complex or a
   complex dfloat context. */
static void
_aeq_random_float_ctx(gr_ctx_t ctx, slong * prec, flint_rand_t state)
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

/* Entries are small Gaussian integers scaled by 1/17: a well-conditioned
   family whose eigenvalues are generically simple.  With probability
   1/4 the matrix is made upper triangular with a repeated diagonal so
   that the non-simple branches of the eigenvector solve are exercised;
   residuals stay small there too, since only backward stability is
   being tested. */
static int
_aeq_random_matrix(gr_mat_t A, slong n, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    slong sz = ctx->sizeof_elem;
    gr_ptr t, im;
    int status = GR_SUCCESS;
    int triangular = (n_randint(state, 4) == 0);

    GR_TMP_INIT2(t, im, ctx);
    status |= gr_i(im, ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            gr_ptr e = GR_MAT_ENTRY(A, i, j, sz);

            if (triangular && j < i)
            {
                status |= gr_zero(e, ctx);
                continue;
            }

            if (triangular && i == j)
            {
                /* repeated diagonal entries */
                status |= gr_set_si(e, (slong) n_randint(state, 3), ctx);
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

/* err <- || A R - R diag(E) ||_max  (right = 1)
   err <- || L A - diag(E) L ||_max  (right = 0) */
static int
_aeq_residual(gr_ptr err, const gr_mat_t A, const gr_mat_t V, gr_srcptr E, int right, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(A, ctx), i, j;
    slong sz = ctx->sizeof_elem;
    gr_mat_t T, S;
    int status = GR_SUCCESS;

    gr_mat_init(T, n, n, ctx);
    gr_mat_init(S, n, n, ctx);

    if (right)
    {
        status |= gr_mat_mul(T, A, V, ctx);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                status |= gr_mul(GR_MAT_ENTRY(S, i, j, sz), GR_MAT_ENTRY(V, i, j, sz), GR_ENTRY(E, j, sz), ctx);
    }
    else
    {
        status |= gr_mat_mul(T, V, A, ctx);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                status |= gr_mul(GR_MAT_ENTRY(S, i, j, sz), GR_MAT_ENTRY(V, i, j, sz), GR_ENTRY(E, i, sz), ctx);
    }

    status |= gr_mat_sub(T, T, S, ctx);
    status |= gr_mat_norm_max(err, T, ctx);

    gr_mat_clear(T, ctx);
    gr_mat_clear(S, ctx);
    return status;
}

TEST_FUNCTION_START(gr_mat_eig_qr, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, L, R;
        gr_ptr E, err, tol, amag, vmag, tr, sum;
        slong n, prec, i;
        slong sz;
        int status = GR_SUCCESS;
        int haveL, haveR;

        _aeq_random_float_ctx(ctx, &prec, state);
        sz = ctx->sizeof_elem;

        n = n_randint(state, 13);
        haveL = n_randint(state, 2);
        haveR = n_randint(state, 2);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(L, n, n, ctx);
        gr_mat_init(R, n, n, ctx);
        E = gr_heap_init_vec(n, ctx);
        GR_TMP_INIT5(err, tol, amag, vmag, tr, ctx);
        GR_TMP_INIT(sum, ctx);

        status |= _aeq_random_matrix(A, n, state, ctx);

        status |= gr_mat_eig_qr(E, haveL ? L : NULL, haveR ? R : NULL,
                                       A, NULL, 0, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_abort();
        }

        /* Backward-stable QR gives residuals of order n * 2^-prec * ||A||;
           allow a generous constant. */
        status |= gr_one(tol, ctx);
        status |= gr_mul_2exp_si(tol, tol, -prec + 16, ctx);
        status |= gr_mul_ui(tol, tol, n + 1, ctx);
        status |= gr_mat_norm_max(amag, A, ctx);

        /* Trace check: sum of eigenvalues = trace(A). */
        status |= gr_mat_trace(tr, A, ctx);
        status |= gr_zero(sum, ctx);
        for (i = 0; i < n; i++)
            status |= gr_add(sum, sum, GR_ENTRY(E, i, sz), ctx);
        status |= gr_sub(sum, sum, tr, ctx);
        status |= gr_abs(err, sum, ctx);
        status |= gr_mul(vmag, tol, amag, ctx);
        status |= gr_mul_ui(vmag, vmag, n + 1, ctx);

        if (status == GR_SUCCESS && gr_le(err, vmag, ctx) == T_FALSE)
        {
            flint_printf("FAIL: eigenvalue sum != trace\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("E = "); _gr_vec_print(E, n, ctx); flint_printf("\n");
            flint_printf("err = "); gr_println(err, ctx);
            flint_printf("tol = "); gr_println(vmag, ctx);
            flint_abort();
        }

        if (haveR)
        {
            status |= _aeq_residual(err, A, R, E, 1, ctx);
            status |= gr_mat_norm_max(vmag, R, ctx);
            status |= gr_mul(vmag, vmag, amag, ctx);
            status |= gr_mul(vmag, vmag, tol, ctx);

            if (status == GR_SUCCESS && gr_le(err, vmag, ctx) == T_FALSE)
            {
                flint_printf("FAIL: right eigenvector residual\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("R = "); gr_mat_print(R, ctx); flint_printf("\n");
                flint_printf("E = "); _gr_vec_print(E, n, ctx); flint_printf("\n");
                flint_printf("err = "); gr_println(err, ctx);
                flint_printf("tol = "); gr_println(vmag, ctx);
                flint_abort();
            }
        }

        if (haveL)
        {
            status |= _aeq_residual(err, A, L, E, 0, ctx);
            status |= gr_mat_norm_max(vmag, L, ctx);
            status |= gr_mul(vmag, vmag, amag, ctx);
            status |= gr_mul(vmag, vmag, tol, ctx);

            if (status == GR_SUCCESS && gr_le(err, vmag, ctx) == T_FALSE)
            {
                flint_printf("FAIL: left eigenvector residual\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("L = "); gr_mat_print(L, ctx); flint_printf("\n");
                flint_printf("E = "); _gr_vec_print(E, n, ctx); flint_printf("\n");
                flint_printf("err = "); gr_println(err, ctx);
                flint_printf("tol = "); gr_println(vmag, ctx);
                flint_abort();
            }
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d in verification\n", status);
            gr_ctx_println(ctx);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(L, ctx);
        gr_mat_clear(R, ctx);
        gr_heap_clear_vec(E, n, ctx);
        GR_TMP_CLEAR5(err, tol, amag, vmag, tr, ctx);
        GR_TMP_CLEAR(sum, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
