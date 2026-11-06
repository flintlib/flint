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
    gr_mat_schur_reorder moves T[ifst,ifst] to position ilst by a
    unitary similarity.  Checks: T stays upper triangular with the
    lower triangle exactly zero, the diagonal is exactly the expected
    cyclic permutation of the original (the swapped entries are assigned
    rather than left to cancellation), V is unitary and V^H T0 V = T.
*/

static void
_schr_random_float_ctx(gr_ctx_t ctx, slong * prec, flint_rand_t state)
{
    *prec = 64 + n_randint(state, 5) * 64;

    if (n_randint(state, 2))
        gr_ctx_init_complex_float_acf(ctx, *prec);
    else
        GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, *prec, 0));
}

TEST_FUNCTION_START(gr_mat_schur_reorder, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t T, T0, V, VH, W, X;
        gr_ptr d0, t, im, err, tol, tmag;
        slong n, prec, i, j, ifst, ilst;
        slong sz;
        slong * perm;
        int status = GR_SUCCESS;

        _schr_random_float_ctx(ctx, &prec, state);
        sz = ctx->sizeof_elem;

        n = 1 + n_randint(state, 12);
        ifst = n_randint(state, n);
        ilst = n_randint(state, n);

        gr_mat_init(T, n, n, ctx);
        gr_mat_init(T0, n, n, ctx);
        gr_mat_init(V, n, n, ctx);
        gr_mat_init(VH, n, n, ctx);
        gr_mat_init(W, n, n, ctx);
        gr_mat_init(X, n, n, ctx);
        d0 = gr_heap_init_vec(n, ctx);
        GR_TMP_INIT5(t, im, err, tol, tmag, ctx);
        perm = flint_malloc(n * sizeof(slong));

        status |= gr_i(im, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = i; j < n; j++)
            {
                gr_ptr e = GR_MAT_ENTRY(T, i, j, sz);
                status |= gr_set_si(e, (slong) n_randint(state, 21) - 10, ctx);
                status |= gr_set_si(t, (slong) n_randint(state, 21) - 10, ctx);
                status |= gr_mul(t, t, im, ctx);
                status |= gr_add(e, e, t, ctx);
                status |= gr_div_ui(e, e, 7, ctx);
            }
        }

        status |= gr_mat_set(T0, T, ctx);
        for (i = 0; i < n; i++)
            status |= gr_set(GR_ENTRY(d0, i, sz), GR_MAT_ENTRY(T, i, i, sz), ctx);
        status |= gr_mat_one(V, ctx);

        status |= gr_mat_schur_reorder(T, V, ifst, ilst, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("T0 = "); gr_mat_print(T0, ctx); flint_printf("\n");
            flint_abort();
        }

        /* exact structural checks */
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                if (gr_is_zero(GR_MAT_ENTRY(T, i, j, sz), ctx) != T_TRUE)
                {
                    flint_printf("FAIL: lower triangle not exactly zero\n");
                    gr_ctx_println(ctx);
                    flint_printf("T = "); gr_mat_print(T, ctx); flint_printf("\n");
                    flint_abort();
                }

        for (i = 0; i < n; i++)
            perm[i] = i;
        if (ifst > ilst)
            for (i = ifst; i > ilst; i--) { slong s = perm[i]; perm[i] = perm[i - 1]; perm[i - 1] = s; }
        else
            for (i = ifst + 1; i <= ilst; i++) { slong s = perm[i]; perm[i] = perm[i - 1]; perm[i - 1] = s; }

        for (i = 0; i < n; i++)
            if (gr_equal(GR_MAT_ENTRY(T, i, i, sz), GR_ENTRY(d0, perm[i], sz), ctx) != T_TRUE)
            {
                flint_printf("FAIL: diagonal not exactly permuted (ifst = %wd, ilst = %wd)\n", ifst, ilst);
                gr_ctx_println(ctx);
                flint_printf("T0 = "); gr_mat_print(T0, ctx); flint_printf("\n");
                flint_printf("T = "); gr_mat_print(T, ctx); flint_printf("\n");
                flint_abort();
            }

        /* approximate checks */
        status |= gr_one(tol, ctx);
        status |= gr_mul_2exp_si(tol, tol, -prec + 12, ctx);
        status |= gr_mul_ui(tol, tol, n + 1, ctx);
        status |= gr_mat_norm_max(tmag, T0, ctx);

        status |= gr_mat_conjugate_transpose(VH, V, ctx);
        status |= gr_mat_mul(W, VH, V, ctx);
        status |= gr_mat_one(X, ctx);
        status |= gr_mat_sub(W, W, X, ctx);
        status |= gr_mat_norm_max(err, W, ctx);

        if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
        {
            flint_printf("FAIL: V not unitary\n");
            gr_ctx_println(ctx);
            flint_printf("err = "); gr_println(err, ctx);
            flint_printf("tol = "); gr_println(tol, ctx);
            flint_abort();
        }

        status |= gr_mat_mul(X, T0, V, ctx);
        status |= gr_mat_mul(W, VH, X, ctx);
        status |= gr_mat_sub(W, W, T, ctx);
        status |= gr_mat_norm_max(err, W, ctx);
        status |= gr_mul(tol, tol, tmag, ctx);

        if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
        {
            flint_printf("FAIL: V^H T0 V != T\n");
            gr_ctx_println(ctx);
            flint_printf("T0 = "); gr_mat_print(T0, ctx); flint_printf("\n");
            flint_printf("T = "); gr_mat_print(T, ctx); flint_printf("\n");
            flint_printf("err = "); gr_println(err, ctx);
            flint_printf("tol = "); gr_println(tol, ctx);
            flint_abort();
        }

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: status = %d in verification\n", status);
            flint_abort();
        }

        gr_mat_clear(T, ctx);
        gr_mat_clear(T0, ctx);
        gr_mat_clear(V, ctx);
        gr_mat_clear(VH, ctx);
        gr_mat_clear(W, ctx);
        gr_mat_clear(X, ctx);
        gr_heap_clear_vec(d0, n, ctx);
        GR_TMP_CLEAR5(t, im, err, tol, tmag, ctx);
        flint_free(perm);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
