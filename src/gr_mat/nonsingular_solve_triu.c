/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

/* todo: efficient column extraction */
int
gr_mat_nonsingular_solve_triu_classical(gr_mat_t X,
        const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong i, j, n, m;
    gr_ptr tmp, inv;
    gr_ptr s;
    int use_division = 0;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);

    n = U->r;
    m = B->c;

    /* silence compiler warning. for whatever reason gcc
       complains that inv may be uninitialized in this function
       though exactly the same code in solve_tril generates
       no warning. */
    inv = NULL;

    if (!unit)
    {
        GR_TMP_INIT_VEC(inv, n, ctx);
        for (i = 0; i < n; i++)
        {
            status = gr_inv(GR_ENTRY(inv, i, sz), GR_MAT_ENTRY(U, i, i, sz), ctx);
            if (status != GR_SUCCESS)
            {
                use_division = 1;
                status = GR_SUCCESS;
                break;
            }
        }
    }

    GR_TMP_INIT(s, ctx);
    tmp = flint_malloc(sz * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            set_shallow(GR_ENTRY(tmp, j, sz), GR_MAT_ENTRY(X, j, i, sz), ctx);

        for (j = n - 1; j >= 0; j--)
        {
            status |= _gr_vec_dot(s, GR_MAT_ENTRY(B, j, i, sz), 1, GR_MAT_ENTRY(U, j, j + 1, sz), GR_ENTRY(tmp, j + 1, sz), n - j - 1, ctx);

            if (!unit)
            {
                if (use_division)
                    status |= gr_div(GR_ENTRY(tmp, j, sz), s, GR_MAT_ENTRY(U, j, j, sz), ctx);
                else
                    status |= gr_mul(GR_ENTRY(tmp, j, sz), s, GR_ENTRY(inv, j, sz), ctx);
            }
            else
                gr_swap(GR_ENTRY(tmp, j, sz), s, ctx);

            if (status != GR_SUCCESS)
            {
                for (j = 0; j < n; j++)
                    set_shallow(GR_MAT_ENTRY(X, j, i, sz), GR_ENTRY(tmp, j, sz), ctx);
                goto cleanup;
            }
        }

        for (j = 0; j < n; j++)
            set_shallow(GR_MAT_ENTRY(X, j, i, sz), GR_ENTRY(tmp, j, sz), ctx);
    }

cleanup:
    if (!unit)
    {
        GR_TMP_CLEAR_VEC(inv, n, ctx);
    }

    flint_free(tmp);
    GR_TMP_CLEAR(s, ctx);

    return status;
}

int
gr_mat_nonsingular_solve_triu_recursive(gr_mat_t X,
        const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    gr_mat_t UA, UB, UD, XX, XY, BX, BY, T;
    slong r, n, m;
    int status = GR_SUCCESS;

    n = U->r;
    m = B->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return GR_SUCCESS;

    /*
    Denoting inv(M) by M^, we have:
    [A B]^ [X]  ==  [A^ (X - B D^ Y)]
    [0 D]  [Y]  ==  [    D^ Y       ]
    */

    gr_mat_window_init(UA, U, 0, 0, r, r, ctx);
    gr_mat_window_init(UB, U, 0, r, r, n, ctx);
    gr_mat_window_init(UD, U, r, r, n, n, ctx);
    gr_mat_window_init(BX, B, 0, 0, r, m, ctx);
    gr_mat_window_init(BY, B, r, 0, n, m, ctx);
    gr_mat_window_init(XX, X, 0, 0, r, m, ctx);
    gr_mat_window_init(XY, X, r, 0, n, m, ctx);

    status |= gr_mat_nonsingular_solve_triu(XY, UD, BY, unit, ctx);

    if (status == GR_SUCCESS)
    {
        /* gr_mat_submul(XX, BX, UB, XY); */
        gr_mat_init(T, UB->r, XY->c, ctx);
        status |= gr_mat_mul(T, UB, XY, ctx);
        status |= gr_mat_sub(XX, BX, T, ctx);
        gr_mat_clear(T, ctx);

        status |= gr_mat_nonsingular_solve_triu(XX, UA, XX, unit, ctx);
    }

    gr_mat_window_clear(UA, ctx);
    gr_mat_window_clear(UB, ctx);
    gr_mat_window_clear(UD, ctx);
    gr_mat_window_clear(BX, ctx);
    gr_mat_window_clear(BY, ctx);
    gr_mat_window_clear(XX, ctx);
    gr_mat_window_clear(XY, ctx);

    return status;
}

int
gr_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U,
                                    const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    /* todo: tune thresholds */
    if (B->r < 10 || B->c < 10)
        return gr_mat_nonsingular_solve_triu_classical(X, U, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive(X, U, B, unit, ctx);
}
