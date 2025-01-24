/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

/* todo: use dot products */
int
_gr_mat_charpoly_danilevsky_inplace(gr_ptr p, gr_mat_t A, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    gr_method_binary_op addmul = GR_BINARY_OP(ctx, ADDMUL);
    slong n, n_input;
    slong i, j, k;
    gr_ptr V, W, T;
    gr_ptr t, b;
    gr_ptr c, h;
    truth_t is_zero;
    slong plen;
    slong sz = ctx->sizeof_elem;
    int have_addmul;
    int status = GR_SUCCESS;

#define MAT_ENTRY(i, j) GR_MAT_ENTRY(A, i, j, sz)
#define POL_ENTRY(i) GR_ENTRY(p, i, sz)

    have_addmul = (addmul != (gr_method_binary_op) gr_generic_addmul);

#define ADDMUL(z, x, y) \
    if (have_addmul) \
    { \
        status |= addmul(z, x, y, ctx); \
    } \
    else \
    { \
        status |= mul(c, x, y, ctx); \
        status |= add(z, z, c, ctx); \
    } \

    n = n_input = A->r;

    if (n == 0)
        return gr_one(POL_ENTRY(0), ctx);

    if (n == 1)
    {
        status |= gr_neg(POL_ENTRY(0), MAT_ENTRY(0, 0), ctx);
        status |= gr_one(POL_ENTRY(1), ctx);
        return status;
    }

    GR_TMP_INIT2(c, h, ctx);
    GR_TMP_INIT_VEC(t, n_input + 1, ctx);
    GR_TMP_INIT_VEC(b, n_input + 1, ctx);
    GR_TMP_INIT_VEC(V, n_input, ctx);
    GR_TMP_INIT_VEC(W, n_input, ctx);
    GR_TMP_INIT_VEC(T, n_input, ctx);

    status |= gr_one(p, ctx);

    i = 1;
    plen = 1;

    while (i < n)
    {
        status |= gr_set(h, MAT_ENTRY(n - i, n - i - 1), ctx);

        while (1)
        {
            is_zero = gr_is_zero(h, ctx);
            if (is_zero == T_FALSE)
                break;
            if (is_zero == T_UNKNOWN)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            k = 1;
            while (k < n - i)
            {
                is_zero = gr_is_zero(MAT_ENTRY(n - i, n - i - k - 1), ctx);
                if (is_zero == T_FALSE)
                    break;
                if (is_zero == T_UNKNOWN)
                {
                    status = GR_UNABLE;
                    goto cleanup;
                }
                k++;
            }

            if (k == n - i)
            {
                status |= gr_one(GR_ENTRY(b, i, sz), ctx);
                for (k = 1; k <= i; k++)
                    status |= gr_neg(GR_ENTRY(b, k - 1, sz), MAT_ENTRY(n - i, n - k), ctx);

                status |= _gr_poly_mul(t, p, plen, b, i + 1, ctx);
                plen += i;
                _gr_vec_swap(p, t, plen, ctx);

                n -= i;
                i = 1;

                if (n == 1)
                {
                    status |= gr_one(GR_ENTRY(b, 1, sz), ctx);
                    status |= gr_neg(b, MAT_ENTRY(0, 0), ctx);
                    status |= _gr_poly_mul(t, p, plen, b, 2, ctx);
                    plen += 1;
                    _gr_vec_swap(p, t, plen, ctx);

                    goto cleanup;
                }
            }
            else
            {
                status |= gr_mat_swap_rows(A, NULL, n - i - k - 1, n - i - 1, ctx);

                for (j = 1; j <= n - i + 1; j++)
                {
                    gr_swap(MAT_ENTRY(j - 1, n - i - k - 1), MAT_ENTRY(j - 1, n - i - 1), ctx);
                }
            }

            status |= gr_set(h, MAT_ENTRY(n - i, n - i - 1), ctx);
        }

        status |= gr_neg(h, h, ctx);
        status |= gr_inv(h, h, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        for (j = 1; j <= n; j++)
        {
            status |= mul(GR_ENTRY(V, j - 1, sz), MAT_ENTRY(n - i, j - 1), h, ctx);
            status |= gr_set(GR_ENTRY(W, j - 1, sz), MAT_ENTRY(n - i, j - 1), ctx);
        }

        status |= gr_neg(h, h, ctx);

        for (j = 1; j <= n - i; j++)
        {
            for (k = 1; k <= n - i - 1; k++)
            {
                ADDMUL(MAT_ENTRY(j - 1, k - 1), MAT_ENTRY(j - 1, n - i - 1), GR_ENTRY(V, k - 1, sz))
            }

            for (k = n - i + 1; k <= n; k++)
            {
                ADDMUL(MAT_ENTRY(j - 1, k - 1), MAT_ENTRY(j - 1, n - i - 1), GR_ENTRY(V, k - 1, sz))
            }

            status |= mul(MAT_ENTRY(j - 1, n - i - 1), MAT_ENTRY(j - 1, n - i - 1), h, ctx);
        }

        /* todo: rewrite as dot */
        for (j = 1; j <= n - i - 1; j++)
        {
            status |= mul(MAT_ENTRY(n - i - 1, j - 1), MAT_ENTRY(n - i - 1, j - 1), GR_ENTRY(W, n - i - 1, sz), ctx);

            for (k = 1; k < n - i; k++)
            {
                ADDMUL(MAT_ENTRY(n - i - 1, j - 1), MAT_ENTRY(k - 1, j - 1), GR_ENTRY(W, k - 1, sz))
            }
        }

        /* todo: rewrite as dot */
        for (j = n - i; j <= n - 1; j++)
        {
            status |= mul(MAT_ENTRY(n - i - 1, j - 1), MAT_ENTRY(n - i - 1, j - 1), GR_ENTRY(W, n - i - 1, sz), ctx);

            for (k = 1; k < n - i; k++)
            {
                ADDMUL(MAT_ENTRY(n - i - 1, j - 1), MAT_ENTRY(k - 1, j - 1), GR_ENTRY(W, k - 1, sz))
            }

            status |= add(MAT_ENTRY(n - i - 1, j - 1), MAT_ENTRY(n - i - 1, j - 1), GR_ENTRY(W, j, sz), ctx);
        }

        status |= mul(MAT_ENTRY(n - i - 1, n - 1),
                        MAT_ENTRY(n - i - 1, n - 1), GR_ENTRY(W, n - i - 1, sz), ctx);

        /* todo: rewrite as dot */
        for (k = 1; k < n - i; k++)
        {
            ADDMUL(MAT_ENTRY(n - i - 1, n - 1), MAT_ENTRY(k - 1, n - 1), GR_ENTRY(W, k - 1, sz))
        }

        i++;
    }

    status |= gr_one(GR_ENTRY(b, n, sz), ctx);
    for (i = 1; i <= n; i++)
        status |= gr_neg(GR_ENTRY(b, i - 1, sz), MAT_ENTRY(0, n - i), ctx);
    status |= _gr_poly_mul(t, p, plen, b, n + 1, ctx);
    plen += n;
    _gr_vec_swap(p, t, plen, ctx);

cleanup:

    GR_TMP_CLEAR2(c, h, ctx);
    GR_TMP_CLEAR_VEC(t, n_input + 1, ctx);
    GR_TMP_CLEAR_VEC(b, n_input + 1, ctx);
    GR_TMP_CLEAR_VEC(V, n_input, ctx);
    GR_TMP_CLEAR_VEC(W, n_input, ctx);
    GR_TMP_CLEAR_VEC(T, n_input, ctx);

    return status;
}

int
_gr_mat_charpoly_danilevsky(gr_ptr p, const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_mat_t T;
    int status = GR_SUCCESS;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_mat_init(T, mat->r, mat->c, ctx);
    status |= gr_mat_set(T, mat, ctx);
    status |= _gr_mat_charpoly_danilevsky_inplace(p, T, ctx);
    gr_mat_clear(T, ctx);
    return status;

}

int
gr_mat_charpoly_danilevsky(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;
    gr_poly_fit_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_danilevsky(cp->coeffs, mat, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
