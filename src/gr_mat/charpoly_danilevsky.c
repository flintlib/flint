/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

int
_gr_mat_charpoly_danilevsky_inplace(gr_ptr p, gr_mat_t A, gr_ctx_t ctx)
{
    slong n, n_input;
    slong i, j, k;
    slong Astride = A->stride;
    gr_ptr V, W, T;
    gr_ptr t, b;
    gr_ptr c, h;
    truth_t is_zero;
    slong plen;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

#define MAT_ENTRY(i, j) GR_MAT_ENTRY(A, i, j, sz)
#define POL_ENTRY(i) GR_ENTRY(p, i, sz)

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

        status |= _gr_vec_set(W, MAT_ENTRY(n - i, 0), n, ctx);
        status |= _gr_vec_mul_scalar(V, MAT_ENTRY(n - i, 0), n, h, ctx);

        status |= gr_neg(h, h, ctx);

        for (j = 1; j <= n - i; j++)
        {
            gr_srcptr c = MAT_ENTRY(j - 1, n - i - 1);
            gr_ptr row = MAT_ENTRY(j - 1, 0);

            status |= _gr_vec_addmul_scalar(row, V, n - i - 1, c, ctx);
            status |= _gr_vec_addmul_scalar(GR_ENTRY(row, n - i, sz), GR_ENTRY(V, n - i, sz), i, c, ctx);
            status |= gr_mul(MAT_ENTRY(j - 1, n - i - 1), MAT_ENTRY(j - 1, n - i - 1), h, ctx);
        }

        /* The following dot products need to write to a temporary, because the
           destination is aliased with an input vector and some dot implementations
           don't support this aliasing. */

        for (j = 1; j <= n - i - 1; j++)
        {
            status |= _gr_vec_dot_strided(t, NULL, 0, MAT_ENTRY(0, j - 1), Astride, W, 1, n - i, ctx);
            gr_swap(MAT_ENTRY(n - i - 1, j - 1), t, ctx);
        }

        for (j = n - i; j <= n - 1; j++)
        {
            status |= _gr_vec_dot_strided(t, GR_ENTRY(W, j, sz), 0,
                MAT_ENTRY(0, j - 1), Astride, W, 1, n - i, ctx);
            gr_swap(MAT_ENTRY(n - i - 1, j - 1), t, ctx);
        }

        status |= _gr_vec_dot_strided(t, NULL, 0, MAT_ENTRY(0, n - 1), Astride, W, 1, n - i, ctx);
        gr_swap(MAT_ENTRY(n - i - 1, n - 1), t, ctx);

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
