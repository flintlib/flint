/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/*
    Reduction to upper Hessenberg form by Householder reflections,
    retaining the reflectors so that the unitary factor can be
    reconstructed.

    On output, the upper Hessenberg part of A contains the reduced
    matrix; the strictly lower part (below the subdiagonal) together
    with tau contains the Householder reflectors.  Specifically, step i
    (for i = n-1, ..., 2) uses the reflector

        v = (A[i,0], ..., A[i,i-2], tau[i])

    of length i.  The caller is responsible for zeroing the strictly
    lower part if a clean Hessenberg matrix is wanted; see
    gr_mat_hessenberg_expand_q for reconstructing the unitary factor.

    tau must have length n.

    This is the analogue of LAPACK's zgehd2, except that the reflectors
    are stored by rows rather than by columns, which lets both
    applications of each reflector be expressed as unit-stride vector
    operations.
*/
int
gr_mat_hessenberg_reduce(gr_mat_t A, gr_ptr tau, gr_ctx_t ctx)
{
    slong i, j, k, n;
    slong sz = ctx->sizeof_elem;
    gr_ptr H, G, F, f, ff;
    gr_ptr v, cv, w;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    status |= _gr_vec_zero(tau, n, ctx);

    if (n <= 2)
        return status;

    GR_TMP_INIT5(H, G, F, f, ff, ctx);

    GR_TMP_INIT_VEC(v, n, ctx);
    GR_TMP_INIT_VEC(cv, n, ctx);
    GR_TMP_INIT_VEC(w, n, ctx);

    for (i = n - 1; i >= 2; i--)
    {
        /* cv = conj(A[i, 0:i]); H = sum_{k<i} |A[i,k]|^2 */
        for (k = 0; k < i; k++)
            status |= gr_conj(GR_ENTRY(cv, k, sz), GR_MAT_ENTRY(A, i, k, sz), ctx);

        status |= _gr_vec_dot(H, NULL, 0, cv, GR_MAT_ENTRY(A, i, 0, sz), i, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        /* The row is already zero: nothing to do. */
        if (gr_is_zero(H, ctx) == T_TRUE)
        {
            status |= gr_zero(GR_ENTRY(tau, i, sz), ctx);
            status |= gr_zero(GR_MAT_ENTRY(A, i, i - 1, sz), ctx);
            continue;
        }

        status |= gr_set(F, GR_MAT_ENTRY(A, i, i - 1, sz), ctx);
        status |= gr_abs(f, F, ctx);
        status |= gr_sqrt(G, H, ctx);

        /* A[i,i-1] = -G */
        status |= gr_neg(GR_MAT_ENTRY(A, i, i - 1, sz), G, ctx);

        if (gr_is_zero(f, ctx) == T_TRUE)
        {
            /* tau[i] = G */
            status |= gr_set(GR_ENTRY(tau, i, sz), G, ctx);
        }
        else
        {
            /* ff = F / f;  tau[i] = F + G*ff;  A[i,i-1] *= ff */
            status |= gr_div(ff, F, f, ctx);
            status |= gr_mul(GR_ENTRY(tau, i, sz), G, ff, ctx);
            status |= gr_add(GR_ENTRY(tau, i, sz), GR_ENTRY(tau, i, sz), F, ctx);
            status |= gr_mul(GR_MAT_ENTRY(A, i, i - 1, sz),
                             GR_MAT_ENTRY(A, i, i - 1, sz), ff, ctx);
        }

        /* H = 1 / sqrt(H + G*f) */
        status |= gr_addmul(H, G, f, ctx);
        status |= gr_sqrt(H, H, ctx);
        status |= gr_inv(H, H, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        status |= gr_mul(GR_ENTRY(tau, i, sz), GR_ENTRY(tau, i, sz), H, ctx);
        status |= _gr_vec_mul_scalar(GR_MAT_ENTRY(A, i, 0, sz),
                                     GR_MAT_ENTRY(A, i, 0, sz), i - 1, H, ctx);

        /* The reflector is v = (A[i,0], ..., A[i,i-2], tau[i]).  Note that
           A[i,i-1] holds the subdiagonal entry of the Hessenberg matrix,
           not a component of v, so we gather v into a temporary. */
        status |= _gr_vec_set(v, GR_MAT_ENTRY(A, i, 0, sz), i - 1, ctx);
        status |= gr_set(GR_ENTRY(v, i - 1, sz), GR_ENTRY(tau, i, sz), ctx);

        for (k = 0; k < i; k++)
            status |= gr_conj(GR_ENTRY(cv, k, sz), GR_ENTRY(v, k, sz), ctx);

        /* Apply from the right: A[0:i, 0:i] -= (A[0:i, 0:i] conj(v)) v^T.
           Rows are contiguous, so each j is one dot and one axpy. */
        for (j = 0; j < i; j++)
        {
            status |= _gr_vec_dot(G, NULL, 0, cv, GR_MAT_ENTRY(A, j, 0, sz), i, ctx);
            status |= _gr_vec_submul_scalar(GR_MAT_ENTRY(A, j, 0, sz), v, i, G, ctx);
        }

        /* Apply from the left: A[0:i, 0:n] -= conj(v) (v^T A[0:i, 0:n]).
           Accumulating w = v^T A by rows keeps everything unit-stride
           (a column-wise formulation would need strided dots). */
        status |= _gr_vec_zero(w, n, ctx);

        for (k = 0; k < i; k++)
            status |= _gr_vec_addmul_scalar(w, GR_MAT_ENTRY(A, k, 0, sz), n,
                                            GR_ENTRY(v, k, sz), ctx);

        for (k = 0; k < i; k++)
            status |= _gr_vec_submul_scalar(GR_MAT_ENTRY(A, k, 0, sz), w, n,
                                            GR_ENTRY(cv, k, sz), ctx);

        if (status != GR_SUCCESS)
            goto cleanup;
    }

cleanup:
    GR_TMP_CLEAR5(H, G, F, f, ff, ctx);
    GR_TMP_CLEAR_VEC(v, n, ctx);
    GR_TMP_CLEAR_VEC(cv, n, ctx);
    GR_TMP_CLEAR_VEC(w, n, ctx);

    return status;
}

/*
    Given the output of gr_mat_hessenberg_reduce (in Q, which must be a
    copy of the reduced matrix, still containing the reflectors) together
    with tau, overwrites Q with the unitary matrix such that

        Q^H A Q = H

    where A is the original matrix and H the Hessenberg form.

    Analogue of LAPACK's zunghr.
*/
int
gr_mat_hessenberg_expand_q(gr_mat_t Q, gr_srcptr tau, gr_ctx_t ctx)
{
    slong i, j, k, n;
    slong sz = ctx->sizeof_elem;
    gr_ptr v, cv, w;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(Q, ctx);

    if (n != gr_mat_ncols(Q, ctx))
        return GR_DOMAIN;

    if (n == 0)
        return GR_SUCCESS;

    if (n == 1)
        return gr_one(GR_MAT_ENTRY(Q, 0, 0, sz), ctx);

    /* The leading 2x2 block is the identity. */
    status |= gr_one(GR_MAT_ENTRY(Q, 0, 0, sz), ctx);
    status |= gr_one(GR_MAT_ENTRY(Q, 1, 1, sz), ctx);
    status |= gr_zero(GR_MAT_ENTRY(Q, 0, 1, sz), ctx);
    status |= gr_zero(GR_MAT_ENTRY(Q, 1, 0, sz), ctx);

    if (n == 2)
        return status;

    GR_TMP_INIT_VEC(v, n, ctx);
    GR_TMP_INIT_VEC(cv, n, ctx);
    GR_TMP_INIT_VEC(w, n, ctx);

    for (i = 2; i < n; i++)
    {
        if (gr_is_zero(GR_ENTRY(tau, i, sz), ctx) != T_TRUE)
        {
            /* v = (Q[i,0], ..., Q[i,i-2], tau[i]) */
            status |= _gr_vec_set(v, GR_MAT_ENTRY(Q, i, 0, sz), i - 1, ctx);
            status |= gr_set(GR_ENTRY(v, i - 1, sz), GR_ENTRY(tau, i, sz), ctx);

            for (k = 0; k < i; k++)
                status |= gr_conj(GR_ENTRY(cv, k, sz), GR_ENTRY(v, k, sz), ctx);

            /* Q[0:i, 0:i] -= conj(v) (v^T Q[0:i, 0:i]) */
            status |= _gr_vec_zero(w, i, ctx);

            for (k = 0; k < i; k++)
                status |= _gr_vec_addmul_scalar(w, GR_MAT_ENTRY(Q, k, 0, sz), i,
                                                GR_ENTRY(v, k, sz), ctx);

            for (k = 0; k < i; k++)
                status |= _gr_vec_submul_scalar(GR_MAT_ENTRY(Q, k, 0, sz), w, i,
                                                GR_ENTRY(cv, k, sz), ctx);
        }

        status |= gr_one(GR_MAT_ENTRY(Q, i, i, sz), ctx);

        for (j = 0; j < i; j++)
        {
            status |= gr_zero(GR_MAT_ENTRY(Q, j, i, sz), ctx);
            status |= gr_zero(GR_MAT_ENTRY(Q, i, j, sz), ctx);
        }
    }

    GR_TMP_CLEAR_VEC(v, n, ctx);
    GR_TMP_CLEAR_VEC(cv, n, ctx);
    GR_TMP_CLEAR_VEC(w, n, ctx);

    return status;
}
