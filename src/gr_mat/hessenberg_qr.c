/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018, 2026 Fredrik Johansson

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
    All scalar work here (Givens parameters, deflation tests, shifts) is
    O(n) per sweep, against O(n^2) for the rotations themselves, so we do
    it in the working ring rather than in a separate real or magnitude
    ring.  Ordering is only ever needed between magnitudes, which
    gr_cmpabs supplies; this keeps the code free of any assumption about
    how a complex ring stores its components.
*/

/*
    A plane rotation with real cosine: c, s and r with
    [c, conj(s); -s, c] [x; y] = [r; 0] and c^2 + |s|^2 = 1, c >= 0.
    With v = sqrt(|x|^2 + |y|^2) and x != 0,

        c = |x| / v,   s = conj(x) y / (|x| v),   r = x v / |x|

    (for x = 0: c = 0, s = y / v, r = v).  A real cosine makes the
    diagonal entries of the 2 x 2 transformation real, which saves a
    quarter of the real multiplications in each rotation applied (as
    gr_mul of a complex number by an exactly real one, or a
    _gr_vec_apply_2x2 that checks for real entries).  gr_abs2 is exact
    on the imaginary part, so |x|^2 and everything computed from it
    alone is a genuine non-negative real element.  r may be NULL.
*/
static int
_gr_givens(gr_ptr c, gr_ptr s, gr_ptr r, gr_srcptr x, gr_srcptr y,
    gr_ptr t1, gr_ptr t2, gr_ptr t3, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_abs2(t1, x, ctx);
    status |= gr_abs2(t2, y, ctx);
    status |= gr_add(t2, t1, t2, ctx);      /* v^2 */

    if (gr_is_zero(t2, ctx) == T_TRUE)
    {
        status |= gr_one(c, ctx);
        status |= gr_zero(s, ctx);
        if (r != NULL)
            status |= gr_zero(r, ctx);
    }
    else if (gr_is_zero(t1, ctx) == T_TRUE)
    {
        status |= gr_sqrt(t2, t2, ctx);
        if (r != NULL)
            status |= gr_set(r, t2, ctx);
        status |= gr_div(s, y, t2, ctx);
        status |= gr_zero(c, ctx);
    }
    else
    {
        status |= gr_sqrt(t3, t1, ctx);     /* |x| */
        status |= gr_sqrt(c, t2, ctx);      /* v */
        status |= gr_mul(t3, t3, c, ctx);
        status |= gr_inv(t3, t3, ctx);      /* w = 1 / (|x| v) */
        if (r != NULL)
        {
            status |= gr_mul(t2, t2, t3, ctx);
            status |= gr_mul(r, x, t2, ctx);  /* x v^2 w */
        }
        status |= gr_conj(s, x, ctx);
        status |= gr_mul(s, s, y, ctx);
        status |= gr_mul(s, s, t3, ctx);    /* conj(x) y w */
        status |= gr_mul(c, t1, t3, ctx);   /* |x|^2 w */
    }

    return status;
}

/*
    One implicitly shifted QR sweep on the active submatrix
    A[n0:n1, n0:n1].  If Q is not NULL, the accumulated transformation is
    applied to its columns.
*/
int
gr_mat_hessenberg_qr_step(gr_mat_t A, gr_mat_t Q, slong n0, slong n1, gr_srcptr shift, gr_ctx_t ctx)
{
    slong j, n;
    slong sz = ctx->sizeof_elem;
    gr_ptr c, s, cc, cs, v, t1, t2, negs, negcs;
    gr_ptr M;       /* ML = M[0:4], MR = M[4:8] */
    gr_ptr ML, MR;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    GR_TMP_INIT5(c, s, cc, cs, v, ctx);
    GR_TMP_INIT2(t1, t2, ctx);
    GR_TMP_INIT2(negs, negcs, ctx);

    /* ML and MR only ever hold shallow references to the rotation
       parameters, which stay alive for as long as they are in use.
       Building them with gr_set would mean eight deep copies (and
       potentially eight allocations) per rotation. */
    M = GR_TMP_ALLOC(8 * sz);
    ML = M;
    MR = GR_ENTRY(M, 4, sz);

/*
    Given a rotation (c, s) with |c|^2 + |s|^2 = 1, the left action on a
    pair of rows is [conj(c), conj(s); -s, c] and the right action on a
    pair of columns is [c, s; -conj(s), conj(c)] (c is real here, but
    conj(c) costs nothing and keeps this valid for any c).
*/
#define SET_ROTATION() \
    do { \
        status |= gr_conj(cc, c, ctx); \
        status |= gr_conj(cs, s, ctx); \
        status |= gr_neg(negs, s, ctx); \
        status |= gr_neg(negcs, cs, ctx); \
        gr_set_shallow(GR_ENTRY(ML, 0, sz), cc, ctx); \
        gr_set_shallow(GR_ENTRY(ML, 1, sz), cs, ctx); \
        gr_set_shallow(GR_ENTRY(ML, 2, sz), negs, ctx); \
        gr_set_shallow(GR_ENTRY(ML, 3, sz), c, ctx); \
        gr_set_shallow(GR_ENTRY(MR, 0, sz), c, ctx); \
        gr_set_shallow(GR_ENTRY(MR, 1, sz), s, ctx); \
        gr_set_shallow(GR_ENTRY(MR, 2, sz), negcs, ctx); \
        gr_set_shallow(GR_ENTRY(MR, 3, sz), cc, ctx); \
    } while (0)

    /* Initial rotation, determined by the shift. */
    status |= gr_sub(v, GR_MAT_ENTRY(A, n0, n0, sz), shift, ctx);
    status |= _gr_givens(c, s, NULL, v, GR_MAT_ENTRY(A, n0 + 1, n0, sz), t1, t2, cc, ctx);

    SET_ROTATION();

    status |= gr_mat_apply_2x2_rows(A, n0, n0 + 1, n0, n, ML, ctx);
    status |= gr_mat_apply_2x2_cols(A, n0, n0 + 1, 0, FLINT_MIN(n1, n0 + 3), MR, ctx);

    if (Q != NULL)
        status |= gr_mat_apply_2x2_cols(Q, n0, n0 + 1, 0, n, MR, ctx);

    /* Chase the bulge down the subdiagonal. */
    for (j = n0; j < n1 - 2; j++)
    {
        status |= _gr_givens(c, s, v, GR_MAT_ENTRY(A, j + 1, j, sz),
            GR_MAT_ENTRY(A, j + 2, j, sz), t1, t2, cc, ctx);
        gr_swap(GR_MAT_ENTRY(A, j + 1, j, sz), v, ctx);
        status |= gr_zero(GR_MAT_ENTRY(A, j + 2, j, sz), ctx);

        SET_ROTATION();

        status |= gr_mat_apply_2x2_rows(A, j + 1, j + 2, j + 1, n, ML, ctx);
        status |= gr_mat_apply_2x2_cols(A, j + 1, j + 2, 0, FLINT_MIN(n1, j + 4), MR, ctx);

        if (Q != NULL)
            status |= gr_mat_apply_2x2_cols(Q, j + 1, j + 2, 0, n, MR, ctx);
    }

#undef SET_ROTATION

    GR_TMP_CLEAR5(c, s, cc, cs, v, ctx);
    GR_TMP_CLEAR2(t1, t2, ctx);
    GR_TMP_CLEAR2(negs, negcs, ctx);
    GR_TMP_FREE(M, 8 * sz);

    return status;
}


/*
    Runs the shifted QR algorithm on the Hessenberg matrix A until it is
    upper triangular (to within tol), accumulating the transformation
    into Q if Q is not NULL.

    tol may be NULL, in which case a tolerance is derived from the
    working precision of ctx.  maxiter <= 0 selects a default.

    Returns GR_UNABLE if the iteration fails to converge.
*/
int
gr_mat_hessenberg_qr(gr_mat_t A, gr_mat_t Q, gr_srcptr tol, slong maxiter, gr_ctx_t ctx)
{
    slong n, k, n0, n1, iter, total_iter, prec;
    slong sz = ctx->sizeof_elem;
    gr_ptr norm, eps, ts, tm, um;
    gr_ptr shift, s, t, a, b;
    int cmp;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n <= 1)
        return GR_SUCCESS;

    if (gr_ctx_get_real_prec(&prec, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    GR_TMP_INIT5(norm, eps, ts, tm, um, ctx);
    GR_TMP_INIT5(shift, s, t, a, b, ctx);

    /* A is Hessenberg, so this is the norm of its nonzero part. */
    status |= gr_mat_norm_frobenius(norm, A, ctx);
    status |= gr_div_ui(norm, norm, n, ctx);

    if (status != GR_SUCCESS)
        goto cleanup;

    if (gr_is_zero(norm, ctx) == T_TRUE)
        goto cleanup;

    if (tol == NULL)
    {
        status |= gr_one(eps, ctx);
        status |= gr_mul_2exp_si(eps, eps, -prec, ctx);
        status |= gr_div_ui(eps, eps, 100 * n, ctx);
    }
    else
    {
        status |= gr_set(eps, tol, ctx);
    }

    if (maxiter <= 0)
        maxiter = 14 * n + prec / 10;

    /* The active submatrix is A[n0:n1, n0:n1]. */
    n0 = 0;
    n1 = n;

    iter = total_iter = 0;

    while (1)
    {
        k = n0;

        while (k + 1 < n1)
        {
            /* ts = max(|A[k,k]|, |A[k+1,k+1]|), kept as a ring element:
               only its magnitude is ever used, so no square root is
               needed.  (The acb code uses |a| + |b|; the difference is a
               factor of at most two, absorbed by the tolerance.) */
            if (gr_abs_ge(GR_MAT_ENTRY(A, k, k, sz),
                          GR_MAT_ENTRY(A, k + 1, k + 1, sz), ctx) == T_TRUE)
                status |= gr_set(ts, GR_MAT_ENTRY(A, k, k, sz), ctx);
            else
                status |= gr_set(ts, GR_MAT_ENTRY(A, k + 1, k + 1, sz), ctx);

            /* if ts < eps * norm, ts = norm */
            status |= gr_mul(tm, eps, norm, ctx);
            if (gr_abs_lt(ts, tm, ctx) == T_TRUE)
                status |= gr_set(ts, norm, ctx);

            /* if |A[k+1,k]| < eps * ts, deflate */
            status |= gr_mul(um, eps, ts, ctx);
            if (gr_abs_lt(GR_MAT_ENTRY(A, k + 1, k, sz), um, ctx) == T_TRUE)
                break;

            k++;
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        if (k + 1 < n1)
        {
            /* Deflation found at position (k+1, k). */
            status |= gr_zero(GR_MAT_ENTRY(A, k + 1, k, sz), ctx);
            n0 = k + 1;
            iter = 0;

            if (n0 + 1 >= n1)
            {
                /* A block of size at most two has converged. */
                n0 = 0;
                n1 = k + 1;

                if (n1 < 2)
                {
                    status = GR_SUCCESS;
                    goto cleanup;
                }
            }
        }
        else
        {
            if (iter % 30 == 10)
            {
                status |= gr_set(shift, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz), ctx);
            }
            else if (iter % 30 == 20)
            {
                status |= gr_abs(shift, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz), ctx);
            }
            else if (iter % 30 == 29)
            {
                status |= gr_set(shift, norm, ctx);
            }
            else
            {
                /* Wilkinson shift: the eigenvalue of the trailing 2x2
                   block closest to A[n1-1,n1-1]. */
                status |= gr_add(t, GR_MAT_ENTRY(A, n1 - 2, n1 - 2, sz),
                                    GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), ctx);
                status |= gr_sub(a, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz),
                                    GR_MAT_ENTRY(A, n1 - 2, n1 - 2, sz), ctx);
                status |= gr_sqr(a, a, ctx);
                status |= gr_mul(b, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz),
                                    GR_MAT_ENTRY(A, n1 - 2, n1 - 1, sz), ctx);
                status |= gr_mul_2exp_si(b, b, 2, ctx);
                status |= gr_add(s, a, b, ctx);

                /* A complex square root handles the negative-real case
                   directly; no manual branch selection is needed. */
                status |= gr_sqrt(s, s, ctx);

                status |= gr_add(a, t, s, ctx);
                status |= gr_sub(b, t, s, ctx);
                status |= gr_mul_2exp_si(a, a, -1, ctx);
                status |= gr_mul_2exp_si(b, b, -1, ctx);

                status |= gr_sub(s, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), a, ctx);
                status |= gr_sub(t, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), b, ctx);

                status |= gr_cmpabs(&cmp, s, t, ctx);

                if (status != GR_SUCCESS)
                    goto cleanup;

                if (cmp > 0)
                    status |= gr_set(shift, b, ctx);
                else
                    status |= gr_set(shift, a, ctx);
            }

            iter++;
            total_iter++;

            status |= gr_mat_hessenberg_qr_step(A, Q, n0, n1, shift, ctx);

            if (status != GR_SUCCESS)
                goto cleanup;

            if (iter > maxiter)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
        }
    }

cleanup:
    GR_TMP_CLEAR5(norm, eps, ts, tm, um, ctx);
    GR_TMP_CLEAR5(shift, s, t, a, b, ctx);

    return status;
}
