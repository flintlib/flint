/*
    Copyright (C) 2026 Fredrik Johansson

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
    Reordering of a complex Schur form.  Over an algebraically closed
    (or numerically complex) ring every block is 1x1, so an adjacent
    swap is a single Givens rotation and no 2x2 case arises -- this is
    the analogue of LAPACK's ztrexc, not dtrexc.

    Aggressive early deflation needs this: after the trailing window is
    reduced to Schur form, eigenvalues whose spike entry is not
    negligible have to be moved out of the way so that the search can
    continue past them.  It is also useful on its own, for grouping
    selected eigenvalues into a leading invariant subspace.
*/

/*
    Swaps the diagonal entries T[k-1,k-1] and T[k,k] of the upper
    triangular matrix T by a unitary similarity, and postmultiplies V by
    the same rotation (V may be NULL).

    The eigenvector of the 2x2 block [[t11, t12], [0, t22]] belonging to
    t22 is (t12, t22 - t11); normalising it gives the first column of
    the rotation, and (-conj(s), conj(c)) completes the basis.
*/
static int
gr_mat_schur_swap_adjacent(gr_mat_t T, gr_mat_t V, slong k, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(T, ctx);
    slong sz = ctx->sizeof_elem;
    gr_ptr c, s, cc, cs, v, t1, t2, t11, t22, negs, negcs;
    gr_ptr M, ML, MR;
    int status = GR_SUCCESS;

    GR_TMP_INIT5(c, s, cc, cs, v, ctx);
    GR_TMP_INIT5(t1, t2, t11, t22, negs, ctx);
    GR_TMP_INIT(negcs, ctx);

    M = GR_TMP_ALLOC(8 * sz);
    ML = M;
    MR = GR_ENTRY(M, 4, sz);

    status |= gr_set(t11, GR_MAT_ENTRY(T, k - 1, k - 1, sz), ctx);
    status |= gr_set(t22, GR_MAT_ENTRY(T, k, k, sz), ctx);

    /* c = T[k-1,k],  s = T[k,k] - T[k-1,k-1] */
    status |= gr_set(c, GR_MAT_ENTRY(T, k - 1, k, sz), ctx);
    status |= gr_sub(s, t22, t11, ctx);

    status |= gr_abs2(t1, c, ctx);
    status |= gr_abs2(t2, s, ctx);
    status |= gr_add(v, t1, t2, ctx);
    status |= gr_sqrt(v, v, ctx);

    if (gr_is_zero(v, ctx) == T_TRUE)
    {
        /* The two eigenvalues are equal and the off-diagonal entry
           vanishes; there is nothing to do. */
        goto cleanup;
    }

    status |= gr_inv(v, v, ctx);
    status |= gr_mul(c, c, v, ctx);
    status |= gr_mul(s, s, v, ctx);

    status |= gr_conj(cc, c, ctx);
    status |= gr_conj(cs, s, ctx);
    status |= gr_neg(negs, s, ctx);
    status |= gr_neg(negcs, cs, ctx);

    /* G = [[c, -conj(s)], [s, conj(c)]];  rows get G^H, columns get G */
    gr_set_shallow(GR_ENTRY(ML, 0, sz), cc, ctx);
    gr_set_shallow(GR_ENTRY(ML, 1, sz), cs, ctx);
    gr_set_shallow(GR_ENTRY(ML, 2, sz), negs, ctx);
    gr_set_shallow(GR_ENTRY(ML, 3, sz), c, ctx);
    gr_set_shallow(GR_ENTRY(MR, 0, sz), c, ctx);
    gr_set_shallow(GR_ENTRY(MR, 1, sz), s, ctx);
    gr_set_shallow(GR_ENTRY(MR, 2, sz), negcs, ctx);
    gr_set_shallow(GR_ENTRY(MR, 3, sz), cc, ctx);

    /* Rows k-1, k are zero left of column k-1; columns k-1, k are zero
       below row k. */
    status |= gr_mat_apply_2x2_rows(T, k - 1, k, k - 1, n, ML, ctx);
    status |= gr_mat_apply_2x2_cols(T, k - 1, k, 0, k + 1, MR, ctx);

    if (V != NULL)
        status |= gr_mat_apply_2x2_cols(V, k - 1, k, 0, gr_mat_nrows(V, ctx), MR, ctx);

    /* Enforce the exact result rather than relying on cancellation. */
    status |= gr_set(GR_MAT_ENTRY(T, k - 1, k - 1, sz), t22, ctx);
    status |= gr_set(GR_MAT_ENTRY(T, k, k, sz), t11, ctx);
    status |= gr_zero(GR_MAT_ENTRY(T, k, k - 1, sz), ctx);

cleanup:
    GR_TMP_CLEAR5(c, s, cc, cs, v, ctx);
    GR_TMP_CLEAR5(t1, t2, t11, t22, negs, ctx);
    GR_TMP_CLEAR(negcs, ctx);
    GR_TMP_FREE(M, 8 * sz);

    return status;
}

/*
    Moves the eigenvalue at diagonal position ifst to position ilst,
    shifting the eigenvalues in between by one.  T must be upper
    triangular; V, if not NULL, is postmultiplied by the accumulated
    unitary factor.
*/
int
gr_mat_schur_reorder(gr_mat_t T, gr_mat_t V, slong ifst, slong ilst, gr_ctx_t ctx)
{
    slong k, n;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(T, ctx);

    if (n != gr_mat_ncols(T, ctx))
        return GR_DOMAIN;

    if (ifst < 0 || ifst >= n || ilst < 0 || ilst >= n)
        return GR_DOMAIN;

    if (V != NULL && gr_mat_ncols(V, ctx) != n)
        return GR_DOMAIN;

    if (ifst > ilst)
    {
        for (k = ifst; k > ilst; k--)
            status |= gr_mat_schur_swap_adjacent(T, V, k, ctx);
    }
    else
    {
        for (k = ifst + 1; k <= ilst; k++)
            status |= gr_mat_schur_swap_adjacent(T, V, k, ctx);
    }

    return status;
}
